// robust/cdt.rs — Constrained Delaunay triangulation on exact 2D points.
//
// Triangulates the interior of one input triangle after robust/arrangement.rs
// has collected every intersection point and constraint segment that lands on
// it (paper §6.4). All predicates are exact (robust/exact), so there are no
// epsilon decisions anywhere: point location, edge flips, and constraint
// recovery all reason on true signs. Delaunay-ness itself is not load-bearing
// for the boolean — validity and constraint preservation are — but Delaunay
// flips avoid the near-degenerate sub-triangles that would otherwise stress
// downstream float code.
//
// Scale note: an arrangement holds one mesh triangle's intersections, so
// point counts are small (tens). The implementation favors obviously-correct
// O(n²) scans (linear point location, fixpoint legalization) over
// sophisticated locality structures.

use super::exact::approx::{incircle_a, orient2d_a};
use super::exact::predicates::{homog2_of, incircle_h, orient2d_h, Homog2, TriLoc};
use super::exact::rational::{rat_to_f64, R2};
use super::exact::Sign;

/// One triangle of the CDT. Vertices are CCW; edge i runs v[i] → v[i+1],
/// adj[i] is the triangle index across that edge (-1 on the hull), con[i]
/// marks a constrained edge that flips must never cross.
#[derive(Clone, Debug)]
struct Tri {
    v: [usize; 3],
    adj: [i32; 3],
    con: [bool; 3],
    alive: bool,
}

struct Cdt {
    /// Homogenized once per triangulation; exact predicates reuse these.
    pts: Vec<Homog2>,
    /// Correctly rounded f64 approximations (relative error ≤ ε), computed
    /// once; the semi-static filters in exact/approx.rs certify most
    /// predicate signs from these alone, escalating to `pts` only on
    /// near-degeneracies.
    apts: Vec<[f64; 2]>,
    tris: Vec<Tri>,
    /// Suspect triangles for queue-based Lawson legalization.
    suspects: Vec<usize>,
    /// For each vertex, the most recently created triangle containing it.
    /// Invariant: always alive — every operation that kills a triangle
    /// containing v pushes (and records) a replacement containing v first.
    /// Constraint recovery uses this to rotate around a vertex instead of
    /// scanning the whole triangulation.
    vert_tri: Vec<u32>,
    /// Rank of each point in exact lexicographic (x, y) coordinate order.
    /// Cocircular-tie flips key on this so the triangulation is a function
    /// of the point coordinates, not of construction history — coincident
    /// coplanar triangles must tile their shared region identically or the
    /// cell complex sees one physical sheet crossing as several
    /// vi-distinct walls with understated winding steps.
    rank: Vec<u32>,
}

/// Constrained Delaunay triangulation of `points` inside the triangle formed
/// by `points[0..3]` (any winding). Every other point must lie inside that
/// triangle or on its boundary, all points must be pairwise distinct, and
/// every constraint (index pair) must be free of interior points and proper
/// crossings with other constraints — robust/arrangement.rs guarantees all
/// three. Returns CCW index triangles exactly covering the input triangle.
pub fn triangulate(points: &[R2], constraints: &[(usize, usize)]) -> Vec<[usize; 3]> {
    assert!(points.len() >= 3, "need the three corner points");
    let hom: Vec<Homog2> = points.iter().map(homog2_of).collect();
    let mut corners = [0usize, 1, 2];
    let orient = orient2d_h(&hom[0], &hom[1], &hom[2]);
    assert!(
        orient != Sign::Zero,
        "degenerate base triangle in CDT input"
    );
    if orient == Sign::Neg {
        corners.swap(1, 2);
    }

    let apts: Vec<[f64; 2]> = points
        .iter()
        .map(|p| [rat_to_f64(&p.x), rat_to_f64(&p.y)])
        .collect();
    let mut by_coord: Vec<usize> = (0..points.len()).collect();
    by_coord.sort_by(|&i, &j| {
        points[i]
            .x
            .cmp(&points[j].x)
            .then_with(|| points[i].y.cmp(&points[j].y))
    });
    let mut rank = vec![0u32; points.len()];
    for (r, &i) in by_coord.iter().enumerate() {
        rank[i] = r as u32;
    }
    let mut cdt = Cdt {
        pts: hom,
        apts,
        tris: vec![Tri {
            v: corners,
            adj: [-1, -1, -1],
            con: [false, false, false],
            alive: true,
        }],
        suspects: Vec::new(),
        vert_tri: vec![0; points.len()],
        rank,
    };

    for p in 3..cdt.pts.len() {
        let first_new = cdt.tris.len();
        cdt.insert_point(p);
        cdt.seed_suspects(first_new);
        cdt.legalize_suspects();
    }
    for &(a, b) in constraints {
        debug_assert_ne!(a, b, "zero-length constraint");
        let first_new = cdt.tris.len();
        cdt.insert_constraint(a, b);
        cdt.seed_suspects(first_new);
        cdt.legalize_suspects();
    }

    cdt.tris.iter().filter(|t| t.alive).map(|t| t.v).collect()
}

impl Cdt {
    /// orient2d with the approx filter first, exact fallback. All CDT sign
    /// tests funnel through here (and `nondelaunay`) so generic-position
    /// queries never touch the bignum tier.
    #[inline]
    fn o2(&self, i: usize, j: usize, k: usize) -> Sign {
        orient2d_a(self.apts[i], self.apts[j], self.apts[k])
            .unwrap_or_else(|| orient2d_h(&self.pts[i], &self.pts[j], &self.pts[k]))
    }

    /// Incircle sign with the approx filter first, exact fallback. `Pos`
    /// means `d` is strictly inside the circumcircle of the (CCW) triangle;
    /// `Zero` is an exact cocircular tie.
    #[inline]
    fn incircle(&self, tri: [usize; 3], d: usize) -> Sign {
        match incircle_a(
            self.apts[tri[0]],
            self.apts[tri[1]],
            self.apts[tri[2]],
            self.apts[d],
        ) {
            Some(s) => s,
            None => incircle_h(
                &self.pts[tri[0]],
                &self.pts[tri[1]],
                &self.pts[tri[2]],
                &self.pts[d],
            ),
        }
    }

    /// Strict Delaunay violation (see [`Cdt::incircle`]).
    #[inline]
    fn nondelaunay(&self, tri: [usize; 3], d: usize) -> bool {
        self.incircle(tri, d) == Sign::Pos
    }

    /// Canonical identity of the diagonal {u, v}: the pair of coordinate
    /// ranks, low first. Ranks order points by exact coordinates, so two
    /// coincident coplanar triangles — which project into the same
    /// dominant-axis plane and therefore hand this CDT the same coordinates
    /// for shared points — agree on this key even though their local point
    /// indices differ.
    #[inline]
    fn diag_key(&self, u: usize, v: usize) -> (u32, u32) {
        let (ru, rv) = (self.rank[u], self.rank[v]);
        (ru.min(rv), ru.max(rv))
    }

    /// `point_in_tri_2d` over the filtered predicate (same TriLoc semantics
    /// as `predicates::point_in_tri_2d_h`; triangle is CCW by construction).
    fn loc_in_tri(&self, p: usize, v: [usize; 3]) -> TriLoc {
        let s0 = self.o2(v[0], v[1], p);
        let s1 = self.o2(v[1], v[2], p);
        let s2 = self.o2(v[2], v[0], p);
        if s0 == Sign::Neg || s1 == Sign::Neg || s2 == Sign::Neg {
            return TriLoc::Outside;
        }
        match (s0 == Sign::Zero, s1 == Sign::Zero, s2 == Sign::Zero) {
            (false, false, false) => TriLoc::Inside,
            (true, false, false) => TriLoc::OnEdge(0),
            (false, true, false) => TriLoc::OnEdge(1),
            (false, false, true) => TriLoc::OnEdge(2),
            (true, false, true) => TriLoc::OnVertex(0),
            (true, true, false) => TriLoc::OnVertex(1),
            (false, true, true) => TriLoc::OnVertex(2),
            (true, true, true) => TriLoc::Outside,
        }
    }

    /// Record `t` as the latest triangle containing each of its vertices.
    #[inline]
    fn record(&mut self, t: usize) {
        for v in self.tris[t].v {
            self.vert_tri[v] = t as u32;
        }
    }

    /// A live triangle containing vertex `a`. The `vert_tri` invariant makes
    /// the recorded entry always alive; the scan is a defensive fallback.
    fn live_tri_with(&self, a: usize) -> usize {
        let cand = self.vert_tri[a] as usize;
        if cand < self.tris.len() && self.tris[cand].alive && self.tris[cand].v.contains(&a) {
            return cand;
        }
        debug_assert!(false, "vert_tri invariant broken for vertex {a}");
        (0..self.tris.len())
            .rev()
            .find(|&i| self.tris[i].alive && self.tris[i].v.contains(&a))
            .expect("vertex not in any live triangle")
    }

    /// Rotate around vertex `a` (both directions from its recorded triangle,
    /// so boundary fans are fully covered) applying `f` to each incident
    /// triangle until it returns Some.
    fn rotate_around<T>(&self, a: usize, mut f: impl FnMut(usize) -> Option<T>) -> Option<T> {
        let seed = self.live_tri_with(a);
        for dirn in 0..2 {
            let mut cur = seed;
            loop {
                if dirn == 0 || cur != seed {
                    if let Some(out) = f(cur) {
                        return Some(out);
                    }
                }
                let t = &self.tris[cur];
                let ia = (0..3)
                    .position(|k| t.v[k] == a)
                    .expect("rotation left the vertex fan");
                // dir 0 crosses the edge leaving `a`; dir 1 the edge entering.
                let e_step = if dirn == 0 { ia } else { (ia + 2) % 3 };
                let n = t.adj[e_step];
                if n < 0 || n as usize == seed {
                    break;
                }
                cur = n as usize;
            }
        }
        None
    }

    /// Index of `t`'s edge shared with triangle `n` (panics if not adjacent —
    /// an internal-consistency violation, not an input condition).
    fn shared_edge(&self, t: usize, n: usize) -> usize {
        (0..3)
            .find(|&e| self.tris[t].adj[e] == n as i32)
            .expect("adjacency tables out of sync")
    }

    /// Point back the neighbor across `(t, e)` (if any) to `new_t`.
    fn rewire(&mut self, t: usize, e: usize, new_t: usize) {
        let n = self.tris[t].adj[e];
        if n >= 0 {
            let ne = self.shared_edge(n as usize, t);
            self.tris[n as usize].adj[ne] = new_t as i32;
        }
    }

    fn insert_point(&mut self, p: usize) {
        let (t, loc) = self.locate(p);
        match loc {
            TriLoc::Inside => self.split_interior(t, p),
            TriLoc::OnEdge(e) => self.split_edge(t, e as usize, p),
            TriLoc::OnVertex(_) => {
                debug_assert!(false, "duplicate point reached CDT: {p}");
            }
            TriLoc::Outside => unreachable!("point outside base triangle"),
        }
    }

    fn locate(&self, p: usize) -> (usize, TriLoc) {
        // Visibility walk from the most recent live triangle: step across any
        // edge whose line strictly separates the point (points are inserted
        // before constraints, so the triangulation is Delaunay here and the
        // walk terminates — Edelsbrunner). The step cap is a safety net; on
        // overrun the exhaustive scan below still answers.
        if let Some(start) = (0..self.tris.len()).rev().find(|&i| self.tris[i].alive) {
            let mut cur = start;
            let mut steps = 0usize;
            let cap = 4 * self.tris.len() + 16;
            'walk: while steps < cap {
                steps += 1;
                let t = &self.tris[cur];
                let loc = self.loc_in_tri(p, t.v);
                if loc != TriLoc::Outside {
                    return (cur, loc);
                }
                for e in 0..3 {
                    // CCW triangle: strictly outside edge e ⇔ orient2d neg.
                    if self.o2(t.v[e], t.v[(e + 1) % 3], p) == Sign::Neg {
                        let n = t.adj[e];
                        if n >= 0 {
                            cur = n as usize;
                            continue 'walk;
                        }
                    }
                }
                break; // no separating edge with a neighbor — fall through
            }
        }
        for (i, t) in self.tris.iter().enumerate() {
            if !t.alive {
                continue;
            }
            if self.loc_in_tri(p, t.v) != TriLoc::Outside {
                return (i, self.loc_in_tri(p, t.v));
            }
        }
        unreachable!("point {p} not located in any live triangle");
    }

    /// Split triangle `t` = (a,b,c) at interior point `p` into (a,b,p),
    /// (b,c,p), (c,a,p).
    fn split_interior(&mut self, t: usize, p: usize) {
        let Tri { v, adj, con, .. } = self.tris[t].clone();
        let base = self.tris.len();
        // New triangle k (k=0,1,2) covers edge k of the old triangle; its
        // inner edges connect to the cyclic neighbors.
        for k in 0..3 {
            let next = base + (k + 1) % 3;
            let prev = base + (k + 2) % 3;
            self.tris.push(Tri {
                v: [v[k], v[(k + 1) % 3], p],
                adj: [adj[k], next as i32, prev as i32],
                con: [con[k], false, false],
                alive: true,
            });
        }
        for k in 0..3 {
            self.rewire(t, k, base + k);
            self.record(base + k);
        }
        self.tris[t].alive = false;
    }

    /// Split edge `e` of triangle `t` at point `p` lying exactly on it —
    /// two new triangles per side of the edge.
    fn split_edge(&mut self, t: usize, e: usize, p: usize) {
        let n = self.tris[t].adj[e];
        self.split_edge_one_side(t, e, p);
        if n >= 0 {
            let ne = self.shared_edge(n as usize, t);
            // split_edge_one_side left dangling adjacency into the dead
            // triangle t; splitting the neighbor rewires everything below.
            self.split_edge_one_side(n as usize, ne, p);
            // The four halves now flank the split point: wire them to each
            // other. The two t-halves still point at n and vice versa.
            self.fix_split_pairs(p);
        }
    }

    /// Replace `t` = (a,b,c) with (a,p,c) and (p,b,c), where p lies on edge
    /// e=(a,b). The two halves keep the old adjacency indices toward the
    /// (possibly split) neighbor across e; `fix_split_pairs` repairs those
    /// when both sides exist.
    fn split_edge_one_side(&mut self, t: usize, e: usize, p: usize) {
        let Tri { v, adj, con, .. } = self.tris[t].clone();
        let (a, b, c) = (v[e], v[(e + 1) % 3], v[(e + 2) % 3]);
        let (adj_ab, adj_bc, adj_ca) = (adj[e], adj[(e + 1) % 3], adj[(e + 2) % 3]);
        let (con_ab, con_bc, con_ca) = (con[e], con[(e + 1) % 3], con[(e + 2) % 3]);
        let t1 = self.tris.len();
        let t2 = t1 + 1;
        // (a, p, c): edges (a,p) toward old neighbor, (p,c) inner, (c,a) old.
        self.tris.push(Tri {
            v: [a, p, c],
            adj: [adj_ab, t2 as i32, adj_ca],
            con: [con_ab, false, con_ca],
            alive: true,
        });
        // (p, b, c): edges (p,b) toward old neighbor, (b,c) old, (c,p) inner.
        self.tris.push(Tri {
            v: [p, b, c],
            adj: [adj_ab, adj_bc, t1 as i32],
            con: [con_ab, con_bc, false],
            alive: true,
        });
        self.rewire(t, (e + 1) % 3, t2);
        self.rewire(t, (e + 2) % 3, t1);
        self.record(t1);
        self.record(t2);
        self.tris[t].alive = false;
    }

    /// After both sides of a split edge have been divided, adjacency across
    /// the split still names the two dead triangles. Match the four live
    /// halves incident to `p` whose cross-split edges contain `p`, pairing
    /// edges with identical endpoint sets.
    fn fix_split_pairs(&mut self, p: usize) {
        // Collect (tri, edge) pairs whose adjacency points at a dead triangle.
        let mut dangling: Vec<(usize, usize)> = Vec::new();
        for i in 0..self.tris.len() {
            if !self.tris[i].alive {
                continue;
            }
            for e in 0..3 {
                let n = self.tris[i].adj[e];
                if n >= 0 && !self.tris[n as usize].alive {
                    debug_assert!(
                        self.tris[i].v[e] == p || self.tris[i].v[(e + 1) % 3] == p,
                        "dangling edge must touch the split point"
                    );
                    dangling.push((i, e));
                }
            }
        }
        for i in 0..dangling.len() {
            for j in (i + 1)..dangling.len() {
                let (ti, ei) = dangling[i];
                let (tj, ej) = dangling[j];
                let vi = [self.tris[ti].v[ei], self.tris[ti].v[(ei + 1) % 3]];
                let vj = [self.tris[tj].v[ej], self.tris[tj].v[(ej + 1) % 3]];
                if vi[0] == vj[1] && vi[1] == vj[0] {
                    self.tris[ti].adj[ei] = tj as i32;
                    self.tris[tj].adj[ej] = ti as i32;
                }
            }
        }
    }

    /// Enqueue every triangle from `first_new` onward (the ones an insert or
    /// constraint pass just created) as legalization suspects.
    fn seed_suspects(&mut self, first_new: usize) {
        for t in first_new..self.tris.len() {
            self.suspects.push(t);
        }
    }

    /// Queue-based Lawson legalization: flip every non-constrained, strictly
    /// non-Delaunay, flippable edge reachable from the suspect set. A flip
    /// enqueues the two replacement triangles, and their neighbors get
    /// re-tested through the shared edges when those triangles are examined,
    /// so the worklist reaches everything the old global fixpoint rescan
    /// reached. Termination is argued at [`Cdt::try_flip`].
    fn legalize_suspects(&mut self) {
        while let Some(t) = self.suspects.pop() {
            if t >= self.tris.len() || !self.tris[t].alive {
                continue;
            }
            for e in 0..3 {
                if self.try_flip(t, e) {
                    // t is dead; its replacements are the last two triangles,
                    // and the flipped neighbor's replacement edges are on them.
                    let n = self.tris.len();
                    self.suspects.push(n - 2);
                    self.suspects.push(n - 1);
                    break;
                }
            }
        }
    }

    /// Flip edge `e` of `t` if it is internal, unconstrained, its quad is
    /// strictly convex, and the flip either repairs a strict Delaunay
    /// violation or resolves an exact cocircular tie toward the canonical
    /// (lowest [`Cdt::diag_key`]) diagonal. Returns whether a flip happened.
    ///
    /// Tie flips make the triangulation of a cocircular quad a function of
    /// its point coordinates instead of its construction history, so
    /// coincident coplanar triangles tile shared regions identically (found
    /// on Thingi10K #36088, where a doubled sheet's two copies picked
    /// opposite diagonals of an exact square and the cell complex saw four
    /// [1,0]-step walls where the geometry has [2,0]-step stacks).
    ///
    /// Termination: a strict flip strictly lowers the lifted-paraboloid
    /// volume; a tie flip keeps it equal (that is what cocircular means for
    /// the lift) while strictly lowering the multiset of diagonal keys. The
    /// pair (lift volume, key multiset) therefore strictly decreases
    /// lexicographically on every flip.
    fn try_flip(&mut self, t: usize, e: usize) -> bool {
        if self.tris[t].con[e] {
            return false;
        }
        let n = self.tris[t].adj[e];
        if n < 0 {
            return false;
        }
        let n = n as usize;
        let ne = self.shared_edge(n, t);
        let (a, b, c) = (
            self.tris[t].v[e],
            self.tris[t].v[(e + 1) % 3],
            self.tris[t].v[(e + 2) % 3],
        );
        let d = self.tris[n].v[(ne + 2) % 3];
        match self.incircle(self.tris[t].v, d) {
            Sign::Neg => return false,
            Sign::Zero => {
                // Cocircular: flip only toward the canonical diagonal.
                if self.diag_key(c, d) >= self.diag_key(a, b) {
                    return false;
                }
            }
            Sign::Pos => {}
        }
        // Strictly convex quad a-d-b-c (new triangles must be CCW)?
        if self.o2(a, d, c) != Sign::Pos || self.o2(d, b, c) != Sign::Pos {
            return false;
        }
        self.flip(t, e, n, ne);
        true
    }

    /// Replace triangles t=(a,b,c) and n=(b,a,d) with (a,d,c) and (d,b,c).
    fn flip(&mut self, t: usize, e: usize, n: usize, ne: usize) {
        let (a, b, c) = (
            self.tris[t].v[e],
            self.tris[t].v[(e + 1) % 3],
            self.tris[t].v[(e + 2) % 3],
        );
        let d = self.tris[n].v[(ne + 2) % 3];
        debug_assert_eq!(self.tris[n].v[ne], b);
        debug_assert_eq!(self.tris[n].v[(ne + 1) % 3], a);

        let adj_bc = self.tris[t].adj[(e + 1) % 3];
        let con_bc = self.tris[t].con[(e + 1) % 3];
        let adj_ca = self.tris[t].adj[(e + 2) % 3];
        let con_ca = self.tris[t].con[(e + 2) % 3];
        let adj_ad = self.tris[n].adj[(ne + 1) % 3];
        let con_ad = self.tris[n].con[(ne + 1) % 3];
        let adj_db = self.tris[n].adj[(ne + 2) % 3];
        let con_db = self.tris[n].con[(ne + 2) % 3];

        let t1 = self.tris.len();
        let t2 = t1 + 1;
        // (a, d, c): edges (a,d) old, (d,c) diagonal, (c,a) old.
        self.tris.push(Tri {
            v: [a, d, c],
            adj: [adj_ad, t2 as i32, adj_ca],
            con: [con_ad, false, con_ca],
            alive: true,
        });
        // (d, b, c): edges (d,b) old, (b,c) old, (c,d) diagonal.
        self.tris.push(Tri {
            v: [d, b, c],
            adj: [adj_db, adj_bc, t1 as i32],
            con: [con_db, con_bc, false],
            alive: true,
        });
        self.rewire(t, (e + 1) % 3, t2);
        self.rewire(t, (e + 2) % 3, t1);
        self.rewire(n, (ne + 1) % 3, t1);
        self.rewire(n, (ne + 2) % 3, t2);
        self.record(t1);
        self.record(t2);
        self.tris[t].alive = false;
        self.tris[n].alive = false;
    }

    /// Recover constraint edge (a,b) by Anglada's cavity retriangulation:
    /// walk the corridor of triangles the open segment pierces, delete them,
    /// and re-triangulate the two pseudo-polygon halves flanking the segment
    /// (each necessarily gains edge (a,b)). Edge-flip recovery — flipping
    /// away crossing edges one at a time — is NOT used because it can cycle:
    /// when a flipped diagonal still crosses the segment, the next search can
    /// select it and flip it straight back (observed on Thingi10K 1075458 −
    /// 91115, where two edges alternated forever and the triangle arena grew
    /// without bound). The corridor walk is bounded by the live triangle
    /// count and the retriangulation recursion by the chain length, so this
    /// terminates unconditionally.
    /// Precondition (from robust/arrangement.rs): no vertex lies strictly
    /// inside segment (a,b) and no other constraint properly crosses it.
    fn insert_constraint(&mut self, a: usize, b: usize) {
        if self.mark_if_present(a, b) {
            return;
        }
        let (crossed, left, right) = self.collect_corridor(a, b);
        self.retriangulate_corridor(a, b, &crossed, &left, &right);
        assert!(
            self.mark_if_present(a, b),
            "cavity retriangulation did not produce constraint edge ({a},{b})"
        );
    }

    /// The corridor of triangles pierced by open segment (a,b): the triangle
    /// list in walk order plus the chains of corridor vertices strictly left
    /// and strictly right of the directed line a→b, each ordered from a to b.
    /// Chains may repeat a vertex (pseudo-polygon pinch); consecutive
    /// duplicates cannot occur.
    fn collect_corridor(&self, a: usize, b: usize) -> (Vec<usize>, Vec<usize>, Vec<usize>) {
        // Entry: the incident triangle of `a` whose opposite edge the segment
        // leaves through. No vertex lies strictly inside (a,b) and (a,b) is
        // not an edge (mark_if_present said so), so the exit is a strict edge
        // crossing.
        let (t0, e0) = self
            .rotate_around(a, |t| {
                let ia = (0..3).position(|k| self.tris[t].v[k] == a).unwrap();
                let e_opp = (ia + 1) % 3;
                let u = self.tris[t].v[e_opp];
                let v = self.tris[t].v[(e_opp + 1) % 3];
                (u != b && v != b && self.strictly_crosses(a, b, u, v)).then_some((t, e_opp))
            })
            .expect("segment leaves its start vertex through a crossed edge");
        // CCW triangle (a, u, v): u is strictly right of a→b, v strictly left.
        let u = self.tris[t0].v[e0];
        let v = self.tris[t0].v[(e0 + 1) % 3];
        let mut crossed = vec![t0];
        let mut left = vec![v];
        let mut right = vec![u];
        let (mut t, mut e) = (t0, e0);
        loop {
            assert!(
                crossed.len() <= self.tris.len(),
                "corridor walk revisited a triangle — triangulation corrupt"
            );
            let n = self.tris[t].adj[e];
            assert!(
                n >= 0,
                "segment crossed the hull — point outside base triangle"
            );
            let n = n as usize;
            let ne = self.shared_edge(n, t);
            crossed.push(n);
            let w = self.tris[n].v[(ne + 2) % 3];
            if w == b {
                return (crossed, left, right);
            }
            match self.o2(a, b, w) {
                Sign::Pos => left.push(w),
                Sign::Neg => right.push(w),
                Sign::Zero => panic!(
                    "vertex {w} lies on constraint segment ({a},{b}) — precondition violated"
                ),
            }
            // The segment exits `n` through whichever remaining edge it
            // strictly crosses (exactly one, since w is strictly off-line).
            e = (1..3)
                .map(|k| (ne + k) % 3)
                .find(|&e2| {
                    let p = self.tris[n].v[e2];
                    let q = self.tris[n].v[(e2 + 1) % 3];
                    self.strictly_crosses(a, b, p, q)
                })
                .expect("segment must exit the corridor triangle it entered");
            t = n;
        }
    }

    /// Is `v` strictly inside the circumcircle of triangle (a,b,c), whatever
    /// its orientation? Collinear (a,b,c) has no circumcircle: false.
    fn in_circumcircle(&self, a: usize, b: usize, c: usize, v: usize) -> bool {
        match self.o2(a, b, c) {
            Sign::Pos => self.nondelaunay([a, b, c], v),
            Sign::Neg => self.nondelaunay([a, c, b], v),
            Sign::Zero => false,
        }
    }

    /// Anglada's pseudo-polygon triangulation: the region bounded by directed
    /// base edge (a,b) and `chain` (the boundary path from a to b, region on
    /// the left of a→b). Picks the Delaunay apex c — the chain vertex whose
    /// circumcircle with the base is empty of chain vertices — emits CCW
    /// triangle (a,b,c), and recurses on the two sub-chains.
    fn triangulate_pseudo(&self, a: usize, b: usize, chain: &[usize], out: &mut Vec<[usize; 3]>) {
        if chain.is_empty() {
            return;
        }
        let mut ci = 0usize;
        for i in 1..chain.len() {
            // A collinear "apex" has no circumcircle and can never be the
            // Delaunay apex; any candidate supersedes it.
            let replace = self.o2(a, b, chain[ci]) == Sign::Zero
                || self.in_circumcircle(a, b, chain[ci], chain[i]);
            if replace {
                ci = i;
            }
        }
        let c = chain[ci];
        assert!(
            self.o2(a, b, c) == Sign::Pos,
            "pseudo-polygon apex must lie strictly left of its base"
        );
        self.triangulate_pseudo(a, c, &chain[..ci], out);
        self.triangulate_pseudo(c, b, &chain[ci + 1..], out);
        out.push([a, b, c]);
    }

    /// Replace the corridor triangles with the CDTs of the two pseudo-polygon
    /// halves and wire the new patch into the surrounding triangulation.
    fn retriangulate_corridor(
        &mut self,
        a: usize,
        b: usize,
        crossed: &[usize],
        left: &[usize],
        right: &[usize],
    ) {
        // Fx hashing (unseeded): both maps below are probe-only scratch —
        // `boundary` is looked up per patch edge and `open` is a pairing
        // buffer drained by key, neither is ever iterated.
        use rustc_hash::FxHashMap as HashMap;
        // Cavity boundary: edges of corridor triangles whose neighbor is
        // outside the corridor. Keyed by undirected endpoints (each boundary
        // edge has exactly one inside face, so keys are unique). Interior
        // edges shared by two corridor triangles vanish; a constrained one
        // (possible only at a pseudo-polygon pinch) must be re-marked after
        // the rebuild.
        let mut boundary: HashMap<(usize, usize), (i32, bool)> = HashMap::default();
        let mut interior_con: Vec<(usize, usize)> = Vec::new();
        for &t in crossed {
            for e in 0..3 {
                let n = self.tris[t].adj[e];
                let u = self.tris[t].v[e];
                let v = self.tris[t].v[(e + 1) % 3];
                let key = (u.min(v), u.max(v));
                if n < 0 || !crossed.contains(&(n as usize)) {
                    boundary.insert(key, (n, self.tris[t].con[e]));
                } else if self.tris[t].con[e] {
                    interior_con.push(key);
                }
            }
        }

        // Left half: region left of a→b, chain ordered a to b. Right half:
        // region left of b→a, so its chain must run b to a.
        let mut new_tris: Vec<[usize; 3]> = Vec::new();
        self.triangulate_pseudo(a, b, left, &mut new_tris);
        let right_rev: Vec<usize> = right.iter().rev().copied().collect();
        self.triangulate_pseudo(b, a, &right_rev, &mut new_tris);

        // Push the patch, then wire adjacency: boundary edges reconnect to
        // the outside; the rest pair up patch-internally (including (a,b)
        // between the two halves).
        let base = self.tris.len();
        for &tv in &new_tris {
            self.tris.push(Tri {
                v: tv,
                adj: [-1, -1, -1],
                con: [false, false, false],
                alive: true,
            });
        }
        let mut open: HashMap<(usize, usize), (usize, usize)> = HashMap::default();
        for t in base..self.tris.len() {
            for e in 0..3 {
                let u = self.tris[t].v[e];
                let v = self.tris[t].v[(e + 1) % 3];
                let key = (u.min(v), u.max(v));
                if let Some(&(bn, bcon)) = boundary.get(&key) {
                    self.tris[t].adj[e] = bn;
                    self.tris[t].con[e] = bcon;
                    if bn >= 0 {
                        let n = bn as usize;
                        let ne = (0..3)
                            .find(|&k| {
                                let nu = self.tris[n].v[k];
                                let nv = self.tris[n].v[(k + 1) % 3];
                                (nu.min(nv), nu.max(nv)) == key
                            })
                            .expect("cavity boundary neighbor lost its edge");
                        self.tris[n].adj[ne] = t as i32;
                    }
                } else if let Some((t2, e2)) = open.remove(&key) {
                    self.tris[t].adj[e] = t2 as i32;
                    self.tris[t2].adj[e2] = t as i32;
                } else {
                    open.insert(key, (t, e));
                }
            }
        }
        assert!(
            open.is_empty(),
            "cavity retriangulation left unmatched patch edges"
        );
        // vert_tri invariant: record replacements before killing the corridor.
        for t in base..self.tris.len() {
            self.record(t);
        }
        for &t in crossed {
            self.tris[t].alive = false;
        }
        for (u, v) in interior_con {
            assert!(
                self.mark_if_present(u, v),
                "constrained pinch edge ({u},{v}) was not reproduced by the rebuild"
            );
        }
    }

    /// If edge (a,b) exists, set its constrained flag on both sides and
    /// return true. Rotation around `a` visits every incident triangle, so
    /// absence there is definitive — no global scan.
    fn mark_if_present(&mut self, a: usize, b: usize) -> bool {
        let found = self.rotate_around(a, |t| {
            (0..3).find_map(|e| {
                let u = self.tris[t].v[e];
                let v = self.tris[t].v[(e + 1) % 3];
                ((u, v) == (a, b) || (u, v) == (b, a)).then_some((t, e))
            })
        });
        match found {
            Some((t, e)) => {
                self.tris[t].con[e] = true;
                let n = self.tris[t].adj[e];
                if n >= 0 {
                    let ne = self.shared_edge(n as usize, t);
                    self.tris[n as usize].con[ne] = true;
                }
                true
            }
            None => false,
        }
    }

    /// Strict proper-crossing test of edge (u,v) against segment (a,b).
    fn strictly_crosses(&self, a: usize, b: usize, u: usize, v: usize) -> bool {
        let su = self.o2(a, b, u);
        let sv = self.o2(a, b, v);
        if su == Sign::Zero || sv == Sign::Zero || su == sv {
            return false;
        }
        let sa = self.o2(u, v, a);
        let sb = self.o2(u, v, b);
        sa != Sign::Zero && sb != Sign::Zero && sa != sb
    }
}

#[cfg(test)]
#[path = "cdt_tests.rs"]
mod tests;
