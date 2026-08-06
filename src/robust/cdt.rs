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

    cdt.tris
        .iter()
        .filter(|t| t.alive)
        .map(|t| t.v)
        .collect()
}

impl Cdt {
    /// orient2d with the approx filter first, exact fallback. All CDT sign
    /// tests funnel through here (and `nondelaunay`) so generic-position
    /// queries never touch BigInt.
    #[inline]
    fn o2(&self, i: usize, j: usize, k: usize) -> Sign {
        orient2d_a(self.apts[i], self.apts[j], self.apts[k])
            .unwrap_or_else(|| orient2d_h(&self.pts[i], &self.pts[j], &self.pts[k]))
    }

    /// Strict-incircle with the approx filter first, exact fallback.
    #[inline]
    fn nondelaunay(&self, tri: [usize; 3], d: usize) -> bool {
        match incircle_a(
            self.apts[tri[0]],
            self.apts[tri[1]],
            self.apts[tri[2]],
            self.apts[d],
        ) {
            Some(s) => s == Sign::Pos,
            None => {
                incircle_h(
                    &self.pts[tri[0]],
                    &self.pts[tri[1]],
                    &self.pts[tri[2]],
                    &self.pts[d],
                ) == Sign::Pos
            }
        }
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
    /// reached. Exact incircle ties (cocircular quads) never flip, which
    /// guarantees termination.
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

    /// Flip edge `e` of `t` if it is internal, unconstrained, strictly
    /// non-Delaunay, and its quad is strictly convex. Returns whether a flip
    /// happened.
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
        // Strict Delaunay violation?
        if !self.nondelaunay(self.tris[t].v, d) {
            return false;
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

    /// Recover constraint edge (a,b) by flipping away crossing edges
    /// (Anglada's algorithm), then mark it constrained on both sides.
    /// Precondition (from robust/arrangement.rs): no vertex lies strictly
    /// inside segment (a,b) and no other constraint properly crosses it.
    fn insert_constraint(&mut self, a: usize, b: usize) {
        let mut guard = 0usize;
        while !self.mark_if_present(a, b) {
            guard += 1;
            assert!(
                guard <= 4 * self.tris.len() * self.tris.len() + 64,
                "constraint recovery failed to converge — precondition violated"
            );
            let (t, e) = self
                .find_flippable_crossing(a, b)
                .expect("a crossing edge with a convex quad must exist");
            let n = self.tris[t].adj[e] as usize;
            let ne = self.shared_edge(n, t);
            self.flip(t, e, n, ne);
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

    /// Is the internal edge `e` of `t` flippable (strictly convex quad)?
    fn flippable(&self, t: usize, e: usize) -> bool {
        let n = self.tris[t].adj[e];
        if n < 0 {
            return false;
        }
        let ne = self.shared_edge(n as usize, t);
        let u = self.tris[t].v[e];
        let v = self.tris[t].v[(e + 1) % 3];
        let c = self.tris[t].v[(e + 2) % 3];
        let d = self.tris[n as usize].v[(ne + 2) % 3];
        self.o2(u, d, c) == Sign::Pos && self.o2(d, v, c) == Sign::Pos
    }

    /// Find an edge strictly crossing segment (a,b) whose quad is strictly
    /// convex (flippable). Anglada's lemma guarantees one exists while any
    /// crossing remains.
    /// Walk the triangles pierced by segment (a,b), returning the first
    /// crossing edge whose quad is strictly convex (flippable). Anglada's
    /// lemma guarantees one exists while any crossing remains; the walk
    /// visits exactly the crossing edges — O(crossings), not O(triangles).
    fn find_flippable_crossing(&self, a: usize, b: usize) -> Option<(usize, usize)> {
        // Entry: the incident triangle of `a` whose opposite edge the
        // segment leaves through. (No vertex lies strictly inside (a,b) —
        // arrangement precondition — so the exit is a strict edge crossing.)
        let entry = self.rotate_around(a, |t| {
            let ia = (0..3).position(|k| self.tris[t].v[k] == a).unwrap();
            let e_opp = (ia + 1) % 3;
            let u = self.tris[t].v[e_opp];
            let v = self.tris[t].v[(e_opp + 1) % 3];
            (u != b && v != b && self.strictly_crosses(a, b, u, v)).then_some((t, e_opp))
        });
        let Some((mut t, mut e)) = entry else {
            return self.find_flippable_crossing_scan(a, b);
        };
        loop {
            debug_assert!(
                !self.tris[t].con[e],
                "constraint crosses an existing constraint — arrangement bug"
            );
            if self.flippable(t, e) {
                return Some((t, e));
            }
            // March into the neighbor; the segment exits it through one of
            // the two remaining edges (or ends at b).
            let n = self.tris[t].adj[e];
            if n < 0 {
                return self.find_flippable_crossing_scan(a, b);
            }
            let n = n as usize;
            let ne = self.shared_edge(n, t);
            let mut next = None;
            for k in 1..3 {
                let e2 = (ne + k) % 3;
                let u = self.tris[n].v[e2];
                let v = self.tris[n].v[(e2 + 1) % 3];
                if u == a || u == b || v == a || v == b {
                    continue;
                }
                if self.strictly_crosses(a, b, u, v) {
                    next = Some(e2);
                    break;
                }
            }
            match next {
                Some(e2) => {
                    t = n;
                    e = e2;
                }
                // Walk exhausted without a flippable edge (collinear corner
                // cases can end it early) — the global scan is the slow but
                // complete answer, and behaves exactly like the pre-walk code.
                None => return self.find_flippable_crossing_scan(a, b),
            }
        }
    }

    /// The original exhaustive search — fallback when the segment walk
    /// terminates without an answer.
    fn find_flippable_crossing_scan(&self, a: usize, b: usize) -> Option<(usize, usize)> {
        for t in 0..self.tris.len() {
            if !self.tris[t].alive {
                continue;
            }
            for e in 0..3 {
                let u = self.tris[t].v[e];
                let v = self.tris[t].v[(e + 1) % 3];
                if u == a || u == b || v == a || v == b {
                    continue;
                }
                if self.strictly_crosses(a, b, u, v) && self.flippable(t, e) {
                    return Some((t, e));
                }
            }
        }
        None
    }
}

/// Delaunay test on raw points: is the fourth point strictly inside the
/// circumcircle of the (CCW) triangle? Kept for the validation in
/// robust/cdt_tests.rs.
#[cfg(test)]
pub(super) fn is_strictly_non_delaunay(pts: &[R2], tri: [usize; 3], d: usize) -> bool {
    incircle_h(
        &homog2_of(&pts[tri[0]]),
        &homog2_of(&pts[tri[1]]),
        &homog2_of(&pts[tri[2]]),
        &homog2_of(&pts[d]),
    ) == Sign::Pos
}

#[cfg(test)]
#[path = "cdt_tests.rs"]
mod tests;
