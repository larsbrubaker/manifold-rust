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

use super::exact::predicates::{
    homog2_of, incircle_h, orient2d_h, point_in_tri_2d_h, Homog2, TriLoc,
};
use super::exact::rational::R2;
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
    /// Homogenized once per triangulation; every predicate reuses these.
    pts: Vec<Homog2>,
    tris: Vec<Tri>,
    /// Suspect triangles for queue-based Lawson legalization.
    suspects: Vec<usize>,
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

    let mut cdt = Cdt {
        pts: hom,
        tris: vec![Tri {
            v: corners,
            adj: [-1, -1, -1],
            con: [false, false, false],
            alive: true,
        }],
        suspects: Vec::new(),
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
                let loc = point_in_tri_2d_h(
                    &self.pts[p],
                    &self.pts[t.v[0]],
                    &self.pts[t.v[1]],
                    &self.pts[t.v[2]],
                );
                if loc != TriLoc::Outside {
                    return (cur, loc);
                }
                for e in 0..3 {
                    // CCW triangle: strictly outside edge e ⇔ orient2d neg.
                    if orient2d_h(
                        &self.pts[t.v[e]],
                        &self.pts[t.v[(e + 1) % 3]],
                        &self.pts[p],
                    ) == Sign::Neg
                    {
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
            let loc = point_in_tri_2d_h(
                &self.pts[p],
                &self.pts[t.v[0]],
                &self.pts[t.v[1]],
                &self.pts[t.v[2]],
            );
            if loc != TriLoc::Outside {
                return (i, loc);
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
        if !is_strictly_non_delaunay_h(&self.pts, self.tris[t].v, d) {
            return false;
        }
        // Strictly convex quad a-d-b-c (new triangles must be CCW)?
        if orient2d_h(&self.pts[a], &self.pts[d], &self.pts[c]) != Sign::Pos
            || orient2d_h(&self.pts[d], &self.pts[b], &self.pts[c]) != Sign::Pos
        {
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
    /// return true.
    fn mark_if_present(&mut self, a: usize, b: usize) -> bool {
        for t in 0..self.tris.len() {
            if !self.tris[t].alive {
                continue;
            }
            for e in 0..3 {
                let u = self.tris[t].v[e];
                let v = self.tris[t].v[(e + 1) % 3];
                if (u, v) == (a, b) || (u, v) == (b, a) {
                    self.tris[t].con[e] = true;
                    let n = self.tris[t].adj[e];
                    if n >= 0 {
                        let ne = self.shared_edge(n as usize, t);
                        self.tris[n as usize].con[ne] = true;
                    }
                    return true;
                }
            }
        }
        false
    }

    /// Find an edge strictly crossing segment (a,b) whose quad is strictly
    /// convex (flippable). Anglada's lemma guarantees one exists while any
    /// crossing remains.
    fn find_flippable_crossing(&self, a: usize, b: usize) -> Option<(usize, usize)> {
        let pa = &self.pts[a];
        let pb = &self.pts[b];
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
                let pu = &self.pts[u];
                let pv = &self.pts[v];
                // Strict crossing test.
                let su = orient2d_h(pa, pb, pu);
                let sv = orient2d_h(pa, pb, pv);
                if su == Sign::Zero || sv == Sign::Zero || su == sv {
                    continue;
                }
                let sa = orient2d_h(pu, pv, pa);
                let sb = orient2d_h(pu, pv, pb);
                if sa == Sign::Zero || sb == Sign::Zero || sa == sb {
                    continue;
                }
                debug_assert!(
                    !self.tris[t].con[e],
                    "constraint crosses an existing constraint — arrangement bug"
                );
                // Flippable (strictly convex quad)?
                let n = self.tris[t].adj[e];
                if n < 0 {
                    continue;
                }
                let ne = self.shared_edge(n as usize, t);
                let c = self.tris[t].v[(e + 2) % 3];
                let d = self.tris[n as usize].v[(ne + 2) % 3];
                if orient2d_h(&self.pts[u], &self.pts[d], &self.pts[c]) == Sign::Pos
                    && orient2d_h(&self.pts[d], &self.pts[v], &self.pts[c]) == Sign::Pos
                {
                    return Some((t, e));
                }
            }
        }
        None
    }
}

/// Delaunay test used by `try_flip`: is the fourth point strictly inside the
/// circumcircle of the (CCW) triangle? Split out so tests can call it too.
pub(super) fn is_strictly_non_delaunay(pts: &[R2], tri: [usize; 3], d: usize) -> bool {
    incircle_h(
        &homog2_of(&pts[tri[0]]),
        &homog2_of(&pts[tri[1]]),
        &homog2_of(&pts[tri[2]]),
        &homog2_of(&pts[d]),
    ) == Sign::Pos
}

/// Internal Delaunay test on the cached homogenized points.
fn is_strictly_non_delaunay_h(pts: &[Homog2], tri: [usize; 3], d: usize) -> bool {
    incircle_h(&pts[tri[0]], &pts[tri[1]], &pts[tri[2]], &pts[d]) == Sign::Pos
}

#[cfg(test)]
#[path = "cdt_tests.rs"]
mod tests;
