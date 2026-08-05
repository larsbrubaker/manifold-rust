// robust/intersection_graph.rs — From two triangle soups to classified-ready
// pieces (paper §6).
//
// Pipeline stage between the narrow phase (robust/tri_tri.rs) and
// classification (robust/classify.rs):
//   1. AABB broad phase over the P×Q triangle pairs, exact narrow phase.
//   2. Distribute each pair's intersection primitives to both triangles.
//   3. For coplanar overlaps, cross-copy each side's other primitives
//      (clipped to the overlap region) so both sides subdivide the shared
//      region identically.
//   4. Global registries force a common subdivision of (a) every original
//      mesh edge and (b) every intersection segment, by feeding all split
//      points to every arrangement that sees the same geometry — exact-key
//      edge matching downstream depends on this.
//   5. Build the per-triangle arrangements (robust/arrangement.rs +
//      robust/cdt.rs) and emit `Piece`s: outward-oriented sub-triangles (or
//      whole untouched triangles) tagged with their origin.
//
// Everything is exact; broad-phase boxes are conservative f64.

use std::collections::{BTreeMap, BTreeSet};

use num_rational::BigRational;
use num_traits::{One, Zero};

use crate::linalg::Vec3;
use crate::types::Box;

use super::arrangement::{self, ArrangementInput};
use super::exact::rational::R3;
use super::tri_tri::{tri_tri_intersect, TriTriIsect};

/// Canonical (lexicographically sorted) exact edge between two points.
pub type EdgeKey = (R3, R3);

pub fn edge_key(a: &R3, b: &R3) -> EdgeKey {
    if a <= b {
        (a.clone(), b.clone())
    } else {
        (b.clone(), a.clone())
    }
}

/// One output fragment: a sub-triangle of an arranged input triangle, or an
/// untouched whole triangle. `v` is wound to match the input mesh's outward
/// orientation.
#[derive(Clone, Debug)]
pub struct Piece {
    /// 0 = first operand (P), 1 = second operand (Q).
    pub mesh: u8,
    /// Index of the originating triangle in its soup.
    pub tri: usize,
    pub v: [R3; 3],
}

/// Everything classification and assembly need.
pub struct IntersectionGraph {
    pub pieces: Vec<Piece>,
    /// Canonical keys of every arrangement constraint edge — the exact
    /// intersection sub-segments the classification rings live on.
    pub isect_edges: BTreeSet<EdgeKey>,
    /// True when any P×Q pair intersected at all.
    pub any_intersections: bool,
}

/// A pair's primitives after distribution: segments (including coplanar
/// boundary edges) and isolated points.
#[derive(Clone, Debug, Default)]
struct TriPrims {
    points: Vec<(R3, usize)>,
    segments: Vec<(R3, R3, usize)>,
}

fn tri_box(t: &[Vec3; 3]) -> Box {
    let mut b = Box::from_points(t[0], t[1]);
    b.union_point(t[2]);
    b
}

fn is_degenerate(t: &[Vec3; 3]) -> bool {
    use super::exact::predicates::tri_normal_r;
    tri_normal_r(
        &R3::from_vec3(t[0]),
        &R3::from_vec3(t[1]),
        &R3::from_vec3(t[2]),
    )
    .is_zero()
}

/// Exact: p collinear with (a,b) and within the closed segment.
fn point_on_segment(p: &R3, a: &R3, b: &R3) -> bool {
    super::exact::predicates::point_on_segment_r(p, a, b)
}

/// Clip segment (a,b) to a convex coplanar polygon (2D test via projection
/// on the polygon's own plane). Returns a positive-length sub-segment or
/// None. Used to cross-copy primitives into coplanar overlap regions.
fn clip_segment_to_polygon(a: &R3, b: &R3, poly: &[R3]) -> Option<(R3, R3)> {
    use super::exact::predicates::{orient2d_r, tri_normal_r};
    use super::exact::rational::R2;
    use super::exact::Sign;
    use super::tri_tri::dominant_axis;

    debug_assert!(poly.len() >= 3);
    let n = tri_normal_r(&poly[0], &poly[1], &poly[2]);
    let axis = dominant_axis(&n);
    let mut pts2: Vec<R2> = poly.iter().map(|p| p.project_drop(axis)).collect();
    if orient2d_r(&pts2[0], &pts2[1], &pts2[2]) == Sign::Neg {
        pts2.reverse();
    }
    let a2 = a.project_drop(axis);
    let b2 = b.project_drop(axis);
    let dir = b2.sub(&a2);

    // Parametric clip of [0,1] against each CCW edge halfplane.
    let mut t0 = BigRational::zero();
    let mut t1 = BigRational::one();
    for i in 0..pts2.len() {
        let e0 = &pts2[i];
        let e1 = &pts2[(i + 1) % pts2.len()];
        let edge = e1.sub(e0);
        // Signed distance numerators of a2 + t*dir against the edge line:
        // f(t) = cross(edge, a2 + t*dir - e0) = fa + t * fd.
        let fa = edge.cross(&a2.sub(e0));
        let fd = edge.cross(&dir);
        if fd.is_zero() {
            if fa < BigRational::zero() {
                return None; // parallel and strictly outside
            }
            continue;
        }
        let t_hit = -&fa / &fd;
        if fd > BigRational::zero() {
            // entering: f grows with t → require t >= t_hit
            if t_hit > t0 {
                t0 = t_hit;
            }
        } else if t_hit < t1 {
            t1 = t_hit;
        }
        if t0 >= t1 {
            return None;
        }
    }
    if t0 >= t1 {
        return None;
    }
    let seg = |t: &BigRational| a.add(&b.sub(a).scale(t));
    Some((seg(&t0), seg(&t1)))
}

/// Build the intersection graph for soups `p` and `q` (each triangle wound
/// outward; degenerate triangles are dropped here, paper §5).
pub fn build_graph(p: &[[Vec3; 3]], q: &[[Vec3; 3]]) -> IntersectionGraph {
    let t_all = crate::timing::start();
    let meshes: [&[[Vec3; 3]]; 2] = [p, q];
    let live: [Vec<bool>; 2] = [
        p.iter().map(|t| !is_degenerate(t)).collect(),
        q.iter().map(|t| !is_degenerate(t)).collect(),
    ];

    // 1. Broad + narrow phase. (O(|P|·|Q|) box pruning; the Collider BVH is
    // a later perf upgrade — correctness identical.)
    let p_boxes: Vec<Box> = p.iter().map(tri_box).collect();
    let q_boxes: Vec<Box> = q.iter().map(tri_box).collect();

    // Per-(mesh, tri) primitive lists; provenance = pair index.
    let mut prims: [Vec<TriPrims>; 2] = [
        vec![TriPrims::default(); p.len()],
        vec![TriPrims::default(); q.len()],
    ];
    // Coplanar overlap regions per pair, for the cross-copy step:
    // (p_tri, q_tri, polygon).
    let mut coplanar_regions: Vec<(usize, usize, Vec<R3>)> = Vec::new();
    let mut any_intersections = false;
    let mut pair_count = 0usize;

    for (pi, pt) in p.iter().enumerate() {
        if !live[0][pi] {
            continue;
        }
        for (qi, qt) in q.iter().enumerate() {
            if !live[1][qi] || !p_boxes[pi].does_overlap_box(&q_boxes[qi]) {
                continue;
            }
            let isect = tri_tri_intersect(*pt, *qt);
            let pair = pair_count;
            match isect {
                TriTriIsect::None => continue,
                TriTriIsect::Point(x) => {
                    prims[0][pi].points.push((x.clone(), pair));
                    prims[1][qi].points.push((x, pair));
                }
                TriTriIsect::Segment(x, y) => {
                    prims[0][pi].segments.push((x.clone(), y.clone(), pair));
                    prims[1][qi].segments.push((x, y, pair));
                }
                TriTriIsect::Coplanar { polygon, .. } => {
                    for i in 0..polygon.len() {
                        let a = polygon[i].clone();
                        let b = polygon[(i + 1) % polygon.len()].clone();
                        prims[0][pi].segments.push((a.clone(), b.clone(), pair));
                        prims[1][qi].segments.push((a, b, pair));
                    }
                    coplanar_regions.push((pi, qi, polygon));
                }
            }
            any_intersections = true;
            pair_count += 1;
        }
    }

    crate::timing::print("robust: pair narrow phase", t_all);
    let t_cross = crate::timing::start();

    // 3. Cross-copy primitives through coplanar overlap regions so both
    // sides see identical geometry inside the shared area. Clip against the
    // region to avoid dragging unrelated geometry across.
    for (pi, qi, poly) in &coplanar_regions {
        let from_p: TriPrims = prims[0][*pi].clone();
        let from_q: TriPrims = prims[1][*qi].clone();
        let copy = |src: &TriPrims, dst: &mut TriPrims| {
            for (a, b, prov) in &src.segments {
                if let Some((ca, cb)) = clip_segment_to_polygon(a, b, poly) {
                    if !dst
                        .segments
                        .iter()
                        .any(|(x, y, pv)| pv == prov && ((x, y) == (&ca, &cb) || (x, y) == (&cb, &ca)))
                    {
                        dst.segments.push((ca, cb, *prov));
                    }
                }
            }
            for (pt, prov) in &src.points {
                if clip_segment_to_polygon(pt, pt, poly).is_some()
                    || point_in_polygon_coplanar(pt, poly)
                {
                    if !dst.points.iter().any(|(x, pv)| pv == prov && x == pt) {
                        dst.points.push((pt.clone(), *prov));
                    }
                }
            }
        };
        copy(&from_p, &mut prims[1][*qi]);
        copy(&from_q, &mut prims[0][*pi]);
    }

    crate::timing::print("robust: coplanar cross-copy", t_cross);
    let t_cand = crate::timing::start();

    // 4a. Candidate points per intersected triangle.
    let mut candidates: [Vec<Option<Vec<R3>>>; 2] = [
        vec![None; p.len()],
        vec![None; q.len()],
    ];
    for m in 0..2 {
        for ti in 0..meshes[m].len() {
            let pr = &prims[m][ti];
            if pr.points.is_empty() && pr.segments.is_empty() {
                continue;
            }
            let input = ArrangementInput {
                points: pr.points.clone(),
                segments: pr.segments.clone(),
            };
            candidates[m][ti] = Some(arrangement::candidate_points(meshes[m][ti], &input));
        }
    }

    crate::timing::print("robust: candidate points", t_cand);
    let t_reg = crate::timing::start();

    // 4b. Original-edge registry: split points on each mesh edge (geometric
    // identity — soups have no reliable connectivity).
    let mut edge_registry: [BTreeMap<EdgeKey, BTreeSet<R3>>; 2] =
        [BTreeMap::new(), BTreeMap::new()];
    for m in 0..2 {
        for ti in 0..meshes[m].len() {
            let Some(cands) = &candidates[m][ti] else { continue };
            let t = meshes[m][ti];
            let corners = [
                R3::from_vec3(t[0]),
                R3::from_vec3(t[1]),
                R3::from_vec3(t[2]),
            ];
            for e in 0..3 {
                let a = &corners[e];
                let b = &corners[(e + 1) % 3];
                let key = edge_key(a, b);
                for pt in cands {
                    if pt != a && pt != b && point_on_segment(pt, a, b) {
                        edge_registry[m].entry(key.clone()).or_default().insert(pt.clone());
                    }
                }
            }
        }
    }

    // 4c. Intersection-segment registry: for every pair segment, gather the
    // split points both sides know about.
    let mut seg_splits: BTreeMap<EdgeKey, BTreeSet<R3>> = BTreeMap::new();
    for m in 0..2 {
        for ti in 0..meshes[m].len() {
            let Some(cands) = &candidates[m][ti] else { continue };
            for (a, b, _prov) in &prims[m][ti].segments {
                let key = edge_key(a, b);
                for pt in cands {
                    if pt != a && pt != b && point_on_segment(pt, a, b) {
                        seg_splits.entry(key.clone()).or_default().insert(pt.clone());
                    }
                }
            }
        }
    }

    crate::timing::print("robust: split registries", t_reg);
    let t_arr = crate::timing::start();

    // 5. Build arrangements and emit pieces.
    let mut pieces: Vec<Piece> = Vec::new();
    let mut isect_edges: BTreeSet<EdgeKey> = BTreeSet::new();

    for m in 0..2 {
        for ti in 0..meshes[m].len() {
            if !live[m][ti] {
                continue;
            }
            let t = meshes[m][ti];
            let corners = [
                R3::from_vec3(t[0]),
                R3::from_vec3(t[1]),
                R3::from_vec3(t[2]),
            ];
            // Boundary split points for this triangle.
            let mut extra: BTreeSet<R3> = BTreeSet::new();
            for e in 0..3 {
                let key = edge_key(&corners[e], &corners[(e + 1) % 3]);
                if let Some(set) = edge_registry[m].get(&key) {
                    extra.extend(set.iter().cloned());
                }
            }
            let pr = &prims[m][ti];
            // Split points along this triangle's intersection segments
            // discovered by the other side.
            for (a, b, _) in &pr.segments {
                if let Some(set) = seg_splits.get(&edge_key(a, b)) {
                    extra.extend(set.iter().cloned());
                }
            }

            if pr.points.is_empty() && pr.segments.is_empty() && extra.is_empty() {
                // Untouched triangle → whole piece.
                pieces.push(Piece {
                    mesh: m as u8,
                    tri: ti,
                    v: corners,
                });
                continue;
            }

            let mut input = ArrangementInput {
                points: pr.points.clone(),
                segments: pr.segments.clone(),
            };
            for pt in extra {
                input.points.push((pt, usize::MAX));
            }
            let arr = arrangement::build(t, &input);
            for (u, w) in arr.constraints.keys() {
                isect_edges.insert(edge_key(&arr.points3[*u], &arr.points3[*w]));
            }
            for st in &arr.tris {
                let (a, b, c) = (st[0], st[1], st[2]);
                let v = if arr.flipped {
                    [
                        arr.points3[a].clone(),
                        arr.points3[c].clone(),
                        arr.points3[b].clone(),
                    ]
                } else {
                    [
                        arr.points3[a].clone(),
                        arr.points3[b].clone(),
                        arr.points3[c].clone(),
                    ]
                };
                pieces.push(Piece {
                    mesh: m as u8,
                    tri: ti,
                    v,
                });
            }
        }
    }

    crate::timing::print("robust: arrangements", t_arr);

    IntersectionGraph {
        pieces,
        isect_edges,
        any_intersections,
    }
}

/// Exact point-in-convex-polygon test for a point on the polygon's plane.
fn point_in_polygon_coplanar(p: &R3, poly: &[R3]) -> bool {
    use super::exact::predicates::{orient2d_r, tri_normal_r};
    use super::exact::Sign;
    use super::tri_tri::dominant_axis;

    let n = tri_normal_r(&poly[0], &poly[1], &poly[2]);
    let axis = dominant_axis(&n);
    let mut pts2: Vec<_> = poly.iter().map(|q| q.project_drop(axis)).collect();
    if orient2d_r(&pts2[0], &pts2[1], &pts2[2]) == Sign::Neg {
        pts2.reverse();
    }
    let p2 = p.project_drop(axis);
    for i in 0..pts2.len() {
        if orient2d_r(&pts2[i], &pts2[(i + 1) % pts2.len()], &p2) == Sign::Neg {
            return false;
        }
    }
    true
}
