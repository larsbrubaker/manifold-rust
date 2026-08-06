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

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

use num_rational::BigRational;
use num_traits::{One, Zero};

use crate::linalg::Vec3;
use crate::types::Box;

use super::arrangement::{self, ArrangementInput};
use super::exact::rational::R3;
use super::tri_tri::{tri_tri_intersect, TriTriIsect};

/// Canonical (sorted) edge between two interned vertex ids. Downstream
/// stages (classify rings, propagate flood fill) key their maps on these
/// integers instead of exact rational point pairs — vertex interning at
/// piece-emission time makes id equality coincide with exact geometric
/// identity.
pub type EdgeKey = (u32, u32);

pub fn edge_key(a: u32, b: u32) -> EdgeKey {
    if a <= b {
        (a, b)
    } else {
        (b, a)
    }
}

/// Canonical (lexicographically sorted) exact edge between two points —
/// local key for the split-point registries built before interning exists.
type GeoEdgeKey = (R3, R3);

fn geo_edge_key(a: &R3, b: &R3) -> GeoEdgeKey {
    if a <= b {
        (a.clone(), b.clone())
    } else {
        (b.clone(), a.clone())
    }
}

/// One output fragment: a sub-triangle of an arranged input triangle, or an
/// untouched whole triangle. `v` is wound to match the input mesh's outward
/// orientation; `vi` are the interned ids of the same three vertices.
#[derive(Clone, Debug)]
pub struct Piece {
    /// 0 = first operand (P), 1 = second operand (Q).
    pub mesh: u8,
    /// Index of the originating triangle in its soup.
    pub tri: usize,
    pub v: [R3; 3],
    pub vi: [u32; 3],
}

/// Everything classification and assembly need.
pub struct IntersectionGraph {
    pub pieces: Vec<Piece>,
    /// Interned unique vertices; `Piece::vi` and `EdgeKey` index into this.
    pub verts: Vec<R3>,
    /// Canonical keys of every arrangement constraint edge — the exact
    /// intersection sub-segments the classification rings live on.
    pub isect_edges: HashSet<EdgeKey>,
    /// True when any P×Q pair intersected at all.
    pub any_intersections: bool,
}

/// Exact-point interner: one id per distinct rational point.
#[derive(Default)]
pub struct VertInterner {
    map: HashMap<R3, u32>,
    pub verts: Vec<R3>,
}

impl VertInterner {
    pub fn intern(&mut self, p: &R3) -> u32 {
        if let Some(&id) = self.map.get(p) {
            return id;
        }
        let id = self.verts.len() as u32;
        self.map.insert(p.clone(), id);
        self.verts.push(p.clone());
        id
    }
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

    // 1. Broad + narrow phase. The broad phase is a BVH (the same Collider
    // the exact engine uses) over Q's triangle boxes, queried with each P
    // triangle's box — O((|P|+|Q|)·log|Q|) instead of the all-pairs box
    // sweep. Candidates are re-sorted to ascending qi per pi, so the pair
    // provenance ids match the exhaustive loop exactly (only genuinely
    // intersecting pairs consume an id, and the exact narrow phase decides
    // those identically regardless of broad-phase method).
    let p_boxes: Vec<Box> = p.iter().map(tri_box).collect();
    let q_boxes: Vec<Box> = q.iter().map(tri_box).collect();

    let scene_box = q_boxes
        .iter()
        .enumerate()
        .filter(|(qi, _)| live[1][*qi])
        .fold(Box::new(), |acc, (_, b)| acc.union_box(b));
    let mut q_order: Vec<usize> = (0..q.len()).filter(|&qi| live[1][qi]).collect();
    q_order.sort_by_key(|&qi| crate::sort::morton_code(q_boxes[qi].center(), &scene_box));
    let leaf_boxes: Vec<Box> = q_order.iter().map(|&qi| q_boxes[qi]).collect();
    let leaf_morton: Vec<u32> = q_order
        .iter()
        .map(|&qi| crate::sort::morton_code(q_boxes[qi].center(), &scene_box))
        .collect();
    let collider = crate::collider::Collider::new(leaf_boxes, leaf_morton);

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

    let mut candidates_q: Vec<usize> = Vec::new();
    for (pi, pt) in p.iter().enumerate() {
        if !live[0][pi] {
            continue;
        }
        candidates_q.clear();
        collider.collisions_one(&p_boxes[pi], pi, |_, leaf| {
            candidates_q.push(q_order[leaf]);
        });
        candidates_q.sort_unstable();
        for &qi in &candidates_q {
            let qt = &q[qi];
            if !p_boxes[pi].does_overlap_box(&q_boxes[qi]) {
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
    let t_self = crate::timing::start();

    // 2b. Self-intersections: cut each mesh along its own P×P / Q×Q contact
    // segments (beyond ordinary adjacency). Without these cuts a piece could
    // straddle a fold of a self-overlapping operand, making "is this piece
    // an interior wall of its own solid" ill-defined; with them, both
    // winding numbers the classification needs are constant per flood-fill
    // component (robust/propagate.rs never crosses constraint edges).
    // Broad phase: per-mesh BVH, same approach as the cross-mesh loop above
    // (candidates re-sorted so provenance ids stay deterministic).
    for m in 0..2 {
        let (tris, boxes) = if m == 0 {
            (p, &p_boxes)
        } else {
            (q, &q_boxes)
        };
        let self_scene = boxes
            .iter()
            .enumerate()
            .filter(|(i, _)| live[m][*i])
            .fold(Box::new(), |acc, (_, b)| acc.union_box(b));
        let mut order: Vec<usize> = (0..tris.len()).filter(|&i| live[m][i]).collect();
        order.sort_by_key(|&i| crate::sort::morton_code(boxes[i].center(), &self_scene));
        let self_collider = crate::collider::Collider::new(
            order.iter().map(|&i| boxes[i]).collect(),
            order
                .iter()
                .map(|&i| crate::sort::morton_code(boxes[i].center(), &self_scene))
                .collect(),
        );
        let mut cands: Vec<usize> = Vec::new();
        for i in 0..tris.len() {
            if !live[m][i] {
                continue;
            }
            cands.clear();
            self_collider.collisions_one(&boxes[i], i, |_, leaf| {
                cands.push(order[leaf]);
            });
            cands.sort_unstable();
            for &j in &cands {
                if j <= i || !boxes[i].does_overlap_box(&boxes[j]) {
                    continue;
                }
                let Some(segs) = real_self_contact(tris[i], tris[j]) else {
                    continue;
                };
                for (x, y) in segs {
                    let pair = pair_count;
                    prims[m][i].segments.push((x.clone(), y.clone(), pair));
                    prims[m][j].segments.push((x, y, pair));
                    pair_count += 1;
                }
            }
        }
    }

    crate::timing::print("robust: self-intersection cuts", t_self);
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
    let mut edge_registry: [BTreeMap<GeoEdgeKey, BTreeSet<R3>>; 2] =
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
                let key = geo_edge_key(a, b);
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
    let mut seg_splits: BTreeMap<GeoEdgeKey, BTreeSet<R3>> = BTreeMap::new();
    for m in 0..2 {
        for ti in 0..meshes[m].len() {
            let Some(cands) = &candidates[m][ti] else { continue };
            for (a, b, _prov) in &prims[m][ti].segments {
                let key = geo_edge_key(a, b);
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
    let mut isect_edges: HashSet<EdgeKey> = HashSet::new();
    let mut interner = VertInterner::default();

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
                let key = geo_edge_key(&corners[e], &corners[(e + 1) % 3]);
                if let Some(set) = edge_registry[m].get(&key) {
                    extra.extend(set.iter().cloned());
                }
            }
            let pr = &prims[m][ti];
            // Split points along this triangle's intersection segments
            // discovered by the other side.
            for (a, b, _) in &pr.segments {
                if let Some(set) = seg_splits.get(&geo_edge_key(a, b)) {
                    extra.extend(set.iter().cloned());
                }
            }

            if pr.points.is_empty() && pr.segments.is_empty() && extra.is_empty() {
                // Untouched triangle → whole piece.
                let vi = [
                    interner.intern(&corners[0]),
                    interner.intern(&corners[1]),
                    interner.intern(&corners[2]),
                ];
                pieces.push(Piece {
                    mesh: m as u8,
                    tri: ti,
                    v: corners,
                    vi,
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
                isect_edges.insert(edge_key(
                    interner.intern(&arr.points3[*u]),
                    interner.intern(&arr.points3[*w]),
                ));
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
                let vi = [
                    interner.intern(&v[0]),
                    interner.intern(&v[1]),
                    interner.intern(&v[2]),
                ];
                pieces.push(Piece {
                    mesh: m as u8,
                    tri: ti,
                    v,
                    vi,
                });
            }
        }
    }

    crate::timing::print("robust: arrangements", t_arr);

    IntersectionGraph {
        pieces,
        verts: interner.verts,
        isect_edges,
        any_intersections,
    }
}

/// Real self-intersection of one triangle pair from the same mesh: the
/// contact of `t1` and `t2` reduced by ordinary mesh adjacency. Shared-vertex
/// point contacts and (sub-)segments of a shared edge are the normal way
/// neighboring triangles of a closed mesh touch and yield `None`; anything
/// else is a genuine self-intersection whose segments must cut the surface,
/// so that every emitted piece lies on a single sheet level of its own mesh
/// (robust/mod.rs classifies own-mesh winding per component).
fn real_self_contact(t1: [Vec3; 3], t2: [Vec3; 3]) -> Option<Vec<(R3, R3)>> {
    use super::exact::filtered::orient3d;
    use super::exact::Sign;

    // Shared vertex positions (exact identity) between the pair.
    let shared: Vec<R3> = t1
        .iter()
        .filter(|v| t2.contains(v))
        .map(|v| R3::from_vec3(*v))
        .collect();

    // Adjacency fast paths — the overwhelming bulk of same-mesh box-overlap
    // pairs are edge- or vertex-neighbors whose only contact is that shared
    // simplex, which never needs a cut. Both shortcuts are exact (filtered
    // predicates escalate to rationals when uncertain).
    if shared.len() == 2 {
        // Edge-adjacent: a non-coplanar neighbor can only meet along the
        // shared edge itself (benign). Coplanar pairs fall through.
        if let Some(&opp) = t2.iter().find(|v| !t1.contains(v)) {
            if orient3d(t1[0], t1[1], t1[2], opp) != Sign::Zero {
                return None;
            }
        }
    } else if shared.len() == 1 {
        // Vertex-adjacent: if t2's two non-shared corners lie strictly on
        // one side of t1's plane, the contact is exactly the shared vertex —
        // an isolated point, no cut.
        let others: Vec<Sign> = t2
            .iter()
            .filter(|v| !t1.contains(v))
            .map(|&v| orient3d(t1[0], t1[1], t1[2], v))
            .collect();
        if others.len() == 2 && others[0] != Sign::Zero && others[0] == others[1] {
            return None;
        }
    }

    match tri_tri_intersect(t1, t2) {
        TriTriIsect::None => None,
        // Isolated point contacts (vertex-on-face, edge-through-edge) have
        // zero area on both sides: they never change which sheet a region
        // is on, so they need no cut.
        TriTriIsect::Point(_) => None,
        TriTriIsect::Segment(x, y) => {
            let benign = shared.len() >= 2
                && point_on_segment(&x, &shared[0], &shared[1])
                && point_on_segment(&y, &shared[0], &shared[1]);
            (!benign).then(|| vec![(x, y)])
        }
        // Positive-area coplanar overlap (a fold or doubled patch): cut both
        // triangles along the overlap region's boundary.
        TriTriIsect::Coplanar { polygon, .. } => Some(
            (0..polygon.len())
                .map(|i| {
                    (
                        polygon[i].clone(),
                        polygon[(i + 1) % polygon.len()].clone(),
                    )
                })
                .collect(),
        ),
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
