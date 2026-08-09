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

use std::collections::BTreeSet;

// Fx hashing instead of SipHash. Every map/set below is probe-only or has an
// order-invariant consumer (documented per site); the hasher is unseeded, so
// even iteration order is stable across runs — output cannot depend on it.
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};

use num_rational::BigRational;
use num_traits::{One, Zero};

use crate::linalg::Vec3;
use crate::types::Box;

use super::arrangement::{self, ArrangementInput};
use super::exact::rational::{r3_eq, R3, R3Key};
use super::exact::Sign;
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
/// Wrapped in [`R3Key`] so the registry hashes structurally instead of
/// comparing rationals per probe.
type GeoEdgeKey = (R3Key, R3Key);

fn geo_edge_key(a: &R3, b: &R3) -> GeoEdgeKey {
    if a <= b {
        (R3Key(a.clone()), R3Key(b.clone()))
    } else {
        (R3Key(b.clone()), R3Key(a.clone()))
    }
}

/// Canonical original-mesh edge keyed by raw coordinate bits — original
/// edges always join exact f64 vertices, so the boundary-split registry
/// never needs rational keys (and untouched triangles probe it for free).
type BitEdgeKey = ([u64; 3], [u64; 3]);

fn bit_edge_key(a: Vec3, b: Vec3) -> BitEdgeKey {
    let (ka, kb) = (f64_key(a), f64_key(b));
    if ka <= kb {
        (ka, kb)
    } else {
        (kb, ka)
    }
}

/// One output fragment: a sub-triangle of an arranged input triangle, or an
/// untouched whole triangle. `v` is wound to match the input mesh's outward
/// orientation; `vi` are the interned ids of the same three vertices.
#[derive(Clone, Copy, Debug)]
pub struct Piece {
    /// 0 = first operand (P), 1 = second operand (Q).
    pub mesh: u8,
    /// Index of the originating triangle in its soup.
    pub tri: usize,
    /// Interned vertex ids (indices into `IntersectionGraph::verts`), wound
    /// to the input mesh's outward orientation. Pieces carry no coordinates
    /// of their own — the shared tables keep untouched triangles free of
    /// rational clones entirely.
    pub vi: [u32; 3],
}

/// Everything classification and assembly need.
pub struct IntersectionGraph {
    pub pieces: Vec<Piece>,
    /// Interned unique vertices; `Piece::vi` and `EdgeKey` index into this.
    pub verts: Vec<R3>,
    /// Correctly rounded f64 approximation per interned vertex (exact for
    /// input vertices) — float filters and output assembly read these
    /// instead of re-rounding rationals.
    pub verts_f64: Vec<Vec3>,
    /// Canonical keys of every arrangement constraint edge — the exact
    /// intersection sub-segments the classification rings live on.
    pub isect_edges: HashSet<EdgeKey>,
    /// True when any P×Q pair intersected at all.
    pub any_intersections: bool,
}

impl IntersectionGraph {
    /// The three exact vertices of a piece.
    pub fn piece_verts(&self, pi: usize) -> [&R3; 3] {
        let vi = self.pieces[pi].vi;
        [
            &self.verts[vi[0] as usize],
            &self.verts[vi[1] as usize],
            &self.verts[vi[2] as usize],
        ]
    }
}

/// Exact-point interner: one id per distinct point, with two disjoint key
/// spaces. f64-representable points (all input vertices, and any constructed
/// point that rounds exactly) key on their coordinate bits — no rational
/// hashing, so untouched input triangles intern for the cost of a HashMap
/// probe. Only genuinely non-representable constructed points use the
/// rational map. `verts_f64` caches the correctly rounded approximation of
/// every id (exact for bit-keyed points), which downstream float filters
/// and output assembly reuse instead of re-rounding.
/// Order invariance: both maps are probe-only (`get`/`entry`, never
/// iterated); ids come from `verts.len()` at insertion time, so they depend
/// only on the sequential call order, not on the hasher.
#[derive(Default)]
pub struct VertInterner {
    map: HashMap<R3Key, u32>,
    fmap: HashMap<[u64; 3], u32>,
    pub verts: Vec<R3>,
    pub verts_f64: Vec<Vec3>,
}

fn f64_key(v: Vec3) -> [u64; 3] {
    // Normalize -0.0 so it shares an id with +0.0 (they are the same
    // rational point).
    let norm = |x: f64| if x == 0.0 { 0.0f64 } else { x }.to_bits();
    [norm(v.x), norm(v.y), norm(v.z)]
}

impl VertInterner {
    /// Intern an exact-f64 point (input mesh vertices): zero rational work
    /// on hits; one `R3::from_vec3` on first sight, for the exact table.
    pub fn intern_f64(&mut self, v: Vec3) -> u32 {
        let key = f64_key(v);
        if let Some(&id) = self.fmap.get(&key) {
            return id;
        }
        let id = self.verts.len() as u32;
        self.fmap.insert(key, id);
        self.verts.push(R3::from_vec3(v));
        self.verts_f64.push(v);
        id
    }

    /// Intern an exact rational point. Representable points route to the
    /// f64 key space so both paths agree on ids.
    pub fn intern(&mut self, p: &R3) -> u32 {
        let rounded = p.to_vec3_rounded();
        if r3_eq(&R3::from_vec3(rounded), p) {
            return self.intern_f64(rounded);
        }
        let next = self.verts.len() as u32;
        match self.map.entry(R3Key(p.clone())) {
            std::collections::hash_map::Entry::Occupied(e) => *e.get(),
            std::collections::hash_map::Entry::Vacant(e) => {
                e.insert(next);
                self.verts.push(p.clone());
                self.verts_f64.push(rounded);
                next
            }
        }
    }
}

/// A pair's primitives after distribution: segments (including coplanar
/// boundary edges) and isolated points.
#[derive(Clone, Debug, Default)]
struct TriPrims {
    points: Vec<(R3, usize)>,
    segments: Vec<(R3, R3, usize)>,
}

pub(super) fn tri_box(t: &[Vec3; 3]) -> Box {
    let mut b = Box::from_points(t[0], t[1]);
    b.union_point(t[2]);
    b
}

pub(super) fn is_degenerate(t: &[Vec3; 3]) -> bool {
    // Certified-nonzero f64 cross first (magnitude-permanent bound, matching
    // exact/approx.rs conventions); only near-degenerate triangles pay for
    // the rational cross.
    const EPS: f64 = f64::EPSILON * 0.5;
    let u = t[1] - t[0];
    let v = t[2] - t[0];
    let n = crate::linalg::cross(u, v);
    let m = |k: usize| t[0][k].abs() + t[1][k].abs() + t[2][k].abs();
    let (mx, my, mz) = (m(0), m(1), m(2));
    if n.x.abs() > 16.0 * EPS * my * mz
        || n.y.abs() > 16.0 * EPS * mz * mx
        || n.z.abs() > 16.0 * EPS * mx * my
    {
        return false;
    }
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

/// Correctly rounded f64 approximation of an exact point (relative error
/// ≤ ε per coordinate) for the semi-static prefilters in exact/approx.rs.
fn approx3(p: &R3) -> [f64; 3] {
    use super::exact::rational::rat_to_f64;
    [rat_to_f64(&p.x), rat_to_f64(&p.y), rat_to_f64(&p.z)]
}

/// Filtered point-on-segment: the approx prefilter rejects the generic case
/// without touching BigInt; only near-incidences run the exact test.
/// Conservative 3D box `[min; 3], [max; 3]` around a segment's exact
/// endpoints from their correctly rounded approximations, inflated past
/// rounding error — a point exactly on the segment is never rejected by
/// testing its approximation against this box. Mirrors the 2D prefilter in
/// robust/arrangement.rs; the registry sweeps below were quadratic in exact
/// comparisons without it.
#[inline]
fn seg_box3(a: [f64; 3], b: [f64; 3]) -> [[f64; 3]; 2] {
    let pad = |x: f64| x.abs() * 1e-15 + f64::MIN_POSITIVE;
    let mut lo = [0.0; 3];
    let mut hi = [0.0; 3];
    for k in 0..3 {
        let (l, h) = (a[k].min(b[k]), a[k].max(b[k]));
        lo[k] = l - pad(l);
        hi[k] = h + pad(h);
    }
    [lo, hi]
}

#[inline]
fn box3_contains(b: &[[f64; 3]; 2], p: [f64; 3]) -> bool {
    (0..3).all(|k| p[k] >= b[0][k] && p[k] <= b[1][k])
}

fn point_on_segment_f(p_approx: [f64; 3], p: &R3, a_approx: [f64; 3], a: &R3, b_approx: [f64; 3], b: &R3) -> bool {
    match super::exact::approx::not_on_segment_a(p_approx, a_approx, b_approx) {
        Some(false) => false,
        _ => point_on_segment(p, a, b),
    }
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
    build_graph_with_token(p, q, None).expect("uncancellable build_graph cannot cancel")
}

/// [`build_graph`] with cooperative cancellation. Returns `None` when the
/// token fires. Checks run per triangle in every phase and inside the
/// arrangement sweeps — heavily self-intersecting inputs spend minutes in
/// per-triangle quadratic loops, and a cancel that only top-level phases
/// notice can overshoot its deadline by that much (Thingi10K #42211 ran
/// 565 s past a 60 s cancel before this plumbing).
pub fn build_graph_with_token(
    p: &[[Vec3; 3]],
    q: &[[Vec3; 3]],
    token: Option<&crate::cancel::CancelToken>,
) -> Option<IntersectionGraph> {
    let cancelled = || crate::cancel::is_cancelled(token);
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
        if cancelled() {
            return None;
        }
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
        // The exact narrow phase per triangle is pure; workers return each
        // triangle's (j, segments) contacts and per-worker stats, and the
        // sequential merge assigns provenance pair ids in (i, j) order —
        // identical to the sequential sweep.
        let mut n_pairs = 0usize;
        let mut n_cut = 0usize;
        let mut stats = SelfCutStats::default();
        let contact_results = crate::par::maybe_par_map_ct(tris.len(), 64, token, |i| {
            let mut local = SelfCutStats::default();
            let mut contacts: Vec<(usize, Vec<(R3, R3)>)> = Vec::new();
            let mut local_pairs = 0usize;
            if live[m][i] {
                let mut cands: Vec<usize> = Vec::new();
                self_collider.collisions_one(&boxes[i], i, |_, leaf| {
                    cands.push(order[leaf]);
                });
                cands.sort_unstable();
                for &j in &cands {
                    if j <= i || !boxes[i].does_overlap_box(&boxes[j]) {
                        continue;
                    }
                    local_pairs += 1;
                    if let Some(segs) = real_self_contact(tris[i], tris[j], &mut local) {
                        contacts.push((j, segs));
                    }
                }
            }
            (contacts, local, local_pairs)
        })?;
        for (i, (contacts, local, local_pairs)) in contact_results.into_iter().enumerate() {
            n_pairs += local_pairs;
            stats.add(&local);
            for (j, segs) in contacts {
                n_cut += 1;
                for (x, y) in segs {
                    let pair = pair_count;
                    prims[m][i].segments.push((x.clone(), y.clone(), pair));
                    prims[m][j].segments.push((x, y, pair));
                    pair_count += 1;
                }
            }
        }
        crate::timing::print_count(
            &format!("robust: self-cut mesh {m}: {n_pairs} box pairs, {n_cut} cutting"),
        );
        crate::timing::print_count(&format!(
            "robust: self-cut mesh {m} tri_tri exits: {}",
            super::tri_tri::stats::snapshot_and_reset()
        ));
        crate::timing::print_count(&format!(
            "robust: self-cut mesh {m} paths: identical {}, edge-benign {}, vert-benign {}, \
             full {} ({:.3}s: none {}, point {}, seg-benign {})",
            stats.identical,
            stats.edge_benign,
            stats.vert_benign,
            stats.full,
            stats.full_secs,
            stats.full_none,
            stats.full_point,
            stats.full_seg_benign,
        ));
    }

    crate::timing::print("robust: self-intersection cuts", t_self);
    let t_cross = crate::timing::start();

    // 3. Cross-copy primitives through coplanar overlap regions so both
    // sides see identical geometry inside the shared area. Clip against the
    // region to avoid dragging unrelated geometry across.
    for (pi, qi, poly) in &coplanar_regions {
        if cancelled() {
            return None;
        }
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

    // 4a. Candidate points per intersected triangle. Pure per triangle, so
    // the map parallelizes under the bit-identical rule: results land in
    // worklist order regardless of schedule.
    let mut candidates: [Vec<Option<Vec<R3>>>; 2] = [
        vec![None; p.len()],
        vec![None; q.len()],
    ];
    let cand_work: Vec<(usize, usize)> = (0..2)
        .flat_map(|m| (0..meshes[m].len()).map(move |ti| (m, ti)))
        .filter(|&(m, ti)| {
            let pr = &prims[m][ti];
            !pr.points.is_empty() || !pr.segments.is_empty()
        })
        .collect();
    let cand_results = crate::par::maybe_par_map_ct(cand_work.len(), 16, token, |i| {
        let (m, ti) = cand_work[i];
        let pr = &prims[m][ti];
        let input = ArrangementInput {
            points: pr.points.clone(),
            segments: pr.segments.clone(),
        };
        arrangement::candidate_points(meshes[m][ti], &input, token)
    })?;
    for (&(m, ti), cands) in cand_work.iter().zip(cand_results) {
        candidates[m][ti] = Some(cands?);
    }

    crate::timing::print("robust: candidate points", t_cand);
    let t_reg = crate::timing::start();

    // 4b. Original-edge registry: split points on each mesh edge (geometric
    // identity — soups have no reliable connectivity). Bit-keyed: original
    // edges join exact f64 vertices.
    // Both registry sweeps are pure per triangle; workers collect local
    // (key, point) hits and the single-threaded merge inserts them in
    // worklist order. Registry values are sets, so content is
    // order-independent and the merge order only preserves determinism of
    // allocation, not meaning.
    let reg_work: Vec<(usize, usize)> = (0..2)
        .flat_map(|m| (0..meshes[m].len()).map(move |ti| (m, ti)))
        .filter(|&(m, ti)| candidates[m][ti].is_some())
        .collect();

    // Registry values are (points, seen) rather than BTreeSet: probes and
    // dedup go through structural R3Key hashing, so the sequential merge
    // never pays ordered rational comparisons, and the consuming `extra`
    // sets see each point exactly once — the same content BTreeSet gave.
    // Order invariance: the registry is only ever `get`/`entry`-probed, and
    // its point lists land in a `BTreeSet<R3>` (`extra`) below, which sorts.
    let mut edge_registry: [HashMap<BitEdgeKey, (Vec<R3>, HashSet<R3Key>)>; 2] =
        [HashMap::default(), HashMap::default()];
    let edge_hits = crate::par::maybe_par_map_ct(reg_work.len(), 16, token, |i| {
        let (m, ti) = reg_work[i];
        let cands = candidates[m][ti].as_ref().expect("filtered to Some");
        let t = meshes[m][ti];
        let corners = [
            R3::from_vec3(t[0]),
            R3::from_vec3(t[1]),
            R3::from_vec3(t[2]),
        ];
        let ca: [[f64; 3]; 3] = [
            [t[0].x, t[0].y, t[0].z],
            [t[1].x, t[1].y, t[1].z],
            [t[2].x, t[2].y, t[2].z],
        ];
        let cands_a: Vec<[f64; 3]> = cands.iter().map(approx3).collect();
        let mut hits: Vec<(BitEdgeKey, R3)> = Vec::new();
        for e in 0..3 {
            let a = &corners[e];
            let b = &corners[(e + 1) % 3];
            let key = bit_edge_key(t[e], t[(e + 1) % 3]);
            let sbox = seg_box3(ca[e], ca[(e + 1) % 3]);
            for (pt, pt_a) in cands.iter().zip(&cands_a) {
                // A point on the edge lies inside its inflated box; the
                // reject skips the exact comparisons for everything else.
                if !box3_contains(&sbox, *pt_a) {
                    continue;
                }
                if !r3_eq(pt, a)
                    && !r3_eq(pt, b)
                    && point_on_segment_f(*pt_a, pt, ca[e], a, ca[(e + 1) % 3], b)
                {
                    hits.push((key, pt.clone()));
                }
            }
        }
        hits
    })?;
    for (&(m, _), hits) in reg_work.iter().zip(&edge_hits) {
        for (key, pt) in hits {
            let e = edge_registry[m].entry(*key).or_default();
            if e.1.insert(R3Key(pt.clone())) {
                e.0.push(pt.clone());
            }
        }
    }

    // 4c. Intersection-segment registry: for every pair segment, gather the
    // split points both sides know about.
    // Hash-keyed with structural R3Key hashing: the map is only ever probed
    // (entry/get), never iterated, and BTreeMap's exact rational comparisons
    // per probe dominated this phase on segment-heavy meshes.
    // Same order-invariance argument as `edge_registry`: probe-only, and the
    // points it hands out are re-sorted through `extra: BTreeSet<R3>`.
    let mut seg_splits: HashMap<(R3Key, R3Key), (Vec<R3>, HashSet<R3Key>)> = HashMap::default();
    let split_hits = crate::par::maybe_par_map_ct(reg_work.len(), 16, token, |i| {
        let (m, ti) = reg_work[i];
        let cands = candidates[m][ti].as_ref().expect("filtered to Some");
        let cands_a: Vec<[f64; 3]> = cands.iter().map(approx3).collect();
        let mut hits: Vec<(GeoEdgeKey, R3)> = Vec::new();
        for (a, b, _prov) in &prims[m][ti].segments {
            let key = geo_edge_key(a, b);
            let (aa, ba) = (approx3(a), approx3(b));
            let sbox = seg_box3(aa, ba);
            for (pt, pt_a) in cands.iter().zip(&cands_a) {
                if !box3_contains(&sbox, *pt_a) {
                    continue;
                }
                if !r3_eq(pt, a) && !r3_eq(pt, b) && point_on_segment_f(*pt_a, pt, aa, a, ba, b) {
                    hits.push((key.clone(), pt.clone()));
                }
            }
        }
        hits
    })?;
    for hits in &split_hits {
        for (key, pt) in hits {
            let e = seg_splits.entry(key.clone()).or_default();
            if e.1.insert(R3Key(pt.clone())) {
                e.0.push(pt.clone());
            }
        }
    }

    crate::timing::print("robust: split registries", t_reg);
    let t_arr = crate::timing::start();

    // 5. Build arrangements and emit pieces. The per-triangle arrangement
    // (registry probes, CDT, crossings) is pure and runs in parallel; the
    // interner is order-sensitive, so interning and piece emission replay
    // the results strictly in worklist order — outputs are bit-identical to
    // the sequential build.
    enum TriResult {
        /// Untouched triangle → whole piece, interned by f64 bits.
        Untouched,
        Arranged(arrangement::Arrangement),
    }
    let arr_work: Vec<(usize, usize)> = (0..2)
        .flat_map(|m| (0..meshes[m].len()).map(move |ti| (m, ti)))
        .filter(|&(m, ti)| live[m][ti])
        .collect();
    let arr_results = crate::par::maybe_par_map_ct(arr_work.len(), 16, token, |i| {
        let (m, ti) = arr_work[i];
        let t = meshes[m][ti];
        let pr = &prims[m][ti];
        // Boundary split points for this triangle (bit-keyed: uncut
        // triangles probe with zero rational work).
        let mut extra: BTreeSet<R3> = BTreeSet::new();
        for e in 0..3 {
            if let Some(set) = edge_registry[m].get(&bit_edge_key(t[e], t[(e + 1) % 3])) {
                extra.extend(set.0.iter().cloned());
            }
        }
        // Split points along this triangle's intersection segments
        // discovered by the other side.
        for (a, b, _) in &pr.segments {
            if let Some(set) = seg_splits.get(&geo_edge_key(a, b)) {
                extra.extend(set.0.iter().cloned());
            }
        }

        if pr.points.is_empty() && pr.segments.is_empty() && extra.is_empty() {
            return Some(TriResult::Untouched);
        }
        let mut input = ArrangementInput {
            points: pr.points.clone(),
            segments: pr.segments.clone(),
        };
        for pt in extra {
            input.points.push((pt, usize::MAX));
        }
        arrangement::build(t, &input, token).map(TriResult::Arranged)
    })?;

    let mut pieces: Vec<Piece> = Vec::new();
    // Membership-only set (never iterated by this crate); order-invariant.
    let mut isect_edges: HashSet<EdgeKey> = HashSet::default();
    let mut interner = VertInterner::default();
    for (&(m, ti), result) in arr_work.iter().zip(arr_results) {
        let t = meshes[m][ti];
        match result? {
            TriResult::Untouched => {
                pieces.push(Piece {
                    mesh: m as u8,
                    tri: ti,
                    vi: [
                        interner.intern_f64(t[0]),
                        interner.intern_f64(t[1]),
                        interner.intern_f64(t[2]),
                    ],
                });
            }
            TriResult::Arranged(arr) => {
                // Intern each arrangement point once; sub-triangles and
                // constraint edges then only shuffle ids.
                let ids: Vec<u32> = arr.points3.iter().map(|p| interner.intern(p)).collect();
                for (u, w) in arr.constraints.keys() {
                    isect_edges.insert(edge_key(ids[*u], ids[*w]));
                }
                for st in &arr.tris {
                    let (a, b, c) = (st[0], st[1], st[2]);
                    let vi = if arr.flipped {
                        [ids[a], ids[c], ids[b]]
                    } else {
                        [ids[a], ids[b], ids[c]]
                    };
                    pieces.push(Piece {
                        mesh: m as u8,
                        tri: ti,
                        vi,
                    });
                }
            }
        }
    }

    crate::timing::print("robust: arrangements", t_arr);
    crate::timing::print_count(&format!(
        "robust: arrangement phases: {}",
        arrangement::stats::snapshot_and_reset()
    ));

    Some(IntersectionGraph {
        pieces,
        verts: interner.verts,
        verts_f64: interner.verts_f64,
        isect_edges,
        any_intersections,
    })
}

/// Real self-intersection of one triangle pair from the same mesh: the
/// contact of `t1` and `t2` reduced by ordinary mesh adjacency. Shared-vertex
/// point contacts and (sub-)segments of a shared edge are the normal way
/// neighboring triangles of a closed mesh touch and yield `None`; anything
/// else is a genuine self-intersection whose segments must cut the surface,
/// so that every emitted piece lies on a single sheet level of its own mesh
/// (robust/mod.rs classifies own-mesh winding per component).
/// orient3d(t[0], t[1], t[2], v): float filter first, exact integer
/// determinant on escalation. This replaced a cached exact-plane structure
/// (TriPlane) — with intpred's division-free fallback, building planes
/// eagerly per triangle cost more than it ever saved.
fn orient3d_plane(t: &[Vec3; 3], v: Vec3) -> Sign {
    if let Some(s) = super::exact::approx::orient3d_a(
        [t[0].x, t[0].y, t[0].z],
        [t[1].x, t[1].y, t[1].z],
        [t[2].x, t[2].y, t[2].z],
        [v.x, v.y, v.z],
    ) {
        return s;
    }
    super::exact::intpred::orient3d_i(
        [t[0].x, t[0].y, t[0].z],
        [t[1].x, t[1].y, t[1].z],
        [t[2].x, t[2].y, t[2].z],
        [v.x, v.y, v.z],
    )
}

/// Axis of the largest |component| of the (f64) triangle normal. Only a
/// projection *choice*: when the exact normal's chosen component happens to
/// be zero, the projected points go collinear, the exact 2D signs come back
/// Zero, and every shortcut below falls through — sound, just unoptimized.
fn dominant_axis_f64(t: [Vec3; 3]) -> usize {
    let n = crate::linalg::cross(t[1] - t[0], t[2] - t[0]);
    let (ax, ay, az) = (n.x.abs(), n.y.abs(), n.z.abs());
    if az >= ax && az >= ay {
        2
    } else if ay >= ax {
        1
    } else {
        0
    }
}

/// `R3::project_drop` for raw f64 points (same cyclic axis convention).
fn project_f64(v: Vec3, axis: usize) -> crate::linalg::Vec2 {
    match axis {
        0 => crate::linalg::Vec2::new(v.y, v.z),
        1 => crate::linalg::Vec2::new(v.z, v.x),
        _ => crate::linalg::Vec2::new(v.x, v.y),
    }
}

/// Per-path counters for the self-cut narrow phase, printed under
/// MANIFOLD_TIMING to show where box-pair time goes (shortcut hits vs full
/// tri_tri calls and their outcomes).
#[derive(Default, Debug)]
pub(super) struct SelfCutStats {
    identical: usize,
    edge_benign: usize,
    vert_benign: usize,
    full: usize,
    full_none: usize,
    full_point: usize,
    full_seg_benign: usize,
    full_secs: f64,
}

impl SelfCutStats {
    /// Merge a worker's counters. Only diagnostics — `full_secs` summed
    /// across workers reports aggregate CPU seconds, not wall time.
    fn add(&mut self, o: &SelfCutStats) {
        self.identical += o.identical;
        self.edge_benign += o.edge_benign;
        self.vert_benign += o.vert_benign;
        self.full += o.full;
        self.full_none += o.full_none;
        self.full_point += o.full_point;
        self.full_seg_benign += o.full_seg_benign;
        self.full_secs += o.full_secs;
    }
}

/// The self-cut narrow phase: `Some(segments)` when `t1` and `t2` genuinely
/// intersect (a crossing or a positive-area coplanar overlap), `None` when
/// their only contact is ordinary adjacency — a shared edge, a shared
/// vertex, an isolated point, or an exactly duplicated triangle. Also used
/// by `soup::has_self_intersections` to decide engine dispatch.
pub(super) fn real_self_contact(
    t1: [Vec3; 3],
    t2: [Vec3; 3],
    stats: &mut SelfCutStats,
) -> Option<Vec<(R3, R3)>> {
    use super::exact::Sign;

    // Shared vertex positions (exact f64 identity) between the pair. Kept in
    // f64: hundreds of thousands of benign pairs pass through here, and the
    // rational form is only needed by the final Segment-benign check.
    //
    // Exactly identical triangles (all three vertices coincide — doubled
    // surfaces, which some scans apply to their whole mesh) need no cut:
    // both emit whole pieces with identical interned ids, so they land on
    // one wall in robust/cells.rs and their winding steps add (a doubled
    // sheet steps by two, a fold cancels to zero). Cutting them instead
    // would drag every such triangle through the full arrangement pipeline
    // along its own boundary, for nothing.

    // Adjacency fast paths — the overwhelming bulk of same-mesh box-overlap
    // pairs are edge- or vertex-neighbors whose only contact is that shared
    // simplex, which never needs a cut. All shortcuts are exact (filtered
    // predicates escalate to rationals when uncertain); flat triangulated
    // regions make the *coplanar* neighbor cases as common as the generic
    // ones, and without their 2D shortcuts every such pair pays for a full
    // rational coplanar-overlap clip.
    // Stack-allocated shared-vertex list: this runs per box pair (hundreds
    // of thousands on dense meshes) and a heap Vec here is measurable.
    let mut shared_f = [Vec3::default(); 3];
    let mut n_shared = 0usize;
    for &v in &t1 {
        if t2.contains(&v) {
            shared_f[n_shared] = v;
            n_shared += 1;
        }
    }
    let shared_f = &shared_f[..n_shared];
    if shared_f.len() == 3 {
        stats.identical += 1;
        return None;
    }
    if shared_f.len() == 2 {
        if let Some(&opp) = t2.iter().find(|v| !t1.contains(v)) {
            // Non-coplanar edge-neighbors only meet along the shared edge.
            if orient3d_plane(&t1, opp) != Sign::Zero {
                stats.edge_benign += 1;
                return None;
            }
            // Coplanar edge-neighbors: benign exactly when the two opposite
            // corners strictly straddle the shared edge's line within the
            // plane (the flat-plate case) — then the closed half-plane
            // intersection is the shared edge itself.
            if let Some(&own) = t1.iter().find(|v| !t2.contains(v)) {
                let axis = dominant_axis_f64(t1);
                let p2 = |v: Vec3| project_f64(v, axis);
                let s_own =
                    super::exact::filtered::orient2d(p2(shared_f[0]), p2(shared_f[1]), p2(own));
                let s_opp =
                    super::exact::filtered::orient2d(p2(shared_f[0]), p2(shared_f[1]), p2(opp));
                if s_own != Sign::Zero && s_opp != Sign::Zero && s_own != s_opp {
                    stats.edge_benign += 1;
                    return None;
                }
            }
        }
    } else if shared_f.len() == 1 {
        // Vertex-adjacent: if t2's two non-shared corners lie strictly on
        // one side of t1's plane, the contact is exactly the shared vertex —
        // an isolated point, no cut.
        let mut others = [(Vec3::default(), Sign::Zero); 3];
        let mut n_others = 0usize;
        for &v in &t2 {
            if !t1.contains(&v) {
                others[n_others] = (v, orient3d_plane(&t1, v));
                n_others += 1;
            }
        }
        let others = &others[..n_others];
        if others.len() == 2 && others[0].1 != Sign::Zero && others[0].1 == others[1].1 {
            stats.vert_benign += 1;
            return None;
        }
        // Fully coplanar vertex-neighbors (triangle fans on flat regions):
        // an edge through the shared vertex that strictly separates the two
        // triangles certifies the contact is exactly that vertex.
        if others.len() == 2 && others[0].1 == Sign::Zero && others[1].1 == Sign::Zero {
            let axis = dominant_axis_f64(t1);
            let p2 = |v: Vec3| project_f64(v, axis);
            let v0 = shared_f[0];
            let mut own = [Vec3::default(); 3];
            let mut n_own = 0usize;
            for &v in &t1 {
                if !t2.contains(&v) {
                    own[n_own] = v;
                    n_own += 1;
                }
            }
            let own = &own[..n_own];
            let other = [others[0].0, others[1].0];
            // Candidate separators: each triangle's two edges through v0,
            // tested against its own third corner vs the other triangle's
            // two corners.
            let separated = |ea: Vec3, third: Vec3, far: [Vec3; 2]| -> bool {
                let s_t = super::exact::filtered::orient2d(p2(v0), p2(ea), p2(third));
                if s_t == Sign::Zero {
                    return false;
                }
                far.iter().all(|&f| {
                    let s = super::exact::filtered::orient2d(p2(v0), p2(ea), p2(f));
                    s != Sign::Zero && s != s_t
                })
            };
            if own.len() == 2
                && (separated(own[0], own[1], other)
                    || separated(own[1], own[0], other)
                    || separated(other[0], other[1], [own[0], own[1]])
                    || separated(other[1], other[0], [own[0], own[1]]))
            {
                stats.vert_benign += 1;
                return None;
            }
        }
    }

    stats.full += 1;
    let t_full = crate::timing::Stopwatch::start();
    let isect = tri_tri_intersect(t1, t2);
    stats.full_secs += t_full.elapsed_secs();
    match isect {
        TriTriIsect::None => {
            stats.full_none += 1;
            None
        }
        // Isolated point contacts (vertex-on-face, edge-through-edge) have
        // zero area on both sides: they never change which sheet a region
        // is on, so they need no cut.
        TriTriIsect::Point(_) => {
            stats.full_point += 1;
            None
        }
        TriTriIsect::Segment(x, y) => {
            let benign = shared_f.len() >= 2 && {
                let s0 = R3::from_vec3(shared_f[0]);
                let s1 = R3::from_vec3(shared_f[1]);
                point_on_segment(&x, &s0, &s1) && point_on_segment(&y, &s0, &s1)
            };
            if benign {
                stats.full_seg_benign += 1;
            }
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
