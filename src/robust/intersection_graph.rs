// robust/intersection_graph.rs — From two triangle soups to classified-ready
// pieces (paper §6).
//
// Pipeline stage between the narrow phase (robust/tri_tri.rs) and
// classification (robust/cells.rs):
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
//
// This file holds the build pipeline only. Its vocabulary types live in
// robust/graph_types.rs (edge keys, `VertInterner`, `Piece`,
// `IntersectionGraph`), its geometric helpers in robust/graph_geom.rs, and
// the same-mesh narrow phase in robust/graph_self_cut.rs — all re-exported
// here so callers keep using `intersection_graph::` paths.

use std::collections::BTreeSet;

// Fx hashing instead of SipHash. Every map/set below is probe-only or has an
// order-invariant consumer (documented per site); the hasher is unseeded, so
// even iteration order is stable across runs — output cannot depend on it.
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};

use crate::linalg::Vec3;
use crate::types::Box;

use super::arrangement::{self, ArrangementInput};
use super::exact::rational::{r3_eq, R3, R3Key};
use super::tri_tri::{tri_tri_intersect, TriTriIsect};

use super::graph_geom::{
    approx3, box3_contains, clip_segment_to_polygon, point_in_polygon_coplanar, point_on_segment_f,
    seg_box3,
};
use super::graph_types::{bit_edge_key, geo_edge_key, BitEdgeKey, GeoEdgeKey};

// `tri_box` / `is_degenerate` / `real_self_contact` / `SelfCutStats` stay
// crate-internal (robust/soup.rs reaches them through this path).
pub(super) use super::graph_geom::{is_degenerate, tri_box};
pub(super) use super::graph_self_cut::{real_self_contact, SelfCutStats};
pub use super::graph_types::{edge_key, EdgeKey, IntersectionGraph, Piece, VertInterner};

/// A pair's primitives after distribution: segments (including coplanar
/// boundary edges) and isolated points.
#[derive(Clone, Debug, Default)]
struct TriPrims {
    points: Vec<(R3, usize)>,
    segments: Vec<(R3, R3, usize)>,
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
    build_graph_with_progress(p, q, token, None)
}

/// [`build_graph_with_token`] that also reports its five phases to
/// `progress` (see [`crate::progress`]). `None` is exactly
/// [`build_graph_with_token`]: no counter is touched and no branch is taken
/// inside any inner loop.
pub fn build_graph_with_progress(
    p: &[[Vec3; 3]],
    q: &[[Vec3; 3]],
    token: Option<&crate::cancel::CancelToken>,
    progress: Option<&crate::progress::ProgressReporter>,
) -> Option<IntersectionGraph> {
    use crate::progress::{begin_phase, maybe_par_map_ct_progress, Phase};
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
    begin_phase(progress, Phase::NarrowPhase, p.len() as u64);
    for (pi, pt) in p.iter().enumerate() {
        if cancelled() {
            return None;
        }
        if let Some(pr) = progress {
            pr.advance(1);
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
    begin_phase(
        progress,
        Phase::SelfIntersections,
        (p.len() + q.len()) as u64,
    );
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
        let contact_results = maybe_par_map_ct_progress(tris.len(), 64, token, progress, |i| {
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
    begin_phase(progress, Phase::CandidatePoints, cand_work.len() as u64);
    let cand_results = maybe_par_map_ct_progress(cand_work.len(), 16, token, progress, |i| {
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
    // Two sweeps (original edges, then intersection segments) over the same
    // worklist, so the phase total counts it twice.
    begin_phase(progress, Phase::Registries, 2 * reg_work.len() as u64);

    // Registry values are (points, seen) rather than BTreeSet: probes and
    // dedup go through structural R3Key hashing, so the sequential merge
    // never pays ordered rational comparisons, and the consuming `extra`
    // sets see each point exactly once — the same content BTreeSet gave.
    // Order invariance: the registry is only ever `get`/`entry`-probed, and
    // its point lists land in a `BTreeSet<R3>` (`extra`) below, which sorts.
    let mut edge_registry: [HashMap<BitEdgeKey, (Vec<R3>, HashSet<R3Key>)>; 2] =
        [HashMap::default(), HashMap::default()];
    let edge_hits = maybe_par_map_ct_progress(reg_work.len(), 16, token, progress, |i| {
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
    let split_hits = maybe_par_map_ct_progress(reg_work.len(), 16, token, progress, |i| {
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
    begin_phase(progress, Phase::Arrangements, arr_work.len() as u64);
    let arr_results = maybe_par_map_ct_progress(arr_work.len(), 16, token, progress, |i| {
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
