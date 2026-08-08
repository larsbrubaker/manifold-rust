// robust/arrangement.rs — Per-triangle 2D arrangement (paper §6.3–6.4).
//
// Each input triangle that intersects triangles of the other operand gets an
// arrangement: the exact intersection primitives that landed on it (points,
// segments, coplanar-overlap polygon boundaries) are projected into the
// triangle's dominant-axis plane, split at their mutual crossings, and handed
// to robust/cdt.rs as constraints. The result subdivides the triangle into
// sub-triangles whose edges preserve every intersection segment, with a
// provenance map from constraint edges back to the primitives that created
// them — the edge→triangles incidence the classification stage
// (robust/classify.rs) walks.
//
// Everything here is exact: projection is coordinate dropping (bijective on
// the triangle's plane), crossings are rational constructions, and identity
// is rational equality — no tolerances anywhere.

use std::collections::BTreeMap;

use num_traits::Signed;

use crate::linalg::Vec3;

use super::cdt;
use super::exact::approx::orient2d_a;
use super::exact::predicates::{homog2_of, line_line_intersect_2d, orient2d_h, point_in_tri_2d, tri_normal_r, Homog2, TriLoc};
use super::exact::rational::{rat_to_f64, R2, R2Key, R3};
use super::exact::Sign;
use super::tri_tri::{dominant_axis, lift_to_plane};

/// Intersection primitives to arrange on one triangle, each tagged with a
/// caller-defined provenance id (the robust pipeline uses the index of the
/// opposing triangle pair that produced it).
#[derive(Clone, Debug, Default)]
pub struct ArrangementInput {
    /// Isolated contact points (vertex touches, edge-through-edge points).
    pub points: Vec<(R3, usize)>,
    /// Intersection segments, including coplanar-overlap polygon edges.
    /// Endpoints must be distinct and lie inside or on the triangle.
    pub segments: Vec<(R3, R3, usize)>,
}

/// The subdivided triangle.
#[derive(Clone, Debug)]
pub struct Arrangement {
    /// Dropped coordinate of the projection (0=x, 1=y, 2=z).
    pub axis: usize,
    /// Exact 3D points; indices 0..3 are the triangle corners in input order.
    pub points3: Vec<R3>,
    /// Their exact 2D projections (same indexing).
    pub points2: Vec<R2>,
    /// Sub-triangles (indices into points*), CCW in projection space.
    pub tris: Vec<[usize; 3]>,
    /// Constraint edges as (min,max) index pairs → provenance ids of every
    /// primitive that generated the edge.
    pub constraints: BTreeMap<(usize, usize), Vec<usize>>,
    /// True when the 2D projection mirrors the triangle's 3D winding (the
    /// dropped normal component is negative): consumers must swap two
    /// indices of each sub-triangle to recover outward orientation.
    pub flipped: bool,
}


/// Aggregate phase timers for build(), printed under MANIFOLD_TIMING by the
/// pipeline after the arrangements stage. Relaxed atomics, nanoseconds.
pub mod stats {
    use std::sync::atomic::{AtomicU64, Ordering::Relaxed};

    pub static SETUP_NS: AtomicU64 = AtomicU64::new(0);
    pub static NORM_NS: AtomicU64 = AtomicU64::new(0);
    pub static CROSS_NS: AtomicU64 = AtomicU64::new(0);
    pub static ONSEG_NS: AtomicU64 = AtomicU64::new(0);
    pub static CDT_NS: AtomicU64 = AtomicU64::new(0);
    pub static CALLS: AtomicU64 = AtomicU64::new(0);
    pub static SEGS: AtomicU64 = AtomicU64::new(0);
    pub static PTS: AtomicU64 = AtomicU64::new(0);

    pub fn snapshot_and_reset() -> String {
        let take = |a: &AtomicU64| a.swap(0, Relaxed) as f64 * 1e-9;
        let taken = |a: &AtomicU64| a.swap(0, Relaxed);
        format!(
            "setup {:.3}s (norm {:.3}s), crossings {:.3}s, on-seg {:.3}s, cdt {:.3}s ({} calls, {} segs, {} pts)",
            take(&SETUP_NS),
            take(&NORM_NS),
            take(&CROSS_NS),
            take(&ONSEG_NS),
            take(&CDT_NS),
            taken(&CALLS),
            taken(&SEGS),
            taken(&PTS),
        )
    }
}

/// Build the arrangement of `input` on triangle `tri`. Primitives must
/// already be clipped to the triangle (tri_tri output is), with rational
/// coordinates on its plane.
///
/// Returns `None` only when `token` is cancelled: a single segment-heavy
/// triangle can spend minutes in the quadratic sweeps below, so the
/// per-phase checks in build_graph alone leave a cancel unanswered for
/// however long the worst triangle takes.
pub fn build(
    tri: [Vec3; 3],
    input: &ArrangementInput,
    token: Option<&crate::cancel::CancelToken>,
) -> Option<Arrangement> {
    use std::sync::atomic::Ordering::Relaxed;
    let t0 = crate::timing::Stopwatch::start();
    stats::CALLS.fetch_add(1, Relaxed);
    stats::SEGS.fetch_add(input.segments.len() as u64, Relaxed);
    stats::PTS.fetch_add(input.points.len() as u64, Relaxed);
    let corners: [R3; 3] = [
        R3::from_vec3(tri[0]),
        R3::from_vec3(tri[1]),
        R3::from_vec3(tri[2]),
    ];
    let normal = tri_normal_r(&corners[0], &corners[1], &corners[2]);
    debug_assert!(!normal.is_zero(), "degenerate triangle in arrangement");
    let axis = dominant_axis(&normal);
    stats::NORM_NS.fetch_add(t0.elapsed_ns(), Relaxed);

    let mut points3: Vec<R3> = Vec::new();
    let mut points2: Vec<R2> = Vec::new();
    // Hash-keyed dedup via R2Key: division-free structural hashing of the
    // canonical rationals (num-rational's own Hash runs a Euclidean
    // recursion per lookup, which dominated this whole function). Indices
    // are assigned in insertion order, so output stays deterministic.
    let mut index: std::collections::HashMap<R2Key, usize> = std::collections::HashMap::new();
    let mut add_point = |p3: R3, points3: &mut Vec<R3>, points2: &mut Vec<R2>| -> usize {
        let p2 = p3.project_drop(axis);
        let next = points3.len();
        match index.entry(R2Key(p2.clone())) {
            std::collections::hash_map::Entry::Occupied(e) => *e.get(),
            std::collections::hash_map::Entry::Vacant(e) => {
                e.insert(next);
                points3.push(p3);
                points2.push(p2);
                next
            }
        }
    };

    for c in &corners {
        add_point(c.clone(), &mut points3, &mut points2);
    }
    debug_assert_eq!(points3.len(), 3, "corner points must be distinct");

    // Register primitive endpoints / isolated points.
    struct Seg {
        a: usize,
        b: usize,
        prov: usize,
    }
    let mut segs: Vec<Seg> = Vec::with_capacity(input.segments.len());
    for (a3, b3, prov) in &input.segments {
        debug_assert_ne!(a3, b3, "zero-length segment primitive");
        let a = add_point(a3.clone(), &mut points3, &mut points2);
        let b = add_point(b3.clone(), &mut points3, &mut points2);
        debug_assert_ne!(a, b);
        segs.push(Seg {
            a,
            b,
            prov: *prov,
        });
    }
    for (p3, _prov) in &input.points {
        add_point(p3.clone(), &mut points3, &mut points2);
    }

    stats::SETUP_NS.fetch_add(t0.elapsed_ns(), Relaxed);
    let t0 = crate::timing::Stopwatch::start();

    // Mutual proper crossings between segments become new points. Points are
    // homogenized once (Homog2) for the exact fallback, and approximated
    // once (correctly rounded f64) for the semi-static filters that certify
    // generic-position signs without any BigInt work.
    let mut homogs: Vec<Homog2> = points2.iter().map(homog2_of).collect();
    let approx = |p: &R2| -> [f64; 2] { [rat_to_f64(&p.x), rat_to_f64(&p.y)] };
    let mut apts: Vec<[f64; 2]> = points2.iter().map(approx).collect();
    let o2 = |apts: &[[f64; 2]], homogs: &[Homog2], i: usize, j: usize, k: usize| -> Sign {
        orient2d_a(apts[i], apts[j], apts[k])
            .unwrap_or_else(|| orient2d_h(&homogs[i], &homogs[j], &homogs[k]))
    };
    for i in 0..segs.len() {
        if crate::cancel::is_cancelled(token) {
            return None;
        }
        for j in (i + 1)..segs.len() {
            let (ia, ib) = (segs[i].a, segs[i].b);
            let (ic, id) = (segs[j].a, segs[j].b);
            let sc = o2(&apts, &homogs, ia, ib, ic);
            let sd = o2(&apts, &homogs, ia, ib, id);
            let sa = o2(&apts, &homogs, ic, id, ia);
            let sb = o2(&apts, &homogs, ic, id, ib);
            // Strict crossing only — endpoint contacts and collinear overlap
            // are handled by the point-on-segment sweep below.
            if sc != Sign::Zero && sd != Sign::Zero && sc != sd
                && sa != Sign::Zero && sb != Sign::Zero && sa != sb
            {
                let x2 = line_line_intersect_2d(
                    &points2[segs[i].a],
                    &points2[segs[i].b],
                    &points2[segs[j].a],
                    &points2[segs[j].b],
                )
                .expect("properly crossing segments are not parallel");
                let x3 = lift_to_plane(&x2, axis, &corners[0], &normal);
                add_point(x3, &mut points3, &mut points2);
            }
        }
    }
    // Points added by crossings need homogs and approximations too.
    for p in points2.iter().skip(homogs.len()) {
        homogs.push(homog2_of(p));
    }
    for p in points2.iter().skip(apts.len()) {
        apts.push(approx(p));
    }

    debug_assert!(points2.iter().all(|p| {
        point_in_tri_2d(p, &points2[0], &points2[1], &points2[2]) != TriLoc::Outside
    }), "arrangement primitive escapes its triangle");

    stats::CROSS_NS.fetch_add(t0.elapsed_ns(), Relaxed);
    let t0 = crate::timing::Stopwatch::start();

    // Subdivide each segment at every registered point lying exactly on it;
    // consecutive point pairs become constraint edges carrying provenance.
    let mut constraints: BTreeMap<(usize, usize), Vec<usize>> = BTreeMap::new();
    for seg in &segs {
        if crate::cancel::is_cancelled(token) {
            return None;
        }
        let (ha, hb) = (&homogs[seg.a], &homogs[seg.b]);
        // Segment direction cleared of denominators: v = (pb−pa)·(Bw·Aw).
        let vx = &hb.0 * &ha.2 - &ha.0 * &hb.2;
        let vy = &hb.1 * &ha.2 - &ha.1 * &hb.2;
        let vv = &vx * &vx + &vy * &vy;
        // Parameters as unreduced fractions: for point P (homog (Px,Py,Pw)),
        // u = (p−pa)·(Pw·Aw) and t ∝ (u·v)/Pw — ordered and range-checked by
        // cross-multiplication with the positive denominators, no canonical
        // rationals anywhere.
        let mut on_seg: Vec<(num_bigint::BigInt, num_bigint::BigInt, usize)> = Vec::new();
        for (idx, hp) in homogs.iter().enumerate() {
            // Approx filter first: almost every point is certifiably off the
            // segment's line; only near-collinear candidates pay for BigInt.
            match orient2d_a(apts[seg.a], apts[seg.b], apts[idx]) {
                Some(_) => continue, // certified nonzero → not collinear
                None => {}
            }
            if orient2d_h(ha, hb, hp) != Sign::Zero {
                continue;
            }
            let ux = &hp.0 * &ha.2 - &ha.0 * &hp.2;
            let uy = &hp.1 * &ha.2 - &ha.1 * &hp.2;
            let uv = &ux * &vx + &uy * &vy;
            // 0 ≤ t  ⇔  0 ≤ u·v;   t ≤ |dir|²  ⇔  (u·v)·Bw ≤ (v·v)·Pw.
            if !uv.is_negative() && &uv * &hb.2 <= &vv * &hp.2 {
                on_seg.push((uv, hp.2.clone(), idx));
            }
        }
        // t_i < t_j  ⇔  uv_i·Pw_j < uv_j·Pw_i (denominators positive).
        on_seg.sort_by(|a, b| (&a.0 * &b.1).cmp(&(&b.0 * &a.1)));
        debug_assert!(on_seg.len() >= 2, "segment lost its own endpoints");
        for w in on_seg.windows(2) {
            let (u, v) = (w[0].2, w[1].2);
            debug_assert_ne!(u, v, "duplicate points on segment");
            let key = (u.min(v), u.max(v));
            let provs = constraints.entry(key).or_default();
            if !provs.contains(&seg.prov) {
                provs.push(seg.prov);
            }
        }
    }

    stats::ONSEG_NS.fetch_add(t0.elapsed_ns(), Relaxed);
    let t0 = crate::timing::Stopwatch::start();

    let constraint_pairs: Vec<(usize, usize)> = constraints.keys().copied().collect();
    let tris = cdt::triangulate(&points2, &constraint_pairs);
    stats::CDT_NS.fetch_add(t0.elapsed_ns(), Relaxed);

    let axis_comp = match axis {
        0 => &normal.x,
        1 => &normal.y,
        _ => &normal.z,
    };
    let flipped = Sign::of_rat(axis_comp) == Sign::Neg;

    Some(Arrangement {
        axis,
        points3,
        points2,
        tris,
        constraints,
        flipped,
    })
}

/// The exact points this input would register on `tri`, without running the
/// CDT: primitive endpoints, isolated points, and pairwise proper segment
/// crossings. robust/intersection_graph.rs uses this to build the global
/// registries that force a common subdivision on both sides of every shared
/// intersection segment and mesh edge before the real arrangements run.
/// Returns `None` only when `token` is cancelled (see [`build`]).
pub fn candidate_points(
    tri: [Vec3; 3],
    input: &ArrangementInput,
    token: Option<&crate::cancel::CancelToken>,
) -> Option<Vec<R3>> {
    let corners: [R3; 3] = [
        R3::from_vec3(tri[0]),
        R3::from_vec3(tri[1]),
        R3::from_vec3(tri[2]),
    ];
    let normal = tri_normal_r(&corners[0], &corners[1], &corners[2]);
    let axis = dominant_axis(&normal);

    let mut out: Vec<R3> = Vec::new();
    // Membership-only set with division-free R2Key hashing; `out` keeps
    // insertion order so determinism is unaffected.
    let mut seen: std::collections::HashSet<R2Key> = std::collections::HashSet::new();
    let add = |p3: R3, out: &mut Vec<R3>, seen: &mut std::collections::HashSet<R2Key>| {
        if seen.insert(R2Key(p3.project_drop(axis))) {
            out.push(p3);
        }
    };
    for (p, _) in &input.points {
        add(p.clone(), &mut out, &mut seen);
    }
    for (a, b, _) in &input.segments {
        add(a.clone(), &mut out, &mut seen);
        add(b.clone(), &mut out, &mut seen);
    }
    // Same filtered predicate stack as build(): correctly rounded f64
    // approximations certify almost every non-crossing pair, homogenized
    // exact signs only on near-degeneracies — the all-rational version made
    // this sweep a dominant pipeline stage on segment-heavy meshes.
    let segs2: Vec<(R2, R2)> = input
        .segments
        .iter()
        .map(|(a, b, _)| (a.project_drop(axis), b.project_drop(axis)))
        .collect();
    let homogs: Vec<(Homog2, Homog2)> = segs2
        .iter()
        .map(|(a, b)| (homog2_of(a), homog2_of(b)))
        .collect();
    let approx = |p: &R2| -> [f64; 2] { [rat_to_f64(&p.x), rat_to_f64(&p.y)] };
    let apts: Vec<([f64; 2], [f64; 2])> = segs2
        .iter()
        .map(|(a, b)| (approx(a), approx(b)))
        .collect();
    let o2 = |a: ([f64; 2], &Homog2), b: ([f64; 2], &Homog2), c: ([f64; 2], &Homog2)| -> Sign {
        orient2d_a(a.0, b.0, c.0).unwrap_or_else(|| orient2d_h(a.1, b.1, c.1))
    };
    for i in 0..segs2.len() {
        if crate::cancel::is_cancelled(token) {
            return None;
        }
        for j in (i + 1)..segs2.len() {
            let (a, b) = (
                (apts[i].0, &homogs[i].0),
                (apts[i].1, &homogs[i].1),
            );
            let (c, d) = (
                (apts[j].0, &homogs[j].0),
                (apts[j].1, &homogs[j].1),
            );
            let sc = o2(a, b, c);
            let sd = o2(a, b, d);
            let sa = o2(c, d, a);
            let sb = o2(c, d, b);
            if sc != Sign::Zero && sd != Sign::Zero && sc != sd
                && sa != Sign::Zero && sb != Sign::Zero && sa != sb
            {
                let x2 = line_line_intersect_2d(&segs2[i].0, &segs2[i].1, &segs2[j].0, &segs2[j].1)
                    .expect("properly crossing segments are not parallel");
                add(
                    lift_to_plane(&x2, axis, &corners[0], &normal),
                    &mut out,
                    &mut seen,
                );
            }
        }
    }
    Some(out)
}

#[cfg(test)]
#[path = "arrangement_tests.rs"]
mod tests;
