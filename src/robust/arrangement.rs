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
use super::exact::predicates::{homog2_of, line_line_intersect_2d, orient2d_h, orient2d_r, point_in_tri_2d, tri_normal_r, Homog2, TriLoc};
use super::exact::rational::{rat_to_f64, R2, R3};
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


/// Build the arrangement of `input` on triangle `tri`. Primitives must
/// already be clipped to the triangle (tri_tri output is), with rational
/// coordinates on its plane.
pub fn build(tri: [Vec3; 3], input: &ArrangementInput) -> Arrangement {
    let corners: [R3; 3] = [
        R3::from_vec3(tri[0]),
        R3::from_vec3(tri[1]),
        R3::from_vec3(tri[2]),
    ];
    let normal = tri_normal_r(&corners[0], &corners[1], &corners[2]);
    debug_assert!(!normal.is_zero(), "degenerate triangle in arrangement");
    let axis = dominant_axis(&normal);

    let mut points3: Vec<R3> = Vec::new();
    let mut points2: Vec<R2> = Vec::new();
    // Hash-keyed dedup: canonical rationals hash linearly, while a BTreeMap's
    // Ord runs two BigInt cross-multiplications per tree level per lookup.
    // Indices are assigned in insertion order, so output stays deterministic.
    let mut index: std::collections::HashMap<R2, usize> = std::collections::HashMap::new();
    let mut add_point = |p3: R3, points3: &mut Vec<R3>, points2: &mut Vec<R2>| -> usize {
        let p2 = p3.project_drop(axis);
        if let Some(&i) = index.get(&p2) {
            return i;
        }
        let i = points3.len();
        index.insert(p2.clone(), i);
        points3.push(p3);
        points2.push(p2);
        i
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

    // Subdivide each segment at every registered point lying exactly on it;
    // consecutive point pairs become constraint edges carrying provenance.
    let mut constraints: BTreeMap<(usize, usize), Vec<usize>> = BTreeMap::new();
    for seg in &segs {
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

    let constraint_pairs: Vec<(usize, usize)> = constraints.keys().copied().collect();
    let tris = cdt::triangulate(&points2, &constraint_pairs);

    let axis_comp = match axis {
        0 => &normal.x,
        1 => &normal.y,
        _ => &normal.z,
    };
    let flipped = Sign::of_rat(axis_comp) == Sign::Neg;

    Arrangement {
        axis,
        points3,
        points2,
        tris,
        constraints,
        flipped,
    }
}

/// The exact points this input would register on `tri`, without running the
/// CDT: primitive endpoints, isolated points, and pairwise proper segment
/// crossings. robust/intersection_graph.rs uses this to build the global
/// registries that force a common subdivision on both sides of every shared
/// intersection segment and mesh edge before the real arrangements run.
pub fn candidate_points(tri: [Vec3; 3], input: &ArrangementInput) -> Vec<R3> {
    let corners: [R3; 3] = [
        R3::from_vec3(tri[0]),
        R3::from_vec3(tri[1]),
        R3::from_vec3(tri[2]),
    ];
    let normal = tri_normal_r(&corners[0], &corners[1], &corners[2]);
    let axis = dominant_axis(&normal);

    let mut out: Vec<R3> = Vec::new();
    // Membership-only set: hashing beats BTreeSet's cross-multiplying Ord,
    // and `out` keeps insertion order so determinism is unaffected.
    let mut seen: std::collections::HashSet<R2> = std::collections::HashSet::new();
    let add = |p3: R3, out: &mut Vec<R3>, seen: &mut std::collections::HashSet<R2>| {
        if seen.insert(p3.project_drop(axis)) {
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
    let segs2: Vec<(R2, R2)> = input
        .segments
        .iter()
        .map(|(a, b, _)| (a.project_drop(axis), b.project_drop(axis)))
        .collect();
    for i in 0..segs2.len() {
        for j in (i + 1)..segs2.len() {
            let (a, b) = (&segs2[i].0, &segs2[i].1);
            let (c, d) = (&segs2[j].0, &segs2[j].1);
            let sc = orient2d_r(a, b, c);
            let sd = orient2d_r(a, b, d);
            let sa = orient2d_r(c, d, a);
            let sb = orient2d_r(c, d, b);
            if sc != Sign::Zero && sd != Sign::Zero && sc != sd
                && sa != Sign::Zero && sb != Sign::Zero && sa != sb
            {
                let x2 = line_line_intersect_2d(a, b, c, d)
                    .expect("properly crossing segments are not parallel");
                add(
                    lift_to_plane(&x2, axis, &corners[0], &normal),
                    &mut out,
                    &mut seen,
                );
            }
        }
    }
    out
}

#[cfg(test)]
#[path = "arrangement_tests.rs"]
mod tests;
