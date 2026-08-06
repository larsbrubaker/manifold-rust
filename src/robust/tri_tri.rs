// robust/tri_tri.rs — Exact triangle-triangle intersection for the robust
// boolean engine.
//
// Narrow phase behind the Collider broad phase: given one triangle from each
// operand mesh, classify their intersection exactly as nothing, a single
// point, a segment, or (for coplanar pairs) a convex overlap polygon. All
// vertex-vs-plane tests go through the filtered predicates in
// robust/exact/filtered.rs; every constructed point is exact rational
// (robust/exact/predicates.rs), so downstream arrangements
// (robust/arrangement.rs) never see rounded coordinates.
//
// Degenerate (zero-area) input triangles are the caller's responsibility to
// drop beforehand (paper §5 pre-processing); this file debug-asserts that.

use num_rational::BigRational;
use num_traits::{Signed, Zero};

use crate::linalg::Vec3;

use super::exact::filtered::orient3d;
use super::exact::predicates::{line_line_intersect_2d, line_plane_intersect, orient2d_r, tri_normal_r};
use super::exact::rational::{R2, R3};
use super::exact::Sign;

/// Exact intersection of two triangles.
#[derive(Clone, Debug, PartialEq)]
pub enum TriTriIsect {
    None,
    /// Single-point contact (vertex-on-face, vertex-on-edge, edge-through-
    /// edge, or interval intersection collapsing to one point).
    Point(R3),
    /// Proper crossing (or edge/vertex contact with positive length).
    Segment(R3, R3),
    /// Coplanar triangles overlapping with positive area. `polygon` is the
    /// convex overlap region (distinct vertices, no three collinear, no
    /// guaranteed winding); `same_orientation` is true when the two
    /// triangles' normals point the same way, false for opposite planes.
    Coplanar {
        polygon: Vec<R3>,
        same_orientation: bool,
    },
}

/// Dominant-axis choice for the paper's bijective drop-one-coordinate
/// projection: the axis of the exactly-largest |normal| component (ties
/// broken toward z, then y). The chosen component is guaranteed nonzero for
/// a non-degenerate triangle.
pub fn dominant_axis(n: &R3) -> usize {
    let ax = n.x.abs();
    let ay = n.y.abs();
    let az = n.z.abs();
    if az >= ax && az >= ay {
        2
    } else if ay >= ax {
        1
    } else {
        0
    }
}

// Re-exported for the existing call sites; the implementation lives with the
// other integer-only constructions in robust/exact/predicates.rs.
pub use super::exact::predicates::lift_to_plane;

/// Exit-path counters for perf analysis, printed under MANIFOLD_TIMING by
/// the self-cut loop. Relaxed atomics; negligible cost on the hot path.
pub mod stats {
    use std::sync::atomic::{AtomicU64, Ordering::Relaxed};

    pub static PLANE_REJECT: AtomicU64 = AtomicU64::new(0);
    pub static COPLANAR: AtomicU64 = AtomicU64::new(0);
    pub static COPLANAR_SAT: AtomicU64 = AtomicU64::new(0);
    pub static SAT_REJECT: AtomicU64 = AtomicU64::new(0);
    pub static INTERVAL: AtomicU64 = AtomicU64::new(0);
    pub static COPLANAR_NS: AtomicU64 = AtomicU64::new(0);
    pub static PLANE_NS: AtomicU64 = AtomicU64::new(0);
    pub static INTERVAL_NS: AtomicU64 = AtomicU64::new(0);

    pub fn snapshot_and_reset() -> String {
        let take = |a: &AtomicU64| a.swap(0, Relaxed);
        format!(
            "plane-reject {} ({:.3}s signs), coplanar {} (sat {}, {:.3}s), sat-reject {}, interval {} ({:.3}s)",
            take(&PLANE_REJECT),
            take(&PLANE_NS) as f64 * 1e-9,
            take(&COPLANAR),
            take(&COPLANAR_SAT),
            take(&COPLANAR_NS) as f64 * 1e-9,
            take(&SAT_REJECT),
            take(&INTERVAL),
            take(&INTERVAL_NS) as f64 * 1e-9,
        )
    }
}

/// Exact intersection of triangles t1 and t2 (each three finite f64
/// vertices). Symmetric: swapping the arguments yields the same set.
pub fn tri_tri_intersect(t1: [Vec3; 3], t2: [Vec3; 3]) -> TriTriIsect {
    use std::sync::atomic::Ordering::Relaxed;
    let t_signs = std::time::Instant::now();
    // Signs of t2's vertices against t1's plane.
    let s2 = [
        orient3d(t1[0], t1[1], t1[2], t2[0]),
        orient3d(t1[0], t1[1], t1[2], t2[1]),
        orient3d(t1[0], t1[1], t1[2], t2[2]),
    ];
    if all_same_strict(&s2) {
        stats::PLANE_REJECT.fetch_add(1, Relaxed);
        stats::PLANE_NS.fetch_add(t_signs.elapsed().as_nanos() as u64, Relaxed);
        return TriTriIsect::None;
    }
    if s2.iter().all(|s| *s == Sign::Zero) {
        stats::COPLANAR.fetch_add(1, Relaxed);
        let out = coplanar_overlap(t1, t2);
        stats::COPLANAR_NS.fetch_add(t_signs.elapsed().as_nanos() as u64, Relaxed);
        return out;
    }
    // Signs of t1's vertices against t2's plane.
    let s1 = [
        orient3d(t2[0], t2[1], t2[2], t1[0]),
        orient3d(t2[0], t2[1], t2[2], t1[1]),
        orient3d(t2[0], t2[1], t2[2], t1[2]),
    ];
    if all_same_strict(&s1) {
        stats::PLANE_REJECT.fetch_add(1, Relaxed);
        stats::PLANE_NS.fetch_add(t_signs.elapsed().as_nanos() as u64, Relaxed);
        return TriTriIsect::None;
    }
    stats::PLANE_NS.fetch_add(t_signs.elapsed().as_nanos() as u64, Relaxed);
    debug_assert!(
        !s1.iter().all(|s| *s == Sign::Zero),
        "t1 coplanar with t2's plane implies t2 coplanar with t1's — handled above"
    );

    // Both triangles straddle each other's plane, but most such box-pair
    // candidates still miss along the common line. A certified
    // separating-axis check on the raw f64 vertices skips the entire
    // rational interval construction for them.
    {
        let f1 = [
            [t1[0].x, t1[0].y, t1[0].z],
            [t1[1].x, t1[1].y, t1[1].z],
            [t1[2].x, t1[2].y, t1[2].z],
        ];
        let f2 = [
            [t2[0].x, t2[0].y, t2[0].z],
            [t2[1].x, t2[1].y, t2[1].z],
            [t2[2].x, t2[2].y, t2[2].z],
        ];
        if super::exact::approx::sat_edge_axes_disjoint(&f1, &f2) {
            stats::SAT_REJECT.fetch_add(1, Relaxed);
            return TriTriIsect::None;
        }
    }
    stats::INTERVAL.fetch_add(1, Relaxed);
    let t_interval = std::time::Instant::now();

    // Both triangles meet the common line L of the two planes. Overlap the
    // two 1- or 2-point intervals along L entirely in scaled integer
    // arithmetic (no rational constructions, no gcds): per-axis power-of-two
    // scaling maps every vertex to an exact integer, and for points on L the
    // scaled-space parameter dir_s·x_s is a positive multiple of the true
    // parameter dir·x (x−y ∥ dir makes the difference d·|A·dir|²·λ for
    // x−y = λ·dir), so ordering — including its orientation — matches the
    // rational computation exactly. Endpoints stay symbolic; only the 1–2
    // points of the final answer are constructed rationally.
    let out = interval_overlap(t1, t2, &s1, &s2);
    stats::INTERVAL_NS.fetch_add(t_interval.elapsed().as_nanos() as u64, Relaxed);
    out
}

/// Symbolic interval endpoint on the common line L: an original vertex
/// exactly on the other plane, or a strictly straddling edge's crossing.
#[derive(Clone, Copy)]
enum EndPt {
    /// (which_tri: 0|1, vertex index)
    Vert(u8, u8),
    /// (which_tri: 0|1, edge start index i — the edge runs i → (i+1)%3)
    Cross(u8, u8),
}

fn interval_overlap(t1: [Vec3; 3], t2: [Vec3; 3], s1: &[Sign; 3], s2: &[Sign; 3]) -> TriTriIsect {
    use num_bigint::BigInt;

    // Scaled integer coordinates; one common scale per axis across BOTH
    // triangles so cross-triangle parameter comparisons share a basis.
    let sx = super::exact::intpred::scaled_big([
        t1[0].x, t1[1].x, t1[2].x, t2[0].x, t2[1].x, t2[2].x,
    ]);
    let sy = super::exact::intpred::scaled_big([
        t1[0].y, t1[1].y, t1[2].y, t2[0].y, t2[1].y, t2[2].y,
    ]);
    let sz = super::exact::intpred::scaled_big([
        t1[0].z, t1[1].z, t1[2].z, t2[0].z, t2[1].z, t2[2].z,
    ]);
    let v = |k: usize| [&sx[k], &sy[k], &sz[k]];
    let sub = |a: [&BigInt; 3], b: [&BigInt; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let cross = |a: &[BigInt; 3], b: &[BigInt; 3]| {
        [
            &a[1] * &b[2] - &a[2] * &b[1],
            &a[2] * &b[0] - &a[0] * &b[2],
            &a[0] * &b[1] - &a[1] * &b[0],
        ]
    };
    let dot = |a: &[BigInt; 3], b: [&BigInt; 3]| &a[0] * b[0] + &a[1] * b[1] + &a[2] * b[2];

    let n1 = cross(&sub(v(1), v(0)), &sub(v(2), v(0)));
    let n2 = cross(&sub(v(4), v(3)), &sub(v(5), v(3)));
    let dir = cross(&n1, &n2);
    debug_assert!(
        dir.iter().any(|c| !c.is_zero()),
        "non-coplanar intersecting planes"
    );

    // Parameters dir·v of all six vertices, and signed heights against the
    // other triangle's plane (whose signs replicate s1/s2 exactly).
    let du: [BigInt; 6] = std::array::from_fn(|k| dot(&dir, v(k)));
    let h = |k: usize, n: &[BigInt; 3], origin: usize| dot(n, v(k)) - dot(n, v(origin));
    let h1: [BigInt; 3] = std::array::from_fn(|i| h(i, &n2, 3));
    let h2: [BigInt; 3] = std::array::from_fn(|i| h(3 + i, &n1, 0));
    #[cfg(debug_assertions)]
    for i in 0..3 {
        debug_assert_eq!(int_sign(&h1[i]), s1[i], "scaled height disagrees with s1");
        debug_assert_eq!(int_sign(&h2[i]), s2[i], "scaled height disagrees with s2");
    }

    // The ≤2 endpoints of one triangle's crossing with the other's plane, as
    // (unreduced parameter fraction, symbolic point), in the same
    // enumeration order as the rational implementation used (vertices in
    // index order, then edges (0,1), (1,2), (2,0)) so ties break alike.
    let endpoints = |which: u8, s: &[Sign; 3], hs: &[BigInt; 3]| -> Vec<(Frac, EndPt)> {
        let base = if which == 0 { 0 } else { 3 };
        let mut pts = Vec::with_capacity(2);
        for i in 0..3 {
            if s[i] == Sign::Zero {
                pts.push((
                    (du[base + i].clone(), BigInt::from(1)),
                    EndPt::Vert(which, i as u8),
                ));
            }
        }
        for i in 0..3 {
            let j = (i + 1) % 3;
            if s[i] != Sign::Zero && s[j] != Sign::Zero && s[i] != s[j] {
                // x = u + h_u/(h_u−h_v)·(v−u) ⇒
                // dir·x = [(h_u−h_v)·du_u + h_u·(du_v−du_u)] / (h_u−h_v).
                let mut den = &hs[i] - &hs[j];
                let mut num = &den * &du[base + i] + &hs[i] * (&du[base + j] - &du[base + i]);
                if den.sign() == num_bigint::Sign::Minus {
                    den = -den;
                    num = -num;
                }
                pts.push(((num, den), EndPt::Cross(which, i as u8)));
            }
        }
        debug_assert!(!pts.is_empty() && pts.len() <= 2);
        pts
    };
    let pts1 = endpoints(0, s1, &h1);
    let pts2 = endpoints(1, s2, &h2);

    // Per-triangle interval, first-encountered point winning ties (matching
    // the old interval_along).
    let minmax = |pts: Vec<(Frac, EndPt)>| -> ((Frac, EndPt), (Frac, EndPt)) {
        let mut lo = pts[0].clone();
        let mut hi = pts[0].clone();
        for p in &pts[1..] {
            if cmp_frac(&p.0, &lo.0) == std::cmp::Ordering::Less {
                lo = p.clone();
            }
            if cmp_frac(&p.0, &hi.0) == std::cmp::Ordering::Greater {
                hi = p.clone();
            }
        }
        (lo, hi)
    };
    let i1 = minmax(pts1);
    let i2 = minmax(pts2);
    let (lo, lo_pt) = if cmp_frac(&i1.0 .0, &i2.0 .0) != std::cmp::Ordering::Less { i1.0 } else { i2.0 };
    let (hi, hi_pt) = if cmp_frac(&i1.1 .0, &i2.1 .0) != std::cmp::Ordering::Greater { i1.1 } else { i2.1 };

    match cmp_frac(&lo, &hi) {
        std::cmp::Ordering::Greater => TriTriIsect::None,
        std::cmp::Ordering::Equal => TriTriIsect::Point(build_endpoint(lo_pt, &t1, &t2)),
        std::cmp::Ordering::Less => TriTriIsect::Segment(
            build_endpoint(lo_pt, &t1, &t2),
            build_endpoint(hi_pt, &t1, &t2),
        ),
    }
}

#[cfg(debug_assertions)]
fn int_sign(v: &num_bigint::BigInt) -> Sign {
    match v.sign() {
        num_bigint::Sign::Minus => Sign::Neg,
        num_bigint::Sign::NoSign => Sign::Zero,
        num_bigint::Sign::Plus => Sign::Pos,
    }
}

/// Materialize a symbolic interval endpoint as the exact rational point the
/// fully rational implementation would have produced.
fn build_endpoint(e: EndPt, t1: &[Vec3; 3], t2: &[Vec3; 3]) -> R3 {
    let tri = |which: u8| if which == 0 { t1 } else { t2 };
    match e {
        EndPt::Vert(w, i) => R3::from_vec3(tri(w)[i as usize]),
        EndPt::Cross(w, i) => {
            let own = tri(w);
            let other = tri(1 - w);
            let a = R3::from_vec3(own[i as usize]);
            let b = R3::from_vec3(own[(i as usize + 1) % 3]);
            let p: [R3; 3] = [
                R3::from_vec3(other[0]),
                R3::from_vec3(other[1]),
                R3::from_vec3(other[2]),
            ];
            line_plane_intersect(&a, &b, &p[0], &p[1], &p[2])
                .expect("strictly straddling edge cannot be parallel to the plane")
        }
    }
}

/// Unreduced fraction with positive denominator.
type Frac = (num_bigint::BigInt, num_bigint::BigInt);

fn cmp_frac(a: &Frac, b: &Frac) -> std::cmp::Ordering {
    // Denominators positive → cross-multiplication preserves order.
    (&a.0 * &b.1).cmp(&(&b.0 * &a.1))
}

fn all_same_strict(s: &[Sign; 3]) -> bool {
    s[0] != Sign::Zero && s[0] == s[1] && s[1] == s[2]
}

// ─── Coplanar overlap ────────────────────────────────────────────────────────

/// Certified separating-edge pre-reject for coplanar pairs, on the raw f64
/// projection dropping coordinate `axis`. Sound for ANY choice of axis: a
/// projection is linear, so a shared 3D point would project into both
/// projected triangles — strict 2D separation therefore proves 3D
/// disjointness even when the projection degenerates the triangles. Signs
/// come from the exact filtered orient2d, so a `true` answer is certain.
fn coplanar_separated_2d(t1: [Vec3; 3], t2: [Vec3; 3], axis: usize) -> bool {
    use super::exact::filtered::orient2d;
    // Same cyclic drop-axis convention as R3::project_drop.
    let proj = |v: Vec3| match axis {
        0 => crate::linalg::Vec2::new(v.y, v.z),
        1 => crate::linalg::Vec2::new(v.z, v.x),
        _ => crate::linalg::Vec2::new(v.x, v.y),
    };
    let p1 = t1.map(proj);
    let p2 = t2.map(proj);
    let separates = |tri: &[crate::linalg::Vec2; 3], other: &[crate::linalg::Vec2; 3]| {
        (0..3).any(|i| {
            let a = tri[i];
            let b = tri[(i + 1) % 3];
            let s_ref = orient2d(a, b, tri[(i + 2) % 3]);
            s_ref != Sign::Zero
                && other.iter().all(|&q| {
                    let s = orient2d(a, b, q);
                    s != Sign::Zero && s != s_ref
                })
        })
    };
    separates(&p1, &p2) || separates(&p2, &p1)
}

/// Intersection of two coplanar triangles: Sutherland–Hodgman clip of t2
/// against t1 in the exact 2D projection, classified by the dimension of the
/// result (empty / point / segment / convex polygon).
fn coplanar_overlap(t1: [Vec3; 3], t2: [Vec3; 3]) -> TriTriIsect {
    {
        let n = crate::linalg::cross(t1[1] - t1[0], t1[2] - t1[0]);
        let (ax, ay, az) = (n.x.abs(), n.y.abs(), n.z.abs());
        let axis = if az >= ax && az >= ay { 2 } else if ay >= ax { 1 } else { 0 };
        if coplanar_separated_2d(t1, t2, axis) {
            stats::COPLANAR_SAT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            return TriTriIsect::None;
        }
    }
    let r1: [R3; 3] = [
        R3::from_vec3(t1[0]),
        R3::from_vec3(t1[1]),
        R3::from_vec3(t1[2]),
    ];
    let r2: [R3; 3] = [
        R3::from_vec3(t2[0]),
        R3::from_vec3(t2[1]),
        R3::from_vec3(t2[2]),
    ];
    let n1 = tri_normal_r(&r1[0], &r1[1], &r1[2]);
    let n2 = tri_normal_r(&r2[0], &r2[1], &r2[2]);
    debug_assert!(!n1.is_zero() && !n2.is_zero(), "degenerate input triangle");
    let same_orientation = match Sign::of_rat(&n1.dot(&n2)) {
        Sign::Pos => true,
        Sign::Neg => false,
        Sign::Zero => unreachable!("coplanar triangles have parallel normals"),
    };

    let axis = dominant_axis(&n1);
    let mut clip: Vec<R2> = r1.iter().map(|p| p.project_drop(axis)).collect();
    // Normalize the clip triangle to CCW in projection space.
    if orient2d_r(&clip[0], &clip[1], &clip[2]) == Sign::Neg {
        clip.swap(1, 2);
    }
    let mut subject: Vec<R2> = r2.iter().map(|p| p.project_drop(axis)).collect();
    if orient2d_r(&subject[0], &subject[1], &subject[2]) == Sign::Neg {
        subject.swap(1, 2);
    }

    // Clip `subject` against each closed halfplane left of the CCW clip edges.
    let mut poly = subject;
    for i in 0..3 {
        if poly.is_empty() {
            break;
        }
        let c0 = clip[i].clone();
        let c1 = clip[(i + 1) % 3].clone();
        let mut out: Vec<R2> = Vec::with_capacity(poly.len() + 2);
        for k in 0..poly.len() {
            let s = &poly[k];
            let e = &poly[(k + 1) % poly.len()];
            let s_in = orient2d_r(&c0, &c1, s) != Sign::Neg;
            let e_in = orient2d_r(&c0, &c1, e) != Sign::Neg;
            match (s_in, e_in) {
                (true, true) => out.push(e.clone()),
                (true, false) => {
                    let x = line_line_intersect_2d(&c0, &c1, s, e)
                        .expect("strictly crossing edge is not parallel to clip line");
                    out.push(x);
                }
                (false, true) => {
                    let x = line_line_intersect_2d(&c0, &c1, s, e)
                        .expect("strictly crossing edge is not parallel to clip line");
                    out.push(x);
                    out.push(e.clone());
                }
                (false, false) => {}
            }
        }
        poly = out;
    }

    // Canonicalize: drop consecutive duplicates (exact equality) and
    // collinear intermediate vertices.
    let poly = canonical_polygon(poly);
    match poly.len() {
        0 => TriTriIsect::None,
        1 => TriTriIsect::Point(lift_to_plane(&poly[0], axis, &r1[0], &n1)),
        2 => TriTriIsect::Segment(
            lift_to_plane(&poly[0], axis, &r1[0], &n1),
            lift_to_plane(&poly[1], axis, &r1[0], &n1),
        ),
        _ => TriTriIsect::Coplanar {
            polygon: poly
                .iter()
                .map(|p| lift_to_plane(p, axis, &r1[0], &n1))
                .collect(),
            same_orientation,
        },
    }
}

/// Remove exact duplicates and collinear intermediate vertices from a closed
/// polygon; a fully collinear result collapses to its two extreme points, a
/// single repeated point to one point.
fn canonical_polygon(poly: Vec<R2>) -> Vec<R2> {
    // Dedup (cyclic).
    let mut pts: Vec<R2> = Vec::with_capacity(poly.len());
    for p in poly {
        if pts.last() != Some(&p) {
            pts.push(p);
        }
    }
    while pts.len() > 1 && pts.first() == pts.last() {
        pts.pop();
    }
    if pts.len() <= 2 {
        return pts;
    }
    // Fully collinear (possible when the overlap is a shared edge segment
    // that SH clipping walked over several vertices): keep the two extremes.
    let all_collinear = (0..pts.len()).all(|i| {
        let a = &pts[i];
        let b = &pts[(i + 1) % pts.len()];
        let c = &pts[(i + 2) % pts.len()];
        orient2d_r(a, b, c) == Sign::Zero
    });
    if all_collinear {
        // Order along the dominant direction of the point spread.
        let dir = pts
            .iter()
            .skip(1)
            .map(|p| p.sub(&pts[0]))
            .find(|d| !d.is_zero())
            .expect("at least two distinct points");
        let param = |p: &R2| p.sub(&pts[0]).dot(&dir);
        let mut lo = pts[0].clone();
        let mut hi = pts[0].clone();
        let (mut lo_t, mut hi_t) = (BigRational::zero(), BigRational::zero());
        for p in &pts {
            let t = param(p);
            if t < lo_t {
                lo_t = t.clone();
                lo = p.clone();
            }
            if t > hi_t {
                hi_t = t;
                hi = p.clone();
            }
        }
        if lo == hi {
            return vec![lo];
        }
        return vec![lo, hi];
    }
    // Drop collinear intermediate vertices.
    let n = pts.len();
    let keep: Vec<R2> = (0..n)
        .filter(|&i| {
            let prev = &pts[(i + n - 1) % n];
            let next = &pts[(i + 1) % n];
            orient2d_r(prev, &pts[i], next) != Sign::Zero
        })
        .map(|i| pts[i].clone())
        .collect();
    keep
}

#[cfg(test)]
#[path = "tri_tri_tests.rs"]
mod tests;
