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
use super::exact::predicates::{dot_point_raw, line_line_intersect_2d, line_plane_intersect, orient2d_r, tri_normal_int, tri_normal_r};
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

/// Exact intersection of triangles t1 and t2 (each three finite f64
/// vertices). Symmetric: swapping the arguments yields the same set.
pub fn tri_tri_intersect(t1: [Vec3; 3], t2: [Vec3; 3]) -> TriTriIsect {
    // Signs of t2's vertices against t1's plane.
    let s2 = [
        orient3d(t1[0], t1[1], t1[2], t2[0]),
        orient3d(t1[0], t1[1], t1[2], t2[1]),
        orient3d(t1[0], t1[1], t1[2], t2[2]),
    ];
    if all_same_strict(&s2) {
        return TriTriIsect::None;
    }
    if s2.iter().all(|s| *s == Sign::Zero) {
        return coplanar_overlap(t1, t2);
    }
    // Signs of t1's vertices against t2's plane.
    let s1 = [
        orient3d(t2[0], t2[1], t2[2], t1[0]),
        orient3d(t2[0], t2[1], t2[2], t1[1]),
        orient3d(t2[0], t2[1], t2[2], t1[2]),
    ];
    if all_same_strict(&s1) {
        return TriTriIsect::None;
    }
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

    // Both triangles meet the common line L of the two planes. Collect each
    // triangle's (1- or 2-point) intersection with the other's plane and
    // overlap the two intervals along L.
    let pts1 = plane_crossings(&r1, &s1, &r2);
    let pts2 = plane_crossings(&r2, &s2, &r1);
    debug_assert!(!pts1.is_empty() && pts1.len() <= 2);
    debug_assert!(!pts2.is_empty() && pts2.len() <= 2);

    // Integer line direction (positive scale of n1×n2) and unreduced interval
    // parameters — ordering along L only needs comparisons, not canonical
    // rationals, so no gcds anywhere in the interval overlap.
    let n1 = tri_normal_int(&r1[0], &r1[1], &r1[2]);
    let n2 = tri_normal_int(&r2[0], &r2[1], &r2[2]);
    let dir = [
        &n1[1] * &n2[2] - &n1[2] * &n2[1],
        &n1[2] * &n2[0] - &n1[0] * &n2[2],
        &n1[0] * &n2[1] - &n1[1] * &n2[0],
    ];
    debug_assert!(
        dir.iter().any(|c| !c.is_zero()),
        "non-coplanar intersecting planes"
    );

    let i1 = interval_along(&dir, pts1);
    let i2 = interval_along(&dir, pts2);
    let (lo, lo_pt) = if cmp_frac(&i1.0 .0, &i2.0 .0) != std::cmp::Ordering::Less { i1.0 } else { i2.0 };
    let (hi, hi_pt) = if cmp_frac(&i1.1 .0, &i2.1 .0) != std::cmp::Ordering::Greater { i1.1 } else { i2.1 };
    match cmp_frac(&lo, &hi) {
        std::cmp::Ordering::Greater => TriTriIsect::None,
        std::cmp::Ordering::Equal => TriTriIsect::Point(lo_pt),
        std::cmp::Ordering::Less => TriTriIsect::Segment(lo_pt, hi_pt),
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

/// The 1 or 2 points where triangle `tri` (vertex signs `s` against the
/// other plane) meets the plane of `other`. Vertices exactly on the plane
/// contribute themselves; strictly straddling edges contribute constructed
/// intersection points.
fn plane_crossings(tri: &[R3; 3], s: &[Sign; 3], other: &[R3; 3]) -> Vec<R3> {
    let mut pts: Vec<R3> = Vec::with_capacity(2);
    for i in 0..3 {
        if s[i] == Sign::Zero {
            pts.push(tri[i].clone());
        }
    }
    for i in 0..3 {
        let j = (i + 1) % 3;
        if s[i] != Sign::Zero && s[j] != Sign::Zero && s[i] != s[j] {
            let x = line_plane_intersect(&tri[i], &tri[j], &other[0], &other[1], &other[2])
                .expect("strictly straddling edge cannot be parallel to the plane");
            pts.push(x);
        }
    }
    pts
}

/// Interval of points along integer direction `dir` as ((min_param,
/// min_point), (max_param, max_point)), parameters as unreduced fractions.
/// A single input point yields a degenerate interval.
#[allow(clippy::type_complexity)]
fn interval_along(
    dir: &[num_bigint::BigInt; 3],
    mut pts: Vec<R3>,
) -> ((Frac, R3), (Frac, R3)) {
    let first = pts.remove(0);
    let p0 = dot_point_raw(dir, &first);
    let mut lo = ((p0.0.clone(), p0.1.clone()), first.clone());
    let mut hi = (p0, first);
    for p in pts {
        let t = dot_point_raw(dir, &p);
        if cmp_frac(&t, &lo.0) == std::cmp::Ordering::Less {
            lo = ((t.0.clone(), t.1.clone()), p.clone());
        }
        if cmp_frac(&t, &hi.0) == std::cmp::Ordering::Greater {
            hi = (t, p);
        }
    }
    (lo, hi)
}

// ─── Coplanar overlap ────────────────────────────────────────────────────────

/// Intersection of two coplanar triangles: Sutherland–Hodgman clip of t2
/// against t1 in the exact 2D projection, classified by the dimension of the
/// result (empty / point / segment / convex polygon).
fn coplanar_overlap(t1: [Vec3; 3], t2: [Vec3; 3]) -> TriTriIsect {
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
