// robust/exact/predicates.rs — Fully exact predicates and geometric
// constructions on rational points. Ground truth for the robust boolean
// engine: robust/exact/filtered.rs escalates here whenever its float
// filters cannot certify a sign, and the intersection code builds all new
// vertices through the constructions at the bottom of this file.
//
// Orientation conventions (documented once, used everywhere in src/robust):
//   orient2d(a,b,c)   = sign of cross(b-a, c-a); Pos = a,b,c counterclockwise.
//   orient3d(a,b,c,d) = sign of dot(cross(b-a, c-a), d-a); Pos = d on the
//                       side of plane(a,b,c) that its CCW normal points to.
//   incircle(a,b,c,d) with a,b,c CCW: Pos = d strictly inside the circle
//                       through a,b,c. (For CW a,b,c the sign flips.)

use num_rational::BigRational;
use num_traits::Zero;

use super::rational::{R2, R3};
use super::Sign;

// ─── Exact predicates ────────────────────────────────────────────────────────

/// Sign of cross(b-a, c-a). Pos ⇔ a,b,c wind counterclockwise.
pub fn orient2d_r(a: &R2, b: &R2, c: &R2) -> Sign {
    Sign::of_rat(&b.sub(a).cross(&c.sub(a)))
}

/// Sign of dot(cross(b-a, c-a), d-a). Pos ⇔ d lies on the CCW-normal side
/// of the plane through a, b, c; Zero ⇔ the four points are coplanar.
pub fn orient3d_r(a: &R3, b: &R3, c: &R3, d: &R3) -> Sign {
    let u = b.sub(a);
    let v = c.sub(a);
    let w = d.sub(a);
    Sign::of_rat(&u.cross(&v).dot(&w))
}

/// Incircle test. With a,b,c counterclockwise: Pos ⇔ d strictly inside the
/// circumcircle of (a,b,c). Computed as the standard 3×3 determinant of
/// coordinates lifted onto the paraboloid, with rows differenced against d.
pub fn incircle_r(a: &R2, b: &R2, c: &R2, d: &R2) -> Sign {
    let ad = a.sub(d);
    let bd = b.sub(d);
    let cd = c.sub(d);
    let a_lift = ad.dot(&ad);
    let b_lift = bd.dot(&bd);
    let c_lift = cd.dot(&cd);
    let det = &a_lift * bd.cross(&cd) + &b_lift * cd.cross(&ad) + &c_lift * ad.cross(&bd);
    Sign::of_rat(&det)
}

/// Unnormalized CCW normal of triangle (a,b,c): cross(b-a, c-a). Zero vector
/// ⇔ the triangle is degenerate.
pub fn tri_normal_r(a: &R3, b: &R3, c: &R3) -> R3 {
    b.sub(a).cross(&c.sub(a))
}

/// Where a point lies relative to a non-degenerate triangle, in 2D.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TriLoc {
    Inside,
    /// On the open edge i, where edge i runs from vertex i to vertex i+1 (mod 3).
    OnEdge(u8),
    /// Coincident with vertex i.
    OnVertex(u8),
    Outside,
}

/// Locate point p relative to triangle (a,b,c). Works for either winding;
/// a degenerate (zero-area) triangle reports every point as Outside, which
/// is consistent with the pipeline dropping degenerate triangles up front.
pub fn point_in_tri_2d(p: &R2, a: &R2, b: &R2, c: &R2) -> TriLoc {
    let orient = orient2d_r(a, b, c);
    if orient == Sign::Zero {
        return TriLoc::Outside;
    }
    // Normalize so the triangle reads as CCW.
    let normalize = |s: Sign| if orient == Sign::Pos { s } else { s.flip() };
    let s0 = normalize(orient2d_r(a, b, p)); // edge 0: a→b
    let s1 = normalize(orient2d_r(b, c, p)); // edge 1: b→c
    let s2 = normalize(orient2d_r(c, a, p)); // edge 2: c→a
    if s0 == Sign::Neg || s1 == Sign::Neg || s2 == Sign::Neg {
        return TriLoc::Outside;
    }
    match (s0 == Sign::Zero, s1 == Sign::Zero, s2 == Sign::Zero) {
        (false, false, false) => TriLoc::Inside,
        (true, false, false) => TriLoc::OnEdge(0),
        (false, true, false) => TriLoc::OnEdge(1),
        (false, false, true) => TriLoc::OnEdge(2),
        (true, false, true) => TriLoc::OnVertex(0),  // a: on edges c→a and a→b
        (true, true, false) => TriLoc::OnVertex(1),  // b
        (false, true, true) => TriLoc::OnVertex(2),  // c
        (true, true, true) => TriLoc::Outside,       // impossible for orient != 0
    }
}

// ─── Exact constructions ─────────────────────────────────────────────────────

/// Intersection of the line through p,q with the plane of triangle (a,b,c).
///
/// Returns None when the segment direction is parallel to the plane (no
/// unique intersection point — the coplanar-overlap machinery handles that
/// case separately). The caller decides whether the parameter t lies inside
/// [0,1]; use the orient3d signs of p and q for that, not float comparisons.
pub fn line_plane_intersect(p: &R3, q: &R3, a: &R3, b: &R3, c: &R3) -> Option<R3> {
    let n = tri_normal_r(a, b, c);
    let dir = q.sub(p);
    let denom = n.dot(&dir);
    if denom.is_zero() {
        return None;
    }
    let t = n.dot(&a.sub(p)) / denom;
    Some(p.add(&dir.scale(&t)))
}

/// Intersection point of the 2D lines through (a,b) and (c,d). None when the
/// lines are parallel (including collinear — overlap is handled by the
/// caller's collinear branch).
pub fn line_line_intersect_2d(a: &R2, b: &R2, c: &R2, d: &R2) -> Option<R2> {
    let ab = b.sub(a);
    let cd = d.sub(c);
    let denom = ab.cross(&cd);
    if denom.is_zero() {
        return None;
    }
    let t = c.sub(a).cross(&cd) / denom;
    Some(a.add(&ab.scale(&t)))
}

/// Parameter of point x on segment (p,q) along the dominant axis of the
/// segment direction — exact, in [0,1] iff x lies within the segment. The
/// caller guarantees x is on the line through p and q and p != q.
pub fn segment_param(p: &R3, q: &R3, x: &R3) -> BigRational {
    let d = q.sub(p);
    let (num, den) = if !d.x.is_zero() {
        (&x.x - &p.x, d.x)
    } else if !d.y.is_zero() {
        (&x.y - &p.y, d.y)
    } else {
        (&x.z - &p.z, d.z)
    };
    debug_assert!(!den.is_zero(), "segment_param requires p != q");
    num / den
}
