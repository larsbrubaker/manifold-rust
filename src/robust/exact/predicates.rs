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

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{Signed, Zero};

use super::rational::{R2, R3};
use super::Sign;

// ─── Exact predicates ────────────────────────────────────────────────────────
//
// Predicate signs are computed in pure BigInt arithmetic: each point is
// homogenized once ((x, y) = (X/W, Y/W) with W > 0 — num-rational keeps
// denominators positive), and the determinant is scaled through by positive
// denominator products, which preserves its sign. This avoids BigRational's
// gcd normalization on every intermediate operation — the dominant cost of
// the original rational formulation (the CDT's incircle calls on constructed
// intersection points made it ~80% of robust-boolean wall time).

/// (X, Y, W): x = X/W, y = Y/W with W > 0.
#[inline]
fn homog2(p: &R2) -> (BigInt, BigInt, BigInt) {
    let (xn, xd) = (p.x.numer(), p.x.denom());
    let (yn, yd) = (p.y.numer(), p.y.denom());
    (xn * yd, yn * xd, xd * yd)
}

/// (X, Y, Z, W): coordinates over one positive common denominator.
#[inline]
fn homog3(p: &R3) -> (BigInt, BigInt, BigInt, BigInt) {
    let (xn, xd) = (p.x.numer(), p.x.denom());
    let (yn, yd) = (p.y.numer(), p.y.denom());
    let (zn, zd) = (p.z.numer(), p.z.denom());
    let yz = yd * zd;
    (xn * &yz, yn * (xd * zd), zn * (xd * yd), xd * yz)
}

#[inline]
fn sign_of_int(v: &BigInt) -> Sign {
    if v.is_positive() {
        Sign::Pos
    } else if v.is_negative() {
        Sign::Neg
    } else {
        Sign::Zero
    }
}

/// Sign of cross(b-a, c-a). Pos ⇔ a,b,c wind counterclockwise.
pub fn orient2d_r(a: &R2, b: &R2, c: &R2) -> Sign {
    let (ax, ay, aw) = homog2(a);
    let (bx, by, bw) = homog2(b);
    let (cx, cy, cw) = homog2(c);
    // det(b-a, c-a) · Wa²WbWc = (BxWa−AxWb)(CyWa−AyWc) − (ByWa−AyWb)(CxWa−AxWc)
    let ux = &bx * &aw - &ax * &bw;
    let uy = &by * &aw - &ay * &bw;
    let vx = &cx * &aw - &ax * &cw;
    let vy = &cy * &aw - &ay * &cw;
    sign_of_int(&(ux * vy - uy * vx))
}

/// Sign of dot(cross(b-a, c-a), d-a). Pos ⇔ d lies on the CCW-normal side
/// of the plane through a, b, c; Zero ⇔ the four points are coplanar.
pub fn orient3d_r(a: &R3, b: &R3, c: &R3, d: &R3) -> Sign {
    let (ax, ay, az, aw) = homog3(a);
    let (bx, by, bz, bw) = homog3(b);
    let (cx, cy, cz, cw) = homog3(c);
    let (dx, dy, dz, dw) = homog3(d);
    // Each difference row is scaled by the positive factor Wa·W_row; the
    // triple product then carries a positive overall scale.
    let ux = &bx * &aw - &ax * &bw;
    let uy = &by * &aw - &ay * &bw;
    let uz = &bz * &aw - &az * &bw;
    let vx = &cx * &aw - &ax * &cw;
    let vy = &cy * &aw - &ay * &cw;
    let vz = &cz * &aw - &az * &cw;
    let wx = &dx * &aw - &ax * &dw;
    let wy = &dy * &aw - &ay * &dw;
    let wz = &dz * &aw - &az * &dw;
    let det = (&uy * &vz - &uz * &vy) * wx
        + (&uz * &vx - &ux * &vz) * wy
        + (&ux * &vy - &uy * &vx) * wz;
    sign_of_int(&det)
}

/// Incircle test. With a,b,c counterclockwise: Pos ⇔ d strictly inside the
/// circumcircle of (a,b,c). Computed as the standard 3×3 determinant of
/// coordinates lifted onto the paraboloid, with rows differenced against d.
pub fn incircle_r(a: &R2, b: &R2, c: &R2, d: &R2) -> Sign {
    let (dx, dy, dw) = homog2(d);
    // Row i (i = a,b,c): (xi−xd, yi−yd) over the positive denominator WiWd.
    // Scaling row i by (WiWd)² keeps the lift column polynomial:
    //   Ui = (XiWd−XdWi)·WiWd,  Vi = (YiWd−YdWi)·WiWd,
    //   Li = (XiWd−XdWi)² + (YiWd−YdWi)².
    let row = |p: &R2| -> (BigInt, BigInt, BigInt) {
        let (px, py, pw) = homog2(p);
        let nx = &px * &dw - &dx * &pw;
        let ny = &py * &dw - &dy * &pw;
        let s = pw * &dw;
        let lift = &nx * &nx + &ny * &ny;
        (nx * &s, ny * &s, lift)
    };
    let (ux, uy, ul) = row(a);
    let (vx, vy, vl) = row(b);
    let (wx, wy, wl) = row(c);
    let det = ul * (&vx * &wy - &vy * &wx)
        + vl * (&wx * &uy - &wy * &ux)
        + wl * (&ux * &vy - &uy * &vx);
    sign_of_int(&det)
}

/// Exact: p collinear with (a,b) and within the closed segment [a,b].
/// Same integer-only strategy as the sign predicates above — the registry
/// sweeps in robust/intersection_graph.rs call this in a tight loop.
pub fn point_on_segment_r(p: &R3, a: &R3, b: &R3) -> bool {
    let (px, py, pz, pw) = homog3(p);
    let (ax, ay, az, aw) = homog3(a);
    let (bx, by, bz, bw) = homog3(b);
    // ap scaled by the positive PwAw; d = b−a scaled by the positive AwBw.
    let apx = &px * &aw - &ax * &pw;
    let apy = &py * &aw - &ay * &pw;
    let apz = &pz * &aw - &az * &pw;
    let dx = &bx * &aw - &ax * &bw;
    let dy = &by * &aw - &ay * &bw;
    let dz = &bz * &aw - &az * &bw;
    // Collinearity: cross(ap, d) = 0 (positive scales cannot zero a component).
    if !(&apy * &dz - &apz * &dy).is_zero()
        || !(&apz * &dx - &apx * &dz).is_zero()
        || !(&apx * &dy - &apy * &dx).is_zero()
    {
        return false;
    }
    // 0 ≤ ap·d and ap·d ≤ d·d, cleared of their (positive) denominators:
    //   ap·d = S1 / (PwAw·AwBw),  d·d = S2 / (AwBw)²
    //   S1 ≥ 0   and   S1·AwBw ≤ S2·PwAw.
    let s1 = &apx * &dx + &apy * &dy + &apz * &dz;
    if s1.is_negative() {
        return false;
    }
    let s2 = &dx * &dx + &dy * &dy + &dz * &dz;
    s1 * (&aw * &bw) <= s2 * (pw * aw)
}

/// (a−o)·u as an unreduced fraction (numerator, positive denominator) —
/// integer-only, no gcd normalization. For sign tests and cross-multiplied
/// comparisons (e.g. the radial ring sort in robust/classify.rs) the
/// unreduced form is exactly as good as the canonical one and much cheaper
/// to produce.
pub fn dot_diff_raw(a: &R3, o: &R3, u: &R3) -> (BigInt, BigInt) {
    let (ax, ay, az, aw) = homog3(a);
    let (ox, oy, oz, ow) = homog3(o);
    let (ux, uy, uz, uw) = homog3(u);
    // (a−o) scaled by the positive AwOw; dot with u adds a Uw denominator.
    let num = (&ax * &ow - &ox * &aw) * ux
        + (&ay * &ow - &oy * &aw) * uy
        + (&az * &ow - &oz * &aw) * uz;
    (num, aw * ow * uw)
}

/// Triangle normal as a denominator-cleared integer vector: cross(b−a, c−a)
/// scaled by the positive Aw²BwCw. Direction (and zero-ness) match
/// `tri_normal_r`; use where only the normal's direction matters.
pub fn tri_normal_int(a: &R3, b: &R3, c: &R3) -> [BigInt; 3] {
    let (ax, ay, az, aw) = homog3(a);
    let (bx, by, bz, bw) = homog3(b);
    let (cx, cy, cz, cw) = homog3(c);
    let ux = &bx * &aw - &ax * &bw;
    let uy = &by * &aw - &ay * &bw;
    let uz = &bz * &aw - &az * &bw;
    let vx = &cx * &aw - &ax * &cw;
    let vy = &cy * &aw - &ay * &cw;
    let vz = &cz * &aw - &az * &cw;
    [
        &uy * &vz - &uz * &vy,
        &uz * &vx - &ux * &vz,
        &ux * &vy - &uy * &vx,
    ]
}

/// d·p as an unreduced fraction (numerator, positive denominator), for an
/// integer direction `d` and rational point `p`. Comparable across points by
/// cross-multiplication — the segment-interval overlap in robust/tri_tri.rs
/// orders plane-crossing points along the intersection line with this.
pub fn dot_point_raw(d: &[BigInt; 3], p: &R3) -> (BigInt, BigInt) {
    let (px, py, pz, pw) = homog3(p);
    (&d[0] * &px + &d[1] * &py + &d[2] * &pz, pw)
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
    // Integer-only formulation (same exact point as the rational one): the
    // plane normal enters both numerator and denominator of t, so any
    // positive common scale on it cancels — compute it as an integer cross
    // of denominator-cleared edge vectors and never normalize intermediates.
    // Each output coordinate reduces exactly once in BigRational::new.
    let (px, py, pz, pw) = homog3(p);
    let (qx, qy, qz, qw) = homog3(q);
    let (ax, ay, az, aw) = homog3(a);
    let (bx, by, bz, bw) = homog3(b);
    let (cx, cy, cz, cw) = homog3(c);

    // (b−a)·AwBw and (c−a)·AwCw; their cross is the normal up to Aw²BwCw > 0.
    let ux = &bx * &aw - &ax * &bw;
    let uy = &by * &aw - &ay * &bw;
    let uz = &bz * &aw - &az * &bw;
    let vx = &cx * &aw - &ax * &cw;
    let vy = &cy * &aw - &ay * &cw;
    let vz = &cz * &aw - &az * &cw;
    let nx = &uy * &vz - &uz * &vy;
    let ny = &uz * &vx - &ux * &vz;
    let nz = &ux * &vy - &uy * &vx;

    // dir = q−p scaled by PwQw; e = a−p scaled by PwAw.
    let dx = &qx * &pw - &px * &qw;
    let dy = &qy * &pw - &py * &qw;
    let dz = &qz * &pw - &pz * &qw;
    let n_dot_d = &nx * &dx + &ny * &dy + &nz * &dz;
    if n_dot_d.is_zero() {
        return None;
    }
    let ex = &ax * &pw - &px * &aw;
    let ey = &ay * &pw - &py * &aw;
    let ez = &az * &pw - &pz * &aw;
    let n_dot_e = &nx * &ex + &ny * &ey + &nz * &ez;

    // t = (n·e)·Qw / (Aw·(n·d));  x_i = (P_i·Qw·T_d + D_i·T_n) / (Pw·Qw·T_d)
    let t_n = &n_dot_e * &qw;
    let t_d = &aw * &n_dot_d;
    let den = &pw * &qw * &t_d;
    let coord = |pi: &BigInt, di: &BigInt| -> BigRational {
        BigRational::new(pi * &qw * &t_d + di * &t_n, den.clone())
    };
    Some(R3::new(coord(&px, &dx), coord(&py, &dy), coord(&pz, &dz)))
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
