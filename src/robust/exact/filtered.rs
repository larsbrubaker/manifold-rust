// robust/exact/filtered.rs — Float entry points for the exact predicates.
//
// Each predicate first evaluates the determinant in plain f64 together with
// a "permanent" (the same expression with every subtraction replaced by
// addition of absolute values). Shewchuk's static error-bound analysis
// ("Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric
// Predicates", 1997) shows the f64 sign is certain whenever
// |det| > errboundA * permanent; otherwise we escalate to the BigRational
// evaluation in robust/exact/predicates.rs. Unlike Shewchuk we skip the
// adaptive intermediate stages — the exact fallback fires rarely enough
// (see the filter-hit-rate test) that simplicity wins.
//
// The classic bounds assume no underflow/overflow, so any permanent that is
// subnormal, zero, or non-finite also escalates to the exact path.

use crate::linalg::{Vec2, Vec3};

use super::predicates::{incircle_r, orient2d_r, orient3d_r};
use super::rational::{R2, R3};
use super::Sign;

const EPS: f64 = f64::EPSILON * 0.5; // 2^-53, Shewchuk's machine epsilon
const CCW_ERRBOUND_A: f64 = (3.0 + 16.0 * EPS) * EPS;
const O3D_ERRBOUND_A: f64 = (7.0 + 56.0 * EPS) * EPS;
const ICC_ERRBOUND_A: f64 = (10.0 + 96.0 * EPS) * EPS;

/// Counters proving the float filters do their job; only compiled into this
/// crate's own test builds (the hot path stays free of atomics otherwise).
#[cfg(test)]
pub mod stats {
    use std::sync::atomic::{AtomicU64, Ordering};

    pub static FAST: AtomicU64 = AtomicU64::new(0);
    pub static EXACT: AtomicU64 = AtomicU64::new(0);

    pub fn reset() {
        FAST.store(0, Ordering::Relaxed);
        EXACT.store(0, Ordering::Relaxed);
    }

    /// (filter-resolved, exact-fallback) counts since the last reset.
    pub fn snapshot() -> (u64, u64) {
        (FAST.load(Ordering::Relaxed), EXACT.load(Ordering::Relaxed))
    }
}

#[inline]
fn note_fast() {
    #[cfg(test)]
    stats::FAST.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
}

#[inline]
fn note_exact() {
    #[cfg(test)]
    stats::EXACT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
}

/// True when the filter comparison itself is trustworthy: a normal, finite
/// permanent. Subnormal permanents can hide underflowed terms whose error is
/// not covered by the static bound; infinite ones mean the f64 det overflowed.
#[inline]
fn permanent_ok(permanent: f64) -> bool {
    permanent >= f64::MIN_POSITIVE && permanent.is_finite()
}

/// Sign of cross(b-a, c-a); Pos ⇔ a,b,c counterclockwise. Exact.
pub fn orient2d(a: Vec2, b: Vec2, c: Vec2) -> Sign {
    let detleft = (a.x - c.x) * (b.y - c.y);
    let detright = (a.y - c.y) * (b.x - c.x);
    let det = detleft - detright;
    let permanent = detleft.abs() + detright.abs();
    if permanent_ok(permanent) && det.abs() > CCW_ERRBOUND_A * permanent {
        note_fast();
        return Sign::of_f64(det);
    }
    note_exact();
    orient2d_r(&R2::from_vec2(a), &R2::from_vec2(b), &R2::from_vec2(c))
}

/// Sign of dot(cross(b-a, c-a), d-a); Pos ⇔ d on the CCW-normal side of
/// plane(a,b,c), Zero ⇔ coplanar. Exact.
pub fn orient3d(a: Vec3, b: Vec3, c: Vec3, d: Vec3) -> Sign {
    let ux = b.x - a.x;
    let uy = b.y - a.y;
    let uz = b.z - a.z;
    let vx = c.x - a.x;
    let vy = c.y - a.y;
    let vz = c.z - a.z;
    let wx = d.x - a.x;
    let wy = d.y - a.y;
    let wz = d.z - a.z;

    let vywz = vy * wz;
    let vzwy = vz * wy;
    let vzwx = vz * wx;
    let vxwz = vx * wz;
    let vxwy = vx * wy;
    let vywx = vy * wx;

    let det = ux * (vywz - vzwy) + uy * (vzwx - vxwz) + uz * (vxwy - vywx);
    let permanent = (vywz.abs() + vzwy.abs()) * ux.abs()
        + (vzwx.abs() + vxwz.abs()) * uy.abs()
        + (vxwy.abs() + vywx.abs()) * uz.abs();
    if permanent_ok(permanent) && det.abs() > O3D_ERRBOUND_A * permanent {
        note_fast();
        return Sign::of_f64(det);
    }
    note_exact();
    orient3d_r(
        &R3::from_vec3(a),
        &R3::from_vec3(b),
        &R3::from_vec3(c),
        &R3::from_vec3(d),
    )
}

/// Incircle test; with a,b,c CCW: Pos ⇔ d strictly inside the circumcircle
/// of (a,b,c). Exact.
pub fn incircle(a: Vec2, b: Vec2, c: Vec2, d: Vec2) -> Sign {
    let adx = a.x - d.x;
    let ady = a.y - d.y;
    let bdx = b.x - d.x;
    let bdy = b.y - d.y;
    let cdx = c.x - d.x;
    let cdy = c.y - d.y;

    let bdxcdy = bdx * cdy;
    let cdxbdy = cdx * bdy;
    let alift = adx * adx + ady * ady;

    let cdxady = cdx * ady;
    let adxcdy = adx * cdy;
    let blift = bdx * bdx + bdy * bdy;

    let adxbdy = adx * bdy;
    let bdxady = bdx * ady;
    let clift = cdx * cdx + cdy * cdy;

    let det = alift * (bdxcdy - cdxbdy) + blift * (cdxady - adxcdy) + clift * (adxbdy - bdxady);
    let permanent = (bdxcdy.abs() + cdxbdy.abs()) * alift
        + (cdxady.abs() + adxcdy.abs()) * blift
        + (adxbdy.abs() + bdxady.abs()) * clift;
    if permanent_ok(permanent) && det.abs() > ICC_ERRBOUND_A * permanent {
        note_fast();
        return Sign::of_f64(det);
    }
    note_exact();
    incircle_r(
        &R2::from_vec2(a),
        &R2::from_vec2(b),
        &R2::from_vec2(c),
        &R2::from_vec2(d),
    )
}
