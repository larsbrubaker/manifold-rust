// robust/exact/approx.rs — Semi-static float filters over approximate
// coordinates of exact points ("indirect predicates" in the sense of
// Attene 2020, the engine room of interactive exact mesh booleans).
//
// Every point the robust pipeline touches is either an input point (exact
// f64) or a constructed rational whose f64 approximation we obtain by one
// correct rounding (robust/exact/rational.rs::rat_to_f64) — so every
// approximate coordinate x̃ satisfies x̃ = x·(1+δ) with |δ| ≤ ε = 2⁻⁵³.
// A predicate can therefore run entirely in f64 on the approximations and
// certify its sign with a magnitude-based error bound; only near-degenerate
// configurations escalate to the exact BigInt evaluations in predicates.rs.
//
// The bounds here are deliberately conservative (they cover both the
// float-evaluation roundoff and the input perturbation, with slack): an
// over-tight constant would be a soundness bug, an over-loose one only costs
// extra escalations exactly where the exact path was needed anyway. Unlike
// the classic Shewchuk filters in filtered.rs, the "permanents" below use
// coordinate magnitudes rather than formed differences — input perturbation
// error scales with |a|+|c|, not |a−c|, under cancellation.

use super::Sign;

const EPS: f64 = f64::EPSILON * 0.5; // 2⁻⁵³

/// Conservative certified-sign helper: `det` was computed in f64 from
/// ε-perturbed inputs; `bound` ≥ the true error. Returns None when the sign
/// cannot be certified (including any non-finite intermediate).
#[inline]
fn certify(det: f64, bound: f64) -> Option<Sign> {
    if !det.is_finite() || !bound.is_finite() {
        return None;
    }
    if det > bound {
        Some(Sign::Pos)
    } else if det < -bound {
        Some(Sign::Neg)
    } else {
        None
    }
}

/// Filtered orient2d over approximate coordinates: sign of
/// cross(b−a, c−a), or None when uncertain.
///
/// Error analysis (conservative): with per-coordinate relative error ≤ ε and
/// ≤ 3 roundings per product term after differencing, the total error is
/// < 8ε·P where P = (|ax|+|bx|+|cx|)·(|ay|+|by|+|cy|) bounds every term's
/// magnitude sum. We use 16ε·P for slack.
#[inline]
pub fn orient2d_a(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> Option<Sign> {
    let det = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
    let px = a[0].abs() + b[0].abs() + c[0].abs();
    let py = a[1].abs() + b[1].abs() + c[1].abs();
    certify(det, 16.0 * EPS * px * py)
}

/// Filtered incircle over approximate coordinates (a,b,c CCW ⇒ Pos = d
/// strictly inside), or None when uncertain. Conservative bound 64ε·P with P
/// the product of per-axis magnitude sums times the lift magnitude sum.
#[inline]
pub fn incircle_a(a: [f64; 2], b: [f64; 2], c: [f64; 2], d: [f64; 2]) -> Option<Sign> {
    let adx = a[0] - d[0];
    let ady = a[1] - d[1];
    let bdx = b[0] - d[0];
    let bdy = b[1] - d[1];
    let cdx = c[0] - d[0];
    let cdy = c[1] - d[1];
    let alift = adx * adx + ady * ady;
    let blift = bdx * bdx + bdy * bdy;
    let clift = cdx * cdx + cdy * cdy;
    let det = alift * (bdx * cdy - cdx * bdy)
        + blift * (cdx * ady - adx * cdy)
        + clift * (adx * bdy - bdx * ady);

    // Magnitude-based permanent: differences bounded by |p|+|d| sums.
    let mx = a[0].abs().max(b[0].abs()).max(c[0].abs()) + d[0].abs();
    let my = a[1].abs().max(b[1].abs()).max(c[1].abs()) + d[1].abs();
    let lift = mx * mx + my * my;
    certify(det, 64.0 * EPS * lift * mx * my)
}

/// Filtered orient3d over approximate coordinates: sign of
/// dot(cross(b−a, c−a), d−a), or None when uncertain. Conservative 32ε·P.
#[inline]
pub fn orient3d_a(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> Option<Sign> {
    let ux = b[0] - a[0];
    let uy = b[1] - a[1];
    let uz = b[2] - a[2];
    let vx = c[0] - a[0];
    let vy = c[1] - a[1];
    let vz = c[2] - a[2];
    let wx = d[0] - a[0];
    let wy = d[1] - a[1];
    let wz = d[2] - a[2];
    let det = (uy * vz - uz * vy) * wx + (uz * vx - ux * vz) * wy + (ux * vy - uy * vx) * wz;
    let px = a[0].abs() + b[0].abs() + c[0].abs() + d[0].abs();
    let py = a[1].abs() + b[1].abs() + c[1].abs() + d[1].abs();
    let pz = a[2].abs() + b[2].abs() + c[2].abs() + d[2].abs();
    // Every det term is a product of one coordinate-difference per axis.
    let p = py * pz * px;
    certify(det, 32.0 * EPS * p)
}

/// Filtered collinear-and-within-segment prefilter for
/// `predicates::point_on_segment_r`: `Some(false)` when p is certifiably NOT
/// on segment [a,b] (the overwhelmingly common case in the registry sweeps),
/// `None` when the exact test must decide. Never returns `Some(true)` — a
/// certified "on" would need exact zero detection.
#[inline]
pub fn not_on_segment_a(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> Option<bool> {
    let apx = p[0] - a[0];
    let apy = p[1] - a[1];
    let apz = p[2] - a[2];
    let dx = b[0] - a[0];
    let dy = b[1] - a[1];
    let dz = b[2] - a[2];
    // Cross-product components with magnitude-permanent bounds per component.
    let px = p[0].abs() + a[0].abs() + b[0].abs();
    let py = p[1].abs() + a[1].abs() + b[1].abs();
    let pz = p[2].abs() + a[2].abs() + b[2].abs();
    let cx = apy * dz - apz * dy;
    let cy = apz * dx - apx * dz;
    let cz = apx * dy - apy * dx;
    if cx.abs() > 16.0 * EPS * py * pz
        || cy.abs() > 16.0 * EPS * pz * px
        || cz.abs() > 16.0 * EPS * px * py
    {
        return Some(false); // certifiably off the line
    }
    // Certifiably outside the parameter range?
    let s1 = apx * dx + apy * dy + apz * dz;
    let s2 = dx * dx + dy * dy + dz * dz;
    let pd = px * px + py * py + pz * pz; // ≥ every dot magnitude here
    let dot_err = 16.0 * EPS * pd;
    if s1 < -dot_err || s1 > s2 + 2.0 * dot_err {
        return Some(false);
    }
    None
}

#[cfg(test)]
mod tests {
    use super::super::predicates::{incircle_r, orient2d_r, orient3d_r, point_on_segment_r};
    use super::super::rational::{rat_to_f64, R2, R3};
    use super::*;
    use crate::linalg::{Vec2, Vec3};
    use num_rational::BigRational;

    struct Lcg(u64);
    impl Lcg {
        fn next_f64(&mut self, scale: f64) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((self.0 >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0) * scale
        }
    }

    fn approx2(p: &R2) -> [f64; 2] {
        [rat_to_f64(&p.x), rat_to_f64(&p.y)]
    }
    fn approx3(p: &R3) -> [f64; 3] {
        [rat_to_f64(&p.x), rat_to_f64(&p.y), rat_to_f64(&p.z)]
    }

    /// A rational point near (x,y) with a huge denominator — its rounding is
    /// a genuine ε-perturbation, exercising the indirect-point case.
    fn wobble2(x: f64, y: f64, k: i64) -> R2 {
        let tiny = BigRational::new(k.into(), num_bigint::BigInt::from(3u8).pow(40));
        R2::new(
            BigRational::from_float(x).unwrap() + &tiny,
            BigRational::from_float(y).unwrap() - &tiny,
        )
    }

    #[test]
    fn certified_signs_agree_with_exact() {
        let mut rng = Lcg(0xFEEDFACE);
        let mut certified = 0usize;
        for i in 0..4000 {
            let pts: Vec<R2> = (0..4)
                .map(|j| wobble2(rng.next_f64(50.0), rng.next_f64(50.0), (i * 4 + j) as i64))
                .collect();
            let ap: Vec<[f64; 2]> = pts.iter().map(approx2).collect();
            if let Some(s) = orient2d_a(ap[0], ap[1], ap[2]) {
                assert_eq!(s, orient2d_r(&pts[0], &pts[1], &pts[2]), "orient2d #{i}");
                certified += 1;
            }
            if let Some(s) = incircle_a(ap[0], ap[1], ap[2], ap[3]) {
                assert_eq!(
                    s,
                    incircle_r(&pts[0], &pts[1], &pts[2], &pts[3]),
                    "incircle #{i}"
                );
            }
        }
        assert!(
            certified > 3900,
            "filter should certify generic input ({certified}/4000)"
        );
    }

    #[test]
    fn orient3d_filter_agrees_with_exact() {
        let mut rng = Lcg(0xDEADBEA7);
        for i in 0..3000 {
            let p: Vec<R3> = (0..4)
                .map(|_| {
                    R3::from_vec3(Vec3::new(
                        rng.next_f64(20.0),
                        rng.next_f64(20.0),
                        rng.next_f64(20.0),
                    ))
                })
                .collect();
            let ap: Vec<[f64; 3]> = p.iter().map(approx3).collect();
            if let Some(s) = orient3d_a(ap[0], ap[1], ap[2], ap[3]) {
                assert_eq!(s, orient3d_r(&p[0], &p[1], &p[2], &p[3]), "orient3d #{i}");
            }
        }
    }

    #[test]
    fn near_degenerate_defers_to_exact() {
        // Points a hair off a line: the filter must return None (uncertain),
        // never a wrong certified sign.
        let a = Vec2::new(12.0, 12.0);
        let b = Vec2::new(24.0, 24.0);
        for i in -8i64..=8 {
            let c = wobble2(48.0, 48.0, i);
            let ra = R2::from_vec2(a);
            let rb = R2::from_vec2(b);
            let ac = approx2(&c);
            match orient2d_a([a.x, a.y], [b.x, b.y], ac) {
                Some(s) => assert_eq!(s, orient2d_r(&ra, &rb, &c), "i={i}"),
                None => {} // escalation is the correct answer here
            }
        }
    }

    #[test]
    fn not_on_segment_prefilter_is_sound() {
        let mut rng = Lcg(0xBAD5EED);
        let a = R3::from_vec3(Vec3::new(0.0, 0.0, 0.0));
        let b = R3::from_vec3(Vec3::new(10.0, 4.0, 2.0));
        let (aa, ab) = (approx3(&a), approx3(&b));
        let mut rejected = 0usize;
        for i in 0..4000 {
            let p = R3::new(
                BigRational::from_float(rng.next_f64(12.0)).unwrap(),
                BigRational::from_float(rng.next_f64(12.0)).unwrap(),
                BigRational::from_float(rng.next_f64(12.0)).unwrap(),
            );
            let exact = point_on_segment_r(&p, &a, &b);
            match not_on_segment_a(approx3(&p), aa, ab) {
                Some(false) => {
                    assert!(!exact, "prefilter wrongly rejected an on-segment point #{i}");
                    rejected += 1;
                }
                Some(true) => unreachable!(),
                None => {}
            }
        }
        assert!(rejected > 3900, "prefilter should reject generic points ({rejected}/4000)");
        // And points genuinely on the segment always defer or agree.
        for k in 1..20 {
            let t = BigRational::new(k.into(), 21.into());
            let on = a.add(&b.sub(&a).scale(&t));
            assert!(point_on_segment_r(&on, &a, &b));
            assert_ne!(
                not_on_segment_a(approx3(&on), aa, ab),
                Some(false),
                "k={k}"
            );
        }
    }
}
