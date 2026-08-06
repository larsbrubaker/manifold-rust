// robust/exact/intpred.rs — Exact integer evaluation of the raw-f64
// predicates, the fast escalation tier between the float filters
// (robust/exact/filtered.rs) and the BigRational ground truth
// (robust/exact/predicates.rs).
//
// Every finite f64 is a dyadic rational m·2^e. Scaling all four inputs of a
// predicate *per coordinate axis* by a common power of two turns them into
// integers without changing the determinant's sign (the scale factors are
// positive and factor out of each column). The determinant is then evaluated
// in i128 when a static bit-length budget proves it cannot overflow, and in
// BigInt otherwise — either way the sign is exact, with no rational
// normalization (gcds) anywhere. On the degenerate-heavy meshes that defeat
// the float filters (doubled surfaces, large flat regions) this tier is
// 10–50× cheaper than the BigRational fallback it replaces.

use num_bigint::BigInt;
use num_traits::Zero;

use super::Sign;

/// Exact dyadic decomposition: returns (m, e) with v == m·2^e and trailing
/// zeros stripped from m (minimizing later shift widths). v must be finite.
#[inline]
fn decomp(v: f64) -> (i64, i32) {
    debug_assert!(v.is_finite());
    if v == 0.0 {
        return (0, 0);
    }
    let bits = v.to_bits();
    let biased = ((bits >> 52) & 0x7ff) as i32;
    let frac = bits & ((1u64 << 52) - 1);
    let (mut m, mut e) = if biased == 0 {
        (frac, -1074) // subnormal: no hidden bit
    } else {
        (frac | (1u64 << 52), biased - 1075)
    };
    let tz = m.trailing_zeros();
    m >>= tz;
    e += tz as i32;
    (if v < 0.0 { -(m as i64) } else { m as i64 }, e)
}

/// The four values of one coordinate axis as exactly scaled i128 integers,
/// or None when any scaled magnitude would exceed `budget` bits (the caller
/// then takes the BigInt path). Zeros scale to zero regardless of exponent.
#[inline]
fn scaled_i128<const N: usize>(vs: [f64; N], budget: u32) -> Option<[i128; N]> {
    let d = vs.map(decomp);
    let emin = d
        .iter()
        .filter(|(m, _)| *m != 0)
        .map(|(_, e)| *e)
        .min()
        .unwrap_or(0);
    let mut out = [0i128; N];
    for i in 0..N {
        let (m, e) = d[i];
        if m == 0 {
            continue;
        }
        let shift = (e - emin) as u32;
        let bits = 64 - m.unsigned_abs().leading_zeros() + shift;
        if bits > budget {
            return None;
        }
        out[i] = (m as i128) << shift;
    }
    Some(out)
}

/// Same scaling with unbounded BigInt magnitudes. Public: the tri-tri
/// interval overlap (robust/tri_tri.rs) reuses it to compare interval
/// endpoints along the common plane-intersection line in pure integers.
pub fn scaled_big<const N: usize>(vs: [f64; N]) -> [BigInt; N] {
    let d = vs.map(decomp);
    let emin = d
        .iter()
        .filter(|(m, _)| *m != 0)
        .map(|(_, e)| *e)
        .min()
        .unwrap_or(0);
    d.map(|(m, e)| {
        if m == 0 {
            BigInt::zero()
        } else {
            BigInt::from(m) << (e - emin) as u32
        }
    })
}

#[inline]
fn sign_i128(v: i128) -> Sign {
    match v.cmp(&0) {
        std::cmp::Ordering::Less => Sign::Neg,
        std::cmp::Ordering::Equal => Sign::Zero,
        std::cmp::Ordering::Greater => Sign::Pos,
    }
}

fn sign_big(v: &BigInt) -> Sign {
    match v.sign() {
        num_bigint::Sign::Minus => Sign::Neg,
        num_bigint::Sign::NoSign => Sign::Zero,
        num_bigint::Sign::Plus => Sign::Pos,
    }
}

/// Exact sign of the orient2d determinant cross(b-a, c-a) on raw f64 points.
///
/// i128 budget: scaled values ≤ 61 bits ⇒ differences ≤ 62, products ≤ 124,
/// final subtraction ≤ 125 bits — safely inside i128.
pub fn orient2d_i(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> Sign {
    let xs = scaled_i128([a[0], b[0], c[0]], 61);
    let ys = scaled_i128([a[1], b[1], c[1]], 61);
    if let (Some([ax, bx, cx]), Some([ay, by, cy])) = (xs, ys) {
        return sign_i128((bx - ax) * (cy - ay) - (by - ay) * (cx - ax));
    }
    let [ax, bx, cx] = scaled_big([a[0], b[0], c[0]]);
    let [ay, by, cy] = scaled_big([a[1], b[1], c[1]]);
    sign_big(&((&bx - &ax) * (&cy - &ay) - (&by - &ay) * (&cx - &ax)))
}

/// Exact sign of the orient3d determinant dot(cross(b-a, c-a), d-a) on raw
/// f64 points.
///
/// i128 budget: scaled values ≤ 40 bits ⇒ differences ≤ 41, 2-products
/// ≤ 82, their difference ≤ 83, 3-products ≤ 124, 3-term sum ≤ 126 bits.
pub fn orient3d_i(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> Sign {
    let xs = scaled_i128([a[0], b[0], c[0], d[0]], 40);
    let ys = scaled_i128([a[1], b[1], c[1], d[1]], 40);
    let zs = scaled_i128([a[2], b[2], c[2], d[2]], 40);
    if let (Some([ax, bx, cx, dx]), Some([ay, by, cy, dy]), Some([az, bz, cz, dz])) = (xs, ys, zs) {
        let (ux, uy, uz) = (bx - ax, by - ay, bz - az);
        let (vx, vy, vz) = (cx - ax, cy - ay, cz - az);
        let (wx, wy, wz) = (dx - ax, dy - ay, dz - az);
        let det = ux * (vy * wz - vz * wy) + uy * (vz * wx - vx * wz) + uz * (vx * wy - vy * wx);
        return sign_i128(det);
    }
    let [ax, bx, cx, dx] = scaled_big([a[0], b[0], c[0], d[0]]);
    let [ay, by, cy, dy] = scaled_big([a[1], b[1], c[1], d[1]]);
    let [az, bz, cz, dz] = scaled_big([a[2], b[2], c[2], d[2]]);
    let (ux, uy, uz) = (&bx - &ax, &by - &ay, &bz - &az);
    let (vx, vy, vz) = (&cx - &ax, &cy - &ay, &cz - &az);
    let (wx, wy, wz) = (&dx - &ax, &dy - &ay, &dz - &az);
    let det = ux * (&vy * &wz - &vz * &wy) + uy * (&vz * &wx - &vx * &wz)
        + uz * (&vx * &wy - &vy * &wx);
    sign_big(&det)
}

#[cfg(test)]
mod tests {
    use super::super::predicates::{orient2d_r, orient3d_r};
    use super::super::rational::{R2, R3};
    use super::*;
    use crate::linalg::{Vec2, Vec3};

    fn o2_ref(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> Sign {
        orient2d_r(
            &R2::from_vec2(Vec2::new(a[0], a[1])),
            &R2::from_vec2(Vec2::new(b[0], b[1])),
            &R2::from_vec2(Vec2::new(c[0], c[1])),
        )
    }

    fn o3_ref(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> Sign {
        orient3d_r(
            &R3::from_vec3(Vec3::new(a[0], a[1], a[2])),
            &R3::from_vec3(Vec3::new(b[0], b[1], b[2])),
            &R3::from_vec3(Vec3::new(c[0], c[1], c[2])),
            &R3::from_vec3(Vec3::new(d[0], d[1], d[2])),
        )
    }

    /// Deterministic LCG so failures reproduce.
    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            self.0
        }
        /// f64 in roughly [-2, 2] with f32-like granularity (the common case:
        /// coordinates originating from f32 MeshGL input).
        fn coord(&mut self) -> f64 {
            (self.next() as i32 as f64) * 2.0f64.powi(-30)
        }
        /// f64 with a wild exponent, exercising the BigInt path.
        fn wild(&mut self) -> f64 {
            let m = self.next() as i32 as f64;
            let e = (self.next() % 500) as i32 - 250;
            let v = m * 2.0f64.powi(e);
            if v.is_finite() { v } else { m }
        }
    }

    #[test]
    fn orient3d_i_matches_rational_generic() {
        let mut rng = Rng(0x5eed);
        for _ in 0..2000 {
            let p: Vec<[f64; 3]> = (0..4)
                .map(|_| [rng.coord(), rng.coord(), rng.coord()])
                .collect();
            assert_eq!(
                orient3d_i(p[0], p[1], p[2], p[3]),
                o3_ref(p[0], p[1], p[2], p[3])
            );
        }
    }

    #[test]
    fn orient3d_i_matches_rational_wild_exponents() {
        let mut rng = Rng(0xbad_cafe);
        for _ in 0..500 {
            let p: Vec<[f64; 3]> = (0..4)
                .map(|_| [rng.wild(), rng.wild(), rng.wild()])
                .collect();
            assert_eq!(
                orient3d_i(p[0], p[1], p[2], p[3]),
                o3_ref(p[0], p[1], p[2], p[3])
            );
        }
    }

    #[test]
    fn orient3d_i_exact_zeros() {
        let mut rng = Rng(0xc0_11a9e);
        for _ in 0..500 {
            // d on the lattice plane through a with directions (b-a), (c-a):
            // small dyadic combination keeps every coordinate exact.
            let a = [rng.coord(), rng.coord(), rng.coord()];
            let b = [rng.coord(), rng.coord(), rng.coord()];
            let c = [rng.coord(), rng.coord(), rng.coord()];
            let s = 0.25;
            let t = 0.5;
            let d = [
                a[0] + s * (b[0] - a[0]) + t * (c[0] - a[0]),
                a[1] + s * (b[1] - a[1]) + t * (c[1] - a[1]),
                a[2] + s * (b[2] - a[2]) + t * (c[2] - a[2]),
            ];
            // The combination is exact only when no rounding occurred; verify
            // against the rational reference either way.
            assert_eq!(orient3d_i(a, b, c, d), o3_ref(a, b, c, d));
        }
    }

    #[test]
    fn orient3d_i_degenerate_inputs() {
        let z = [0.0, -0.0, 0.0];
        let a = [1.0, 2.0, 3.0];
        assert_eq!(orient3d_i(z, z, a, a), Sign::Zero);
        assert_eq!(orient3d_i(a, a, a, a), Sign::Zero);
        let sub = [f64::MIN_POSITIVE / 4.0, 1e300, -1e-300];
        assert_eq!(orient3d_i(z, a, sub, sub), Sign::Zero);
        assert_eq!(
            orient3d_i(z, a, sub, [1.0, 1.0, 1.0]),
            o3_ref(z, a, sub, [1.0, 1.0, 1.0])
        );
    }

    #[test]
    fn orient2d_i_matches_rational() {
        let mut rng = Rng(0x2d2d);
        for _ in 0..2000 {
            let p: Vec<[f64; 2]> = (0..3).map(|_| [rng.coord(), rng.coord()]).collect();
            assert_eq!(orient2d_i(p[0], p[1], p[2]), o2_ref(p[0], p[1], p[2]));
        }
        for _ in 0..500 {
            let p: Vec<[f64; 2]> = (0..3).map(|_| [rng.wild(), rng.wild()]).collect();
            assert_eq!(orient2d_i(p[0], p[1], p[2]), o2_ref(p[0], p[1], p[2]));
        }
        // Exact collinear.
        assert_eq!(
            orient2d_i([0.0, 0.0], [1.0, 1.0], [0.5, 0.5]),
            Sign::Zero
        );
    }
}
