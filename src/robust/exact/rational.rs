// robust/exact/rational.rs — Exact rational points and correctly rounded
// rational→f64 conversion for the robust boolean engine.
//
// Input mesh vertices are finite f64 and convert to BigRational exactly
// (`rat`, `R3::from_vec3`). Constructed points — plane/segment intersections
// built in robust/exact/predicates.rs — stay rational through the whole
// pipeline; only output assembly rounds them back, via `rat_to_f64`, which
// rounds to the *nearest* f64 (ties to even, subnormals and overflow
// handled). That single-rounding guarantee is what lets the robust engine's
// output vertices agree with the exact engine's to the last ulp on
// intersection points, and bit-for-bit on pass-through input vertices.

use num_bigint::BigUint;
use num_rational::BigRational;
use num_traits::{Signed, ToPrimitive, Zero};

use crate::linalg::{Vec2, Vec3};

/// Exact conversion of a finite f64. Every finite f64 is a dyadic rational,
/// so this never loses information.
///
/// Panics on NaN/infinity — mesh import rejects non-finite vertices
/// (`Error::NonFiniteVertex`) long before the robust engine runs, so a
/// non-finite value here is an internal logic error, not bad user input.
#[inline]
pub fn rat(v: f64) -> BigRational {
    BigRational::from_float(v).expect("robust engine: coordinate must be finite")
}

/// 2^e as f64, exact for the full representable range -1074..=1023
/// (subnormal powers of two included). Built by bit manipulation so no
/// intermediate rounding can occur.
#[inline]
fn pow2(e: i64) -> f64 {
    debug_assert!((-1074..=1023).contains(&e), "pow2 exponent out of range: {e}");
    if e >= -1022 {
        f64::from_bits(((e + 1023) as u64) << 52)
    } else {
        f64::from_bits(1u64 << (e + 1074))
    }
}

/// Round a rational to the nearest f64, ties to even — the correctly rounded
/// result, identical to rounding the exact real value once. Values beyond
/// f64 range become ±infinity; values below half the smallest subnormal
/// become (signed) zero.
pub fn rat_to_f64(r: &BigRational) -> f64 {
    if r.is_zero() {
        return 0.0;
    }
    let neg = r.is_negative();
    let n: &BigUint = r.numer().magnitude();
    let d: &BigUint = r.denom().magnitude();

    // Exact floor exponent e: 2^e <= n/d < 2^(e+1).
    let mut e = n.bits() as i64 - d.bits() as i64;
    let ge = if e >= 0 {
        *n >= (d << e as usize)
    } else {
        (n << (-e) as usize) >= *d
    };
    if !ge {
        e -= 1;
    }
    if e > 1023 {
        return if neg { f64::NEG_INFINITY } else { f64::INFINITY };
    }

    // Position of the result's least significant bit. Normal numbers carry
    // 53 bits ending at e-52; subnormals are cut off at 2^-1074.
    let lsb = (e - 52).max(-1074);

    // m = round_nearest_even((n/d) / 2^lsb), computed with one exact
    // integer division plus a remainder comparison.
    let (num, den) = if lsb >= 0 {
        (n.clone(), d << lsb as usize)
    } else {
        (n << (-lsb) as usize, d.clone())
    };
    let q = &num / &den;
    let rem = &num - &q * &den;
    let mut m = q;
    let twice_rem = &rem << 1usize;
    match twice_rem.cmp(&den) {
        std::cmp::Ordering::Greater => m += 1u32,
        std::cmp::Ordering::Equal => {
            if m.bit(0) {
                m += 1u32;
            }
        }
        std::cmp::Ordering::Less => {}
    }

    if m.is_zero() {
        return if neg { -0.0 } else { 0.0 };
    }
    // Rounding up may have crossed into the next binade (m = 2^53) or past
    // the largest finite value (2^1024 -> infinity).
    if m.bits() as i64 - 1 + lsb > 1023 {
        return if neg { f64::NEG_INFINITY } else { f64::INFINITY };
    }
    // m <= 2^53, so both the u64 and the f64 conversion are exact, and
    // m * 2^lsb is representable by construction — the multiply is exact.
    let val = m.to_u64().expect("mantissa fits in u64") as f64 * pow2(lsb);
    if neg {
        -val
    } else {
        val
    }
}

// ─── R2 — exact 2D point/vector ─────────────────────────────────────────────

/// Exact 2D point (or vector). Derived `Ord` is lexicographic, which the
/// arrangement code uses for exact point dedup in BTree maps.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct R2 {
    pub x: BigRational,
    pub y: BigRational,
}

impl R2 {
    #[inline]
    pub fn new(x: BigRational, y: BigRational) -> Self {
        Self { x, y }
    }

    #[inline]
    pub fn from_vec2(v: Vec2) -> Self {
        Self::new(rat(v.x), rat(v.y))
    }

    pub fn to_vec2_rounded(&self) -> Vec2 {
        Vec2::new(rat_to_f64(&self.x), rat_to_f64(&self.y))
    }

    pub fn sub(&self, o: &R2) -> R2 {
        R2::new(&self.x - &o.x, &self.y - &o.y)
    }

    pub fn add(&self, o: &R2) -> R2 {
        R2::new(&self.x + &o.x, &self.y + &o.y)
    }

    pub fn scale(&self, s: &BigRational) -> R2 {
        R2::new(&self.x * s, &self.y * s)
    }

    pub fn dot(&self, o: &R2) -> BigRational {
        &self.x * &o.x + &self.y * &o.y
    }

    /// 2D cross product (z of the 3D cross of the embedded vectors).
    pub fn cross(&self, o: &R2) -> BigRational {
        &self.x * &o.y - &self.y * &o.x
    }

    pub fn is_zero(&self) -> bool {
        self.x.is_zero() && self.y.is_zero()
    }
}

// ─── R3 — exact 3D point/vector ─────────────────────────────────────────────

/// Exact 3D point (or vector). Derived `Ord` is lexicographic (x, y, z) for
/// exact vertex welding in output assembly.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct R3 {
    pub x: BigRational,
    pub y: BigRational,
    pub z: BigRational,
}

impl R3 {
    #[inline]
    pub fn new(x: BigRational, y: BigRational, z: BigRational) -> Self {
        Self { x, y, z }
    }

    #[inline]
    pub fn from_vec3(v: Vec3) -> Self {
        Self::new(rat(v.x), rat(v.y), rat(v.z))
    }

    pub fn to_vec3_rounded(&self) -> Vec3 {
        Vec3::new(rat_to_f64(&self.x), rat_to_f64(&self.y), rat_to_f64(&self.z))
    }

    pub fn sub(&self, o: &R3) -> R3 {
        R3::new(&self.x - &o.x, &self.y - &o.y, &self.z - &o.z)
    }

    pub fn add(&self, o: &R3) -> R3 {
        R3::new(&self.x + &o.x, &self.y + &o.y, &self.z + &o.z)
    }

    pub fn scale(&self, s: &BigRational) -> R3 {
        R3::new(&self.x * s, &self.y * s, &self.z * s)
    }

    pub fn dot(&self, o: &R3) -> BigRational {
        &self.x * &o.x + &self.y * &o.y + &self.z * &o.z
    }

    pub fn cross(&self, o: &R3) -> R3 {
        R3::new(
            &self.y * &o.z - &self.z * &o.y,
            &self.z * &o.x - &self.x * &o.z,
            &self.x * &o.y - &self.y * &o.x,
        )
    }

    pub fn is_zero(&self) -> bool {
        self.x.is_zero() && self.y.is_zero() && self.z.is_zero()
    }

    /// Drop the coordinate at `axis` (0=x, 1=y, 2=z), keeping the other two
    /// in cyclic order — the paper's bijective dominant-axis projection used
    /// to embed per-triangle 2D arrangements.
    pub fn project_drop(&self, axis: usize) -> R2 {
        match axis {
            0 => R2::new(self.y.clone(), self.z.clone()),
            1 => R2::new(self.z.clone(), self.x.clone()),
            2 => R2::new(self.x.clone(), self.y.clone()),
            _ => unreachable!("axis must be 0, 1, or 2"),
        }
    }
}
