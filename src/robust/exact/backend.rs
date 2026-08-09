// robust/exact/backend.rs — The single seam between the robust engine and its
// arbitrary-precision arithmetic library.
//
// Every module in src/robust that needs unbounded integers or rationals goes
// through the aliases and re-exports here; nothing else in the crate names
// `num_bigint` or `num_rational` directly. Swapping the backend (e.g. to a
// small-value-inline bignum such as dashu) is therefore a change to this file
// plus the small number of backend-coupled call sites inventoried below.
//
// The types are aliases, not newtypes, on purpose: the call sites use the
// backend's operator impls, `Ratio`-style constructors and `num_traits`
// implementations directly, and wrapping them would either hide those or
// require re-implementing several hundred trait impls with no behavioral gain.
//
// ─── Phase-2 checklist: backend-coupled hot spots ───────────────────────────
//
// These are the places whose *correctness* (not just compilation) depends on
// how the backend represents and normalizes values. A backend swap must
// re-verify each one.
//
//  1. `rational.rs::hash_rat` / `rat_fields_eq`, backing `R2Key` / `R3Key`.
//     Structural hashing that reaches into the backend's representation:
//       - `Rational::numer()` / `denom()` -> `&Int`
//       - `Int::sign()` compared against `IntSign::Minus` (sign is hashed
//         separately because the digit iterator yields magnitude only)
//       - `Int::iter_u64_digits()` — little-endian 64-bit limbs of the
//         magnitude, with no leading zero limb, so equal magnitudes always
//         produce identical digit sequences.
//     A replacement backend must expose an equivalent canonical limb/word
//     view (or the hash must be rewritten, e.g. over a byte serialization).
//     Note the digit iterator's *length* is not hashed; the `0xfeed`
//     separator between numerator and denominator is what keeps the stream
//     unambiguous. Preserve that if the limb width changes.
//
//  2. Canonicality invariant behind (1). `R2Key`/`R3Key` equality is
//     *field* equality (numer and denom compared independently), which
//     equals value equality only if every stored rational is fully reduced
//     with a positive denominator. Today that holds because all stored
//     rationals come from `Rational::new` (gcd-reducing), `from_float`, or
//     the arithmetic operators — all of which canonicalize. `new_raw`
//     (no reduction) has no remaining call site in src/; if one is
//     reintroduced, those values must never become keys. A replacement
//     backend must auto-reduce in its equivalent of `Ratio::new` and must
//     normalize the sign onto the numerator, or (1) breaks silently.
//
//  3. `predicates.rs` construction sites that rely on "exactly one gcd":
//     `line_plane_intersect`, `line_line_intersect_2d`, `lift_to_plane` and
//     `segment_param` each build every output coordinate with a single
//     `Rational::new`, so the reduction cost is paid exactly once per
//     coordinate. The gcd is a measured hot spot; re-measure it on a swap.
//
//  4. `intpred.rs` `Int` fallback tier: `scaled_big` (exact dyadic scaling
//     via `Int::from(i64) << u32`) and `sign_big` (`Int::sign()` mapped to
//     the engine's `Sign`). This is the unbounded tier below i64/i128 and
//     is where big-integer performance shows up on degenerate meshes.
//
//  5. `rational.rs::rat_to_f64` — correctly rounded (nearest, ties to even)
//     rational -> f64. It uses, on the *unsigned magnitude* type:
//       `Rational::numer()/denom()` -> `Int`, `Int::magnitude()` -> `&Uint`,
//       `Uint::bits()`, shifts by usize, `Ord`, `/`, `-`, `*`,
//       `Uint::bit(0)`, `+= u32`, `Uint::to_u64()`.
//     It must NOT be replaced by any backend-provided `to_f64()`:
//     num-rational's `ToPrimitive for Ratio` is not correctly rounded, and
//     other backends differ too. Keep this hand-rolled routine on a swap.
//     (`rat_to_f64` is the only rational->f64 rounding path in the engine.)
//
//  6. `predicates.rs::Homog2(pub Int, pub Int, pub Int)` — a homogeneous 2D
//     point holding backend integers directly, and the `[Int; 3]` normals
//     from `tri_normal_int` / `dot_point_raw` / `dot_diff_raw`, plus
//     `tri_tri.rs`'s `type Frac = (Int, Int)`. These are public-ish
//     API surfaces whose signatures change with the backend type.
//
//  7. `Rational::from_float` (`rational.rs::rat`) must be exact for every
//     finite f64 and must return `None`/error for NaN and infinity — the
//     engine relies on the panic path never triggering for finite input.

/// Unbounded signed integer.
pub type Int = num_bigint::BigInt;

/// Unsigned magnitude companion of [`Int`] (used by `rat_to_f64`).
pub type Uint = num_bigint::BigUint;

/// Exact rational, always stored fully reduced with the sign on the numerator.
pub type Rational = num_rational::BigRational;

/// Sign tag of an [`Int`], as returned by `Int::sign()`. Named `IntSign` to
/// avoid colliding with the engine's own `super::Sign` predicate result.
pub use num_bigint::Sign as IntSign;

/// The `num_traits` items the exact-arithmetic call sites need, re-exported so
/// backend-dependent trait paths live in one place too.
pub use num_traits::{One, Signed, ToPrimitive, Zero};
