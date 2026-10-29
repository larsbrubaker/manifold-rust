//! Pure Rust port of the [Manifold](https://github.com/elalish/manifold) 3D
//! geometry library: guaranteed-manifold boolean operations (union,
//! intersection, difference) on triangle meshes, plus constructors, convex
//! hull, Minkowski sum/difference, SDF meshing, smoothing and 2D
//! cross-sections.
//!
//! The port follows the C++ reference's algorithms (v3.5.4) and passes its
//! ported test suite; that suite, with exact expected values, is the bar, and
//! a change is judged on accuracy, robustness, performance and usefulness to
//! callers rather than on matching the C++ more closely. Where the output
//! differs from the C++ on purpose, a comment at the code says how and why;
//! the differences still open are listed in `PORTING_PLAN.md`. The optional
//! `parallel` feature adds rayon parallelism that is restricted to
//! determinism-preserving sites, so its results remain bit-identical to the
//! sequential build.
//!
//! # Example
//!
//! ```
//! use manifold_rust::manifold::Manifold;
//! use manifold_rust::linalg::Vec3;
//!
//! let cube = Manifold::cube(Vec3::new(1.0, 1.0, 1.0), true);
//! let sphere = Manifold::sphere(0.6, 32);
//!
//! let difference = cube.difference(&sphere);
//! assert!(difference.volume() > 0.0);
//! assert_eq!(difference.status(), manifold_rust::types::Error::NoError);
//! ```
//!
//! The main entry point is [`manifold::Manifold`]; meshes move in and out via
//! [`types::MeshGL`] / [`types::MeshGL64`]. For other precision / index
//! combinations, such as f64 vertices with u32 indices, use the generic
//! [`types::MeshGLP`] with [`manifold::Manifold::from_mesh_glp`] /
//! [`manifold::Manifold::get_mesh_glp`]. Lower-level modules mirror the C++
//! source layout and are public primarily for testing and advanced use — their
//! APIs may change before 1.0.

// Clippy lints this crate does not follow, for the whole crate. Other crates
// pin this crate's output bit for bit (manifold-sharp's oracle tests,
// MatterCAD's golden scenes), so a lint whose fix can move a floating-point
// result, or that would churn hundreds of ported lines for style alone, is
// allowed here rather than fixed:
// - `excessive_precision`, `approx_constant`: numeric literals stay as the C++
//   reference writes them; a std constant can differ from the literal in the
//   last bits.
// - `manual_clamp`: `x.min(hi).max(lo)` maps NaN to a bound where `clamp`
//   passes NaN through (and panics on lo > hi).
// - `neg_multiply`, `assign_op_pattern`, `eq_op`: `x * -1.0`, `a = a + b` on
//   the vector types (whose `+=` is a separate impl) and `x - x` (the NaN and
//   infinity checks in `math.rs`) are written as computed.
// - `needless_range_loop`, `too_many_arguments`, `type_complexity`,
//   `field_reassign_with_default`: index loops, signatures and step-by-step
//   struct setup mirror the C++ so the two can be read side by side.
#![allow(
    clippy::excessive_precision,
    clippy::approx_constant,
    clippy::manual_clamp,
    clippy::neg_multiply,
    clippy::assign_op_pattern,
    clippy::eq_op,
    clippy::needless_range_loop,
    clippy::too_many_arguments,
    clippy::type_complexity,
    clippy::field_reassign_with_default
)]

pub mod boolean3;
pub mod boolean_result;
pub mod boolean_stage_progress;
pub mod cancel;
pub mod collider;
pub mod constructors;
pub mod convex_dilation;
pub mod convex_erosion;
mod convex_patches;
pub mod cross_section;
pub mod csg_tree;
pub mod disjoint_sets;
pub mod edge_op;
pub mod face_op;
pub mod impl_mesh;
mod impl_shapes;
mod impl_transform;
pub mod interp_tri;
pub mod linalg;
pub mod manifold;
pub mod math;
pub mod minkowski;
pub mod par;
pub mod polygon;
pub mod progress;
pub mod properties;
pub mod quickhull;
pub mod robust;
pub mod sdf;
pub mod smoothing;
pub mod sort;
pub mod subdivision;
pub mod svd;
pub mod timing;
pub mod tree2d;
pub mod types;
