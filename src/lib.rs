// manifold-rust: Pure Rust port of the Manifold 3D geometry library
// https://github.com/elalish/manifold
//
// Porting in progress — see PORTING_PLAN.md for current status.

pub mod linalg;
pub mod types;
pub mod polygon;
pub mod impl_mesh;
pub mod sort;
pub mod constructors;
pub mod face_op;
pub mod edge_op;
pub mod properties;
pub mod svd;
pub mod smoothing;
pub mod collider;
pub mod tree2d;
pub mod boolean3;
pub mod boolean_result;
pub mod csg_tree;
pub mod manifold;
pub mod cross_section;
pub mod subdivision;
pub mod sdf;
pub mod minkowski;
pub mod quickhull;
pub mod disjoint_sets;
pub mod math;
pub mod interp_tri;
