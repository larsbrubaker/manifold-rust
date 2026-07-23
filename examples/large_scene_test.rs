// Rust port of cpp-reference/manifold/extras/large_scene_test.cpp.
//
// Unions an n x n x n grid of unit spheres (minus the origin one) and prints
// the same "nTri = N, time = S sec" line as the C++ driver. C++ builds a lazy
// CSG tree via repeated `scene.Boolean(sphere, Add)` and evaluates once at
// NumTri(); on evaluation the nested Adds collapse into one n-ary union that
// BatchUnion processes. We build that collapsed form directly with
// CsgNode::op_n — same evaluation path (BatchUnion + Compose + heap-ordered
// BatchBoolean), without recursing through ~n^3 nested tree levels (the C++
// collapse uses an explicit stack; our collect_children recurses per level and
// a left-deep 8000-node chain would risk stack overflow — flat n-ary is the
// post-collapse equivalent).
//
// Run with: cargo run --release --example large_scene_test [n]

use std::sync::Arc;
use std::time::Instant;

use manifold_rust::csg_tree::{CsgLeafNode, CsgNode};
use manifold_rust::linalg::{mat4_to_mat3x4, translation_matrix, Vec3};
use manifold_rust::manifold::Manifold;

fn main() {
    let n: i32 = std::env::args()
        .nth(1)
        .and_then(|a| a.parse().ok())
        .unwrap_or(20);

    println!("n = {}", n);

    let start = Instant::now();
    // One template sphere; each grid instance is a leaf carrying a lazy
    // translation, exactly like C++ Sphere(1).Translate(...) before evaluation.
    let sphere = Manifold::sphere(1.0, 0);
    // Share one mesh across all leaves, as C++ shares the Impl via shared_ptr.
    let sphere_impl = Arc::new(sphere.as_impl().clone());

    let mut leaves: Vec<CsgNode> = Vec::new();
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                if i == 0 && j == 0 && k == 0 {
                    continue;
                }
                let leaf = CsgLeafNode {
                    p_impl: Arc::clone(&sphere_impl),
                    transform: mat4_to_mat3x4(translation_matrix(Vec3::new(
                        i as f64, j as f64, k as f64,
                    ))),
                };
                leaves.push(CsgNode::leaf_node(leaf));
            }
        }
    }
    let scene = CsgNode::op_n(manifold_rust::types::OpType::Add, leaves).evaluate();
    let n_tri = scene.num_tri();
    let elapsed = start.elapsed().as_secs_f64();
    std::hint::black_box(&scene);
    println!("nTri = {}, time = {} sec", n_tri, elapsed);
}
