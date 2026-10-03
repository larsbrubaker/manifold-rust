// Tests for how composing disjoint meshes keeps each input's runs — the
// MeshGL run table (run_index / run_original_id / run_transform) that
// boolean3::compose_meshes builds for Manifold::compose, the disjoint-union
// fast path in boolean3::boolean_with_token and csg_tree's batch_union.
//
// C++ Compose (csg_tree.cpp:384-410) offsets node i's meshIDs by
// i * meshIDCounter before IncrementMeshIDs, so instanced copies of one mesh
// stay separate runs, each with its own transform, ordered node by node. The
// expected values below were captured from the v3.5.2 C++ reference
// (cpp-reference/manifold built as a static library, Manifold::Compose /
// operator+ / GetMeshGL on the same inputs).

use super::*;
use crate::csg_tree::{CsgLeafNode, CsgNode};
use crate::linalg::Mat3x4;

const IDENTITY: [f32; 12] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0];

fn translated(x: f32, y: f32, z: f32) -> [f32; 12] {
    let mut t = IDENTITY;
    t[9] = x;
    t[10] = y;
    t[11] = z;
    t
}

/// Asserts the run table matches C++: one run per input node, in node order,
/// each `tris_per_run` triangles, all carrying `original_id`, with the given
/// per-run transforms.
fn assert_runs(m: &Manifold, original_id: i32, tris_per_run: u32, transforms: &[[f32; 12]]) {
    let gl = m.get_mesh_gl(-1);
    let n = transforms.len();
    let expected_index: Vec<u32> = (0..=n as u32).map(|r| 3 * tris_per_run * r).collect();
    assert_eq!(gl.run_index, expected_index, "run_index");
    assert_eq!(
        gl.run_original_id,
        vec![original_id as u32; n],
        "run_original_id"
    );
    let expected_transform: Vec<f32> = transforms.iter().flatten().copied().collect();
    assert_eq!(gl.run_transform, expected_transform, "run_transform");
}

/// C++: Manifold::Compose({cube, cube.Translate({3, 0, 0})}) gives two runs
/// of the cube's original ID, identity then the translation. Composing the
/// same mesh twice must not merge the copies into one run.
#[test]
fn test_compose_keeps_each_instanced_copy_as_its_own_run() {
    let cube = Manifold::cube(Vec3::splat(1.0), false);
    let composed = Manifold::compose(&[cube.clone(), cube.translate(Vec3::new(3.0, 0.0, 0.0))]);
    assert_runs(
        &composed,
        cube.original_id(),
        12,
        &[IDENTITY, translated(3.0, 0.0, 0.0)],
    );
}

/// C++: Compose of three copies keeps node order — the (0, 5, 0) copy first
/// even though it sorts last geometrically and none of the copies' mesh IDs
/// differ.
#[test]
fn test_compose_orders_instanced_copies_node_by_node() {
    let cube = Manifold::cube(Vec3::splat(1.0), false);
    let composed = Manifold::compose(&[
        cube.translate(Vec3::new(0.0, 5.0, 0.0)),
        cube.clone(),
        cube.translate(Vec3::new(3.0, 0.0, 0.0)),
    ]);
    assert_runs(
        &composed,
        cube.original_id(),
        12,
        &[
            translated(0.0, 5.0, 0.0),
            IDENTITY,
            translated(3.0, 0.0, 0.0),
        ],
    );
}

/// C++: cube + cube.Translate({3, 0, 0}) — disjoint, so C++ BatchUnion
/// Composes them; Rust takes the disjoint-Add fast path in
/// boolean_with_token. Same two runs as Compose.
#[test]
fn test_disjoint_union_of_instanced_copies_keeps_both_runs() {
    let cube = Manifold::cube(Vec3::splat(1.0), false);
    let sum = cube.union(&cube.translate(Vec3::new(3.0, 0.0, 0.0)));
    assert_runs(
        &sum,
        cube.original_id(),
        12,
        &[IDENTITY, translated(3.0, 0.0, 0.0)],
    );
}

/// The CSG tree's batch_union composes disjoint leaves that share one
/// impl and differ only in their lazy leaf transform: each leaf keeps its
/// own run, in leaf order.
#[test]
fn test_csg_batch_union_keeps_each_instanced_leaf_as_its_own_run() {
    let cube = Manifold::cube(Vec3::splat(1.0), false);
    let leaf = |x: f64| {
        let t = Mat3x4::from_cols(
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(x, 0.0, 0.0),
        );
        CsgNode::leaf_node(CsgLeafNode::with_transform(cube.as_impl().clone(), t))
    };
    let tree = CsgNode::op_n(OpType::Add, vec![leaf(0.0), leaf(3.0), leaf(6.0)]);
    let sum = Manifold::from_impl(tree.evaluate());
    assert_runs(
        &sum,
        cube.original_id(),
        12,
        &[
            IDENTITY,
            translated(3.0, 0.0, 0.0),
            translated(6.0, 0.0, 0.0),
        ],
    );
}
