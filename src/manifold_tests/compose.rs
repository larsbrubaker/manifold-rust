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
// operator+ / GetMeshGL on the same inputs). The last test pins a full-MeshGL
// hash of a CSG-tree union whose pairs run in parallel under the `parallel`
// feature (csg_tree::batch_boolean), so the run table cannot follow the order
// in which those pairs reserve mesh IDs.

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

/// FNV-1a 64 over every MeshGL64 field, each array prefixed by its length
/// (u64 LE): num_prop, vert_properties, tri_verts, merge_from_vert,
/// merge_to_vert, run_index, run_original_id (minus `base_id`, as u32 LE, so
/// the hash does not depend on where the process-wide ID counter stood),
/// run_transform, face_id, halfedge_tangent, run_flags, then tolerance.
/// Floats hash as their IEEE bits, integers as LE bytes of their MeshGL64
/// width.
fn mesh_gl64_hash(gl: &crate::types::MeshGL64, base_id: u32) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    let mut eat = |bytes: &[u8]| {
        for &b in bytes {
            h = (h ^ u64::from(b)).wrapping_mul(0x0100_0000_01b3);
        }
    };
    let len = |n: usize| (n as u64).to_le_bytes();
    eat(&gl.num_prop.to_le_bytes());
    eat(&len(gl.vert_properties.len()));
    gl.vert_properties
        .iter()
        .for_each(|x| eat(&x.to_bits().to_le_bytes()));
    for ints in [
        &gl.tri_verts,
        &gl.merge_from_vert,
        &gl.merge_to_vert,
        &gl.run_index,
    ] {
        eat(&len(ints.len()));
        ints.iter().for_each(|x| eat(&x.to_le_bytes()));
    }
    eat(&len(gl.run_original_id.len()));
    for id in &gl.run_original_id {
        eat(&id.wrapping_sub(base_id).to_le_bytes());
    }
    eat(&len(gl.run_transform.len()));
    gl.run_transform
        .iter()
        .for_each(|x| eat(&x.to_bits().to_le_bytes()));
    eat(&len(gl.face_id.len()));
    gl.face_id.iter().for_each(|x| eat(&x.to_le_bytes()));
    eat(&len(gl.halfedge_tangent.len()));
    gl.halfedge_tangent
        .iter()
        .for_each(|x| eat(&x.to_bits().to_le_bytes()));
    eat(&len(gl.run_flags.len()));
    eat(&gl.run_flags);
    eat(&gl.tolerance.to_bits().to_le_bytes());
    h
}

/// Four spheres sharing one original ID — two of 8 segments at y = 0, two
/// of 256 at y = 2.2 — each rotated 45 degrees about z through its leaf's
/// lazy transform, so the leaves' bounding boxes overlap (batch_union cannot
/// compose them) while the meshes are disjoint (each pair takes the
/// disjoint-Add fast path). Under `parallel` the pairs of a batch_boolean
/// round reserve mesh IDs in whatever order they finish; the run table must
/// not follow that order.
fn rotated_instanced_sphere_union() -> u64 {
    let id = Manifold::reserve_ids(1);
    let sphere = |segments: i32| {
        let mut gl = Manifold::sphere(1.0, segments).get_mesh_gl64(-1);
        gl.run_original_id = vec![id];
        gl.run_index.clear();
        gl.run_transform.clear();
        gl.run_flags.clear();
        gl.face_id.clear();
        Manifold::from_mesh_gl64(&gl)
    };
    let (lo, hi) = (sphere(8), sphere(256));
    let r = std::f64::consts::FRAC_1_SQRT_2;
    let leaf = |m: &Manifold, x: f64, y: f64| {
        let t = Mat3x4::from_cols(
            Vec3::new(r, r, 0.0),
            Vec3::new(-r, r, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(x, y, 0.0),
        );
        CsgNode::leaf_node(CsgLeafNode::with_transform(m.as_impl().clone(), t))
    };
    let leaves = vec![
        leaf(&lo, 0.0, 0.0),
        leaf(&lo, 2.2, 0.0),
        leaf(&hi, 0.0, 2.2),
        leaf(&hi, 2.2, 2.2),
    ];
    let sum = Manifold::from_impl(CsgNode::op_n(OpType::Add, leaves).evaluate());
    let gl = sum.get_mesh_gl64(-1);
    assert_eq!(gl.run_original_id, vec![id; 4], "one run per sphere");
    mesh_gl64_hash(&gl, id)
}

const ROTATED_INSTANCED_SPHERES_HASH: u64 = 0xe55d_84ec_b9f6_2e40;

/// The full MeshGL of the union above is pinned, and with `parallel` it is
/// the same on 1 and on 8 threads, repeated, because compose_meshes ranks
/// runs by input order rather than by the absolute mesh IDs.
#[test]
fn test_batch_union_of_rotated_instanced_leaves_is_independent_of_scheduling() {
    #[cfg(feature = "parallel")]
    for threads in [1, 8] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        for rep in 0..8 {
            assert_eq!(
                pool.install(rotated_instanced_sphere_union),
                ROTATED_INSTANCED_SPHERES_HASH,
                "{threads} threads, repetition {rep}"
            );
        }
    }
    assert_eq!(
        rotated_instanced_sphere_union(),
        ROTATED_INSTANCED_SPHERES_HASH
    );
}
