// Copyright 2026 Lars Brubaker
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// Tests for the CSG tree (src/csg_tree.rs): leaf transforms, union /
// intersection / difference evaluation, and the batch reduction rounds whose
// output must not depend on the thread count.

use super::*;
use crate::linalg::{mat4_to_mat3x4, translation_matrix, Vec3};

#[test]
fn test_csg_tree_union_disjoint() {
    let a = ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(Vec3::new(
        0.0, 0.0, 0.0,
    ))));
    let b = ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(Vec3::new(
        3.0, 0.0, 0.0,
    ))));
    let tree = CsgNode::op(OpType::Add, CsgNode::leaf(a), CsgNode::leaf(b));
    let result = tree.evaluate();
    assert_eq!(result.num_tri(), 24);
}

#[test]
fn test_csg_tree_union_overlapping() {
    let a = ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(Vec3::new(
        0.0, 0.0, 0.0,
    ))));
    let b = ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(Vec3::new(
        0.5, 0.0, 0.0,
    ))));
    let tree = CsgNode::op(OpType::Add, CsgNode::leaf(a), CsgNode::leaf(b));
    let result = tree.evaluate();
    assert!(
        result.num_tri() > 0,
        "Overlapping union should produce non-empty mesh"
    );
}

#[test]
fn test_csg_tree_intersection() {
    let a = ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(Vec3::new(
        0.0, 0.0, 0.0,
    ))));
    let b = ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(Vec3::new(
        0.5, 0.0, 0.0,
    ))));
    let tree = CsgNode::op(OpType::Intersect, CsgNode::leaf(a), CsgNode::leaf(b));
    let result = tree.evaluate();
    assert!(
        result.num_tri() > 0,
        "Overlapping intersection should produce non-empty mesh"
    );
}

#[test]
fn test_csg_tree_subtract() {
    let a = ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(Vec3::new(
        0.0, 0.0, 0.0,
    ))));
    let b = ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(Vec3::new(
        0.5, 0.0, 0.0,
    ))));
    let tree = CsgNode::op(OpType::Subtract, CsgNode::leaf(a), CsgNode::leaf(b));
    let result = tree.evaluate();
    assert!(
        result.num_tri() > 0,
        "Subtraction should produce non-empty mesh"
    );
}

#[test]
fn test_batch_boolean_three_cubes() {
    let a = CsgLeafNode::new(ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(
        Vec3::new(0.0, 0.0, 0.0),
    ))));
    let b = CsgLeafNode::new(ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(
        Vec3::new(0.5, 0.0, 0.0),
    ))));
    let c = CsgLeafNode::new(ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(
        Vec3::new(1.0, 0.0, 0.0),
    ))));
    let mut children = vec![a, b, c];
    let result = batch_boolean(OpType::Add, &mut children, None);
    let mesh = result.get_impl();
    assert!(
        mesh.num_tri() > 0,
        "BatchBoolean of 3 overlapping cubes should produce non-empty mesh"
    );
}

#[test]
fn test_batch_union_disjoint() {
    let a = CsgLeafNode::new(ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(
        Vec3::new(0.0, 0.0, 0.0),
    ))));
    let b = CsgLeafNode::new(ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(
        Vec3::new(3.0, 0.0, 0.0),
    ))));
    let c = CsgLeafNode::new(ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(
        Vec3::new(6.0, 0.0, 0.0),
    ))));
    let mut children = vec![a, b, c];
    let result = batch_union(&mut children, None);
    let mesh = result.get_impl();
    // Three disjoint cubes should compose without boolean, giving 36 tris
    assert_eq!(
        mesh.num_tri(),
        36,
        "BatchUnion of 3 disjoint cubes should have 36 tris"
    );
}

#[test]
fn test_csg_n_ary_union() {
    // N-ary union of 4 disjoint cubes
    let nodes: Vec<CsgNode> = (0..4)
        .map(|i| {
            CsgNode::leaf(ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(
                Vec3::new(i as f64 * 3.0, 0.0, 0.0),
            ))))
        })
        .collect();
    let tree = CsgNode::op_n(OpType::Add, nodes);
    let result = tree.evaluate();
    assert_eq!(
        result.num_tri(),
        48,
        "N-ary union of 4 disjoint cubes should have 48 tris"
    );
}

#[test]
fn test_lazy_leaf_transform_applied_on_evaluate() {
    // Regression: get_impl discarded ManifoldImpl::transform's return value
    // (it is not in-place), so lazily-transformed leaves evaluated at the
    // origin. Two disjoint cubes — one translated via the *leaf* transform,
    // not baked into the mesh — must union to 24 tris, not collapse to 12.
    let cube = ManifoldImpl::cube(&Mat3x4::identity());
    let a = CsgLeafNode::new(cube.clone());
    let b = CsgLeafNode::new(cube)
        .apply_transform(mat4_to_mat3x4(translation_matrix(Vec3::new(3.0, 0.0, 0.0))));
    let bbox = b.get_impl().bbox;
    assert!(
        bbox.min.x >= 2.9 && bbox.max.x <= 4.1,
        "lazy transform not applied by get_impl: bbox.x = [{}, {}]",
        bbox.min.x,
        bbox.max.x
    );
    let tree = CsgNode::op(OpType::Add, CsgNode::leaf_node(a), CsgNode::leaf_node(b));
    assert_eq!(tree.evaluate().num_tri(), 24);
}

#[test]
fn test_tree_transforms() {
    // Test that transforms compose correctly through the tree
    let a = ManifoldImpl::cube(&Mat3x4::identity());
    let leaf = CsgLeafNode::new(a);
    let translated =
        leaf.apply_transform(mat4_to_mat3x4(translation_matrix(Vec3::new(5.0, 0.0, 0.0))));
    let bbox = translated.get_bounding_box();
    assert!(
        bbox.min.x > 4.0,
        "Translated bbox min.x should be > 4.0, got {}",
        bbox.min.x
    );
    assert!(
        bbox.max.x < 6.5,
        "Translated bbox max.x should be < 6.5, got {}",
        bbox.max.x
    );
}

/// Pins a union of 64 overlapping cubes to the sequential rounds' output,
/// and with `parallel` on 1 and 8 threads. Coordinates are exact, so the
/// hash holds on every platform.
#[test]
fn test_batch_union_rounds_keep_the_sequential_output() {
    let union = || {
        let leaves = (0..64)
            .map(|i| {
                let at = Vec3::new(f64::from(i % 8) * 0.5, f64::from(i / 8) * 0.5, 0.0);
                CsgNode::leaf(ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(at))))
            })
            .collect();
        let gl =
            crate::manifold::Manifold::from_impl(CsgNode::op_n(OpType::Add, leaves).evaluate())
                .get_mesh_gl64(-1);
        let mut h: u64 = 0xcbf2_9ce4_8422_2325;
        let bytes = gl
            .vert_properties
            .iter()
            .flat_map(|x| x.to_bits().to_le_bytes());
        let ints = gl.tri_verts.iter().chain(&gl.run_index);
        for b in bytes.chain(ints.flat_map(|x| x.to_le_bytes())) {
            h = (h ^ u64::from(b)).wrapping_mul(0x0100_0000_01b3);
        }
        h
    };
    #[cfg(feature = "parallel")]
    for threads in [1, 8] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        assert_eq!(
            pool.install(union),
            0x0045_f3cc_e459_28c2,
            "{threads} threads"
        );
    }
    assert_eq!(union(), 0x0045_f3cc_e459_28c2);
}

/// Everything a batch result exposes about its mesh IDs: the relation, with
/// each mesh ID replaced by its rank, and the whole `MeshGL64`. Ranks rather
/// than values because other tests reserve IDs from the same process-wide
/// counter while this one runs; the order of the IDs is what decides the run
/// and triangle order of the output.
fn mesh_id_fingerprint(m: &ManifoldImpl) -> String {
    let rel = &m.mesh_relation;
    let rank: std::collections::HashMap<i32, usize> = rel
        .mesh_id_transform
        .keys()
        .enumerate()
        .map(|(r, &id)| (id, r))
        .collect();
    let tri_ref: Vec<_> = rel
        .tri_ref
        .iter()
        .map(|t| {
            (
                rank.get(&t.mesh_id),
                t.original_id,
                t.face_id,
                t.coplanar_id,
            )
        })
        .collect();
    let relations: Vec<_> = rel.mesh_id_transform.values().collect();
    let gl = crate::manifold::Manifold::from_impl(m.clone()).get_mesh_gl64(-1);
    format!("{}\n{relations:?}\n{tri_ref:?}\n{gl:?}", rel.original_id)
}

/// Batch rounds hand out mesh IDs in pair order at any thread count. Each
/// boolean reserves IDs from a process-wide counter, so pairs running side by
/// side reserved them in whatever order they reached that point.
///
/// Eight instances of one sphere share its mesh ID, so all runs of the result
/// share an original ID and `MeshGL` orders them by mesh ID alone. Leaves 2j
/// and 2j+1 overlap, so the first round runs four real booleans; the pairs
/// are 10 apart, so every later round composes two disjoint operands, which
/// keeps both operands' IDs and sorts them by value. The output's run order
/// is then the order in which the rounds reserved their IDs.
#[test]
fn test_batch_rounds_assign_mesh_ids_in_pair_order() {
    let sphere = crate::manifold::Manifold::sphere(1.0, 48).as_impl().clone();
    let batch = || {
        let mut leaves: Vec<CsgLeafNode> = (0..8)
            .map(|i| {
                let at = Vec3::new(f64::from(i / 2) * 10.0 + f64::from(i % 2) * 0.5, 0.0, 0.0);
                CsgLeafNode::new(sphere.clone())
                    .apply_transform(mat4_to_mat3x4(translation_matrix(at)))
            })
            .collect();
        let result = batch_boolean(OpType::Add, &mut leaves, None).get_impl();
        // The runs' x offsets name the leaves in output order.
        let gl = crate::manifold::Manifold::from_impl(result.clone()).get_mesh_gl64(-1);
        let runs: Vec<f64> = gl.run_transform.chunks(12).map(|t| t[9]).collect();
        (runs, mesh_id_fingerprint(&result))
    };
    // The order the sequential build gives, pinned so that every build and
    // thread count is held to it.
    let order = [10.5, 10.0, 0.5, 0.0, 30.5, 30.0, 20.5, 20.0];
    let (runs, expected) = batch();
    assert_eq!(runs, order);
    #[cfg(not(feature = "parallel"))]
    assert!(batch().1 == expected, "a second evaluation");
    #[cfg(feature = "parallel")]
    for threads in [1, 2, 3, 4, 8] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        for rep in 0..25 {
            let (runs, fingerprint) = pool.install(batch);
            assert_eq!(runs, order, "{threads} threads, rep {rep}");
            assert!(fingerprint == expected, "{threads} threads, rep {rep}");
        }
    }
}
