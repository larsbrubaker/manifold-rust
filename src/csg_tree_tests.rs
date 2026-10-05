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

// Tests for csg_tree.rs: CSG tree evaluation (union / intersection /
// subtract, n-ary flattening, batch reduction order) and the progress
// parameter of `CsgNode::evaluate_with_token_and_progress`.

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
    let result = batch_boolean(OpType::Add, &mut children, None, None, None);
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
    let result = batch_union(&mut children, None, None, None);
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

/// Five overlapping cubes in a row: no two are disjoint along the chain, so
/// `batch_union` cannot compose them away and the tree reduces through real
/// pairwise booleans (enough of them to take `batch_boolean`'s heap path).
fn overlapping_cube_row() -> CsgNode {
    let leaves: Vec<CsgNode> = (0..5)
        .map(|i| {
            CsgNode::leaf(ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(
                Vec3::new(i as f64 * 0.5, i as f64 * 0.1, 0.0),
            ))))
        })
        .collect();
    CsgNode::op_n(OpType::Add, leaves)
}

/// Collects every `(phase, fraction)` callback a run emits.
fn recording_reporter() -> (
    crate::progress::ProgressReporter,
    std::sync::Arc<std::sync::Mutex<Vec<(crate::progress::Phase, Option<f64>)>>>,
) {
    let events = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
    let sink = std::sync::Arc::clone(&events);
    let reporter = crate::progress::ProgressReporter::new(move |phase, fraction| {
        sink.lock().unwrap().push((phase, fraction));
    });
    (reporter, events)
}

#[test]
fn evaluate_with_progress_reports_and_ends_on_a_closed_phase() {
    let (reporter, events) = recording_reporter();
    let result = overlapping_cube_row().evaluate_with_token_and_progress(None, Some(&reporter));
    assert!(result.num_tri() > 0);

    let events = events.lock().unwrap().clone();
    assert!(!events.is_empty(), "an n-ary boolean must report progress");
    // Every fraction a consumer sees is a valid bar position.
    for (_, f) in &events {
        if let Some(f) = f {
            assert!((0.0..=1.0).contains(f), "fraction {f} out of range");
        }
    }
    // The stream ends with the final pairwise boolean's closing report: for
    // a determinate phase that is exactly 1.0, for an indeterminate one
    // (`ExactBoolean`, `Winding`, `Assemble`) it is `None`.
    let (_, last) = *events.last().unwrap();
    assert!(
        last.is_none() || last == Some(1.0),
        "the last report must close a phase, got {last:?}"
    );
}

#[test]
fn evaluate_with_progress_matches_evaluate_with_token() {
    let (reporter, events) = recording_reporter();
    let plain = overlapping_cube_row().evaluate_with_token(None);
    let watched = overlapping_cube_row().evaluate_with_token_and_progress(None, Some(&reporter));
    let unwatched = overlapping_cube_row().evaluate_with_token_and_progress(None, None);
    assert!(!events.lock().unwrap().is_empty());

    for other in [&watched, &unwatched] {
        assert_eq!(other.status, plain.status);
        assert_eq!(other.num_vert(), plain.num_vert());
        assert_eq!(other.num_tri(), plain.num_tri());
        let bits = |m: &ManifoldImpl| -> Vec<u64> {
            m.vert_pos
                .iter()
                .flat_map(|p| [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()])
                .collect()
        };
        assert_eq!(bits(other), bits(&plain), "vertex positions differ");
        assert_eq!(
            other
                .halfedge
                .iter()
                .map(|h| h.start_vert)
                .collect::<Vec<_>>(),
            plain
                .halfedge
                .iter()
                .map(|h| h.start_vert)
                .collect::<Vec<_>>(),
            "topology differs"
        );
    }
}
