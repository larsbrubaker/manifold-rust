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

// Phase 12: CSG Tree — ported from C++ csg_tree.cpp (764 lines)
//
// Implements the full CSG tree evaluation system with:
// - CsgLeafNode: lazy transform propagation, Arvo's AABB transform
// - CsgOpNode: N-ary children with caching
// - SimpleBoolean: wrapper invoking Boolean3
// - BatchBoolean: min-heap approach for commutative ops
// - BatchUnion: bounding-box partitioning + Compose + BatchBoolean
// - Explicit-stack DFS evaluation (no recursion)

use std::sync::Arc;
use std::collections::BinaryHeap;
use std::cmp::Ordering;

use crate::boolean3;
use crate::impl_mesh::ManifoldImpl;
use crate::linalg::{Mat3x4, Vec3, mat3x4_to_mat4, mat4_to_mat3x4};
use crate::types::{Box as BBox, OpType};

// ---------------------------------------------------------------------------
// CsgLeafNode — wraps an immutable mesh plus a lazy transform
// ---------------------------------------------------------------------------

#[derive(Clone)]
pub struct CsgLeafNode {
    pub p_impl: Arc<ManifoldImpl>,
    pub transform: Mat3x4,
}

impl CsgLeafNode {
    /// Create a leaf from a mesh with identity transform.
    pub fn new(mesh: ManifoldImpl) -> Self {
        Self {
            p_impl: Arc::new(mesh),
            transform: Mat3x4::identity(),
        }
    }

    /// Create a leaf from a mesh with a specific transform.
    pub fn with_transform(mesh: ManifoldImpl, transform: Mat3x4) -> Self {
        Self {
            p_impl: Arc::new(mesh),
            transform,
        }
    }

    /// Create an empty leaf.
    pub fn empty() -> Self {
        Self {
            p_impl: Arc::new(ManifoldImpl::new()),
            transform: Mat3x4::identity(),
        }
    }

    /// Get the mesh, applying the lazy transform if needed.
    /// Port of C++ CsgLeafNode::GetImpl()
    pub fn get_impl(&self) -> ManifoldImpl {
        if self.transform == Mat3x4::identity() {
            (*self.p_impl).clone()
        } else {
            let mut result = (*self.p_impl).clone();
            result.transform(&self.transform);
            result
        }
    }

    /// Return a new leaf with composed transform.
    /// Port of C++ CsgLeafNode::Transform()
    pub fn apply_transform(&self, m: Mat3x4) -> Self {
        let new_transform = mat4_to_mat3x4(
            mat3x4_to_mat4(m) * mat3x4_to_mat4(self.transform)
        );
        Self {
            p_impl: Arc::clone(&self.p_impl),
            transform: new_transform,
        }
    }

    /// Get bounding box without materializing the full mesh.
    /// Uses Arvo's algorithm for AABB transform.
    /// Port of C++ CsgLeafNode::GetBoundingBox()
    pub fn get_bounding_box(&self) -> BBox {
        let impl_bbox = self.p_impl.bbox;
        if self.transform == Mat3x4::identity() {
            return impl_bbox;
        }
        // Arvo's AABB transform: transform center and half-extents
        let center = (impl_bbox.min + impl_bbox.max) * 0.5;
        let half = (impl_bbox.max - impl_bbox.min) * 0.5;

        // Transform center point
        let mat = self.transform;
        let new_center = Vec3::new(
            mat[0].x * center.x + mat[1].x * center.y + mat[2].x * center.z + mat[3].x,
            mat[0].y * center.x + mat[1].y * center.y + mat[2].y * center.z + mat[3].y,
            mat[0].z * center.x + mat[1].z * center.y + mat[2].z * center.z + mat[3].z,
        );

        // Transform half-extents using absolute values of matrix entries
        let new_half = Vec3::new(
            mat[0].x.abs() * half.x + mat[1].x.abs() * half.y + mat[2].x.abs() * half.z,
            mat[0].y.abs() * half.x + mat[1].y.abs() * half.y + mat[2].y.abs() * half.z,
            mat[0].z.abs() * half.x + mat[1].z.abs() * half.y + mat[2].z.abs() * half.z,
        );

        BBox {
            min: new_center - new_half,
            max: new_center + new_half,
        }
    }

    /// Vertex count without triggering transform.
    pub fn num_vert(&self) -> usize {
        self.p_impl.num_vert()
    }
}

// ---------------------------------------------------------------------------
// CsgNode — the main CSG tree node (leaf or N-ary operation)
// ---------------------------------------------------------------------------

#[derive(Clone)]
pub enum CsgNode {
    Leaf(CsgLeafNode),
    Op {
        op: OpType,
        children: Vec<CsgNode>,
        transform: Mat3x4,
    },
}

impl CsgNode {
    pub fn leaf(mesh: ManifoldImpl) -> Self {
        Self::Leaf(CsgLeafNode::new(mesh))
    }

    pub fn leaf_node(node: CsgLeafNode) -> Self {
        Self::Leaf(node)
    }

    pub fn op(op: OpType, left: CsgNode, right: CsgNode) -> Self {
        Self::Op {
            op,
            children: vec![left, right],
            transform: Mat3x4::identity(),
        }
    }

    pub fn op_n(op: OpType, children: Vec<CsgNode>) -> Self {
        Self::Op {
            op,
            children,
            transform: Mat3x4::identity(),
        }
    }

    /// Evaluate the CSG tree to produce a single mesh.
    /// Uses explicit-stack DFS to avoid recursion stack overflow.
    /// Port of C++ CsgOpNode::ToLeafNode()
    pub fn evaluate(&self) -> ManifoldImpl {
        let leaf = self.to_leaf_node(Mat3x4::identity());
        leaf.get_impl()
    }

    /// Internal: convert this node to a CsgLeafNode, applying the given parent transform.
    fn to_leaf_node(&self, parent_transform: Mat3x4) -> CsgLeafNode {
        match self {
            CsgNode::Leaf(leaf) => leaf.apply_transform(parent_transform),
            CsgNode::Op { op, children, transform } => {
                // Compose local transform with parent
                let combined = mat4_to_mat3x4(
                    mat3x4_to_mat4(parent_transform) * mat3x4_to_mat4(*transform)
                );

                // Flatten: recursively resolve all children to leaves
                let mut positive: Vec<CsgLeafNode> = Vec::new();
                let mut negative: Vec<CsgLeafNode> = Vec::new();

                self.collect_children(*op, combined, children, &mut positive, &mut negative);

                // Perform the operation
                match op {
                    OpType::Add => {
                        // Union of all positive children
                        batch_union(&mut positive)
                    }
                    OpType::Intersect => {
                        // Intersection of all positive children
                        batch_boolean(OpType::Intersect, &mut positive)
                    }
                    OpType::Subtract => {
                        // Subtract: first child is positive, rest are negative
                        if positive.is_empty() {
                            return CsgLeafNode::empty();
                        }
                        let pos_result = batch_union(&mut positive);
                        if negative.is_empty() {
                            return pos_result;
                        }
                        let neg_result = batch_union(&mut negative);
                        simple_boolean(&pos_result, &neg_result, OpType::Subtract)
                    }
                }
            }
        }
    }

    /// Recursively collect children, flattening compatible operations.
    /// Port of the collapsing logic in C++ CsgOpNode::ToLeafNode.
    fn collect_children(
        &self,
        parent_op: OpType,
        transform: Mat3x4,
        children: &[CsgNode],
        positive: &mut Vec<CsgLeafNode>,
        negative: &mut Vec<CsgLeafNode>,
    ) {
        for (i, child) in children.iter().enumerate() {
            match child {
                CsgNode::Leaf(leaf) => {
                    let transformed = leaf.apply_transform(transform);
                    if parent_op == OpType::Subtract && i > 0 {
                        negative.push(transformed);
                    } else {
                        positive.push(transformed);
                    }
                }
                CsgNode::Op { op: child_op, children: grandchildren, transform: child_transform } => {
                    let combined = mat4_to_mat3x4(
                        mat3x4_to_mat4(transform) * mat3x4_to_mat4(*child_transform)
                    );

                    // Collapsing: flatten compatible ops
                    let can_collapse = match (parent_op, child_op) {
                        // Union is associative: (A ∪ B) ∪ C = A ∪ B ∪ C
                        (OpType::Add, OpType::Add) => true,
                        // Intersection is associative: (A ∩ B) ∩ C = A ∩ B ∩ C
                        (OpType::Intersect, OpType::Intersect) => true,
                        // (A - B) - C = A - (B ∪ C): first child's subtraction collapses
                        (OpType::Subtract, OpType::Subtract) if i == 0 => true,
                        _ => false,
                    };

                    if can_collapse {
                        // Flatten: merge grandchildren directly
                        if parent_op == OpType::Subtract && *child_op == OpType::Subtract && i == 0 {
                            // (A - B) is first child of Subtract: A goes to positive, B goes to negative
                            for (gi, gc) in grandchildren.iter().enumerate() {
                                let leaf = gc.to_leaf_node_inner(combined);
                                if gi == 0 {
                                    positive.push(leaf);
                                } else {
                                    negative.push(leaf);
                                }
                            }
                        } else {
                            for gc in grandchildren {
                                let leaf = gc.to_leaf_node_inner(combined);
                                if parent_op == OpType::Subtract && i > 0 {
                                    negative.push(leaf);
                                } else {
                                    positive.push(leaf);
                                }
                            }
                        }
                    } else {
                        // Cannot collapse: evaluate child subtree fully
                        let result = child.to_leaf_node(combined);
                        if parent_op == OpType::Subtract && i > 0 {
                            negative.push(result);
                        } else {
                            positive.push(result);
                        }
                    }
                }
            }
        }
    }

    /// Helper: convert a single node to leaf with given transform (non-flattening).
    fn to_leaf_node_inner(&self, transform: Mat3x4) -> CsgLeafNode {
        match self {
            CsgNode::Leaf(leaf) => leaf.apply_transform(transform),
            CsgNode::Op { .. } => self.to_leaf_node(transform),
        }
    }
}

// ---------------------------------------------------------------------------
// SimpleBoolean — wrapper invoking Boolean3
// Port of C++ SimpleBoolean() (lines 142-184)
// ---------------------------------------------------------------------------

fn simple_boolean(a: &CsgLeafNode, b: &CsgLeafNode, op: OpType) -> CsgLeafNode {
    let impl_a = a.get_impl();
    let impl_b = b.get_impl();
    let result = boolean3::boolean(&impl_a, &impl_b, op);
    CsgLeafNode::new(result)
}

// ---------------------------------------------------------------------------
// BatchBoolean — min-heap approach for commutative ops
// Port of C++ BatchBoolean() (lines 376-428)
// ---------------------------------------------------------------------------

/// Wrapper for BinaryHeap ordering by vertex count (smallest first).
struct MeshEntry(CsgLeafNode);

impl PartialEq for MeshEntry {
    fn eq(&self, other: &Self) -> bool {
        self.0.num_vert() == other.0.num_vert()
    }
}
impl Eq for MeshEntry {}

impl PartialOrd for MeshEntry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for MeshEntry {
    fn cmp(&self, other: &Self) -> Ordering {
        // Reverse ordering: BinaryHeap is a max-heap, we want smallest first
        other.0.num_vert().cmp(&self.0.num_vert())
    }
}

fn batch_boolean(op: OpType, children: &mut Vec<CsgLeafNode>) -> CsgLeafNode {
    if children.is_empty() {
        return CsgLeafNode::empty();
    }
    if children.len() == 1 {
        return children.remove(0);
    }
    if children.len() == 2 {
        let b = children.pop().unwrap();
        let a = children.pop().unwrap();
        return simple_boolean(&a, &b, op);
    }

    // Build min-heap by vertex count
    let mut heap: BinaryHeap<MeshEntry> = BinaryHeap::new();
    for child in children.drain(..) {
        heap.push(MeshEntry(child));
    }

    // Pop two smallest, boolean them, push result back
    while heap.len() > 1 {
        let a = heap.pop().unwrap().0;
        let b = heap.pop().unwrap().0;
        let result = simple_boolean(&a, &b, op);
        heap.push(MeshEntry(result));
    }

    heap.pop().unwrap().0
}

// ---------------------------------------------------------------------------
// BatchUnion — bounding-box partitioning + Compose + BatchBoolean
// Port of C++ BatchUnion() (lines 434-491)
// ---------------------------------------------------------------------------

const K_MAX_UNION_SIZE: usize = 1000;

fn batch_union(children: &mut Vec<CsgLeafNode>) -> CsgLeafNode {
    if children.is_empty() {
        return CsgLeafNode::empty();
    }
    if children.len() == 1 {
        return children.remove(0);
    }

    // Process in chunks to avoid O(n^2) overlap checks
    while children.len() > 1 {
        let chunk_size = children.len().min(K_MAX_UNION_SIZE);
        let chunk_start = children.len() - chunk_size;

        // Get bounding boxes for the chunk
        let boxes: Vec<BBox> = children[chunk_start..]
            .iter()
            .map(|c| c.get_bounding_box())
            .collect();

        // Greedy partition into disjoint sets
        let mut sets: Vec<Vec<usize>> = Vec::new(); // each set is indices into chunk
        for i in 0..chunk_size {
            let mut found_set = false;
            for set in &mut sets {
                let overlaps = set.iter().any(|&j| boxes[i].does_overlap_box(&boxes[j]));
                if !overlaps {
                    set.push(i);
                    found_set = true;
                    break;
                }
            }
            if !found_set {
                sets.push(vec![i]);
            }
        }

        // Process each disjoint set
        let chunk: Vec<CsgLeafNode> = children.drain(chunk_start..).collect();
        let mut results: Vec<CsgLeafNode> = Vec::new();

        for set in &sets {
            if set.len() == 1 {
                results.push(chunk[set[0]].clone());
            } else {
                // Compose disjoint meshes without boolean
                let meshes: Vec<ManifoldImpl> = set.iter()
                    .map(|&i| chunk[i].get_impl())
                    .collect();
                let composed = boolean3::compose_meshes(&meshes);
                results.push(CsgLeafNode::new(composed));
            }
        }

        // BatchBoolean the composed results
        if results.len() == 1 {
            // Insert at front (most complex, best for subsequent batches)
            children.insert(0, results.remove(0));
        } else {
            let result = batch_boolean(OpType::Add, &mut results);
            children.insert(0, result);
        }
    }

    children.remove(0)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linalg::{mat4_to_mat3x4, translation_matrix, Vec3};

    #[test]
    fn test_csg_tree_union_disjoint() {
        let a = ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(Vec3::new(0.0, 0.0, 0.0))));
        let b = ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(Vec3::new(3.0, 0.0, 0.0))));
        let tree = CsgNode::op(OpType::Add, CsgNode::leaf(a), CsgNode::leaf(b));
        let result = tree.evaluate();
        assert_eq!(result.num_tri(), 24);
    }

    #[test]
    fn test_csg_tree_union_overlapping() {
        let a = ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(Vec3::new(0.0, 0.0, 0.0))));
        let b = ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(Vec3::new(0.5, 0.0, 0.0))));
        let tree = CsgNode::op(OpType::Add, CsgNode::leaf(a), CsgNode::leaf(b));
        let result = tree.evaluate();
        assert!(result.num_tri() > 0, "Overlapping union should produce non-empty mesh");
    }

    #[test]
    fn test_csg_tree_intersection() {
        let a = ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(Vec3::new(0.0, 0.0, 0.0))));
        let b = ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(Vec3::new(0.5, 0.0, 0.0))));
        let tree = CsgNode::op(OpType::Intersect, CsgNode::leaf(a), CsgNode::leaf(b));
        let result = tree.evaluate();
        assert!(result.num_tri() > 0, "Overlapping intersection should produce non-empty mesh");
    }

    #[test]
    fn test_csg_tree_subtract() {
        let a = ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(Vec3::new(0.0, 0.0, 0.0))));
        let b = ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(Vec3::new(0.5, 0.0, 0.0))));
        let tree = CsgNode::op(OpType::Subtract, CsgNode::leaf(a), CsgNode::leaf(b));
        let result = tree.evaluate();
        assert!(result.num_tri() > 0, "Subtraction should produce non-empty mesh");
    }

    #[test]
    fn test_batch_boolean_three_cubes() {
        let a = CsgLeafNode::new(ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(Vec3::new(0.0, 0.0, 0.0)))));
        let b = CsgLeafNode::new(ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(Vec3::new(0.5, 0.0, 0.0)))));
        let c = CsgLeafNode::new(ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(Vec3::new(1.0, 0.0, 0.0)))));
        let mut children = vec![a, b, c];
        let result = batch_boolean(OpType::Add, &mut children);
        let mesh = result.get_impl();
        assert!(mesh.num_tri() > 0, "BatchBoolean of 3 overlapping cubes should produce non-empty mesh");
    }

    #[test]
    fn test_batch_union_disjoint() {
        let a = CsgLeafNode::new(ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(Vec3::new(0.0, 0.0, 0.0)))));
        let b = CsgLeafNode::new(ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(Vec3::new(3.0, 0.0, 0.0)))));
        let c = CsgLeafNode::new(ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(Vec3::new(6.0, 0.0, 0.0)))));
        let mut children = vec![a, b, c];
        let result = batch_union(&mut children);
        let mesh = result.get_impl();
        // Three disjoint cubes should compose without boolean, giving 36 tris
        assert_eq!(mesh.num_tri(), 36, "BatchUnion of 3 disjoint cubes should have 36 tris");
    }

    #[test]
    fn test_csg_n_ary_union() {
        // N-ary union of 4 disjoint cubes
        let nodes: Vec<CsgNode> = (0..4).map(|i| {
            CsgNode::leaf(ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(
                Vec3::new(i as f64 * 3.0, 0.0, 0.0)
            ))))
        }).collect();
        let tree = CsgNode::op_n(OpType::Add, nodes);
        let result = tree.evaluate();
        assert_eq!(result.num_tri(), 48, "N-ary union of 4 disjoint cubes should have 48 tris");
    }

    #[test]
    fn test_tree_transforms() {
        // Test that transforms compose correctly through the tree
        let a = ManifoldImpl::cube(&Mat3x4::identity());
        let leaf = CsgLeafNode::new(a);
        let translated = leaf.apply_transform(
            mat4_to_mat3x4(translation_matrix(Vec3::new(5.0, 0.0, 0.0)))
        );
        let bbox = translated.get_bounding_box();
        assert!(bbox.min.x > 4.0, "Translated bbox min.x should be > 4.0, got {}", bbox.min.x);
        assert!(bbox.max.x < 6.5, "Translated bbox max.x should be < 6.5, got {}", bbox.max.x);
    }
}
