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

// Phase 12: CSG Tree
//
// C++ sources: src/csg_tree.cpp (764 lines), src/csg_tree.h
//
// STATUS: MINIMAL — basic tree evaluation only.
// Missing: BatchBoolean, BatchUnion, Compose, lazy transform propagation,
// mesh caching, explicit-stack flattening. These require full boolean3.

use crate::boolean3;
use crate::impl_mesh::ManifoldImpl;
use crate::types::OpType;

#[derive(Clone)]
pub enum CsgNode {
    Leaf(ManifoldImpl),
    Op {
        op: OpType,
        left: Box<CsgNode>,
        right: Box<CsgNode>,
    },
}

impl CsgNode {
    pub fn leaf(mesh: ManifoldImpl) -> Self {
        Self::Leaf(mesh)
    }

    pub fn op(op: OpType, left: CsgNode, right: CsgNode) -> Self {
        Self::Op {
            op,
            left: Box::new(left),
            right: Box::new(right),
        }
    }

    pub fn evaluate(&self) -> ManifoldImpl {
        match self {
            CsgNode::Leaf(mesh) => mesh.clone(),
            CsgNode::Op { op, left, right } => {
                let a = left.evaluate();
                let b = right.evaluate();
                boolean3::boolean(&a, &b, *op)
            }
        }
    }
}

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
}
