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

// Phase 17: Minkowski Sum/Difference
//
// C++ source: src/minkowski.cpp (175 lines)
//
// STATUS: NOT YET IMPLEMENTED
// Requires: QuickHull (done), Boolean3 (not done), CSG BatchBoolean (not done).
// The C++ computes actual Minkowski sums using per-face vertex sums -> hull ->
// BatchBoolean. This placeholder returns bounding-box approximations.

use crate::impl_mesh::ManifoldImpl;
use crate::linalg::{mat4_to_mat3x4, scaling_matrix, translation_matrix, Vec3};

/// Compute the Minkowski sum of two meshes.
///
/// PLACEHOLDER: Returns bounding-box approximation. Will be replaced with
/// proper implementation once Boolean3 and CSG BatchBoolean are ported.
pub fn minkowski_sum(a: &ManifoldImpl, b: &ManifoldImpl) -> ManifoldImpl {
    let min = a.bbox.min + b.bbox.min;
    let max = a.bbox.max + b.bbox.max;
    let t = mat4_to_mat3x4(translation_matrix(min) * scaling_matrix(max - min));
    ManifoldImpl::cube(&t)
}

/// Compute the Minkowski difference of two meshes.
///
/// PLACEHOLDER: Returns bounding-box approximation.
pub fn minkowski_difference(a: &ManifoldImpl, b: &ManifoldImpl) -> ManifoldImpl {
    let min = a.bbox.min - b.bbox.max;
    let max = a.bbox.max - b.bbox.min;
    let t = mat4_to_mat3x4(translation_matrix(min) * scaling_matrix(max - min));
    ManifoldImpl::cube(&t)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linalg::Mat3x4;

    #[test]
    fn test_minkowski_sum_bbox() {
        let a = ManifoldImpl::cube(&Mat3x4::identity());
        let b = ManifoldImpl::cube(&Mat3x4::identity());
        let sum = minkowski_sum(&a, &b);
        assert!((sum.bbox.max.x - 2.0).abs() < 1e-10);
    }
}
