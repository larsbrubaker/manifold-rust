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

// Phase 16: SDF Mesh Generation
//
// C++ source: src/sdf.cpp (538 lines)
//
// STATUS: NOT YET IMPLEMENTED
// The C++ uses Marching Tetrahedra on a body-centered cubic (BCC) grid with
// ITP root-finding for precise surface location. This placeholder generates
// axis-aligned cubes at grid cells where the SDF is positive — it does NOT
// produce correct isosurface geometry.

use crate::boolean3::compose_meshes;
use crate::impl_mesh::ManifoldImpl;
use crate::linalg::{mat4_to_mat3x4, scaling_matrix, translation_matrix, Vec3};
use crate::types::Box as BBox;

/// Generate a mesh from a signed distance function.
///
/// PLACEHOLDER: Currently generates axis-aligned cubes where SDF > 0.
/// Will be replaced with proper Marching Tetrahedra on BCC grid.
pub fn level_set<F: Fn(Vec3) -> f64>(sdf: F, bounds: BBox, edge_length: f64) -> ManifoldImpl {
    if edge_length <= 0.0 {
        return ManifoldImpl::new();
    }

    let size = bounds.size();
    let nx = (size.x / edge_length).ceil().max(1.0) as i32;
    let ny = (size.y / edge_length).ceil().max(1.0) as i32;
    let nz = (size.z / edge_length).ceil().max(1.0) as i32;

    let mut cells = Vec::new();
    for ix in 0..nx {
        for iy in 0..ny {
            for iz in 0..nz {
                let min = Vec3::new(
                    bounds.min.x + ix as f64 * edge_length,
                    bounds.min.y + iy as f64 * edge_length,
                    bounds.min.z + iz as f64 * edge_length,
                );
                let center = min + Vec3::splat(edge_length * 0.5);
                if sdf(center) > 0.0 {
                    let t = mat4_to_mat3x4(translation_matrix(min) * scaling_matrix(Vec3::splat(edge_length)));
                    cells.push(ManifoldImpl::cube(&t));
                }
            }
        }
    }

    compose_meshes(&cells)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_level_set_single_positive_cell() {
        let bounds = BBox::from_points(Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 1.0, 1.0));
        let mesh = level_set(|_p| 1.0, bounds, 1.0);
        assert_eq!(mesh.num_tri(), 12);
    }
}
