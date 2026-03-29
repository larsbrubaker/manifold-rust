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

// Phase 11: Boolean Operations (Core)
//
// C++ sources: src/boolean3.cpp (531 lines), src/boolean_result.cpp (889 lines)
//
// STATUS: NOT YET IMPLEMENTED
// The compose_meshes() helper is genuine — it concatenates disjoint meshes.
// The boolean() function is NOT yet implemented — it requires the full
// edge-face intersection algorithm from boolean3.cpp.

use crate::impl_mesh::ManifoldImpl;
use crate::linalg::IVec3;
use crate::types::OpType;

fn extract_tri_vert(mesh: &ManifoldImpl) -> Vec<IVec3> {
    (0..mesh.num_tri())
        .map(|tri| {
            IVec3::new(
                mesh.halfedge[3 * tri].start_vert,
                mesh.halfedge[3 * tri + 1].start_vert,
                mesh.halfedge[3 * tri + 2].start_vert,
            )
        })
        .collect()
}

fn extract_tri_prop(mesh: &ManifoldImpl) -> Vec<IVec3> {
    (0..mesh.num_tri())
        .map(|tri| {
            IVec3::new(
                mesh.halfedge[3 * tri].prop_vert,
                mesh.halfedge[3 * tri + 1].prop_vert,
                mesh.halfedge[3 * tri + 2].prop_vert,
            )
        })
        .collect()
}

fn property_row(mesh: &ManifoldImpl, row: usize, width: usize) -> Vec<f64> {
    if mesh.num_prop == 0 {
        vec![0.0; width]
    } else {
        let mut out = vec![0.0; width];
        let src = &mesh.properties[row * mesh.num_prop..(row + 1) * mesh.num_prop];
        out[..src.len()].copy_from_slice(src);
        out
    }
}

/// Concatenate multiple disjoint meshes into one. This is a genuine utility
/// used by both boolean operations and CSG compose. It does NOT perform any
/// boolean intersection — the meshes must be non-overlapping for correct results.
pub fn compose_meshes(meshes: &[ManifoldImpl]) -> ManifoldImpl {
    if meshes.is_empty() {
        return ManifoldImpl::new();
    }
    if meshes.len() == 1 {
        return meshes[0].clone();
    }

    let num_prop = meshes.iter().map(|m| m.num_prop).max().unwrap_or(0);
    let mut vert_pos = Vec::new();
    let mut properties = Vec::new();
    let mut tri_vert = Vec::new();
    let mut tri_prop = Vec::new();
    let mut vert_offset = 0i32;
    let mut prop_offset = 0i32;

    for mesh in meshes {
        vert_pos.extend_from_slice(&mesh.vert_pos);

        let old_tri_vert = extract_tri_vert(mesh);
        let old_tri_prop = extract_tri_prop(mesh);
        tri_vert.extend(old_tri_vert.into_iter().map(|t| {
            IVec3::new(t.x + vert_offset, t.y + vert_offset, t.z + vert_offset)
        }));
        tri_prop.extend(old_tri_prop.into_iter().map(|t| {
            IVec3::new(t.x + prop_offset, t.y + prop_offset, t.z + prop_offset)
        }));

        if num_prop > 0 {
            let prop_rows = mesh.num_prop_vert();
            for row in 0..prop_rows {
                properties.extend(property_row(mesh, row, num_prop));
            }
            prop_offset += prop_rows as i32;
        } else {
            prop_offset += mesh.num_prop_vert() as i32;
        }
        vert_offset += mesh.num_vert() as i32;
    }

    let mut out = ManifoldImpl::new();
    out.vert_pos = vert_pos;
    out.num_prop = num_prop;
    out.properties = properties;
    out.create_halfedges(&tri_prop, &tri_vert);
    out.initialize_original();
    out.calculate_bbox();
    out.set_epsilon(-1.0, false);
    out.sort_geometry();
    out.set_normals_and_coplanar();
    out
}

/// Perform a 3D boolean operation on two manifold meshes.
///
/// NOT YET IMPLEMENTED — requires full port of boolean3.cpp edge-face
/// intersection detection and boolean_result.cpp face assembly.
///
/// Currently only handles the trivial non-overlapping case for Add.
/// All other cases panic to prevent silent incorrect results.
pub fn boolean(mesh_a: &ManifoldImpl, mesh_b: &ManifoldImpl, op: OpType) -> ManifoldImpl {
    if mesh_a.is_empty() {
        return match op {
            OpType::Add | OpType::Intersect => mesh_b.clone(),
            OpType::Subtract => ManifoldImpl::new(),
        };
    }
    if mesh_b.is_empty() {
        return match op {
            OpType::Add | OpType::Subtract => mesh_a.clone(),
            OpType::Intersect => ManifoldImpl::new(),
        };
    }

    match op {
        OpType::Add => {
            if !mesh_a.bbox.does_overlap_box(&mesh_b.bbox) {
                // Non-overlapping union is just composition
                compose_meshes(&[mesh_a.clone(), mesh_b.clone()])
            } else {
                unimplemented!(
                    "Boolean Add with overlapping meshes requires full boolean3.cpp port"
                );
            }
        }
        OpType::Intersect => {
            if !mesh_a.bbox.does_overlap_box(&mesh_b.bbox) {
                ManifoldImpl::new()
            } else {
                unimplemented!(
                    "Boolean Intersect requires full boolean3.cpp port"
                );
            }
        }
        OpType::Subtract => {
            if !mesh_a.bbox.does_overlap_box(&mesh_b.bbox) {
                mesh_a.clone()
            } else {
                unimplemented!(
                    "Boolean Subtract requires full boolean3.cpp port"
                );
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linalg::{mat4_to_mat3x4, translation_matrix, Vec3};

    #[test]
    fn test_compose_meshes_disjoint_cubes() {
        let a = ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(Vec3::new(0.0, 0.0, 0.0))));
        let b = ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(Vec3::new(3.0, 0.0, 0.0))));
        let c = compose_meshes(&[a, b]);
        assert_eq!(c.num_tri(), 24);
        assert_eq!(c.num_vert(), 16);
    }

    #[test]
    fn test_boolean_add_disjoint() {
        let a = ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(Vec3::new(0.0, 0.0, 0.0))));
        let b = ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(Vec3::new(3.0, 0.0, 0.0))));
        let c = boolean(&a, &b, OpType::Add);
        assert_eq!(c.num_tri(), 24);
        assert_eq!(c.num_vert(), 16);
    }

    #[test]
    fn test_boolean_intersect_disjoint_empty() {
        let a = ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(Vec3::new(0.0, 0.0, 0.0))));
        let b = ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(Vec3::new(3.0, 0.0, 0.0))));
        let c = boolean(&a, &b, OpType::Intersect);
        assert!(c.is_empty());
    }

    #[test]
    fn test_boolean_subtract_disjoint() {
        let a = ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(Vec3::new(0.0, 0.0, 0.0))));
        let b = ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(Vec3::new(3.0, 0.0, 0.0))));
        let c = boolean(&a, &b, OpType::Subtract);
        assert_eq!(c.num_tri(), 12);
    }
}
