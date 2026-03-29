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

// Phase 15: Subdivision
//
// C++ source: src/subdivision.cpp (811 lines)
//
// STATUS: BASIC MIDPOINT SUBDIVISION ONLY
// Missing: Partition class with cached triangulations, curvature-based
// refinement (refine_to_tolerance), edge-length-based refinement
// (refine_to_length), tangent-based Bezier vertex placement, property
// interpolation, quad subdivision support.
//
// The current midpoint subdivision is geometrically correct for uniform
// refinement but does NOT match the C++ output for smooth/curved surfaces.

use std::collections::HashMap;

use crate::impl_mesh::ManifoldImpl;
use crate::linalg::{IVec3, Vec3};

fn tri_indices(mesh: &ManifoldImpl, tri: usize) -> [i32; 3] {
    [
        mesh.halfedge[3 * tri].start_vert,
        mesh.halfedge[3 * tri + 1].start_vert,
        mesh.halfedge[3 * tri + 2].start_vert,
    ]
}

/// Simple midpoint subdivision: each triangle is split into 4 by inserting
/// edge midpoints. This is a placeholder for the full C++ Subdivide() which
/// uses tangent-based Bezier interpolation for smooth surfaces.
pub fn subdivide_impl(mesh: &ManifoldImpl, levels: usize) -> ManifoldImpl {
    if levels == 0 || mesh.is_empty() {
        return mesh.clone();
    }

    let mut current = mesh.clone();
    for _ in 0..levels {
        let mut verts = current.vert_pos.clone();
        let mut edge_mid: HashMap<(i32, i32), i32> = HashMap::new();
        let mut tris = Vec::new();

        for tri in 0..current.num_tri() {
            let [a, b, c] = tri_indices(&current, tri);
            let ab = midpoint(&mut verts, &mut edge_mid, a, b);
            let bc = midpoint(&mut verts, &mut edge_mid, b, c);
            let ca = midpoint(&mut verts, &mut edge_mid, c, a);

            tris.push(IVec3::new(a, ab, ca));
            tris.push(IVec3::new(ab, b, bc));
            tris.push(IVec3::new(ca, bc, c));
            tris.push(IVec3::new(ab, bc, ca));
        }

        let mut next = ManifoldImpl::new();
        next.vert_pos = verts;
        next.create_halfedges(&tris, &[]);
        next.initialize_original();
        next.calculate_bbox();
        next.set_epsilon(-1.0, false);
        next.sort_geometry();
        next.set_normals_and_coplanar();
        current = next;
    }

    current
}

fn midpoint(
    verts: &mut Vec<Vec3>,
    edge_mid: &mut HashMap<(i32, i32), i32>,
    a: i32,
    b: i32,
) -> i32 {
    let key = if a < b { (a, b) } else { (b, a) };
    if let Some(&idx) = edge_mid.get(&key) {
        return idx;
    }
    let pa = verts[a as usize];
    let pb = verts[b as usize];
    let idx = verts.len() as i32;
    verts.push((pa + pb) * 0.5);
    edge_mid.insert(key, idx);
    idx
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linalg::Mat3x4;

    #[test]
    fn test_subdivide_cube_once() {
        let cube = ManifoldImpl::cube(&Mat3x4::identity());
        let sub = subdivide_impl(&cube, 1);
        assert_eq!(sub.num_tri(), cube.num_tri() * 4);
        assert!(sub.num_vert() > cube.num_vert());
    }
}
