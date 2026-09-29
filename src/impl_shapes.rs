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

// impl_shapes.rs — built-in polyhedron constructors for ManifoldImpl
//
// Ports the Impl(Shape, mat3x4) constructor from src/impl.cpp: the
// tetrahedron, cube and octahedron vertex/triangle tables and the shared
// build pipeline (create_halfedges -> initialize_original -> bbox/epsilon ->
// sort_geometry -> set_normals_and_coplanar). Split out of impl_mesh.rs,
// which defines the ManifoldImpl struct and the pipeline stages; this file
// adds a second `impl ManifoldImpl` block. constructors.rs builds the other
// primitives on top of these.

use crate::impl_mesh::ManifoldImpl;
use crate::linalg::{IVec3, Mat3x4, Vec3};

impl ManifoldImpl {
    // -----------------------------------------------------------------------
    // Shape constructors
    // -----------------------------------------------------------------------

    pub fn tetrahedron(transform: &Mat3x4) -> Self {
        let vert_pos_raw: Vec<[f64; 3]> = vec![
            [-1.0, -1.0, 1.0],
            [-1.0, 1.0, -1.0],
            [1.0, -1.0, -1.0],
            [1.0, 1.0, 1.0],
        ];
        let tri_verts: Vec<IVec3> = vec![
            IVec3::new(2, 0, 1),
            IVec3::new(0, 3, 1),
            IVec3::new(2, 3, 0),
            IVec3::new(3, 2, 1),
        ];
        Self::from_shape(vert_pos_raw, tri_verts, transform)
    }

    pub fn cube(transform: &Mat3x4) -> Self {
        let vert_pos_raw: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 1.0, 0.0],
            [0.0, 1.0, 1.0],
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 0.0],
            [1.0, 1.0, 1.0],
        ];
        let tri_verts: Vec<IVec3> = vec![
            IVec3::new(1, 0, 4),
            IVec3::new(2, 4, 0),
            IVec3::new(1, 3, 0),
            IVec3::new(3, 1, 5),
            IVec3::new(3, 2, 0),
            IVec3::new(3, 7, 2),
            IVec3::new(5, 4, 6),
            IVec3::new(5, 1, 4),
            IVec3::new(6, 4, 2),
            IVec3::new(7, 6, 2),
            IVec3::new(7, 3, 5),
            IVec3::new(7, 5, 6),
        ];
        Self::from_shape(vert_pos_raw, tri_verts, transform)
    }

    pub fn octahedron(transform: &Mat3x4) -> Self {
        let vert_pos_raw: Vec<[f64; 3]> = vec![
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
        ];
        let tri_verts: Vec<IVec3> = vec![
            IVec3::new(0, 2, 4),
            IVec3::new(1, 5, 3),
            IVec3::new(2, 1, 4),
            IVec3::new(3, 5, 0),
            IVec3::new(1, 3, 4),
            IVec3::new(0, 5, 2),
            IVec3::new(3, 0, 4),
            IVec3::new(2, 5, 1),
        ];
        Self::from_shape(vert_pos_raw, tri_verts, transform)
    }

    fn from_shape(vert_pos_raw: Vec<[f64; 3]>, tri_verts: Vec<IVec3>, transform: &Mat3x4) -> Self {
        use crate::linalg::Vec4;
        let mut m = Self::new();
        m.vert_pos = vert_pos_raw
            .iter()
            .map(|v| {
                let p = Vec3::new(v[0], v[1], v[2]);
                // Apply transform: m * vec4(p, 1)
                *transform * Vec4::new(p.x, p.y, p.z, 1.0)
            })
            .collect();

        m.create_halfedges(&tri_verts, &[]);
        m.initialize_original();
        m.calculate_bbox();
        m.set_epsilon(-1.0, false);
        m.sort_geometry();
        m.set_normals_and_coplanar();
        m
    }
}
