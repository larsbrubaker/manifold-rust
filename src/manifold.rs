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

use crate::boolean3;
use crate::cross_section::CrossSection;
use crate::constructors;
use crate::math;
use crate::impl_mesh::ManifoldImpl;
use crate::linalg::{mat4_to_mat3x4, normalize, rotation_matrix, rotation_quat_axis_angle, scaling_matrix, translation_matrix, IVec3, Mat3x4, Vec2, Vec3};
use crate::minkowski;
use crate::quickhull;
use crate::sdf;
use crate::subdivision;
use crate::types::{Error, MeshGL, MeshGL64, OpType, Polygons};

#[derive(Clone)]
pub struct Manifold {
    imp: ManifoldImpl,
}

impl Default for Manifold {
    fn default() -> Self {
        Self::new()
    }
}

impl Manifold {
    pub fn new() -> Self {
        Self { imp: ManifoldImpl::new() }
    }

    pub fn empty() -> Self {
        Self::new()
    }

    pub fn from_impl(imp: ManifoldImpl) -> Self {
        Self { imp }
    }

    pub fn as_impl(&self) -> &ManifoldImpl {
        &self.imp
    }

    pub fn tetrahedron() -> Self {
        Self::from_impl(ManifoldImpl::tetrahedron(&Mat3x4::identity()))
    }

    pub fn cube(size: Vec3, center: bool) -> Self {
        let translation = if center {
            translation_matrix(-size * 0.5)
        } else {
            translation_matrix(Vec3::new(0.0, 0.0, 0.0))
        };
        let transform = mat4_to_mat3x4(translation * scaling_matrix(size));
        Self::from_impl(ManifoldImpl::cube(&transform))
    }

    pub fn cylinder(height: f64, radius_low: f64, radius_high: f64, circular_segments: i32) -> Self {
        Self::from_impl(constructors::cylinder(
            height,
            radius_low,
            radius_high,
            circular_segments,
            false,
        ))
    }

    pub fn sphere(radius: f64, circular_segments: i32) -> Self {
        const K_HALF_PI: f64 = std::f64::consts::FRAC_PI_2;

        // Build unit octahedron and subdivide
        // C++: n = (circularSegments + 3) / 4, then Subdivide with n-1 edge cuts.
        // C++ Subdivide splits each edge into n parts in one pass (n^2 tris per original).
        // Our subdivide_impl does recursive midpoint splits (2^levels parts per edge).
        // So we need levels = ceil(log2(n)) to approximate the C++ behavior.
        let identity = mat4_to_mat3x4(scaling_matrix(Vec3::splat(1.0)));
        let base = ManifoldImpl::octahedron(&identity);
        let n = if circular_segments > 0 { (circular_segments + 3) / 4 } else { 1 };
        let levels = if n <= 1 { 0 } else { (n as f64).log2().ceil() as usize };
        let mut mesh = subdivision::subdivide_impl(&base, levels);

        // Map subdivided octahedron vertices onto the sphere surface
        // (matches C++: v = cos(π/2 * (1 - v)); v = radius * normalize(v))
        for v in mesh.vert_pos.iter_mut() {
            let mapped = Vec3::new(
                math::cos(K_HALF_PI * (1.0 - v.x)),
                math::cos(K_HALF_PI * (1.0 - v.y)),
                math::cos(K_HALF_PI * (1.0 - v.z)),
            );
            let n = normalize(mapped);
            *v = if n.x.is_nan() { Vec3::splat(0.0) } else { Vec3::new(n.x * radius, n.y * radius, n.z * radius) };
        }

        // Rebuild mesh metadata after vertex positions changed
        mesh.calculate_bbox();
        mesh.set_epsilon(-1.0, false);
        mesh.sort_geometry();
        mesh.set_normals_and_coplanar();

        Self::from_impl(mesh)
    }

    pub fn extrude(cross_section: &Polygons, height: f64, n_divisions: i32, twist_degrees: f64, scale_top: Vec2) -> Self {
        Self::from_impl(constructors::extrude(cross_section, height, n_divisions, twist_degrees, scale_top))
    }

    pub fn revolve(cross_section: &Polygons, circular_segments: i32, revolve_degrees: f64) -> Self {
        Self::from_impl(constructors::revolve(cross_section, circular_segments, revolve_degrees))
    }

    pub fn from_mesh_gl(mesh: &MeshGL) -> Self {
        let mut imp = ManifoldImpl::new();
        let num_prop = mesh.num_prop as usize;
        imp.num_prop = num_prop.saturating_sub(3);

        imp.vert_pos = (0..mesh.num_vert())
            .map(|i| {
                let p = mesh.get_vert_pos(i);
                Vec3::new(p[0] as f64, p[1] as f64, p[2] as f64)
            })
            .collect();

        if imp.num_prop > 0 {
            imp.properties = (0..mesh.num_vert())
                .flat_map(|i| {
                    let offset = i * num_prop;
                    mesh.vert_properties[offset + 3..offset + num_prop]
                        .iter()
                        .map(|&v| v as f64)
                        .collect::<Vec<_>>()
                })
                .collect();
        }

        let tri: Vec<IVec3> = (0..mesh.num_tri())
            .map(|i| {
                let t = mesh.get_tri_verts(i);
                IVec3::new(t[0] as i32, t[1] as i32, t[2] as i32)
            })
            .collect();
        imp.create_halfedges(&tri, &[]);
        imp.initialize_original();
        imp.calculate_bbox();
        imp.set_epsilon(mesh.tolerance as f64, false);
        imp.sort_geometry();
        imp.set_normals_and_coplanar();
        Self { imp }
    }

    pub fn get_mesh_gl(&self, normal_idx: i32) -> MeshGL {
        let _ = normal_idx;
        let mut out = MeshGL::default();
        let extra_prop = self.imp.num_prop;
        let num_prop = 3 + extra_prop as u32;
        out.num_prop = num_prop;
        let prop_rows = if self.imp.num_prop == 0 { self.imp.num_vert() } else { self.imp.num_prop_vert() };
        let mut prop_pos = vec![Vec3::new(0.0, 0.0, 0.0); prop_rows];
        let mut prop_seen = vec![false; prop_rows];
        for edge in &self.imp.halfedge {
            if edge.prop_vert >= 0 && edge.start_vert >= 0 {
                let p = edge.prop_vert as usize;
                if !prop_seen[p] {
                    prop_seen[p] = true;
                    prop_pos[p] = self.imp.vert_pos[edge.start_vert as usize];
                }
            }
        }

        for row in 0..prop_rows {
            let pos = prop_pos[row];
            out.vert_properties.extend([pos.x as f32, pos.y as f32, pos.z as f32]);
            if self.imp.num_prop > 0 {
                let base = row * self.imp.num_prop;
                for p in 0..self.imp.num_prop {
                    out.vert_properties.push(self.imp.properties[base + p] as f32);
                }
            }
        }

        for tri in 0..self.imp.num_tri() {
            out.tri_verts.extend([
                self.imp.halfedge[3 * tri].prop_vert as u32,
                self.imp.halfedge[3 * tri + 1].prop_vert as u32,
                self.imp.halfedge[3 * tri + 2].prop_vert as u32,
            ]);
        }

        if !self.imp.halfedge_tangent.is_empty() {
            for t in &self.imp.halfedge_tangent {
                out.halfedge_tangent.extend([t.x as f32, t.y as f32, t.z as f32, t.w as f32]);
            }
        }
        out.tolerance = self.imp.tolerance as f32;
        out.run_index = vec![0, out.tri_verts.len() as u32];
        out.run_original_id = vec![self.imp.mesh_relation.original_id.max(0) as u32];
        out
    }

    pub fn get_mesh_gl64(&self, normal_idx: i32) -> MeshGL64 {
        let mesh = self.get_mesh_gl(normal_idx);
        MeshGL64 {
            num_prop: mesh.num_prop as u64,
            vert_properties: mesh.vert_properties.into_iter().map(|v| v as f64).collect(),
            tri_verts: mesh.tri_verts.into_iter().map(|v| v as u64).collect(),
            merge_from_vert: mesh.merge_from_vert.into_iter().map(|v| v as u64).collect(),
            merge_to_vert: mesh.merge_to_vert.into_iter().map(|v| v as u64).collect(),
            run_index: mesh.run_index.into_iter().map(|v| v as u64).collect(),
            run_original_id: mesh.run_original_id,
            run_transform: mesh.run_transform.into_iter().map(|v| v as f64).collect(),
            face_id: mesh.face_id.into_iter().map(|v| v as u64).collect(),
            halfedge_tangent: mesh.halfedge_tangent.into_iter().map(|v| v as f64).collect(),
            tolerance: mesh.tolerance as f64,
        }
    }

    pub fn num_vert(&self) -> usize { self.imp.num_vert() }
    pub fn num_tri(&self) -> usize { self.imp.num_tri() }
    pub fn num_edge(&self) -> usize { self.imp.num_edge() }
    pub fn num_prop(&self) -> usize { self.imp.num_prop }
    pub fn num_prop_vert(&self) -> usize { self.imp.num_prop_vert() }
    pub fn is_empty(&self) -> bool { self.imp.is_empty() }
    pub fn status(&self) -> Error { self.imp.status }
    pub fn volume(&self) -> f64 { self.imp.get_property(crate::properties::Property::Volume).abs() }
    pub fn surface_area(&self) -> f64 { self.imp.get_property(crate::properties::Property::SurfaceArea) }
    pub fn matches_tri_normals(&self) -> bool { self.imp.matches_tri_normals() }
    pub fn num_degenerate_tris(&self) -> i32 { self.imp.num_degenerate_tris() }
    pub fn get_tolerance(&self) -> f64 { self.imp.epsilon }

    /// Port of C++ Manifold::Genus()
    pub fn genus(&self) -> i32 {
        let chi = self.num_vert() as i32 - self.imp.num_edge() as i32 + self.num_tri() as i32;
        1 - chi / 2
    }

    /// Port of C++ Manifold::OriginalID()
    pub fn original_id(&self) -> i32 {
        self.imp.mesh_relation.original_id
    }

    /// Port of C++ Manifold::AsOriginal()
    /// Removes all mesh relations and recreates as an original mesh.
    pub fn as_original(&self) -> Self {
        let mut out = self.imp.clone();
        out.initialize_original();
        out.set_normals_and_coplanar();
        Self::from_impl(out)
    }

    /// Port of C++ Manifold::ReserveIDs()
    pub fn reserve_ids(n: u32) -> u32 {
        crate::impl_mesh::reserve_ids(n)
    }

    /// Port of C++ Manifold::SetTolerance()
    pub fn set_tolerance(&self, tolerance: f64) -> Self {
        let mut out = self.imp.clone();
        if tolerance >= out.epsilon {
            out.epsilon = tolerance;
            // When tolerance increases, simplify the mesh
            crate::edge_op::simplify_topology(&mut out, 0);
            out.sort_geometry();
            out.set_normals_and_coplanar();
        }
        Self::from_impl(out)
    }

    /// Port of C++ Manifold::Simplify()
    pub fn simplify(&self, tolerance: f64) -> Self {
        let mut out = self.imp.clone();
        let old_tolerance = out.epsilon;
        let mut tol = tolerance;
        if tol == 0.0 {
            tol = old_tolerance;
        }
        if tol > old_tolerance {
            out.epsilon = tol;
            out.set_normals_and_coplanar();
        }
        crate::edge_op::simplify_topology(&mut out, 0);
        out.sort_geometry();
        out.epsilon = old_tolerance;
        Self::from_impl(out)
    }

    /// Port of C++ Manifold::WarpBatch()
    pub fn warp_batch<F: Fn(&mut [Vec3])>(&self, warp_fn: F) -> Self {
        if self.is_empty() {
            return self.clone();
        }
        let mut out = self.imp.clone();
        warp_fn(&mut out.vert_pos);
        out.calculate_bbox();
        out.sort_geometry();
        out.set_normals_and_coplanar();
        Self::from_impl(out)
    }

    pub fn translate(&self, v: Vec3) -> Self {
        let t = mat4_to_mat3x4(translation_matrix(v));
        Self::from_impl(self.imp.transform(&t))
    }

    /// Rotate by Euler angles in degrees: first about X, then Y, then Z.
    pub fn rotate(&self, x_degrees: f64, y_degrees: f64, z_degrees: f64) -> Self {
        let to_rad = std::f64::consts::PI / 180.0;
        let qx = rotation_quat_axis_angle(Vec3::new(1.0, 0.0, 0.0), x_degrees * to_rad);
        let qy = rotation_quat_axis_angle(Vec3::new(0.0, 1.0, 0.0), y_degrees * to_rad);
        let qz = rotation_quat_axis_angle(Vec3::new(0.0, 0.0, 1.0), z_degrees * to_rad);
        use crate::linalg::qmul;
        let q = qmul(qz, qmul(qy, qx));
        let t = mat4_to_mat3x4(rotation_matrix(q));
        Self::from_impl(self.imp.transform(&t))
    }

    pub fn scale(&self, v: Vec3) -> Self {
        let t = mat4_to_mat3x4(scaling_matrix(v));
        Self::from_impl(self.imp.transform(&t))
    }

    pub fn transform(&self, m: &Mat3x4) -> Self {
        Self::from_impl(self.imp.transform(m))
    }

    /// Mirror this manifold over the plane defined by the given normal vector.
    /// If the normal is zero-length, returns an empty manifold.
    pub fn mirror(&self, normal: Vec3) -> Self {
        let len = (normal.x * normal.x + normal.y * normal.y + normal.z * normal.z).sqrt();
        if len == 0.0 {
            return Self::empty();
        }
        let n = Vec3::new(normal.x / len, normal.y / len, normal.z / len);
        // Householder reflection: M = I - 2*n*n^T
        let m = Mat3x4::from_cols(
            Vec3::new(1.0 - 2.0 * n.x * n.x, -2.0 * n.x * n.y, -2.0 * n.x * n.z),
            Vec3::new(-2.0 * n.y * n.x, 1.0 - 2.0 * n.y * n.y, -2.0 * n.y * n.z),
            Vec3::new(-2.0 * n.z * n.x, -2.0 * n.z * n.y, 1.0 - 2.0 * n.z * n.z),
            Vec3::new(0.0, 0.0, 0.0),
        );
        Self::from_impl(self.imp.transform(&m))
    }

    /// Warp the mesh by applying a function to each vertex position.
    /// Does not check for self-intersection.
    pub fn warp<F: Fn(&mut Vec3)>(&self, warp_fn: F) -> Self {
        if self.is_empty() {
            return self.clone();
        }
        let mut out = self.imp.clone();
        for v in out.vert_pos.iter_mut() {
            warp_fn(v);
        }
        out.calculate_bbox();
        out.sort_geometry();
        out.set_normals_and_coplanar();
        Self::from_impl(out)
    }

    /// Get the bounding box of this manifold.
    pub fn bounding_box(&self) -> crate::types::Box {
        crate::types::Box {
            min: self.imp.bbox.min,
            max: self.imp.bbox.max,
        }
    }

    /// Split this manifold into two using a cutter manifold.
    /// Returns (intersection, difference).
    pub fn split(&self, cutter: &Self) -> (Self, Self) {
        let intersection = self.intersection(cutter);
        let difference = self.difference(cutter);
        (intersection, difference)
    }

    /// Split this manifold by a plane defined by a normal and offset from origin.
    /// Returns (in direction of normal, opposite direction).
    pub fn split_by_plane(&self, normal: Vec3, origin_offset: f64) -> (Self, Self) {
        if self.is_empty() {
            return (Self::empty(), Self::empty());
        }
        let halfspace = Self::halfspace(&self.imp.bbox, normal, origin_offset);
        self.split(&halfspace)
    }

    /// Trim this manifold by a half-space, keeping only the part in the direction
    /// of the normal vector.
    pub fn trim_by_plane(&self, normal: Vec3, origin_offset: f64) -> Self {
        if self.is_empty() {
            return Self::empty();
        }
        let halfspace = Self::halfspace(&self.imp.bbox, normal, origin_offset);
        self.intersection(&halfspace)
    }

    /// Apply batch boolean operations on a list of manifolds.
    pub fn batch_boolean(manifolds: &[Self], op: OpType) -> Self {
        if manifolds.is_empty() {
            return Self::empty();
        }
        let mut result = manifolds[0].clone();
        for m in &manifolds[1..] {
            result = result.boolean(m, op);
        }
        result
    }

    /// Internal helper: create a halfspace (large cube) for plane splitting.
    fn halfspace(bbox: &crate::types::Box, normal: Vec3, origin_offset: f64) -> Self {
        let n = normalize(normal);
        let cutter = Self::cube(Vec3::splat(2.0), true).translate(Vec3::new(1.0, 0.0, 0.0));
        let center = bbox.center();
        let size_len = (bbox.size().x * bbox.size().x + bbox.size().y * bbox.size().y + bbox.size().z * bbox.size().z).sqrt();
        let dist = ((center.x - n.x * origin_offset).powi(2)
            + (center.y - n.y * origin_offset).powi(2)
            + (center.z - n.z * origin_offset).powi(2)).sqrt()
            + 0.5 * size_len;
        let cutter = cutter.scale(Vec3::splat(dist)).translate(Vec3::new(origin_offset, 0.0, 0.0));
        let y_deg = -math::asin(n.z).to_degrees();
        let z_deg = math::atan2(n.y, n.x).to_degrees();
        cutter.rotate(0.0, y_deg, z_deg)
    }

    pub fn boolean(&self, other: &Self, op: OpType) -> Self {
        Self::from_impl(boolean3::boolean(&self.imp, &other.imp, op))
    }

    pub fn union(&self, other: &Self) -> Self {
        self.boolean(other, OpType::Add)
    }

    pub fn difference(&self, other: &Self) -> Self {
        self.boolean(other, OpType::Subtract)
    }

    pub fn intersection(&self, other: &Self) -> Self {
        self.boolean(other, OpType::Intersect)
    }

    pub fn calculate_curvature(&self, gaussian_idx: i32, mean_idx: i32) -> Self {
        let mut out = self.imp.clone();
        out.calculate_curvature(gaussian_idx, mean_idx);
        Self::from_impl(out)
    }

    /// Set per-vertex properties via a callback function.
    /// `num_prop` is the new number of properties per vertex (>= 0).
    /// `prop_func` receives `(new_prop_slice, position, old_prop_slice)` for each vertex.
    /// The callback writes into `new_prop_slice` (length `num_prop`).
    pub fn set_properties<F>(&self, num_prop: usize, prop_func: F) -> Self
    where
        F: Fn(&mut [f64], Vec3, &[f64]),
    {
        let mut out = self.imp.clone();
        let old_num_prop = out.num_prop;
        let old_properties = out.properties.clone();

        if num_prop == 0 {
            out.properties.clear();
        } else {
            let num_prop_vert = out.num_prop_vert();
            out.properties = vec![0.0; num_prop * num_prop_vert];
            let num_tri = out.num_tri();
            for tri in 0..num_tri {
                for i in 0..3 {
                    let edge = out.halfedge[3 * tri + i];
                    let vert = edge.start_vert as usize;
                    let prop_vert = edge.prop_vert as usize;
                    let pos = out.vert_pos[vert];
                    let old_slice = if old_num_prop > 0 && prop_vert * old_num_prop < old_properties.len() {
                        &old_properties[old_num_prop * prop_vert..old_num_prop * prop_vert + old_num_prop]
                    } else {
                        &[]
                    };
                    prop_func(
                        &mut out.properties[num_prop * prop_vert..num_prop * prop_vert + num_prop],
                        pos,
                        old_slice,
                    );
                }
            }
        }

        out.num_prop = num_prop;
        Self::from_impl(out)
    }

    /// Port of C++ Manifold::CalculateNormals()
    /// Fills in vertex properties for normals. Edges sharper than
    /// min_sharp_angle (degrees) get separate normals on each side.
    pub fn calculate_normals(&self, normal_idx: usize, min_sharp_angle: f64) -> Self {
        let mut out = self.imp.clone();
        out.set_normals(normal_idx as i32, min_sharp_angle);
        Self::from_impl(out)
    }

    pub fn smooth_out(&self, min_sharp_angle: f64, min_smoothness: f64) -> Self {
        let mut out = self.imp.clone();
        let sharpened = out.sharpen_edges(min_sharp_angle, min_smoothness);
        out.create_tangents(sharpened);
        Self::from_impl(out)
    }

    pub fn smooth_by_normals(&self, normal_idx: usize) -> Self {
        let mut out = self.imp.clone();
        out.create_tangents_from_normals(normal_idx);
        Self::from_impl(out)
    }

    /// Port of C++ Manifold::Refine(int n)
    /// Splits every edge into n pieces, sub-triangulating each face.
    pub fn refine(&self, n: i32) -> Self {
        if n <= 1 || self.imp.is_empty() {
            return self.clone();
        }
        let mut out = self.imp.clone();
        let _vert_bary = out.subdivide(&|_vec, _t0, _t1| n - 1, false);
        out.calculate_bbox();
        out.set_epsilon(-1.0, false);
        out.sort_geometry();
        out.set_normals_and_coplanar();
        Self::from_impl(out)
    }

    /// Port of C++ Manifold::RefineToLength(double length)
    pub fn refine_to_length(&self, length: f64) -> Self {
        if length <= 0.0 || self.imp.is_empty() {
            return self.clone();
        }
        let mut out = self.imp.clone();
        let _vert_bary = out.subdivide(
            &|edge_vec, _t0, _t1| {
                let edge_len = (edge_vec.x * edge_vec.x + edge_vec.y * edge_vec.y
                    + edge_vec.z * edge_vec.z)
                    .sqrt();
                ((edge_len / length).ceil() as i32 - 1).max(0)
            },
            true,
        );
        out.calculate_bbox();
        out.set_epsilon(-1.0, false);
        out.sort_geometry();
        out.set_normals_and_coplanar();
        Self::from_impl(out)
    }

    /// Port of C++ Manifold::RefineToTolerance(double tolerance)
    pub fn refine_to_tolerance(&self, tolerance: f64) -> Self {
        if tolerance <= 0.0 || self.imp.is_empty() {
            return self.clone();
        }
        // Similar to C++ — uses tolerance with tangent-based divisions
        let mut out = self.imp.clone();
        let _vert_bary = out.subdivide(
            &|edge_vec, tangent0, tangent1| {
                use crate::linalg::Vec4;
                let edge_len = (edge_vec.x * edge_vec.x + edge_vec.y * edge_vec.y
                    + edge_vec.z * edge_vec.z)
                    .sqrt();
                if edge_len == 0.0 {
                    return 0;
                }
                // Approximate curvature from tangent deviation
                let t0_len = (tangent0.x * tangent0.x + tangent0.y * tangent0.y
                    + tangent0.z * tangent0.z)
                    .sqrt();
                let t1_len = (tangent1.x * tangent1.x + tangent1.y * tangent1.y
                    + tangent1.z * tangent1.z)
                    .sqrt();
                if t0_len == 0.0 && t1_len == 0.0 {
                    return 0; // flat edge
                }
                let max_sag = (t0_len + t1_len) * 0.5;
                ((max_sag / tolerance).sqrt().ceil() as i32 - 1).max(0)
            },
            true,
        );
        out.calculate_bbox();
        out.set_epsilon(-1.0, false);
        out.sort_geometry();
        out.set_normals_and_coplanar();
        Self::from_impl(out)
    }

    pub fn compose(parts: &[Self]) -> Self {
        let impls: Vec<_> = parts.iter().map(|m| m.imp.clone()).collect();
        Self::from_impl(boolean3::compose_meshes(&impls))
    }

    pub fn decompose(&self) -> Vec<Self> {
        vec![self.clone()]
    }

    pub fn level_set<F: Fn(Vec3) -> f64>(sdf_fn: F, bounds: crate::types::Box, edge_length: f64) -> Self {
        Self::from_impl(sdf::level_set(sdf_fn, bounds, edge_length, 0.0, -1.0))
    }

    pub fn hull(points: &[Vec3]) -> Self {
        Self::from_impl(quickhull::convex_hull(points))
    }

    /// Compute the convex hull of this manifold's vertices.
    pub fn convex_hull(&self) -> Self {
        Self::from_impl(quickhull::convex_hull(&self.imp.vert_pos))
    }

    /// Compute the minimum gap between two manifolds within search_length.
    /// Returns search_length if no closer points found within that range.
    pub fn min_gap(&self, other: &Self, search_length: f64) -> f64 {
        self.imp.min_gap(&other.imp, search_length)
    }

    /// Compute the convex hull of multiple manifolds' combined vertices.
    pub fn hull_manifolds(manifolds: &[Self]) -> Self {
        let all_verts: Vec<Vec3> = manifolds
            .iter()
            .flat_map(|m| m.imp.vert_pos.iter().copied())
            .collect();
        Self::from_impl(quickhull::convex_hull(&all_verts))
    }

    pub fn minkowski_sum(&self, other: &Self) -> Self {
        Self::from_impl(minkowski::minkowski_sum(&self.imp, &other.imp))
    }

    pub fn minkowski_difference(&self, other: &Self) -> Self {
        Self::from_impl(minkowski::minkowski_difference(&self.imp, &other.imp))
    }

    pub fn cross_section_square(size: f64) -> CrossSection {
        CrossSection::square(size)
    }

    pub fn cross_section_circle(radius: f64, segments: i32) -> CrossSection {
        CrossSection::circle(radius, segments)
    }
}

// Operator overloads: + for union, - for difference, ^ for intersection
// Matches C++ operator+(Manifold), operator-(Manifold), operator^(Manifold)

impl std::ops::Add for Manifold {
    type Output = Self;
    fn add(self, rhs: Self) -> Self { self.union(&rhs) }
}

impl std::ops::Add<&Manifold> for Manifold {
    type Output = Self;
    fn add(self, rhs: &Self) -> Self { self.union(rhs) }
}

impl std::ops::Add<&Manifold> for &Manifold {
    type Output = Manifold;
    fn add(self, rhs: &Manifold) -> Manifold { self.union(rhs) }
}

impl std::ops::AddAssign for Manifold {
    fn add_assign(&mut self, rhs: Self) { *self = self.union(&rhs); }
}

impl std::ops::AddAssign<&Manifold> for Manifold {
    fn add_assign(&mut self, rhs: &Self) { *self = self.union(rhs); }
}

impl std::ops::Sub for Manifold {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self { self.difference(&rhs) }
}

impl std::ops::Sub<&Manifold> for Manifold {
    type Output = Self;
    fn sub(self, rhs: &Self) -> Self { self.difference(rhs) }
}

impl std::ops::Sub<&Manifold> for &Manifold {
    type Output = Manifold;
    fn sub(self, rhs: &Manifold) -> Manifold { self.difference(rhs) }
}

impl std::ops::SubAssign for Manifold {
    fn sub_assign(&mut self, rhs: Self) { *self = self.difference(&rhs); }
}

impl std::ops::SubAssign<&Manifold> for Manifold {
    fn sub_assign(&mut self, rhs: &Self) { *self = self.difference(rhs); }
}

impl std::ops::BitXor for Manifold {
    type Output = Self;
    fn bitxor(self, rhs: Self) -> Self { self.intersection(&rhs) }
}

impl std::ops::BitXor<&Manifold> for Manifold {
    type Output = Self;
    fn bitxor(self, rhs: &Self) -> Self { self.intersection(rhs) }
}

impl std::ops::BitXor<&Manifold> for &Manifold {
    type Output = Manifold;
    fn bitxor(self, rhs: &Manifold) -> Manifold { self.intersection(rhs) }
}

impl std::ops::BitXorAssign for Manifold {
    fn bitxor_assign(&mut self, rhs: Self) { *self = self.intersection(&rhs); }
}

impl std::ops::BitXorAssign<&Manifold> for Manifold {
    fn bitxor_assign(&mut self, rhs: &Self) { *self = self.intersection(rhs); }
}

#[cfg(test)]
#[path = "manifold_tests/mod.rs"]
mod tests;
