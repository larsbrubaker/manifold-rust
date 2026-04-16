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
use crate::types::{Error, MeshGL, MeshGL64, OpType, Polygons, RayHit, Rect};

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

    /// Create an empty manifold with a specific error status.
    pub fn make_empty(status: crate::types::Error) -> Self {
        let mut imp = ManifoldImpl::new();
        imp.make_empty(status);
        Self { imp }
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
        if size.x < 0.0 || size.y < 0.0 || size.z < 0.0
            || crate::linalg::length(size) == 0.0
        {
            return Self::make_empty(crate::types::Error::InvalidConstruction);
        }
        let translation = if center {
            translation_matrix(-size * 0.5)
        } else {
            translation_matrix(Vec3::new(0.0, 0.0, 0.0))
        };
        let transform = mat4_to_mat3x4(translation * scaling_matrix(size));
        Self::from_impl(ManifoldImpl::cube(&transform))
    }

    pub fn cylinder(height: f64, radius_low: f64, radius_high: f64, circular_segments: i32) -> Self {
        Self::cylinder_centered(height, radius_low, radius_high, circular_segments, false)
    }

    pub fn cylinder_centered(height: f64, radius_low: f64, radius_high: f64, circular_segments: i32, center: bool) -> Self {
        if height <= 0.0 || radius_low < 0.0 {
            return Self::make_empty(crate::types::Error::InvalidConstruction);
        }
        if radius_low == 0.0 && radius_high <= 0.0 {
            return Self::make_empty(crate::types::Error::InvalidConstruction);
        }
        Self::from_impl(constructors::cylinder(
            height,
            radius_low,
            radius_high,
            circular_segments,
            center,
        ))
    }

    pub fn sphere(radius: f64, circular_segments: i32) -> Self {
        if radius <= 0.0 {
            return Self::make_empty(crate::types::Error::InvalidConstruction);
        }
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
        if cross_section.is_empty() || height <= 0.0 {
            return Self::make_empty(crate::types::Error::InvalidConstruction);
        }
        Self::from_impl(constructors::extrude(cross_section, height, n_divisions, twist_degrees, scale_top))
    }

    pub fn revolve(cross_section: &Polygons, circular_segments: i32, revolve_degrees: f64) -> Self {
        if cross_section.is_empty() {
            return Self::make_empty(crate::types::Error::InvalidConstruction);
        }
        Self::from_impl(constructors::revolve(cross_section, circular_segments, revolve_degrees))
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
    pub fn get_tolerance(&self) -> f64 { self.imp.tolerance }

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
        if self.is_empty() { return self.clone(); }
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
        if self.is_empty() { return self.clone(); }
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
        if self.is_empty() { return self.clone(); }
        let mut out = self.imp.clone();
        // C++ uses tolerance_ (not epsilon_) throughout Simplify()
        let old_tolerance = out.tolerance;
        let mut tol = tolerance;
        if tol == 0.0 {
            tol = old_tolerance;
        }
        if tol > old_tolerance {
            out.tolerance = tol;
            out.set_normals_and_coplanar();
        }
        crate::edge_op::simplify_topology(&mut out, 0);
        out.sort_geometry();
        out.tolerance = old_tolerance;
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
        out.mesh_relation.original_id = -1;
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
        out.mesh_relation.original_id = -1;
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

    /// Slice this manifold at the given Z height, returning the cross-section
    /// as a CrossSection. Mirrors C++ `Manifold::Slice`.
    pub fn slice(&self, height: f64) -> CrossSection {
        if self.is_empty() {
            return CrossSection::new(vec![]);
        }
        let polys = self.imp.slice(height);
        CrossSection::new(polys)
    }

    /// Project this manifold onto the XY plane, returning the silhouette
    /// as a CrossSection. Mirrors C++ `Manifold::Project`.
    pub fn project(&self) -> CrossSection {
        if self.is_empty() {
            return CrossSection::new(vec![]);
        }
        let polys = self.imp.project();
        CrossSection::from_polygons_fill(polys)
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
        if self.is_empty() { return self.clone(); }
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
        if self.is_empty() { return self.clone(); }
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
        if self.is_empty() { return self.clone(); }
        let mut out = self.imp.clone();
        out.set_normals(normal_idx as i32, min_sharp_angle);
        Self::from_impl(out)
    }

    /// Port of C++ Manifold::Smooth(MeshGL, sharpenedEdges)
    /// Constructs a smooth version of the input mesh by creating tangents.
    /// The actual triangle resolution is unchanged; use Refine() to
    /// interpolate to a higher-resolution curve.
    pub fn smooth(mesh_gl: &MeshGL, sharpened_edges: &[crate::types::Smoothness]) -> Self {
        use crate::types::Smoothness;

        // Assign sequential faceIDs if not present
        let mut mesh_tmp = mesh_gl.clone();
        let num_tri = mesh_tmp.num_tri();
        mesh_tmp.face_id.resize(num_tri, 0);
        for i in 0..num_tri {
            mesh_tmp.face_id[i] = i as u32;
        }

        let mut m = Self::from_mesh_gl(&mesh_tmp);
        if m.is_empty() {
            return m;
        }

        // UpdateSharpenedEdges + CreateTangents
        let sharpened: Vec<Smoothness> = sharpened_edges.to_vec();
        let updated = m.imp.update_sharpened_edges(&sharpened);
        m.imp.create_tangents(updated);

        // Restore original faceIDs
        let num_tri_impl = m.imp.num_tri();
        for i in 0..num_tri_impl {
            if i < m.imp.mesh_relation.tri_ref.len() {
                let face_id = m.imp.mesh_relation.tri_ref[i].face_id;
                if mesh_gl.face_id.len() == num_tri && face_id >= 0 && (face_id as usize) < num_tri {
                    m.imp.mesh_relation.tri_ref[i].face_id = mesh_gl.face_id[face_id as usize] as i32;
                } else {
                    m.imp.mesh_relation.tri_ref[i].face_id = -1;
                }
            }
        }

        m
    }

    pub fn smooth_out(&self, min_sharp_angle: f64, min_smoothness: f64) -> Self {
        if self.is_empty() { return self.clone(); }
        let mut out = self.imp.clone();
        if min_smoothness == 0.0 {
            // C++ path: use normals-based tangents, then restore original properties
            let saved_num_prop = out.num_prop;
            let saved_properties = out.properties.clone();
            let saved_halfedge = out.halfedge.clone();
            out.set_normals(0, min_sharp_angle);
            out.create_tangents_from_normals(0);
            // Restore original properties (removing temporary normals)
            out.num_prop = saved_num_prop;
            out.properties = saved_properties;
            out.halfedge = saved_halfedge;
        } else {
            let sharpened = out.sharpen_edges(min_sharp_angle, min_smoothness);
            out.create_tangents(sharpened);
        }
        Self::from_impl(out)
    }

    pub fn smooth_by_normals(&self, normal_idx: usize) -> Self {
        if self.is_empty() { return self.clone(); }
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
        if !out.valid_tangents() {
            out.make_empty(crate::types::Error::InvalidTangents);
            return Self::from_impl(out);
        }
        let old = out.clone();
        let had_tangents = out.halfedge_tangent.len() == out.halfedge.len();
        let vert_bary = out.subdivide(&|_vec, _t0, _t1| n - 1, false);
        if had_tangents && !vert_bary.is_empty() {
            crate::interp_tri::interp_tri(&mut out.vert_pos, &vert_bary, &old);
        }
        out.halfedge_tangent.clear();
        out.calculate_bbox();
        out.set_epsilon(-1.0, false);
        out.sort_geometry();
        if had_tangents {
            out.set_normals_and_coplanar();
        } else {
            crate::face_op::calculate_vert_normals(&mut out);
        }
        out.mesh_relation.original_id = -1;
        Self::from_impl(out)
    }

    /// Port of C++ Manifold::RefineToLength(double length)
    pub fn refine_to_length(&self, length: f64) -> Self {
        let length = length.abs();
        if length == 0.0 || self.imp.is_empty() {
            return self.clone();
        }
        let mut out = self.imp.clone();
        if !out.valid_tangents() {
            out.make_empty(crate::types::Error::InvalidTangents);
            return Self::from_impl(out);
        }
        let old = out.clone();
        let had_tangents = out.halfedge_tangent.len() == out.halfedge.len();
        let vert_bary = out.subdivide(
            &|edge_vec, _t0, _t1| {
                let edge_len = (edge_vec.x * edge_vec.x + edge_vec.y * edge_vec.y
                    + edge_vec.z * edge_vec.z)
                    .sqrt();
                // C++: static_cast<int>(la::length(edge) / length) — truncation
                (edge_len / length) as i32
            },
            false,
        );
        if had_tangents && !vert_bary.is_empty() {
            crate::interp_tri::interp_tri(&mut out.vert_pos, &vert_bary, &old);
        }
        out.halfedge_tangent.clear();
        out.calculate_bbox();
        out.set_epsilon(-1.0, false);
        out.sort_geometry();
        if had_tangents {
            out.set_normals_and_coplanar();
        } else {
            crate::face_op::calculate_vert_normals(&mut out);
        }
        out.mesh_relation.original_id = -1;
        Self::from_impl(out)
    }

    /// Port of C++ Manifold::RefineToTolerance(double tolerance)
    pub fn refine_to_tolerance(&self, tolerance: f64) -> Self {
        let tolerance = tolerance.abs();
        if tolerance == 0.0 || self.imp.is_empty() {
            return self.clone();
        }
        // C++ only refines when tangents are present
        let mut out = self.imp.clone();
        let had_tangents = out.halfedge_tangent.len() == out.halfedge.len();
        if !had_tangents {
            return self.clone();
        }
        if !out.valid_tangents() {
            out.make_empty(crate::types::Error::InvalidTangents);
            return Self::from_impl(out);
        }
        let old = out.clone();
        let vert_bary = out.subdivide(
            &|edge_vec, tangent0, tangent1| {
                let edge_len = (edge_vec.x * edge_vec.x + edge_vec.y * edge_vec.y
                    + edge_vec.z * edge_vec.z)
                    .sqrt();
                if edge_len == 0.0 {
                    return 0;
                }
                let edge_norm = Vec3::new(
                    edge_vec.x / edge_len,
                    edge_vec.y / edge_len,
                    edge_vec.z / edge_len,
                );
                let t_start = Vec3::new(tangent0.x, tangent0.y, tangent0.z);
                let t_end = Vec3::new(tangent1.x, tangent1.y, tangent1.z);
                // Perpendicular to edge
                let dot_s = edge_norm.x * t_start.x + edge_norm.y * t_start.y
                    + edge_norm.z * t_start.z;
                let start = Vec3::new(
                    t_start.x - edge_norm.x * dot_s,
                    t_start.y - edge_norm.y * dot_s,
                    t_start.z - edge_norm.z * dot_s,
                );
                let dot_e = edge_norm.x * t_end.x + edge_norm.y * t_end.y
                    + edge_norm.z * t_end.z;
                let end = Vec3::new(
                    t_end.x - edge_norm.x * dot_e,
                    t_end.y - edge_norm.y * dot_e,
                    t_end.z - edge_norm.z * dot_e,
                );
                // Circular arc result plus heuristic term for non-circular curves
                let len_start = (start.x * start.x + start.y * start.y
                    + start.z * start.z)
                    .sqrt();
                let len_end =
                    (end.x * end.x + end.y * end.y + end.z * end.z).sqrt();
                let diff = Vec3::new(
                    start.x - end.x,
                    start.y - end.y,
                    start.z - end.z,
                );
                let len_diff =
                    (diff.x * diff.x + diff.y * diff.y + diff.z * diff.z)
                        .sqrt();
                let d = 0.5 * (len_start + len_end) + len_diff;
                (3.0 * d / (4.0 * tolerance)).sqrt() as i32
            },
            true,
        );
        if had_tangents && !vert_bary.is_empty() {
            crate::interp_tri::interp_tri(&mut out.vert_pos, &vert_bary, &old);
        }
        out.halfedge_tangent.clear();
        out.calculate_bbox();
        out.set_epsilon(-1.0, false);
        out.sort_geometry();
        if had_tangents {
            out.set_normals_and_coplanar();
        } else {
            crate::face_op::calculate_vert_normals(&mut out);
        }
        out.mesh_relation.original_id = -1;
        Self::from_impl(out)
    }

    pub fn compose(parts: &[Self]) -> Self {
        let impls: Vec<_> = parts.iter().map(|m| m.imp.clone()).collect();
        Self::from_impl(boolean3::compose_meshes(&impls))
    }

    pub fn decompose(&self) -> Vec<Self> {
        use crate::disjoint_sets::DisjointSets;

        let num_vert = self.imp.num_vert();
        if num_vert == 0 {
            // Propagate error status: errored manifolds decompose to [self]
            if self.imp.status != Error::NoError {
                return vec![self.clone()];
            }
            return vec![];
        }

        let uf = DisjointSets::new(num_vert as u32);
        for he in &self.imp.halfedge {
            if he.is_forward() {
                uf.unite(he.start_vert as u32, he.end_vert as u32);
            }
        }

        let mut component_indices = vec![0i32; num_vert];
        let num_components = uf.connected_components(&mut component_indices);

        if num_components <= 1 {
            return vec![self.clone()];
        }

        let num_tri = self.imp.num_tri();
        let mut meshes = Vec::new();

        for comp in 0..num_components {
            let mut imp = ManifoldImpl::new();
            imp.tolerance = self.imp.tolerance;

            // Collect vertices belonging to this component
            let vert_new2old: Vec<i32> = (0..num_vert as i32)
                .filter(|&v| component_indices[v as usize] == comp)
                .collect();
            let n_vert = vert_new2old.len();
            if n_vert == 0 { continue; }

            imp.vert_pos = vert_new2old.iter().map(|&v| self.imp.vert_pos[v as usize]).collect();
            if !self.imp.vert_normal.is_empty() {
                imp.vert_normal = vert_new2old.iter()
                    .map(|&v| self.imp.vert_normal[v as usize]).collect();
            }

            // Collect faces belonging to this component
            let face_new2old: Vec<usize> = (0..num_tri)
                .filter(|&f| {
                    let sv = self.imp.halfedge[3 * f].start_vert;
                    sv >= 0 && component_indices[sv as usize] == comp
                })
                .collect();

            if face_new2old.is_empty() { continue; }

            // Copy full data from original, then gather_faces will filter
            imp.halfedge = self.imp.halfedge.clone();
            imp.face_normal = self.imp.face_normal.clone();
            imp.halfedge_tangent = self.imp.halfedge_tangent.clone();
            imp.num_prop = self.imp.num_prop;
            imp.properties = self.imp.properties.clone();
            imp.mesh_relation = self.imp.mesh_relation.clone();

            crate::sort::gather_faces(&mut imp, &face_new2old);
            crate::sort::reindex_verts(&mut imp, &vert_new2old, self.imp.num_vert());
            imp.calculate_bbox();
            imp.sort_geometry();

            meshes.push(Self::from_impl(imp));
        }

        meshes
    }

    pub fn level_set<F: Fn(Vec3) -> f64>(sdf_fn: F, bounds: crate::types::Box, edge_length: f64) -> Self {
        Self::from_impl(sdf::level_set(sdf_fn, bounds, edge_length, 0.0, -1.0))
    }

    pub fn level_set_with_level<F: Fn(Vec3) -> f64>(sdf_fn: F, bounds: crate::types::Box, edge_length: f64, level: f64) -> Self {
        Self::from_impl(sdf::level_set(sdf_fn, bounds, edge_length, level, -1.0))
    }

    pub fn hull(points: &[Vec3]) -> Self {
        Self::from_impl(quickhull::convex_hull(points))
    }

    /// Compute the convex hull of this manifold's vertices.
    pub fn convex_hull(&self) -> Self {
        if self.is_empty() { return self.clone(); }
        Self::from_impl(quickhull::convex_hull(&self.imp.vert_pos))
    }

    /// Compute the minimum gap between two manifolds within search_length.
    /// Returns search_length if no closer points found within that range.
    pub fn min_gap(&self, other: &Self, search_length: f64) -> f64 {
        self.imp.min_gap(&other.imp, search_length)
    }

    /// Cast a ray segment from `origin` to `endpoint`, returning all triangle
    /// intersections sorted by parametric distance along the segment.
    /// Mirrors C++ `Manifold::RayCast(vec3, vec3)`.
    pub fn ray_cast(&self, origin: Vec3, endpoint: Vec3) -> Vec<RayHit> {
        crate::boolean3::ray_cast(&self.imp, origin, endpoint)
    }

    /// Compute the convex hull of multiple manifolds' combined vertices.
    /// If any manifold is errored, propagates its error status.
    pub fn hull_manifolds(manifolds: &[Self]) -> Self {
        // Propagate any error from inputs
        for m in manifolds {
            if m.imp.status != Error::NoError {
                return m.clone();
            }
        }
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

#[path = "manifold_meshgl.rs"]
mod meshgl;

#[cfg(test)]
#[path = "manifold_tests/mod.rs"]
mod tests;
