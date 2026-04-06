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
        // C++: n = (circularSegments + 3) / 4, then subdivide n-1 times
        let identity = mat4_to_mat3x4(scaling_matrix(Vec3::splat(1.0)));
        let base = ManifoldImpl::octahedron(&identity);
        let n = if circular_segments > 0 { (circular_segments + 3) / 4 } else { 1 };
        let levels = (n - 1).max(0) as usize;
        let mut mesh = subdivision::subdivide_impl(&base, levels);

        // Map subdivided octahedron vertices onto the sphere surface
        // (matches C++: v = cos(π/2 * (1 - v)); v = radius * normalize(v))
        for v in mesh.vert_pos.iter_mut() {
            let mapped = Vec3::new(
                (K_HALF_PI * (1.0 - v.x)).cos(),
                (K_HALF_PI * (1.0 - v.y)).cos(),
                (K_HALF_PI * (1.0 - v.z)).cos(),
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
    pub fn is_empty(&self) -> bool { self.imp.is_empty() }
    pub fn status(&self) -> Error { self.imp.status }
    pub fn volume(&self) -> f64 { self.imp.get_property(crate::properties::Property::Volume).abs() }
    pub fn surface_area(&self) -> f64 { self.imp.get_property(crate::properties::Property::SurfaceArea) }
    pub fn matches_tri_normals(&self) -> bool { self.imp.matches_tri_normals() }
    pub fn num_degenerate_tris(&self) -> i32 { self.imp.num_degenerate_tris() }
    pub fn genus(&self) -> i32 {
        let chi = self.num_vert() as i32 - self.imp.num_edge() as i32 + self.num_tri() as i32;
        1 - chi / 2
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

    pub fn calculate_normals(&self, normal_idx: usize) -> Self {
        let mut out = self.imp.clone();
        let prop_rows = out.num_prop_vert();
        if out.num_prop < normal_idx + 3 {
            let mut new_props = vec![0.0; prop_rows * (normal_idx + 3)];
            for row in 0..prop_rows {
                for p in 0..out.num_prop {
                    new_props[row * (normal_idx + 3) + p] = out.properties[row * out.num_prop + p];
                }
            }
            out.properties = new_props;
            out.num_prop = normal_idx + 3;
        }

        let mut prop_seen = vec![false; prop_rows];
        for edge in &out.halfedge {
            let prop = edge.prop_vert as usize;
            if !prop_seen[prop] {
                prop_seen[prop] = true;
                let normal = out.vert_normal[edge.start_vert as usize];
                let base = prop * out.num_prop + normal_idx;
                out.properties[base] = normal.x;
                out.properties[base + 1] = normal.y;
                out.properties[base + 2] = normal.z;
            }
        }
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

    pub fn refine(&self, n: i32) -> Self {
        if n <= 0 {
            return self.clone();
        }
        Self::from_impl(subdivision::subdivide_impl(&self.imp, n as usize))
    }

    pub fn refine_to_length(&self, length: f64) -> Self {
        if length <= 0.0 || self.imp.is_empty() {
            return self.clone();
        }
        let scale = self.imp.bbox.size().x.max(self.imp.bbox.size().y).max(self.imp.bbox.size().z);
        let levels = ((scale / length).max(1.0).log2().ceil() as usize).min(5);
        Self::from_impl(subdivision::subdivide_impl(&self.imp, levels))
    }

    pub fn refine_to_tolerance(&self, tolerance: f64) -> Self {
        self.refine_to_length(tolerance * 2.0)
    }

    pub fn compose(parts: &[Self]) -> Self {
        let impls: Vec<_> = parts.iter().map(|m| m.imp.clone()).collect();
        Self::from_impl(boolean3::compose_meshes(&impls))
    }

    pub fn decompose(&self) -> Vec<Self> {
        vec![self.clone()]
    }

    pub fn level_set<F: Fn(Vec3) -> f64>(sdf_fn: F, bounds: crate::types::Box, edge_length: f64) -> Self {
        Self::from_impl(sdf::level_set(sdf_fn, bounds, edge_length))
    }

    pub fn hull(points: &[Vec3]) -> Self {
        Self::from_impl(quickhull::convex_hull(points))
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_manifold_cube_counts() {
        let m = Manifold::cube(Vec3::new(1.0, 1.0, 1.0), false);
        assert_eq!(m.num_vert(), 8);
        assert_eq!(m.num_tri(), 12);
    }

    #[test]
    fn test_manifold_transform_translate() {
        let m = Manifold::cube(Vec3::new(1.0, 1.0, 1.0), false).translate(Vec3::new(2.0, 0.0, 0.0));
        let out = m.get_mesh_gl(0);
        let p = out.get_vert_pos(0);
        assert!(p[0] >= 2.0);
    }

    #[test]
    fn test_manifold_union_disjoint() {
        let a = Manifold::cube(Vec3::new(1.0, 1.0, 1.0), false);
        let b = Manifold::cube(Vec3::new(1.0, 1.0, 1.0), false).translate(Vec3::new(3.0, 0.0, 0.0));
        let c = a.union(&b);
        assert_eq!(c.num_tri(), 24);
    }

    #[test]
    fn test_mesh_gl_roundtrip_basic() {
        let m = Manifold::cube(Vec3::new(1.0, 1.0, 1.0), false);
        let mesh = m.get_mesh_gl(0);
        let rebuilt = Manifold::from_mesh_gl(&mesh);
        assert_eq!(rebuilt.num_tri(), m.num_tri());
    }

    #[test]
    fn test_calculate_curvature_keeps_mesh() {
        let m = Manifold::cube(Vec3::new(1.0, 1.0, 1.0), false).calculate_curvature(0, 1);
        let mesh = m.get_mesh_gl(0);
        assert!(mesh.num_prop >= 5);
    }

    #[test]
    fn test_refine_increases_triangles() {
        let m = Manifold::cube(Vec3::new(1.0, 1.0, 1.0), false);
        let r = m.refine(1);
        assert_eq!(r.num_tri(), m.num_tri() * 4);
    }

    #[test]
    fn test_hull_tetrahedron() {
        let hull = Manifold::hull(&[
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ]);
        assert_eq!(hull.num_vert(), 4);
        assert_eq!(hull.num_tri(), 4);
    }

    #[test]
    fn test_hull_cube() {
        let hull = Manifold::hull(&[
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(1.0, 1.0, 0.0),
            Vec3::new(1.0, 0.0, 1.0),
            Vec3::new(0.0, 1.0, 1.0),
            Vec3::new(1.0, 1.0, 1.0),
        ]);
        assert_eq!(hull.num_vert(), 8);
        assert_eq!(hull.num_tri(), 12);
    }

    // -----------------------------------------------------------------------
    // C++ parity tests — ported from cpp-reference/manifold/test/manifold_test.cpp
    // -----------------------------------------------------------------------

    /// C++ TEST(Manifold, Sphere) — n=25, 4*n segments → n*n*8 triangles
    #[test]
    fn test_cpp_sphere_tri_count() {
        let n = 25;
        let sphere = Manifold::sphere(1.0, 4 * n);
        assert_eq!(sphere.num_tri(), (n * n * 8) as usize);
    }

    /// C++ TEST(Manifold, Cylinder) — 10000 segments
    #[test]
    fn test_cpp_cylinder_tri_count() {
        let n = 10000i32;
        let cyl = Manifold::cylinder(2.0, 2.0, 2.0, n);
        assert_eq!(cyl.num_tri(), (4 * n - 4) as usize);
    }

    /// C++ TEST(Manifold, Revolve3) — revolve a circle to make a sphere
    #[test]
    fn test_cpp_revolve3() {
        let circle = crate::cross_section::CrossSection::circle(1.0, 32);
        let sphere = Manifold::revolve(&circle.to_polygons(), 32, 360.0);
        let k_pi = std::f64::consts::PI;
        assert!((sphere.volume() - 4.0 / 3.0 * k_pi).abs() < 0.1);
        assert!((sphere.surface_area() - 4.0 * k_pi).abs() < 0.15);
    }

    /// C++ TEST(Manifold, Transform)
    #[test]
    fn test_cpp_transform() {
        let cube = Manifold::cube(Vec3::splat(1.0), false);
        let translated = cube.translate(Vec3::new(1.0, 2.0, 3.0));
        assert_eq!(translated.num_vert(), 8);
        assert_eq!(translated.num_tri(), 12);
        assert!((translated.volume() - 1.0).abs() < 1e-10);
    }

    /// C++ TEST(Manifold, MirrorUnion) — cube union with its mirror
    #[test]
    fn test_cpp_mirror_union() {
        let cube = Manifold::cube(Vec3::new(5.0, 5.0, 5.0), false)
            .translate(Vec3::new(0.0, 0.0, -3.0));
        let mirrored = cube.scale(Vec3::new(1.0, 1.0, -1.0));
        let result = cube.union(&mirrored);
        assert_eq!(result.genus(), 0);
        assert!((result.volume() - 5.0 * 5.0 * 6.0).abs() < 1e-5);
    }

    #[test]
    fn test_sphere_is_round() {
        let m = Manifold::sphere(1.0, 24);
        // Unit sphere: volume should be ~4π/3 ≈ 4.189
        let vol = m.volume();
        let expected = 4.0 * std::f64::consts::PI / 3.0;
        assert!(
            (vol - expected).abs() < 0.15,
            "Sphere volume should be ~{:.3}, got {:.3}",
            expected,
            vol
        );
        // All vertices should be approximately at radius 1.0
        let mesh = m.get_mesh_gl(0);
        let num_prop = mesh.num_prop as usize;
        let vert_count = if num_prop > 0 { mesh.vert_properties.len() / num_prop } else { 0 };
        for i in 0..vert_count {
            let x = mesh.vert_properties[i * num_prop] as f64;
            let y = mesh.vert_properties[i * num_prop + 1] as f64;
            let z = mesh.vert_properties[i * num_prop + 2] as f64;
            let r = (x * x + y * y + z * z).sqrt();
            assert!(
                (r - 1.0).abs() < 0.01,
                "Vertex {} at ({:.3},{:.3},{:.3}) has radius {:.4}, expected ~1.0",
                i, x, y, z, r
            );
        }
    }
}
