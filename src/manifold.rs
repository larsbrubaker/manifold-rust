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

// manifold.rs — the public `Manifold` handle: a thin, immutable wrapper over
// `ManifoldImpl` (impl_mesh.rs) that ports the C++ `Manifold` API.
//
// This file owns the struct, its constructors and accessors, the geometric
// queries, simplification, affine transforms and warps, property editing,
// compose / decompose, and slice / project. The rest of the API lives in
// child modules so they can reach the private `imp` field:
// manifold_boolean.rs (booleans, plane splits, operator overloads),
// manifold_robust.rs (robust-engine repair entry points),
// manifold_meshgl.rs (MeshGL import/export entry points and the import
// body), manifold_meshgl_export.rs (the export body), manifold_shape.rs
// (primitive and extrusion constructors) and manifold_smooth.rs (normals,
// smoothing, refinement).

use crate::boolean3;
use crate::cross_section::CrossSection;
use crate::impl_mesh::ManifoldImpl;
use crate::linalg::{mat4_to_mat3x4, scaling_matrix, translation_matrix, Mat3, Mat3x4, Vec3};
use crate::types::{Error, RayHit};

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
        Self {
            imp: ManifoldImpl::new(),
        }
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

    pub fn into_impl(self) -> ManifoldImpl {
        self.imp
    }

    pub fn num_vert(&self) -> usize {
        self.imp.num_vert()
    }
    pub fn num_tri(&self) -> usize {
        self.imp.num_tri()
    }
    pub fn num_edge(&self) -> usize {
        self.imp.num_edge()
    }
    pub fn num_prop(&self) -> usize {
        self.imp.num_prop
    }
    pub fn num_prop_vert(&self) -> usize {
        self.imp.num_prop_vert()
    }
    pub fn is_empty(&self) -> bool {
        self.imp.is_empty()
    }
    pub fn status(&self) -> Error {
        self.imp.status
    }
    pub fn volume(&self) -> f64 {
        self.imp
            .get_property(crate::properties::Property::Volume)
            .abs()
    }
    pub fn surface_area(&self) -> f64 {
        self.imp
            .get_property(crate::properties::Property::SurfaceArea)
    }
    pub fn matches_tri_normals(&self) -> bool {
        self.imp.matches_tri_normals()
    }
    pub fn num_degenerate_tris(&self) -> i32 {
        self.imp.num_degenerate_tris()
    }
    pub fn get_tolerance(&self) -> f64 {
        self.imp.tolerance
    }
    pub fn get_epsilon(&self) -> f64 {
        self.imp.epsilon
    }

    /// Port of C++ Manifold::Genus()
    pub fn genus(&self) -> i32 {
        let chi = self.num_vert() as i32 - self.imp.num_edge() as i32 + self.num_tri() as i32;
        1 - chi / 2
    }

    /// Port of C++ Manifold::OriginalID()
    pub fn original_id(&self) -> i32 {
        self.imp.mesh_relation.original_id
    }

    /// Soup impls (non-manifold geometry imported via
    /// [`Manifold::from_mesh_gl_robust`]) support only robust-engine
    /// booleans, transforms, bounding-box queries, hulls, and mesh export.
    /// Pairing-dependent operations call this and return an empty manifold
    /// with [`Error::NotManifold`] instead of walking incomplete halfedges.
    pub(crate) fn require_paired(&self) -> Option<Self> {
        if self.imp.is_soup {
            Some(Self::make_empty(Error::NotManifold))
        } else {
            None
        }
    }

    /// Port of C++ Manifold::AsOriginal()
    /// Removes all mesh relations (original IDs, face IDs, transforms) and
    /// recreates this as an original mesh with a new ID; exported face IDs
    /// then fall back to coplanar groups. To keep caller-assigned face IDs,
    /// import a MeshGL with `face_id` set instead (see
    /// [`Self::from_mesh_gl`]).
    pub fn as_original(&self) -> Self {
        if let Some(e) = self.require_paired() {
            return e;
        }
        if self.is_empty() {
            return self.clone();
        }
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
        if let Some(e) = self.require_paired() {
            return e;
        }
        if self.is_empty() {
            return self.clone();
        }
        let mut out = self.imp.clone();
        // Matches C++ SetTolerance: operate on the `tolerance` field (which
        // drives coplanar grouping in mark_coplanar), not `epsilon`. When
        // raising tolerance, recompute coplanar groups *first* so the
        // colinear-edge collapse in simplify_topology can see the newly
        // co-planar regions, then simplify, then re-sort.
        if tolerance > out.tolerance {
            out.tolerance = tolerance;
            out.set_normals_and_coplanar();
            crate::edge_op::simplify_topology(&mut out, 0);
            out.sort_geometry();
            // Collapsed edges move geometry; the cloned verdict is stale.
            out.invalidate_self_intersects();
        } else {
            // Reducing tolerance: keep it at least epsilon.
            out.tolerance = out.epsilon.max(tolerance);
        }
        Self::from_impl(out)
    }

    /// Port of C++ Manifold::Simplify()
    pub fn simplify(&self, tolerance: f64) -> Self {
        if let Some(e) = self.require_paired() {
            return e;
        }
        if self.is_empty() {
            return self.clone();
        }
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
        // Collapsed edges move geometry; the cloned verdict is stale.
        out.invalidate_self_intersects();
        out.tolerance = old_tolerance;
        Self::from_impl(out)
    }

    /// Port of C++ Manifold::WarpBatch()
    pub fn warp_batch<F: Fn(&mut [Vec3])>(&self, warp_fn: F) -> Self {
        if let Some(e) = self.require_paired() {
            return e;
        }
        if self.is_empty() {
            return self.clone();
        }
        let mut out = self.imp.clone();
        warp_fn(&mut out.vert_pos);
        // Arbitrary deformation: nothing about the old verdict survives.
        out.invalidate_self_intersects();
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
    ///
    /// Port of C++ `CsgNode::Rotate` (csg_tree.cpp): degree-based `sind`/`cosd`
    /// axis matrices composed `rZ*rY*rX`. Quaternions (sin/cos of half-angles in
    /// radians) differ from this by ~1 ULP, which is enough to flip
    /// symbolic-perturbation ties in the boolean kernel on almost-coplanar
    /// inputs (`test_cpp_almost_coplanar` pins this).
    pub fn rotate(&self, x_degrees: f64, y_degrees: f64, z_degrees: f64) -> Self {
        use crate::types::{cosd, sind};
        let (sx, cx) = (sind(x_degrees), cosd(x_degrees));
        let (sy, cy) = (sind(y_degrees), cosd(y_degrees));
        let (sz, cz) = (sind(z_degrees), cosd(z_degrees));
        let rx = Mat3::from_cols(
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, cx, sx),
            Vec3::new(0.0, -sx, cx),
        );
        let ry = Mat3::from_cols(
            Vec3::new(cy, 0.0, -sy),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(sy, 0.0, cy),
        );
        let rz = Mat3::from_cols(
            Vec3::new(cz, sz, 0.0),
            Vec3::new(-sz, cz, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        );
        let m = rz * ry * rx;
        let t = Mat3x4::from_cols(m.x, m.y, m.z, Vec3::new(0.0, 0.0, 0.0));
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
        // Per C++ #1659: propagate an errored input even on the degenerate
        // (zero-length normal) path, which otherwise returns an empty manifold.
        if self.imp.status != Error::NoError {
            return self.clone();
        }
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
        if let Some(e) = self.require_paired() {
            return e;
        }
        if self.is_empty() {
            return self.clone();
        }
        let mut out = self.imp.clone();
        for v in out.vert_pos.iter_mut() {
            warp_fn(v);
        }
        // Arbitrary deformation: nothing about the old verdict survives.
        out.invalidate_self_intersects();
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

    /// Slice this manifold at the given Z height, returning the cross-section
    /// as a CrossSection. C++ `Manifold::Slice` returns raw `Polygons`; this
    /// is C++ `CrossSection(m.Slice(height))`, the Positive-union Polygons
    /// constructor every C++ caller wraps them in.
    ///
    /// The return type is this port's API shape, shipped long ago and
    /// mirrored by manifold-sharp, so changing it is a public-API decision,
    /// not a numerical fix. The numbers match: the contours are bit for bit
    /// those of the C++ `CrossSection(Polygons)` constructor applied to the
    /// C++ raw polygons (`cross_section::ctor_tests::test_slice_and_project_wrap_like_cpp`).
    /// The raw polygons stay reachable in-crate as `as_impl().slice(h)` /
    /// `as_impl().project()`, which the ported `TEST(Smooth, Fillet)` uses.
    pub fn slice(&self, height: f64) -> CrossSection {
        if self.imp.is_soup || self.is_empty() {
            return CrossSection::default();
        }
        CrossSection::new(self.imp.slice(height))
    }

    /// Project this manifold onto the XY plane, returning the silhouette
    /// as a CrossSection. C++ `Manifold::Project` returns raw, often
    /// self-overlapping `Polygons`; this is C++ `CrossSection(m.Project())`,
    /// the Positive-union Polygons constructor its docs recommend. The return
    /// type is an API shape, as for [`Self::slice`]; the numbers match.
    pub fn project(&self) -> CrossSection {
        if self.imp.is_soup || self.is_empty() {
            return CrossSection::default();
        }
        CrossSection::new(self.imp.project())
    }

    pub fn calculate_curvature(&self, gaussian_idx: i32, mean_idx: i32) -> Self {
        if let Some(e) = self.require_paired() {
            return e;
        }
        if self.is_empty() {
            return self.clone();
        }
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
        if let Some(e) = self.require_paired() {
            return e;
        }
        if self.is_empty() {
            return self.clone();
        }
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
                    let old_slice =
                        if old_num_prop > 0 && prop_vert * old_num_prop < old_properties.len() {
                            &old_properties
                                [old_num_prop * prop_vert..old_num_prop * prop_vert + old_num_prop]
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

    pub fn compose(parts: &[Self]) -> Self {
        let impls: Vec<_> = parts.iter().map(|m| m.imp.clone()).collect();
        Self::from_impl(boolean3::compose_meshes(&impls))
    }

    pub fn decompose(&self) -> Vec<Self> {
        use crate::disjoint_sets::DisjointSets;

        if let Some(e) = self.require_paired() {
            return vec![e];
        }
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

        // Each component is gathered from its own vertices, faces and property
        // vertices, in their original order, so the cost does not grow with
        // the number of components. `vert_local` and `face_local` hold every
        // vertex's and face's index within its own component.
        let source = &self.imp;
        let num_tri = source.num_tri();
        let num_prop = source.num_prop;
        let mut vert_local = vec![0i32; num_vert];
        let mut verts: Vec<Vec<usize>> = vec![Vec::new(); num_components as usize];
        for (vert, &comp) in component_indices.iter().enumerate() {
            vert_local[vert] = verts[comp as usize].len() as i32;
            verts[comp as usize].push(vert);
        }
        let mut face_local = vec![0i32; num_tri];
        let mut faces: Vec<Vec<usize>> = vec![Vec::new(); num_components as usize];
        for face in 0..num_tri {
            let start = source.halfedge[3 * face].start_vert;
            if start >= 0 {
                let comp = component_indices[start as usize] as usize;
                face_local[face] = faces[comp].len() as i32;
                faces[comp].push(face);
            }
        }
        let mut prop_local = vec![0i32; source.num_prop_vert()];

        let mut meshes = Vec::new();
        for (verts, faces) in verts.iter().zip(&faces) {
            if faces.is_empty() {
                continue;
            }
            let edges = || faces.iter().flat_map(|&face| 3 * face..3 * face + 3);

            let mut imp = ManifoldImpl::new();
            imp.tolerance = source.tolerance;
            imp.vert_pos = verts.iter().map(|&v| source.vert_pos[v]).collect();
            if !source.vert_normal.is_empty() {
                imp.vert_normal = verts.iter().map(|&v| source.vert_normal[v]).collect();
            }
            if source.face_normal.len() == num_tri {
                imp.face_normal = faces.iter().map(|&f| source.face_normal[f]).collect();
            }
            if !source.halfedge_tangent.is_empty() {
                imp.halfedge_tangent = edges().map(|e| source.halfedge_tangent[e]).collect();
            }
            imp.mesh_relation.original_id = source.mesh_relation.original_id;
            imp.mesh_relation.mesh_id_transform = source.mesh_relation.mesh_id_transform.clone();
            if source.mesh_relation.tri_ref.len() == num_tri {
                imp.mesh_relation.tri_ref = faces
                    .iter()
                    .map(|&f| source.mesh_relation.tri_ref[f])
                    .collect();
            }
            imp.halfedge = edges().map(|e| source.halfedge[e]).collect();

            imp.num_prop = num_prop;
            if num_prop > 0 {
                let mut props: Vec<usize> = imp
                    .halfedge
                    .iter()
                    .filter(|edge| edge.prop_vert >= 0)
                    .map(|edge| edge.prop_vert as usize)
                    .collect();
                props.sort_unstable();
                props.dedup();
                for (new_prop, &old_prop) in props.iter().enumerate() {
                    prop_local[old_prop] = new_prop as i32;
                    imp.properties.extend_from_slice(
                        &source.properties[old_prop * num_prop..(old_prop + 1) * num_prop],
                    );
                }
            }
            for edge in &mut imp.halfedge {
                edge.start_vert = vert_local[edge.start_vert as usize];
                edge.end_vert = vert_local[edge.end_vert as usize];
                if edge.paired_halfedge >= 0 {
                    let paired = edge.paired_halfedge as usize;
                    edge.paired_halfedge = 3 * face_local[paired / 3] + (paired % 3) as i32;
                }
                if num_prop == 0 {
                    edge.prop_vert = edge.start_vert;
                } else if edge.prop_vert >= 0 {
                    edge.prop_vert = prop_local[edge.prop_vert as usize];
                }
            }

            imp.calculate_bbox();
            imp.sort_geometry();

            meshes.push(Self::from_impl(imp));
        }

        meshes
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
}

#[path = "manifold_boolean.rs"]
mod boolean;

#[path = "manifold_robust.rs"]
mod robust;

#[path = "manifold_meshgl.rs"]
mod meshgl;
pub(crate) use meshgl::AssembledRefs;

#[path = "manifold_meshgl_export.rs"]
mod meshgl_export;

#[path = "manifold_shape.rs"]
mod shape;

#[path = "manifold_smooth.rs"]
mod smooth;

#[cfg(test)]
#[path = "manifold_tests/mod.rs"]
mod tests;
