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

// impl_transform.rs — affine transforms of ManifoldImpl
//
// Ports Manifold::Impl::Transform (src/impl.cpp) together with the pieces it
// drives: TransformTangents, FlipTris (src/mesh_fixes.h) for
// negative-determinant transforms, the per-meshID eager transform of recorded
// normals in property slots 0..2 (C++ #1718), and the collider refit that
// keeps the cached BVH valid without a rebuild. Split out of impl_mesh.rs,
// which defines the ManifoldImpl struct; this file adds a second `impl
// ManifoldImpl` block. The collider refit relies on collider.rs
// (is_axis_aligned / transform / update_boxes) and sort.rs (face boxes).

use crate::impl_mesh::ManifoldImpl;
use crate::linalg::{normalize, Mat3x4, Vec3};
use crate::types::Error;

/// Safe normalize: returns zero vector if input is zero or non-finite.
fn safe_normalize(v: Vec3) -> Vec3 {
    let n = normalize(v);
    if n.x.is_finite() {
        n
    } else {
        Vec3::new(0.0, 0.0, 0.0)
    }
}

impl ManifoldImpl {
    /// Eager-transform slot 0..2 of `properties` for propVerts whose meshID
    /// carries hasNormals. Used by both `transform` and `compose` so world-frame
    /// normals stay in sync with vert_pos / face_normal across any sequence of
    /// transforms (including mixed-input Boolean/Compose outputs where some
    /// meshIDs carry normals and others don't). Per C++ #1718.
    ///
    /// `properties` is laid out as `properties[(offset + prop) * stride + i]`,
    /// so callers can target an in-place properties_ vector (offset=0) or a
    /// per-node slice of a combined array (offset=propVertIndices,
    /// stride=numPropOut). Re-normalizes as it transforms so non-orthogonal
    /// transforms (scale) and upstream barycentric interpolation don't leave
    /// non-unit values that compound downstream.
    pub fn eager_transform_prop_normals(
        halfedge: &[crate::types::Halfedge],
        mesh_relation: &crate::types::MeshRelationD,
        normal_transform: crate::linalg::Mat3,
        properties: &mut [f64],
        num_prop_vert: usize,
        stride: usize,
        offset: usize,
    ) {
        // OR semantics (any meshID has normals), unlike all_have_normals():
        // mixed inputs still need the per-meshID iteration below to rotate the
        // with-normals subset.
        if !mesh_relation
            .mesh_id_transform
            .values()
            .any(|m| m.has_normals)
        {
            return;
        }
        let tri_has_normals = |tri: usize| -> bool {
            let mid = mesh_relation.tri_ref[tri].mesh_id;
            mesh_relation
                .mesh_id_transform
                .get(&mid)
                .map(|m| m.has_normals)
                .unwrap_or(false)
        };
        let mut visited = vec![false; num_prop_vert];
        for e in 0..halfedge.len() {
            if !tri_has_normals(e / 3) {
                continue;
            }
            let prop = halfedge[e].prop_vert;
            if prop < 0 {
                continue;
            }
            let prop = prop as usize;
            if visited[prop] {
                continue;
            }
            visited[prop] = true;
            let base = (offset + prop) * stride;
            let n = Vec3::new(properties[base], properties[base + 1], properties[base + 2]);
            let nt = safe_normalize(normal_transform * n);
            properties[base] = nt.x;
            properties[base + 1] = nt.y;
            properties[base + 2] = nt.z;
        }
    }

    // -----------------------------------------------------------------------
    // Transform
    // -----------------------------------------------------------------------

    /// Apply affine transform, returning a new ManifoldImpl.
    pub fn transform(&self, t: &Mat3x4) -> Self {
        use crate::linalg::{Mat3, Vec4};
        let identity = Mat3x4::identity();
        if t == &identity {
            // Clone self — this is a simplified version (full version uses Collider)
            return self.shallow_clone();
        }

        let mut result = Self::new();
        if self.status != Error::NoError {
            result.status = self.status;
            return result;
        }

        result.mesh_relation = self.mesh_relation.clone();
        // Scale epsilon by spectral norm of transform, matching C++:
        // result.epsilon_ *= SpectralNorm(mat3(transform_));
        let m3_for_norm = Mat3::from_cols(
            Vec3::new(t.x.x, t.x.y, t.x.z),
            Vec3::new(t.y.x, t.y.y, t.y.z),
            Vec3::new(t.z.x, t.z.y, t.z.z),
        );
        result.epsilon = self.epsilon * crate::svd::spectral_norm(m3_for_norm);
        result.tolerance = self.tolerance;
        result.num_prop = self.num_prop;
        result.properties = self.properties.clone();
        result.bbox = self.bbox;
        result.halfedge = self.halfedge.clone();
        result.mesh_relation.original_id = -1;
        // Soup impls stay soups across transforms; every step below already
        // guards paired_halfedge < 0.
        result.is_soup = self.is_soup;
        // The self-intersection cache is deliberately *not* carried across:
        // transformed positions are rounded to f64, so an extreme scale can
        // collapse distinct vertices onto each other and create coincident
        // surface that the source mesh did not have. Re-running the detector
        // costs microseconds; propagating a stale `false` costs correctness.

        // Update mesh transforms
        for (_, rel) in result.mesh_relation.mesh_id_transform.iter_mut() {
            // rel.transform = t * Mat4(rel.transform) — combine transforms
            rel.transform = mat3x4_mul_mat3x4(t, &rel.transform);
        }

        // Transform vertex positions
        result.vert_pos = self
            .vert_pos
            .iter()
            .map(|&v| *t * Vec4::new(v.x, v.y, v.z, 1.0))
            .collect();

        // Transform normals (using inverse-transpose of 3x3 part)
        let m3 = Mat3::from_cols(
            Vec3::new(t.x.x, t.x.y, t.x.z),
            Vec3::new(t.y.x, t.y.y, t.y.z),
            Vec3::new(t.z.x, t.z.y, t.z.z),
        );
        let normal_t = m3.inverse().transpose();

        result.face_normal = self
            .face_normal
            .iter()
            .map(|&n| safe_normalize(normal_t * n))
            .collect();
        result.vert_normal = self
            .vert_normal
            .iter()
            .map(|&n| safe_normalize(normal_t * n))
            .collect();

        // Per #1718: the properties clone above doesn't go through the vertPos /
        // faceNormal transform, so eager-transform slot 0..2 per-meshID to keep
        // recorded world-frame normals in sync. tri_ref / hasNormals flags are
        // identical in self and result; iterate by prop vert (winding flip
        // below only reorders halfedges, not prop assignments).
        if self.num_prop >= 3 {
            Self::eager_transform_prop_normals(
                &self.halfedge,
                &self.mesh_relation,
                normal_t,
                &mut result.properties,
                self.num_prop_vert(),
                self.num_prop,
                0,
            );
        }

        let invert = m3.determinant() < 0.0;

        // Transform tangents — C++ TransformTangents
        // Must happen BEFORE FlipTris (matches C++ order)
        if !self.halfedge_tangent.is_empty() {
            result.halfedge_tangent =
                vec![Vec4::new(0.0, 0.0, 0.0, 0.0); self.halfedge_tangent.len()];
            for edge_out in 0..self.halfedge_tangent.len() {
                let edge_in = if invert {
                    let tri = edge_out / 3;
                    let vert = 2 - (edge_out - 3 * tri);
                    let flipped = 3 * tri + vert;
                    self.halfedge[flipped].paired_halfedge as usize
                } else {
                    edge_out
                };
                let old_t = self.halfedge_tangent[edge_in];
                let xyz = m3 * Vec3::new(old_t.x, old_t.y, old_t.z);
                result.halfedge_tangent[edge_out] = Vec4::new(xyz.x, xyz.y, xyz.z, old_t.w);
            }
        }

        if invert {
            // Flip triangle winding — matches C++ FlipTris
            for tri in 0..result.num_tri() {
                // Props belong to corners (start verts), not halfedges: after
                // the flip the new start verts are old starts (0, 2, 1), so
                // the props must follow that order. Pinned v3.5.2 FlipTris let
                // props travel with the swapped halfedges; upstream 422ab6fc
                // (issue #1781) fixed it — see docs/CPP_DIVERGENCES.md.
                let props = [
                    result.halfedge[3 * tri].prop_vert,
                    result.halfedge[3 * tri + 2].prop_vert,
                    result.halfedge[3 * tri + 1].prop_vert,
                ];
                // Swap first and third halfedge within tri
                result.halfedge.swap(3 * tri, 3 * tri + 2);
                // For each halfedge: swap startVert/endVert and remap pairedHalfedge
                for i in 0..3 {
                    let idx = 3 * tri + i;
                    let h = &mut result.halfedge[idx];
                    std::mem::swap(&mut h.start_vert, &mut h.end_vert);
                    h.prop_vert = props[i];
                    // FlipHalfedge: within the paired tri, mirror the edge index
                    let paired = h.paired_halfedge;
                    if paired >= 0 {
                        let p = paired as usize;
                        let p_tri = p / 3;
                        let p_vert = 2 - (p - 3 * p_tri);
                        h.paired_halfedge = (3 * p_tri + p_vert) as i32;
                    }
                }
            }
        }

        result.calculate_bbox();
        result.set_epsilon(result.epsilon, false);

        // Keep the cached collider valid without a full rebuild, mirroring C++
        // Impl::Transform: an axis-aligned transform maps the existing tree's
        // boxes directly; otherwise recompute leaf boxes on the transformed
        // mesh and refit the same tree topology.
        // Soup impls never went through sort_geometry, which is where the
        // collider is built, so they carry a zero-leaf tree. There is no
        // topology to refit in that case — cloning the empty tree and
        // refitting it against real face boxes indexes out of bounds, and
        // update_boxes' debug_assert is compiled out in release. Leave the
        // default collider, exactly as the untransformed soup carried.
        if !result.is_empty() && self.collider.num_leaves() == self.num_tri() {
            if crate::collider::Collider::is_axis_aligned(t) {
                result.collider = self.collider.clone();
                result.collider.transform(t);
            } else {
                result.collider = self.collider.clone();
                let (face_box, _face_morton) = crate::sort::get_face_box_morton(&result);
                result.collider.update_boxes(face_box);
            }
        }
        result
    }

    /// Field-by-field copy used by the identity-transform fast path; the
    /// collider is copied as-is (C++ copies collider_ with the Impl).
    fn shallow_clone(&self) -> Self {
        ManifoldImpl {
            bbox: self.bbox,
            epsilon: self.epsilon,
            tolerance: self.tolerance,
            num_prop: self.num_prop,
            status: self.status,
            vert_pos: self.vert_pos.clone(),
            halfedge: self.halfedge.clone(),
            properties: self.properties.clone(),
            vert_normal: self.vert_normal.clone(),
            face_normal: self.face_normal.clone(),
            halfedge_tangent: self.halfedge_tangent.clone(),
            mesh_relation: self.mesh_relation.clone(),
            collider: self.collider.clone(),
            is_soup: self.is_soup,
            self_intersects: self.self_intersects.clone(),
        }
    }
}

// ---------------------------------------------------------------------------
// Transform helpers
// ---------------------------------------------------------------------------

/// Multiply two Mat3x4 transforms as affine matrices (t1 * t2 = (t1 * to_mat4(t2)).
/// Result is t1 applied after t2.
fn mat3x4_mul_mat3x4(t1: &Mat3x4, t2: &Mat3x4) -> Mat3x4 {
    use crate::linalg::Vec4;
    // Column vectors of t2 (as Vec4 with w=0 for rotation cols, w=1 for translation)
    let c0 = *t1 * Vec4::new(t2.x.x, t2.x.y, t2.x.z, 0.0);
    let c1 = *t1 * Vec4::new(t2.y.x, t2.y.y, t2.y.z, 0.0);
    let c2 = *t1 * Vec4::new(t2.z.x, t2.z.y, t2.z.z, 0.0);
    let c3 = *t1 * Vec4::new(t2.w.x, t2.w.y, t2.w.z, 1.0);
    Mat3x4 {
        x: c0,
        y: c1,
        z: c2,
        w: c3,
    }
}

// ---------------------------------------------------------------------------
// Mat3 for normal transform (we need inverse + transpose from linalg)
// ---------------------------------------------------------------------------

// These are already in linalg.rs but we need to use them here.
// The Mat3 methods: inverse(), transpose(), determinant()
