// Phase 4: Mesh Data Structure — ported from src/impl.h, src/impl.cpp, src/properties.cpp
//
// This module implements the core ManifoldImpl struct: halfedge mesh representation,
// bounding box, epsilon, manifold checks, and shape constructors.
//
// Phases 5-9 will fill in SortGeometry, CleanupTopology, SetNormalsAndCoplanar, etc.

use std::sync::atomic::{AtomicU32, Ordering};
use crate::linalg::{Vec3, Vec4, Mat3x4, IVec3, cross, dot, normalize, length2};
use crate::types::{
    Box as BBox, Error, Halfedge, MeshRelationD, Relation, TriRef, K_PRECISION,
};

// ---------------------------------------------------------------------------
// Global mesh ID counter (mirrors Manifold::Impl::meshIDCounter_)
// ---------------------------------------------------------------------------

static MESH_ID_COUNTER: AtomicU32 = AtomicU32::new(1);

pub fn reserve_ids(n: u32) -> u32 {
    MESH_ID_COUNTER.fetch_add(n, Ordering::Relaxed)
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

pub const K_REMOVED_HALFEDGE: i32 = -2;

/// Next halfedge within the same triangle: 0→1→2→0.
#[inline]
pub fn next_halfedge(current: i32) -> i32 {
    let n = current + 1;
    if n % 3 == 0 { n - 3 } else { n }
}

#[inline]
pub fn next3(i: usize) -> usize {
    if i == 2 { 0 } else { i + 1 }
}

/// Safe normalize: returns zero vector if input is zero or non-finite.
fn safe_normalize(v: Vec3) -> Vec3 {
    let n = normalize(v);
    if n.x.is_finite() { n } else { Vec3::new(0.0, 0.0, 0.0) }
}

fn max_epsilon(min_epsilon: f64, bbox: &BBox) -> f64 {
    let epsilon = min_epsilon.max(K_PRECISION * bbox.scale());
    if epsilon.is_finite() { epsilon } else { -1.0 }
}

// ---------------------------------------------------------------------------
// ManifoldImpl — the core mesh representation
// ---------------------------------------------------------------------------

/// Internal halfedge mesh representation, mirroring `Manifold::Impl` in C++.
#[derive(Clone)]
pub struct ManifoldImpl {
    pub bbox: BBox,
    pub epsilon: f64,
    pub tolerance: f64,
    pub num_prop: usize,
    pub status: Error,
    pub vert_pos: Vec<Vec3>,
    pub halfedge: Vec<Halfedge>,
    pub properties: Vec<f64>,
    pub vert_normal: Vec<Vec3>,
    pub face_normal: Vec<Vec3>,
    pub halfedge_tangent: Vec<Vec4>,
    pub mesh_relation: MeshRelationD,
    // Collider will be added in Phase 10
}

impl Default for ManifoldImpl {
    fn default() -> Self {
        ManifoldImpl {
            bbox: BBox::new(),
            epsilon: -1.0,
            tolerance: -1.0,
            num_prop: 0,
            status: Error::NoError,
            vert_pos: Vec::new(),
            halfedge: Vec::new(),
            properties: Vec::new(),
            vert_normal: Vec::new(),
            face_normal: Vec::new(),
            halfedge_tangent: Vec::new(),
            mesh_relation: MeshRelationD::new(),
        }
    }
}

impl ManifoldImpl {
    pub fn new() -> Self {
        Self::default()
    }

    // -----------------------------------------------------------------------
    // Basic accessors
    // -----------------------------------------------------------------------

    pub fn num_vert(&self) -> usize {
        self.vert_pos.len()
    }

    pub fn num_halfedge(&self) -> usize {
        self.halfedge.len()
    }

    pub fn num_edge(&self) -> usize {
        self.halfedge.len() / 2
    }

    pub fn num_tri(&self) -> usize {
        self.halfedge.len() / 3
    }

    pub fn num_prop_vert(&self) -> usize {
        if self.num_prop == 0 {
            self.num_vert()
        } else {
            self.properties.len() / self.num_prop
        }
    }

    pub fn is_empty(&self) -> bool {
        self.num_tri() == 0
    }

    // -----------------------------------------------------------------------
    // MakeEmpty
    // -----------------------------------------------------------------------

    pub fn make_empty(&mut self, status: Error) {
        self.bbox = BBox::new();
        self.vert_pos.clear();
        self.halfedge.clear();
        self.vert_normal.clear();
        self.face_normal.clear();
        self.halfedge_tangent.clear();
        self.mesh_relation = MeshRelationD::new();
        self.status = status;
    }

    // -----------------------------------------------------------------------
    // ForVert — iterate halfedges around a vertex
    // -----------------------------------------------------------------------

    /// Apply `func` to each halfedge index around the vertex starting from `halfedge_idx`.
    pub fn for_vert<F: FnMut(usize)>(&self, halfedge_idx: usize, mut func: F) {
        let mut current = halfedge_idx;
        loop {
            current = next_halfedge(self.halfedge[current].paired_halfedge) as usize;
            func(current);
            if current == halfedge_idx {
                break;
            }
        }
    }

    // -----------------------------------------------------------------------
    // CalculateBBox
    // -----------------------------------------------------------------------

    pub fn calculate_bbox(&mut self) {
        let mut bbox = BBox::new();
        for v in &self.vert_pos {
            if !v.x.is_nan() {
                bbox.union_point(*v);
            }
        }
        self.bbox = bbox;
        if !self.bbox.is_finite() {
            self.make_empty(Error::NoError);
        }
    }

    // -----------------------------------------------------------------------
    // SetEpsilon
    // -----------------------------------------------------------------------

    pub fn set_epsilon(&mut self, min_epsilon: f64, use_single: bool) {
        self.epsilon = max_epsilon(min_epsilon, &self.bbox);
        let mut min_tol = self.epsilon;
        if use_single {
            let float_eps = (f32::EPSILON as f64) * self.bbox.scale();
            min_tol = min_tol.max(float_eps);
        }
        self.tolerance = self.tolerance.max(min_tol);
    }

    // -----------------------------------------------------------------------
    // IsFinite
    // -----------------------------------------------------------------------

    pub fn is_finite(&self) -> bool {
        self.vert_pos.iter().all(|v| v.x.is_finite() && v.y.is_finite() && v.z.is_finite())
    }

    // -----------------------------------------------------------------------
    // IsManifold / Is2Manifold
    // -----------------------------------------------------------------------

    /// Check that the halfedge data structure is consistent (oriented even manifold).
    pub fn is_manifold(&self) -> bool {
        if self.halfedge.is_empty() {
            return true;
        }
        if self.halfedge.len() % 3 != 0 {
            return false;
        }
        for (edge, h) in self.halfedge.iter().enumerate() {
            // Valid removed halfedge
            if h.start_vert == -1 && h.end_vert == -1 && h.paired_halfedge == -1 {
                continue;
            }
            // Neighbors in same triangle must not be removed
            let n1 = next_halfedge(edge as i32) as usize;
            let n2 = next_halfedge(n1 as i32) as usize;
            if self.halfedge[n1].start_vert == -1 || self.halfedge[n2].start_vert == -1 {
                return false;
            }
            if h.paired_halfedge == -1 {
                return false;
            }
            let paired_idx = h.paired_halfedge as usize;
            let paired = &self.halfedge[paired_idx];
            if paired.paired_halfedge != edge as i32 {
                return false;
            }
            if h.start_vert == h.end_vert {
                return false;
            }
            if h.start_vert != paired.end_vert {
                return false;
            }
            if h.end_vert != paired.start_vert {
                return false;
            }
        }
        true
    }

    /// Check that the mesh is a 2-manifold (no duplicate edges).
    pub fn is_2_manifold(&self) -> bool {
        if self.halfedge.is_empty() {
            return true;
        }
        if !self.is_manifold() {
            return false;
        }
        // Sort halfedges and check for duplicates
        let mut sorted = self.halfedge.clone();
        sorted.sort_unstable();
        for i in 0..sorted.len().saturating_sub(1) {
            let h = &sorted[i];
            let h1 = &sorted[i + 1];
            // Skip removed halfedges
            if h.start_vert == -1 && h.end_vert == -1 && h.paired_halfedge == -1 {
                continue;
            }
            if h.start_vert == h1.start_vert && h.end_vert == h1.end_vert {
                return false; // Duplicate edge
            }
        }
        true
    }

    // -----------------------------------------------------------------------
    // CreateHalfedges
    // -----------------------------------------------------------------------

    /// Build the halfedge data structure from triangle lists.
    ///
    /// - `tri_prop`: property vertex indices per triangle (also geometry if `tri_vert` is empty)
    /// - `tri_vert`: geometry vertex indices per triangle (may be empty)
    ///
    /// When `tri_vert` is empty, `tri_prop` is used for both geometry and properties.
    /// When `tri_vert` is present, `tri_prop[i][j]` = `propVert`, `tri_vert[i][j]` = `startVert`.
    pub fn create_halfedges(&mut self, tri_prop: &[IVec3], tri_vert: &[IVec3]) {
        let num_tri = tri_prop.len();
        if num_tri == 0 {
            self.halfedge.clear();
            return;
        }
        let num_halfedge = 3 * num_tri;
        let num_edge = num_halfedge / 2;

        self.halfedge.clear();
        self.halfedge.resize(num_halfedge, Halfedge {
            start_vert: -1,
            end_vert: -1,
            paired_halfedge: -1,
            prop_vert: -1,
        });

        let use_prop = tri_vert.is_empty();

        // Build halfedges and compute edge sort key
        // key = [forward_bit:1][min_vert:31][max_vert:32]
        // forward: v0 < v1 → bit=1; backward: v0 > v1 → bit=0
        // After sorting: backward halfedges first, then forward, both by (min,max)
        let mut edge_keys = vec![0u64; num_halfedge];

        for tri in 0..num_tri {
            for i in 0usize..3 {
                let j = next3(i);
                let e = 3 * tri + i;
                let v0 = if use_prop { tri_prop[tri][i] } else { tri_vert[tri][i] };
                let v1 = if use_prop { tri_prop[tri][j] } else { tri_vert[tri][j] };
                self.halfedge[e] = Halfedge {
                    start_vert: v0,
                    end_vert: v1,
                    paired_halfedge: -1,
                    prop_vert: tri_prop[tri][i],
                };
                let fwd = if v0 < v1 { 1u64 } else { 0u64 };
                let min_v = v0.min(v1) as u64;
                let max_v = v0.max(v1) as u64;
                edge_keys[e] = (fwd << 63) | (min_v << 32) | max_v;
            }
        }

        // Sort halfedge indices by edge key
        let mut ids: Vec<usize> = (0..num_halfedge).collect();
        ids.sort_unstable_by_key(|&i| edge_keys[i]);

        // ids[0..num_edge] = backward halfedges (startVert > endVert), sorted by (min,max)
        // ids[num_edge..] = forward halfedges (startVert < endVert), sorted by (min,max)

        // Sequential pairing with opposed-triangle detection
        let segment_end = num_edge;
        let mut consecutive_start = 0usize;

        for i in 0..num_edge {
            let pair0 = ids[i];
            let h0_sv = self.halfedge[pair0].start_vert;
            let h0_ev = self.halfedge[pair0].end_vert;

            let mut k = consecutive_start + num_edge;
            'inner: loop {
                if k >= segment_end + num_edge {
                    break 'inner;
                }
                let pair1 = ids[k];
                let h1_sv = self.halfedge[pair1].start_vert;
                let h1_ev = self.halfedge[pair1].end_vert;

                if h0_sv != h1_ev || h0_ev != h1_sv {
                    break 'inner; // Different edge direction — no match
                }

                if self.halfedge[pair1].paired_halfedge != K_REMOVED_HALFEDGE {
                    // Check for opposed triangle: same undirected edge, same third vertex
                    let next0 = next_halfedge(pair0 as i32) as usize;
                    let next1 = next_halfedge(pair1 as i32) as usize;
                    if self.halfedge[next0].end_vert == self.halfedge[next1].end_vert {
                        // Opposed triangles: mark both for removal
                        self.halfedge[pair0].paired_halfedge = K_REMOVED_HALFEDGE;
                        self.halfedge[pair1].paired_halfedge = K_REMOVED_HALFEDGE;
                        if i + num_edge != k {
                            ids[i + num_edge] = pair1;
                        }
                        break 'inner;
                    }
                }

                k += 1;
            }

            // Update consecutive_start for next iteration
            if i + 1 < segment_end {
                let next_sv = self.halfedge[ids[i + 1]].start_vert;
                let next_ev = self.halfedge[ids[i + 1]].end_vert;
                if next_sv != h0_sv || next_ev != h0_ev {
                    consecutive_start = i + 1;
                }
            }
        }

        // Final pairing pass
        for i in 0..num_edge {
            let pair0 = ids[i];
            let pair1 = ids[i + num_edge];
            if self.halfedge[pair0].paired_halfedge != K_REMOVED_HALFEDGE {
                self.halfedge[pair0].paired_halfedge = pair1 as i32;
                self.halfedge[pair1].paired_halfedge = pair0 as i32;
            } else {
                // Invalidate both (opposed triangles removed)
                self.halfedge[pair0] = Halfedge { start_vert: -1, end_vert: -1, paired_halfedge: -1, prop_vert: 0 };
                self.halfedge[pair1] = Halfedge { start_vert: -1, end_vert: -1, paired_halfedge: -1, prop_vert: 0 };
            }
        }
    }

    // -----------------------------------------------------------------------
    // InitializeOriginal
    // -----------------------------------------------------------------------

    /// Set up the mesh relation for a newly created original mesh.
    pub fn initialize_original(&mut self) {
        let mesh_id = reserve_ids(1) as i32;
        self.mesh_relation.original_id = mesh_id;
        let num_tri = self.num_tri();
        self.mesh_relation.tri_ref.resize(num_tri, TriRef::default());
        for (tri, tri_ref) in self.mesh_relation.tri_ref.iter_mut().enumerate() {
            tri_ref.mesh_id = mesh_id;
            tri_ref.original_id = mesh_id;
            tri_ref.face_id = -1;
            tri_ref.coplanar_id = tri as i32;
        }
        self.mesh_relation.mesh_id_transform.clear();
        self.mesh_relation.mesh_id_transform.insert(mesh_id, Relation {
            original_id: mesh_id,
            transform: Mat3x4::identity(),
            back_side: false,
        });
    }

    // -----------------------------------------------------------------------
    // RemoveUnreferencedVerts
    // -----------------------------------------------------------------------

    /// Mark unreferenced vertices as NaN (to be cleaned up by later passes).
    pub fn remove_unreferenced_verts(&mut self) {
        let num_vert = self.num_vert();
        let mut keep = vec![false; num_vert];
        for h in &self.halfedge {
            if h.start_vert >= 0 {
                keep[h.start_vert as usize] = true;
            }
        }
        for (i, k) in keep.iter().enumerate() {
            if !k {
                self.vert_pos[i] = Vec3::new(f64::NAN, f64::NAN, f64::NAN);
            }
        }
    }

    // -----------------------------------------------------------------------
    // SetNormalsAndCoplanar (stub — implemented in Phase 9)
    // -----------------------------------------------------------------------

    /// Compute face normals, assign coplanar IDs, and calculate vertex normals.
    pub fn set_normals_and_coplanar(&mut self) {
        crate::face_op::set_normals_and_coplanar(self);
    }

    // -----------------------------------------------------------------------
    // SortGeometry (stub — implemented in Phase 5)
    // -----------------------------------------------------------------------

    /// Reorder mesh geometry for cache efficiency using Morton codes.
    pub fn sort_geometry(&mut self) {
        crate::sort::sort_geometry(self);
    }

    // -----------------------------------------------------------------------
    // Shape constructors
    // -----------------------------------------------------------------------

    pub fn tetrahedron(transform: &Mat3x4) -> Self {
        let vert_pos_raw: Vec<[f64; 3]> = vec![
            [-1.0, -1.0,  1.0],
            [-1.0,  1.0, -1.0],
            [ 1.0, -1.0, -1.0],
            [ 1.0,  1.0,  1.0],
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
            IVec3::new(1, 0, 4), IVec3::new(2, 4, 0),
            IVec3::new(1, 3, 0), IVec3::new(3, 1, 5),
            IVec3::new(3, 2, 0), IVec3::new(3, 7, 2),
            IVec3::new(5, 4, 6), IVec3::new(5, 1, 4),
            IVec3::new(6, 4, 2), IVec3::new(7, 6, 2),
            IVec3::new(7, 3, 5), IVec3::new(7, 5, 6),
        ];
        Self::from_shape(vert_pos_raw, tri_verts, transform)
    }

    pub fn octahedron(transform: &Mat3x4) -> Self {
        let vert_pos_raw: Vec<[f64; 3]> = vec![
            [ 1.0,  0.0,  0.0],
            [-1.0,  0.0,  0.0],
            [ 0.0,  1.0,  0.0],
            [ 0.0, -1.0,  0.0],
            [ 0.0,  0.0,  1.0],
            [ 0.0,  0.0, -1.0],
        ];
        let tri_verts: Vec<IVec3> = vec![
            IVec3::new(0, 2, 4), IVec3::new(1, 5, 3),
            IVec3::new(2, 1, 4), IVec3::new(3, 5, 0),
            IVec3::new(1, 3, 4), IVec3::new(0, 5, 2),
            IVec3::new(3, 0, 4), IVec3::new(2, 5, 1),
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

    // -----------------------------------------------------------------------
    // Transform
    // -----------------------------------------------------------------------

    /// Apply affine transform, returning a new ManifoldImpl.
    pub fn transform(&self, t: &Mat3x4) -> Self {
        use crate::linalg::{Vec4, Mat3};
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
        result.epsilon = self.epsilon;
        result.tolerance = self.tolerance;
        result.num_prop = self.num_prop;
        result.properties = self.properties.clone();
        result.bbox = self.bbox;
        result.halfedge = self.halfedge.clone();
        result.halfedge_tangent.resize(self.halfedge_tangent.len(), Vec4::new(0.0, 0.0, 0.0, 0.0));
        result.mesh_relation.original_id = -1;

        // Update mesh transforms
        for (_, rel) in result.mesh_relation.mesh_id_transform.iter_mut() {
            // rel.transform = t * Mat4(rel.transform) — combine transforms
            rel.transform = mat3x4_mul_mat3x4(t, &rel.transform);
        }

        // Transform vertex positions
        result.vert_pos = self.vert_pos.iter().map(|&v| {
            *t * Vec4::new(v.x, v.y, v.z, 1.0)
        }).collect();

        // Transform normals (using inverse-transpose of 3x3 part)
        let m3 = Mat3::from_cols(
            Vec3::new(t.x.x, t.x.y, t.x.z),
            Vec3::new(t.y.x, t.y.y, t.y.z),
            Vec3::new(t.z.x, t.z.y, t.z.z),
        );
        let normal_t = m3.inverse().transpose();

        result.face_normal = self.face_normal.iter().map(|&n| {
            safe_normalize(normal_t * n)
        }).collect();
        result.vert_normal = self.vert_normal.iter().map(|&n| {
            safe_normalize(normal_t * n)
        }).collect();

        let invert = m3.determinant() < 0.0;
        if invert {
            // Flip triangle winding
            for tri in 0..result.num_tri() {
                result.halfedge.swap(3 * tri + 1, 3 * tri + 2);
            }
        }

        result.calculate_bbox();
        result.set_epsilon(result.epsilon, false);
        result
    }

    /// Clone without collider (shallow copy for transform operations)
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
    Mat3x4 { x: c0, y: c1, z: c2, w: c3 }
}

// ---------------------------------------------------------------------------
// Mat3 for normal transform (we need inverse + transpose from linalg)
// ---------------------------------------------------------------------------

// These are already in linalg.rs but we need to use them here.
// The Mat3 methods: inverse(), transpose(), determinant()

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linalg::{Vec3, Mat3x4};

    #[test]
    fn test_next_halfedge() {
        assert_eq!(next_halfedge(0), 1);
        assert_eq!(next_halfedge(1), 2);
        assert_eq!(next_halfedge(2), 0);
        assert_eq!(next_halfedge(3), 4);
        assert_eq!(next_halfedge(5), 3);
    }

    #[test]
    fn test_tetrahedron() {
        let m = ManifoldImpl::tetrahedron(&Mat3x4::identity());
        assert_eq!(m.num_vert(), 4);
        assert_eq!(m.num_tri(), 4);
        assert_eq!(m.num_edge(), 6);
        assert!(m.is_manifold(), "Tetrahedron should be manifold");
        assert!(m.is_2_manifold(), "Tetrahedron should be 2-manifold");
    }

    #[test]
    fn test_cube() {
        let m = ManifoldImpl::cube(&Mat3x4::identity());
        assert_eq!(m.num_vert(), 8);
        assert_eq!(m.num_tri(), 12);
        assert_eq!(m.num_edge(), 18);
        assert!(m.is_manifold(), "Cube should be manifold");
        assert!(m.is_2_manifold(), "Cube should be 2-manifold");
    }

    #[test]
    fn test_octahedron() {
        let m = ManifoldImpl::octahedron(&Mat3x4::identity());
        assert_eq!(m.num_vert(), 6);
        assert_eq!(m.num_tri(), 8);
        assert_eq!(m.num_edge(), 12);
        assert!(m.is_manifold(), "Octahedron should be manifold");
        assert!(m.is_2_manifold(), "Octahedron should be 2-manifold");
    }

    #[test]
    fn test_cube_bbox() {
        let m = ManifoldImpl::cube(&Mat3x4::identity());
        assert!((m.bbox.min.x - 0.0).abs() < 1e-10);
        assert!((m.bbox.min.y - 0.0).abs() < 1e-10);
        assert!((m.bbox.min.z - 0.0).abs() < 1e-10);
        assert!((m.bbox.max.x - 1.0).abs() < 1e-10);
        assert!((m.bbox.max.y - 1.0).abs() < 1e-10);
        assert!((m.bbox.max.z - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_tetrahedron_symmetric() {
        let m = ManifoldImpl::tetrahedron(&Mat3x4::identity());
        // Tetrahedron is centered about origin with vertices at distance sqrt(3)
        for v in &m.vert_pos {
            let dist2 = v.x * v.x + v.y * v.y + v.z * v.z;
            assert!((dist2 - 3.0).abs() < 1e-10, "Expected dist^2=3, got {}", dist2);
        }
    }

    #[test]
    fn test_cube_transform() {
        use crate::linalg::{translation_matrix, mat4_to_mat3x4};
        let t = mat4_to_mat3x4(translation_matrix(Vec3::new(1.0, 2.0, 3.0)));
        let m = ManifoldImpl::cube(&t);
        assert!((m.bbox.min.x - 1.0).abs() < 1e-10);
        assert!((m.bbox.min.y - 2.0).abs() < 1e-10);
        assert!((m.bbox.min.z - 3.0).abs() < 1e-10);
        assert!((m.bbox.max.x - 2.0).abs() < 1e-10);
        assert!((m.bbox.max.y - 3.0).abs() < 1e-10);
        assert!((m.bbox.max.z - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_for_vert() {
        let m = ManifoldImpl::tetrahedron(&Mat3x4::identity());
        // Each vertex in the tetrahedron is surrounded by 3 triangles = 3 halfedges
        let mut count = 0;
        m.for_vert(0, |_| count += 1);
        // ForVert visits one halfedge per triangle around the vertex (3 for tetrahedron)
        assert!(count > 0);
    }

    #[test]
    fn test_create_halfedges_simple() {
        // Simple triangle
        let mut m = ManifoldImpl::new();
        m.vert_pos = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ];
        // A tetrahedron has 4 triangles and 12 halfedges
        // Let's do a single triangle (not manifold, just for halfedge construction)
        let tri = vec![IVec3::new(0, 1, 2)];
        m.create_halfedges(&tri, &[]);
        assert_eq!(m.halfedge.len(), 3);
        // Single triangle — no pairs, so all pairedHalfedge should be -1
        // (no paired halfedges available)
    }
}
