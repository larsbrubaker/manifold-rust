// Phase 4: Mesh Data Structure — ported from src/impl.h, src/impl.cpp, src/properties.cpp
//
// This module implements the core ManifoldImpl struct: halfedge mesh representation,
// bounding box, epsilon, manifold checks, halfedge construction, and mesh
// relation bookkeeping. Further `impl ManifoldImpl` blocks live in
// impl_transform.rs (affine transforms, winding flip, collider refit) and
// impl_shapes.rs (tetrahedron / cube / octahedron constructors).
//
// Phases 5-9 will fill in SortGeometry, CleanupTopology, SetNormalsAndCoplanar, etc.

use crate::linalg::{IVec3, Mat3x4, Vec3, Vec4};
use crate::types::{Box as BBox, Error, Halfedge, MeshRelationD, Relation, TriRef, K_PRECISION};
use std::sync::atomic::{AtomicU32, Ordering};

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
    if n % 3 == 0 {
        n - 3
    } else {
        n
    }
}

#[inline]
pub fn next3(i: usize) -> usize {
    if i == 2 {
        0
    } else {
        i + 1
    }
}

fn max_epsilon(min_epsilon: f64, bbox: &BBox) -> f64 {
    let epsilon = min_epsilon.max(K_PRECISION * bbox.scale());
    if epsilon.is_finite() {
        epsilon
    } else {
        -1.0
    }
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
    /// Cached face BVH, built by sort_geometry and updated (not rebuilt) on
    /// transform — mirrors C++ Impl::collider_. Query sites (boolean kernels,
    /// ray cast, face merging, self-intersection) use this instead of
    /// rebuilding the tree per query.
    pub collider: crate::collider::Collider,
    /// True when this impl carries geometrically closed but topologically
    /// non-manifold "triangle soup" imported via `from_mesh_gl_robust`:
    /// halfedge pairing is incomplete (`paired_halfedge == -1` permitted).
    /// Only the robust boolean engine, transforms, bbox, and MeshGL export
    /// accept soup impls; pairing-dependent operations return an empty
    /// result with `Error::NotManifold`. Always false on the strict import
    /// path, so existing behavior is unchanged.
    pub is_soup: bool,
    /// Lazily-resolved verdict of `robust::soup::has_self_intersections` —
    /// whether this impl's own triangles genuinely intersect rather than
    /// merely sharing edges/vertices. Computed at most once per impl (the
    /// scan is a full BVH self-query) and consulted by `Auto` boolean
    /// dispatch, which must route geometrically self-intersecting operands
    /// to the robust engine even when their connectivity is manifold.
    ///
    /// `OnceLock`-backed so it can be filled through a shared `&ManifoldImpl`
    /// from rayon workers under the `parallel` feature. Crate-private
    /// deliberately: `Clone` copies the settled value, so **any** code that
    /// clones an impl and then edits its geometry in place must call
    /// [`ManifoldImpl::invalidate_self_intersects`]. Rebuilds that go through
    /// `create_halfedges` or `make_empty` are covered automatically.
    pub(crate) self_intersects: crate::robust::soup::SelfIntersectCache,
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
            collider: crate::collider::Collider::default(),
            is_soup: false,
            self_intersects: Default::default(),
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
        self.collider = crate::collider::Collider::default();
        self.status = status;
        self.is_soup = false;
        // Geometry is gone; any cached self-intersection verdict is stale.
        self.invalidate_self_intersects();
    }

    /// Drop any cached self-intersection verdict. Must be called by every
    /// operation that edits `vert_pos` or `halfedge` in place on an impl it
    /// cloned (the cache is copied by `Clone`); rebuilds through
    /// `create_halfedges` and `make_empty` do it themselves.
    pub fn invalidate_self_intersects(&mut self) {
        self.self_intersects = Default::default();
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
        self.vert_pos
            .iter()
            .all(|v| v.x.is_finite() && v.y.is_finite() && v.z.is_finite())
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
        // The triangle set is being (re)built, so any earlier verdict about
        // self-intersection no longer describes this geometry.
        self.invalidate_self_intersects();
        let num_tri = tri_prop.len();
        if num_tri == 0 {
            self.halfedge.clear();
            return;
        }
        let num_halfedge = 3 * num_tri;
        let num_edge = num_halfedge / 2;

        self.halfedge.clear();
        self.halfedge.resize(
            num_halfedge,
            Halfedge {
                start_vert: -1,
                end_vert: -1,
                paired_halfedge: -1,
                prop_vert: -1,
            },
        );

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
                let v0 = if use_prop {
                    tri_prop[tri][i]
                } else {
                    tri_vert[tri][i]
                };
                let v1 = if use_prop {
                    tri_prop[tri][j]
                } else {
                    tri_vert[tri][j]
                };
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

        // Sort halfedge indices by edge key. C++ CreateHalfedges uses a
        // STABLE sort here (impl.cpp), and the #1687 fix ensures its parallel
        // stable_sort matches std::stable_sort. When two halfedges share an
        // edge key (duplicate directed edges in degenerate/intermediate
        // meshes) the tie must break on original halfedge-index order, so we
        // use a stable sort to stay bit-identical to C++.
        let mut ids: Vec<usize> = (0..num_halfedge).collect();
        ids.sort_by_key(|&i| edge_keys[i]);

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
                        // Opposed triangles: mark both for removal.
                        // Reorder ids so the remaining valid forward halfedge (at i+num_edge)
                        // moves to position k, and pair1 (the opposed one) goes to i+num_edge.
                        // This matches C++ which does: ids[k] = ids[i+numEdge]; ids[i+numEdge] = pair1;
                        self.halfedge[pair0].paired_halfedge = K_REMOVED_HALFEDGE;
                        self.halfedge[pair1].paired_halfedge = K_REMOVED_HALFEDGE;
                        if i + num_edge != k {
                            ids.swap(k, i + num_edge);
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
                self.halfedge[pair0] = Halfedge {
                    start_vert: -1,
                    end_vert: -1,
                    paired_halfedge: -1,
                    prop_vert: 0,
                };
                self.halfedge[pair1] = Halfedge {
                    start_vert: -1,
                    end_vert: -1,
                    paired_halfedge: -1,
                    prop_vert: 0,
                };
            }
        }
    }

    // -----------------------------------------------------------------------
    // InitializeOriginal
    // -----------------------------------------------------------------------

    /// Set up the mesh relation for a newly created original mesh.
    pub fn initialize_original(&mut self) {
        // Per C++ #1718: preserve the AND-across-old-Relations hasNormals state
        // so AsOriginal keeps the recording when it builds a fresh Relation.
        // Primitives start with an empty map → all_have_normals() is false.
        let had_normals = self.all_have_normals();
        let mesh_id = reserve_ids(1) as i32;
        self.mesh_relation.original_id = mesh_id;
        let num_tri = self.num_tri();
        self.mesh_relation
            .tri_ref
            .resize(num_tri, TriRef::default());
        for (tri, tri_ref) in self.mesh_relation.tri_ref.iter_mut().enumerate() {
            tri_ref.mesh_id = mesh_id;
            tri_ref.original_id = mesh_id;
            tri_ref.face_id = -1;
            tri_ref.coplanar_id = tri as i32;
        }
        self.mesh_relation.mesh_id_transform.clear();
        self.mesh_relation.mesh_id_transform.insert(
            mesh_id,
            Relation {
                original_id: mesh_id,
                transform: Mat3x4::identity(),
                back_side: false,
                has_normals: had_normals,
            },
        );
    }

    /// True only when every meshID carries normals at slot 0..2 — the condition
    /// under which `get_mesh_gl(-1)` can safely auto-substitute that slot. A
    /// mixed Boolean output (some meshIDs with normals, some without) returns
    /// false; the output MeshGL's per-run bit 1 still marks the with-normals
    /// runs individually. AND semantics across meshIDs. Per C++ #1718.
    pub fn all_have_normals(&self) -> bool {
        let map = &self.mesh_relation.mesh_id_transform;
        !map.is_empty() && map.values().all(|m| m.has_normals)
    }

    /// True iff the meshID owning `tri` has hasNormals set. False when the
    /// meshID isn't in mesh_id_transform (treat as no-normals). Per C++ #1718.
    pub fn tri_has_normals(&self, tri: usize) -> bool {
        let mesh_id = self.mesh_relation.tri_ref[tri].mesh_id;
        self.mesh_relation
            .mesh_id_transform
            .get(&mesh_id)
            .map(|m| m.has_normals)
            .unwrap_or(false)
    }

    // -----------------------------------------------------------------------
    // IncrementMeshIDs — port of C++ Manifold::Impl::IncrementMeshIDs()
    // -----------------------------------------------------------------------

    /// Allocates fresh unique mesh IDs and remaps all triRef.meshID values.
    /// This ensures boolean results don't collide with source mesh IDs.
    pub fn increment_mesh_ids(&mut self) {
        use std::collections::{BTreeMap, HashMap};

        // Build old -> new ID mapping. Iteration order determines which old
        // ID gets which fresh ID, so it must be sorted like C++ std::map.
        let old_transforms: BTreeMap<i32, Relation> =
            std::mem::take(&mut self.mesh_relation.mesh_id_transform);
        let num_mesh_ids = old_transforms.len() as u32;
        if num_mesh_ids == 0 {
            return;
        }
        let mut next_mesh_id = reserve_ids(num_mesh_ids) as i32;
        let mut old2new: HashMap<i32, i32> = HashMap::new();
        for (old_id, relation) in old_transforms {
            old2new.insert(old_id, next_mesh_id);
            self.mesh_relation
                .mesh_id_transform
                .insert(next_mesh_id, relation);
            next_mesh_id += 1;
        }

        // Update all triRef.meshID
        for tri_ref in &mut self.mesh_relation.tri_ref {
            if let Some(&new_id) = old2new.get(&tri_ref.mesh_id) {
                tri_ref.mesh_id = new_id;
            }
        }
    }

    // -----------------------------------------------------------------------
    // DedupePropVerts — port of C++ Manifold::Impl::DedupePropVerts()
    // -----------------------------------------------------------------------

    /// Deduplicates property vertices that share identical property values
    /// across paired halfedges within the same mesh.
    pub fn dedupe_prop_verts(&mut self) {
        let num_prop = self.num_prop;
        if num_prop == 0 {
            return;
        }

        let n_edges = self.halfedge.len();
        // Collect (prop0, prop1) pairs for edges where properties match
        let mut vert2vert: Vec<(i32, i32)> = vec![(-1, -1); n_edges];
        for edge_idx in 0..n_edges {
            let edge = self.halfedge[edge_idx];
            if edge.paired_halfedge < 0 {
                continue;
            }
            let edge_face = edge_idx / 3;
            let pair_face = edge.paired_halfedge as usize / 3;

            if self.mesh_relation.tri_ref[edge_face].mesh_id
                != self.mesh_relation.tri_ref[pair_face].mesh_id
            {
                continue;
            }

            let prop0 = self.halfedge[edge_idx].prop_vert;
            let prop1 = self.halfedge[next_halfedge(edge.paired_halfedge) as usize].prop_vert;
            if prop0 < 0 || prop1 < 0 {
                continue;
            }

            let mut prop_equal = true;
            for p in 0..num_prop {
                let idx0 = num_prop * prop0 as usize + p;
                let idx1 = num_prop * prop1 as usize + p;
                if idx0 >= self.properties.len() || idx1 >= self.properties.len() {
                    prop_equal = false;
                    break;
                }
                if self.properties[idx0] != self.properties[idx1] {
                    prop_equal = false;
                    break;
                }
            }
            if prop_equal {
                vert2vert[edge_idx] = (prop0, prop1);
            }
        }

        // Use union-find to merge equivalent property vertices
        let num_prop_vert = self.num_prop_vert();
        let ds = crate::disjoint_sets::DisjointSets::new(num_prop_vert as u32);
        for &(a, b) in &vert2vert {
            if a >= 0 && b >= 0 {
                ds.unite(a as u32, b as u32);
            }
        }
        let mut vert_labels = Vec::new();
        let num_labels = ds.connected_components(&mut vert_labels);

        // Build label -> canonical vert mapping
        let mut label2vert = vec![0i32; num_labels as usize];
        for v in 0..num_prop_vert {
            label2vert[vert_labels[v] as usize] = v as i32;
        }

        // Remap all prop_vert indices
        for edge in &mut self.halfedge {
            if edge.prop_vert >= 0 && (edge.prop_vert as usize) < num_prop_vert {
                edge.prop_vert = label2vert[vert_labels[edge.prop_vert as usize] as usize];
            }
        }
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
}

// ---------------------------------------------------------------------------
#[cfg(test)]
#[path = "impl_mesh_tests.rs"]
mod tests;
