// MeshGL conversion methods for Manifold
// Extracted from manifold.rs for file size management

use crate::impl_mesh::ManifoldImpl;
use crate::linalg::{IVec3, Mat3x4, Vec3};
use crate::types::{MeshGL, MeshGL64};
use super::Manifold;

impl Manifold {
    pub fn from_mesh_gl(mesh: &MeshGL) -> Self {
    let num_vert = mesh.num_vert() as u32;
    let num_tri = mesh.num_tri();

    // Validation checks matching C++ Impl::Impl(const MeshGLP&)
    if num_vert == 0 && num_tri == 0 {
        return Self::make_empty(crate::types::Error::NoError);
    }
    if num_vert < 4 || num_tri < 4 {
        return Self::make_empty(crate::types::Error::NotManifold);
    }
    if mesh.num_prop < 3 {
        return Self::make_empty(crate::types::Error::MissingPositionProperties);
    }
    if mesh.merge_from_vert.len() != mesh.merge_to_vert.len() {
        return Self::make_empty(crate::types::Error::MergeVectorsDifferentLengths);
    }
    if !mesh.run_transform.is_empty()
        && 12 * mesh.run_original_id.len() != mesh.run_transform.len()
    {
        return Self::make_empty(crate::types::Error::TransformWrongLength);
    }
    if !mesh.run_original_id.is_empty()
        && !mesh.run_index.is_empty()
        && mesh.run_original_id.len() + 1 != mesh.run_index.len()
        && mesh.run_original_id.len() != mesh.run_index.len()
    {
        return Self::make_empty(crate::types::Error::RunIndexWrongLength);
    }
    if !mesh.face_id.is_empty() && mesh.face_id.len() != num_tri {
        return Self::make_empty(crate::types::Error::FaceIdWrongLength);
    }
    if !mesh.vert_properties.iter().all(|x| x.is_finite()) {
        return Self::make_empty(crate::types::Error::NonFiniteVertex);
    }
    if !mesh.run_transform.iter().all(|x| x.is_finite()) {
        return Self::make_empty(crate::types::Error::InvalidConstruction);
    }
    if !mesh.halfedge_tangent.iter().all(|x| x.is_finite()) {
        return Self::make_empty(crate::types::Error::InvalidConstruction);
    }

    // Check merge indices are in bounds
    for i in 0..mesh.merge_from_vert.len() {
        if mesh.merge_from_vert[i] >= num_vert || mesh.merge_to_vert[i] >= num_vert {
            return Self::make_empty(crate::types::Error::MergeIndexOutOfBounds);
        }
    }

    // Check tri_verts are in bounds
    for &v in &mesh.tri_verts {
        if v >= num_vert {
            return Self::make_empty(crate::types::Error::VertexOutOfBounds);
        }
    }

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

    // Build prop2vert mapping from merge vectors
    let has_merges = !mesh.merge_from_vert.is_empty();
    let needs_prop_map = imp.num_prop > 0 && has_merges;
    let prop2vert: Vec<i32> = if has_merges {
        let mut p2v: Vec<i32> = (0..num_vert as i32).collect();
        for i in 0..mesh.merge_from_vert.len() {
            p2v[mesh.merge_from_vert[i] as usize] = mesh.merge_to_vert[i] as i32;
        }
        p2v
    } else {
        vec![]
    };

    // Set up mesh relations from runOriginalID (matches C++ MeshGL constructor)
    let mut run_index: Vec<usize> = if mesh.run_index.is_empty() {
        vec![0, 3 * num_tri]
    } else {
        let mut ri: Vec<usize> = mesh.run_index.iter().map(|&v| v as usize).collect();
        if ri.len() == mesh.run_original_id.len() {
            ri.push(3 * num_tri);
        } else if ri.len() == 1 {
            ri.push(3 * num_tri);
        }
        ri
    };
    let mut run_original_id: Vec<u32> = mesh.run_original_id.clone();
    let num_runs = run_original_id.len().max(1);
    let start_id = crate::impl_mesh::reserve_ids(num_runs as u32) as i32;
    if run_original_id.is_empty() {
        run_original_id.push(start_id as u32);
    }

    // Build tri_ref for all input tris
    let mut all_tri_ref: Vec<crate::types::TriRef> = vec![crate::types::TriRef::default(); num_tri];
    for (i, &orig_id) in run_original_id.iter().enumerate() {
        let mesh_id = start_id + i as i32;
        let run_start = if i < run_index.len() { run_index[i] / 3 } else { num_tri };
        let run_end = if i + 1 < run_index.len() { run_index[i + 1] / 3 } else { num_tri };
        for tri in run_start..run_end {
            if tri < all_tri_ref.len() {
                all_tri_ref[tri].mesh_id = mesh_id;
                all_tri_ref[tri].original_id = orig_id as i32;
                all_tri_ref[tri].face_id = if !mesh.face_id.is_empty() && tri < mesh.face_id.len() {
                    mesh.face_id[tri] as i32
                } else {
                    -1
                };
                all_tri_ref[tri].coplanar_id = tri as i32;
            }
        }
        let transform = if mesh.run_transform.len() >= (i + 1) * 12 {
            let m = &mesh.run_transform[i * 12..i * 12 + 12];
            Mat3x4::from_cols(
                Vec3::new(m[0] as f64, m[1] as f64, m[2] as f64),
                Vec3::new(m[3] as f64, m[4] as f64, m[5] as f64),
                Vec3::new(m[6] as f64, m[7] as f64, m[8] as f64),
                Vec3::new(m[9] as f64, m[10] as f64, m[11] as f64),
            )
        } else {
            Mat3x4::identity()
        };
        imp.mesh_relation.mesh_id_transform.insert(mesh_id, crate::types::Relation {
            original_id: orig_id as i32,
            transform,
            back_side: false,
        });
    }

    // Build triangles, filtering out degenerates from merges
    let mut tri_prop: Vec<IVec3> = Vec::with_capacity(num_tri);
    let mut tri_vert: Vec<IVec3> = Vec::new();
    if needs_prop_map {
        tri_vert.reserve(num_tri);
    }
    imp.mesh_relation.tri_ref.clear();
    imp.mesh_relation.tri_ref.reserve(num_tri);

    for i in 0..num_tri {
        let t = mesh.get_tri_verts(i);
        let tri_p = IVec3::new(t[0] as i32, t[1] as i32, t[2] as i32);
        let tri_v = if prop2vert.is_empty() {
            tri_p
        } else {
            IVec3::new(
                prop2vert[t[0] as usize],
                prop2vert[t[1] as usize],
                prop2vert[t[2] as usize],
            )
        };

        // Skip degenerate triangles (where merged verts collapse)
        if tri_v.x != tri_v.y && tri_v.y != tri_v.z && tri_v.z != tri_v.x {
            if needs_prop_map {
                tri_prop.push(tri_p);
                tri_vert.push(tri_v);
            } else {
                tri_prop.push(tri_v);
            }
            imp.mesh_relation.tri_ref.push(all_tri_ref[i]);
        }
    }

    imp.create_halfedges(&tri_prop, &tri_vert);

    // Import halfedge tangents from MeshGL (flat f32 array, 4 per halfedge)
    if !mesh.halfedge_tangent.is_empty() {
        let n_tangents = mesh.halfedge_tangent.len() / 4;
        imp.halfedge_tangent.resize(
            n_tangents,
            crate::linalg::Vec4::new(0.0, 0.0, 0.0, 0.0),
        );
        for i in 0..n_tangents {
            imp.halfedge_tangent[i] = crate::linalg::Vec4::new(
                mesh.halfedge_tangent[4 * i] as f64,
                mesh.halfedge_tangent[4 * i + 1] as f64,
                mesh.halfedge_tangent[4 * i + 2] as f64,
                mesh.halfedge_tangent[4 * i + 3] as f64,
            );
        }
    }

    if !imp.is_manifold() {
        return Self::make_empty(crate::types::Error::NotManifold);
    }

    // A Manifold created from input mesh is never an original
    imp.mesh_relation.original_id = -1;

    // C++ pipeline: DedupePropVerts, SetNormalsAndCoplanar,
    // RemoveDegenerates, RemoveUnreferencedVerts, SortGeometry
    // Note: CleanupTopology omitted — it would fix opposite-face meshes but
    // conflicts with is_manifold check ordering.
    imp.dedupe_prop_verts();
    imp.calculate_bbox();
    imp.set_epsilon(mesh.tolerance as f64, false);
    imp.set_normals_and_coplanar();
    crate::edge_op::remove_degenerates(&mut imp, 0);
    imp.remove_unreferenced_verts();
    imp.sort_geometry();
    Self { imp }
}

pub fn get_mesh_gl(&self, normal_idx: i32) -> MeshGL {
    let _ = normal_idx;
    let mut out = MeshGL::default();
    let num_tri = self.imp.num_tri();
    let extra_prop = self.imp.num_prop;
    out.num_prop = 3 + extra_prop as u32;
    out.tolerance = self.imp.tolerance as f32;

    if !self.imp.halfedge_tangent.is_empty() {
        for t in &self.imp.halfedge_tangent {
            out.halfedge_tangent.extend([t.x as f32, t.y as f32, t.z as f32, t.w as f32]);
        }
    }

    // Sort triangles by (originalID, meshID) for run grouping
    let is_original = self.imp.mesh_relation.original_id >= 0;
    let tri_ref = &self.imp.mesh_relation.tri_ref;
    let mut tri_new2old: Vec<usize> = (0..num_tri).collect();
    if !is_original && !tri_ref.is_empty() {
        tri_new2old.sort_by(|&a, &b| {
            let ra = &tri_ref[a];
            let rb = &tri_ref[b];
            ra.original_id.cmp(&rb.original_id)
                .then(ra.mesh_id.cmp(&rb.mesh_id))
        });
    }

    out.tri_verts.resize(3 * num_tri, 0);
    out.face_id.resize(num_tri, 0);
    let mut mesh_id_transform = self.imp.mesh_relation.mesh_id_transform.clone();
    let mut last_mesh_id = -1i32;

    for new_tri in 0..num_tri {
        let old_tri = tri_new2old[new_tri];
        if old_tri < tri_ref.len() {
            let r = &tri_ref[old_tri];
            out.face_id[new_tri] = (if r.face_id >= 0 { r.face_id } else { r.coplanar_id }) as u32;
        }
        for i in 0..3 {
            let he = &self.imp.halfedge[3 * old_tri + i];
            // First pass: set to geometric vertex (startVert). When numProp > 0,
            // a second pass below will replace these with property-vertex indices.
            out.tri_verts[3 * new_tri + i] = he.start_vert as u32;
        }
        let mesh_id = if old_tri < tri_ref.len() { tri_ref[old_tri].mesh_id } else { 0 };
        if mesh_id != last_mesh_id {
            let rel = mesh_id_transform.remove(&mesh_id).unwrap_or_default();
            out.run_index.push(3 * new_tri as u32);
            out.run_original_id.push(rel.original_id.max(0) as u32);
            for col in 0..4 {
                for row in 0..3 {
                    out.run_transform.push(rel.transform[col][row] as f32);
                }
            }
            last_mesh_id = mesh_id;
        }
    }
    // Add runs for originals that contributed no tris
    for (_id, rel) in &mesh_id_transform {
        out.run_index.push(3 * num_tri as u32);
        out.run_original_id.push(rel.original_id.max(0) as u32);
        for col in 0..4 {
            for row in 0..3 {
                out.run_transform.push(rel.transform[col][row] as f32);
            }
        }
    }
    out.run_index.push(3 * num_tri as u32);

    let num_geom_vert = self.imp.num_vert();
    let num_prop = self.imp.num_prop;
    let total_prop = 3 + num_prop; // xyz + extra

    if num_prop == 0 {
        // No extra properties: positions only, indexed by geometric vertex
        out.vert_properties.resize(3 * num_geom_vert, 0.0);
        for i in 0..num_geom_vert {
            let v = self.imp.vert_pos[i];
            out.vert_properties[3 * i] = v.x as f32;
            out.vert_properties[3 * i + 1] = v.y as f32;
            out.vert_properties[3 * i + 2] = v.z as f32;
        }
        return out;
    }

    // When properties exist: deduplicate (start_vert, prop_vert) pairs, matching
    // C++ GetMeshGLImpl. Each unique (vert, prop) pair gets its own slot in
    // vert_properties. Merge vectors record which slots share the same geometry.
    //
    // C++: vertPropPair[vert] = list of {prop, idx}; vert2idx[vert] = first idx
    let mut vert2idx = vec![-1i32; num_geom_vert];
    let mut vert_prop_pairs: Vec<Vec<(i32, u32)>> = vec![vec![]; num_geom_vert];

    for new_tri in 0..num_tri {
        let old_tri = tri_new2old[new_tri];
        for i in 0..3 {
            let he = &self.imp.halfedge[3 * old_tri + i];
            if he.start_vert < 0 { continue; }
            let vert = he.start_vert as usize;
            let prop = he.prop_vert;

            // Look for existing (vert, prop) pair
            let pairs = &vert_prop_pairs[vert];
            let mut found_idx: Option<u32> = None;
            for &(p, idx) in pairs {
                if p == prop {
                    found_idx = Some(idx);
                    break;
                }
            }

            let idx = if let Some(idx) = found_idx {
                idx
            } else {
                let idx = (out.vert_properties.len() / total_prop) as u32;
                // Write position
                let pos = self.imp.vert_pos[vert];
                out.vert_properties.push(pos.x as f32);
                out.vert_properties.push(pos.y as f32);
                out.vert_properties.push(pos.z as f32);
                // Write extra properties (zeros if prop_vert is invalid)
                if prop >= 0 {
                    let base = prop as usize * num_prop;
                    for p in 0..num_prop {
                        out.vert_properties.push(self.imp.properties[base + p] as f32);
                    }
                } else {
                    for _ in 0..num_prop { out.vert_properties.push(0.0); }
                }
                vert_prop_pairs[vert].push((prop, idx));

                // First slot for this geometric vertex is the canonical merge target.
                // Additional slots get merge entries so from_mesh_gl knows they
                // are coincident with the first slot.
                if vert2idx[vert] == -1 {
                    vert2idx[vert] = idx as i32;
                } else {
                    out.merge_from_vert.push(idx);
                    out.merge_to_vert.push(vert2idx[vert] as u32);
                }
                idx
            };

            out.tri_verts[3 * new_tri + i] = idx;
        }
    }
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
} // impl Manifold
