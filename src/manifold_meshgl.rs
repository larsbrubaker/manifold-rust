// MeshGL conversion methods for Manifold
// Extracted from manifold.rs for file size management

use crate::impl_mesh::ManifoldImpl;
use crate::linalg::{IVec3, Mat3x4, Vec3};
use crate::types::{MeshGL, MeshGL64};
use super::Manifold;

impl Manifold {
    pub fn from_mesh_gl(mesh: &MeshGL) -> Self {
    let mut imp = ManifoldImpl::new();
    let num_prop = mesh.num_prop as usize;
    imp.num_prop = num_prop.saturating_sub(3);
    let num_tri = mesh.num_tri();

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

    let tri: Vec<IVec3> = (0..num_tri)
        .map(|i| {
            let t = mesh.get_tri_verts(i);
            IVec3::new(t[0] as i32, t[1] as i32, t[2] as i32)
        })
        .collect();
    imp.create_halfedges(&tri, &[]);

    // Set up mesh relations from runOriginalID (matches C++ MeshGL constructor)
    let run_index: Vec<usize> = if mesh.run_index.is_empty() {
        vec![0, 3 * num_tri]
    } else {
        mesh.run_index.iter().map(|&v| v as usize).collect()
    };
    let mut run_original_id: Vec<u32> = mesh.run_original_id.clone();
    let num_runs = run_original_id.len().max(1);
    let start_id = crate::impl_mesh::reserve_ids(num_runs as u32) as i32;
    if run_original_id.is_empty() {
        run_original_id.push(start_id as u32);
    }

    imp.mesh_relation.tri_ref.resize(num_tri, crate::types::TriRef::default());
    for (i, &orig_id) in run_original_id.iter().enumerate() {
        let mesh_id = start_id + i as i32;
        let run_start = if i < run_index.len() { run_index[i] / 3 } else { num_tri };
        let run_end = if i + 1 < run_index.len() { run_index[i + 1] / 3 } else { num_tri };
        for tri in run_start..run_end {
            if tri < imp.mesh_relation.tri_ref.len() {
                imp.mesh_relation.tri_ref[tri].mesh_id = mesh_id;
                imp.mesh_relation.tri_ref[tri].original_id = orig_id as i32;
                imp.mesh_relation.tri_ref[tri].face_id = if !mesh.face_id.is_empty() && tri < mesh.face_id.len() {
                    mesh.face_id[tri] as i32
                } else {
                    -1
                };
                imp.mesh_relation.tri_ref[tri].coplanar_id = tri as i32;
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
    // A Manifold created from input mesh is never an original
    imp.mesh_relation.original_id = -1;

    imp.calculate_bbox();
    imp.set_epsilon(mesh.tolerance as f64, false);
    imp.sort_geometry();
    imp.set_normals_and_coplanar();
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

    // Vertex properties
    let prop_rows = if self.imp.num_prop == 0 { self.imp.num_vert() } else { self.imp.num_prop_vert() };
    let mut prop_pos = vec![Vec3::splat(0.0); prop_rows];
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
