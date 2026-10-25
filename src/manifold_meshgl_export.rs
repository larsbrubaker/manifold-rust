// manifold_meshgl_export.rs — MeshGLP export: the generic body behind
// `Manifold::get_mesh_gl`, `get_mesh_gl64` and `get_mesh_glp`.
//
// Split out of manifold_meshgl.rs (which keeps the import path and the public
// entry points) to keep both files under the size limit. It is a child module
// of manifold.rs so it can read the private `imp` field. The body ports the
// C++ `GetMeshGLImpl<Precision, I>` template: the f32 instantiation narrows
// and floors the exported tolerance at `f32::EPSILON * bbox.scale()` exactly
// where the C++ float instantiation does, and the f64 instantiation copies the
// kernel's f64 data out unchanged.

use super::Manifold;
use crate::linalg::Vec3;
use crate::types::{MeshGLP, MeshIndex, MeshPrecision};

impl Manifold {
    /// Generic export shared by every `get_mesh_gl*` entry point; mirrors
    /// C++ `GetMeshGLImpl<Precision, I>`.
    pub(super) fn get_mesh_gl_impl<P: MeshPrecision, I: MeshIndex>(
        &self,
        normal_idx: i32,
    ) -> MeshGLP<P, I> {
        // Per #1718: GetMeshGL(-1) auto-substitutes slot 0 when CalculateNormals
        // recorded world-frame normals on every meshID. A non-negative idx is the
        // legacy interface; >= 0 means "interpret this extra-prop slot as normals
        // and normalize it on export (transforming legacy per-mesh-frame runs)".
        let mut normal_idx = normal_idx;
        if normal_idx < 0 && self.imp.all_have_normals() {
            normal_idx = 0;
        }
        let extra_prop = self.imp.num_prop;
        let update_normals = normal_idx >= 0 && (normal_idx as usize) + 3 <= extra_prop;
        let mut out = MeshGLP::<P, I>::default();
        let num_tri = self.imp.num_tri();
        out.num_prop = I::from_usize(3 + extra_prop);
        // C++: for float output only, floor tolerance at float_epsilon *
        // bBox.Scale() to avoid catastrophic cancellation in RelatedGL checks.
        // The f64 output exports the kernel tolerance unchanged.
        let mut tolerance = self.imp.tolerance;
        if P::IS_SINGLE {
            tolerance = tolerance.max((f32::EPSILON as f64) * self.imp.bbox.scale());
        }
        out.tolerance = P::from_f64(tolerance);

        if !self.imp.halfedge_tangent.is_empty() {
            for t in &self.imp.halfedge_tangent {
                out.halfedge_tangent.extend([
                    P::from_f64(t.x),
                    P::from_f64(t.y),
                    P::from_f64(t.z),
                    P::from_f64(t.w),
                ]);
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
                ra.original_id
                    .cmp(&rb.original_id)
                    .then(ra.mesh_id.cmp(&rb.mesh_id))
            });
        }

        out.tri_verts.resize(3 * num_tri, I::default());
        out.face_id.resize(num_tri, I::default());
        let mut mesh_id_transform = self.imp.mesh_relation.mesh_id_transform.clone();
        let mut last_mesh_id = -1i32;
        // run index that each output triangle belongs to (for per-vertex normal
        // export below). -1 until the first run is pushed.
        let mut tri_run = vec![0usize; num_tri];
        let mut current_run: i32 = -1;

        // run_flags layout (#1718): bit 0 = backside, bit 1 = hasNormals (slot 0..2
        // of the extra properties is world-frame normals; consumers skip
        // re-applying run_transform to it).
        let run_flags_for = |rel: &crate::types::Relation| -> u8 {
            (if rel.back_side { 1u8 } else { 0 }) | (if rel.has_normals { 2u8 } else { 0 })
        };

        for new_tri in 0..num_tri {
            let old_tri = tri_new2old[new_tri];
            if old_tri < tri_ref.len() {
                let r = &tri_ref[old_tri];
                out.face_id[new_tri] = I::from_i32(if r.face_id >= 0 {
                    r.face_id
                } else {
                    r.coplanar_id
                });
            }
            for i in 0..3 {
                let he = &self.imp.halfedge[3 * old_tri + i];
                // First pass: set to geometric vertex (startVert). When numProp > 0,
                // a second pass below will replace these with property-vertex indices.
                out.tri_verts[3 * new_tri + i] = I::from_i32(he.start_vert);
            }
            let mesh_id = if old_tri < tri_ref.len() {
                tri_ref[old_tri].mesh_id
            } else {
                0
            };
            if mesh_id != last_mesh_id {
                let rel = mesh_id_transform.remove(&mesh_id).unwrap_or_default();
                out.run_index.push(I::from_usize(3 * new_tri));
                out.run_original_id.push(rel.original_id.max(0) as u32);
                out.run_flags.push(run_flags_for(&rel));
                // C++: runTransform only emitted for non-original manifolds
                if !is_original {
                    for col in 0..4 {
                        for row in 0..3 {
                            out.run_transform.push(P::from_f64(rel.transform[col][row]));
                        }
                    }
                }
                last_mesh_id = mesh_id;
                current_run += 1;
            }
            tri_run[new_tri] = current_run.max(0) as usize;
        }
        // Add runs for originals that contributed no tris
        for (_id, rel) in &mesh_id_transform {
            out.run_index.push(I::from_usize(3 * num_tri));
            out.run_original_id.push(rel.original_id.max(0) as u32);
            out.run_flags.push(run_flags_for(rel));
            if !is_original {
                for col in 0..4 {
                    for row in 0..3 {
                        out.run_transform.push(P::from_f64(rel.transform[col][row]));
                    }
                }
            }
        }
        out.run_index.push(I::from_usize(3 * num_tri));

        let num_geom_vert = self.imp.num_vert();
        let num_prop = self.imp.num_prop;
        let total_prop = 3 + num_prop; // xyz + extra

        if num_prop == 0 {
            // No extra properties: positions only, indexed by geometric vertex
            out.vert_properties.resize(3 * num_geom_vert, P::default());
            for i in 0..num_geom_vert {
                let v = self.imp.vert_pos[i];
                out.vert_properties[3 * i] = P::from_f64(v.x);
                out.vert_properties[3 * i + 1] = P::from_f64(v.y);
                out.vert_properties[3 * i + 2] = P::from_f64(v.z);
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
                if he.start_vert < 0 {
                    continue;
                }
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
                    out.vert_properties.push(P::from_f64(pos.x));
                    out.vert_properties.push(P::from_f64(pos.y));
                    out.vert_properties.push(P::from_f64(pos.z));
                    // Write extra properties (zeros if prop_vert is invalid)
                    if prop >= 0 {
                        let base = prop as usize * num_prop;
                        for p in 0..num_prop {
                            out.vert_properties
                                .push(P::from_f64(self.imp.properties[base + p]));
                        }
                    } else {
                        for _ in 0..num_prop {
                            out.vert_properties.push(P::default());
                        }
                    }

                    // Per #1718: normalize the requested normal slot on export. Runs
                    // that already carry world-frame normals (hasNormals bit) just
                    // get normalized; legacy runs without the bit are interpreted as
                    // per-mesh-frame and rotated to world via the inverse-frame
                    // transform first.
                    if update_normals {
                        let ni = normal_idx as usize;
                        let off = out.vert_properties.len() - num_prop + ni;
                        // Read back the just-written values: in the f32
                        // instantiation this round-trips through f32 exactly as
                        // the C++ float template does before normalizing.
                        let mut n = Vec3::new(
                            out.vert_properties[off].to_f64(),
                            out.vert_properties[off + 1].to_f64(),
                            out.vert_properties[off + 2].to_f64(),
                        );
                        let run = tri_run[new_tri];
                        let run_has_n = !is_original && (out.run_flags[run] & 2) != 0;
                        if !is_original && !run_has_n {
                            n = normal_transform_for_run(
                                &out.run_transform,
                                run,
                                out.run_flags[run],
                            ) * n;
                        }
                        n = crate::smoothing::safe_normalize(n);
                        out.vert_properties[off] = P::from_f64(n.x);
                        out.vert_properties[off + 1] = P::from_f64(n.y);
                        out.vert_properties[off + 2] = P::from_f64(n.z);
                    }

                    vert_prop_pairs[vert].push((prop, idx));

                    // First slot for this geometric vertex is the canonical merge target.
                    // Additional slots get merge entries so from_mesh_gl knows they
                    // are coincident with the first slot.
                    if vert2idx[vert] == -1 {
                        vert2idx[vert] = idx as i32;
                    } else {
                        out.merge_from_vert.push(I::from_usize(idx as usize));
                        out.merge_to_vert.push(I::from_i32(vert2idx[vert]));
                    }
                    idx
                };

                out.tri_verts[3 * new_tri + i] = I::from_usize(idx as usize);
            }
        }
        out
    }
} // impl Manifold

/// Per-run normal transform for the legacy export path (slot interpreted as
/// per-mesh-frame normals on a run without the hasNormals bit). Reconstructs
/// `NormalTransform(runTransform) * (backside ? -1 : 1)` from the flat
/// column-major run_transform buffer. Matches C++ GetMeshGLImpl. (#1718)
fn normal_transform_for_run<P: MeshPrecision>(
    run_transform: &[P],
    run: usize,
    flags: u8,
) -> crate::linalg::Mat3 {
    use crate::linalg::{Mat3, Vec3};
    let b = 12 * run;
    // run_transform[b + 3*col + row] = transform[col][row]; cols 0..2 = 3x3 part.
    let col = |j: usize| {
        Vec3::new(
            run_transform[b + 3 * j].to_f64(),
            run_transform[b + 3 * j + 1].to_f64(),
            run_transform[b + 3 * j + 2].to_f64(),
        )
    };
    let m3 = Mat3::from_cols(col(0), col(1), col(2));
    let sign = if (flags & 1) != 0 { -1.0 } else { 1.0 };
    // NormalTransform(M) = inverse(transpose(M)).
    m3.transpose().inverse() * sign
}
