// face_op.rs — Phase 7a: Face normals, coplanarity, vertex normals
//
// Ports src/face_op.cpp, the face-normal and coplanarity portions of
// src/impl.cpp (SetNormalsAndCoplanar, CalculateVertNormals), and the
// GetAxisAlignedProjection utility from src/shared.h.

use crate::linalg::{Vec2, Vec3, cross, dot, normalize, length2};
use crate::types::next_halfedge;
use crate::impl_mesh::ManifoldImpl;

// -----------------------------------------------------------------------
// Proj2x3 — 2-row, 3-column projection matrix (drops one axis)
//
// Used to project 3D mesh positions onto a 2D plane for CCW tests and
// triangulation. Mirrors `mat2x3` in the C++ linalg.h library.
// -----------------------------------------------------------------------

/// A 2×3 projection matrix: maps Vec3 → Vec2 via dot products with two rows.
#[derive(Clone, Copy, Debug)]
pub struct Proj2x3 {
    pub row0: Vec3,
    pub row1: Vec3,
}

impl Proj2x3 {
    /// Apply the projection: `[dot(row0, v), dot(row1, v)]`.
    #[inline]
    pub fn apply(&self, v: Vec3) -> Vec2 {
        Vec2::new(dot(self.row0, v), dot(self.row1, v))
    }
}

// -----------------------------------------------------------------------
// GetAxisAlignedProjection
// -----------------------------------------------------------------------

/// Returns a projection matrix that drops the largest-magnitude axis of
/// `normal`, producing a 2D view aligned with the face plane.
///
/// Mirrors `GetAxisAlignedProjection` in `src/shared.h`.
pub fn get_axis_aligned_projection(normal: Vec3) -> Proj2x3 {
    let abs = Vec3::new(normal.x.abs(), normal.y.abs(), normal.z.abs());

    // mat3x2 columns (each col is a Vec3); transposed to get mat2x3 rows.
    let (row0, row1, xyz_max) = if abs.z > abs.x && abs.z > abs.y {
        // Drop Z, keep X and Y
        (Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0), normal.z)
    } else if abs.y > abs.x {
        // Drop Y, keep Z and X
        (Vec3::new(0.0, 0.0, 1.0), Vec3::new(1.0, 0.0, 0.0), normal.y)
    } else {
        // Drop X, keep Y and Z
        (Vec3::new(0.0, 1.0, 0.0), Vec3::new(0.0, 0.0, 1.0), normal.x)
    };

    // If the dominant axis is negative, flip the first row so that the
    // projected winding order is consistent.
    if xyz_max < 0.0 {
        Proj2x3 { row0: Vec3::new(-row0.x, -row0.y, -row0.z), row1 }
    } else {
        Proj2x3 { row0, row1 }
    }
}

// -----------------------------------------------------------------------
// SetNormalsAndCoplanar
// -----------------------------------------------------------------------

/// Computes face normals and flood-fills coplanar face groups, then
/// calls `calculate_vert_normals` to compute per-vertex normals.
///
/// Mirrors `Manifold::Impl::SetNormalsAndCoplanar()` in `src/impl.cpp`.
pub fn set_normals_and_coplanar(mesh: &mut ManifoldImpl) {
    let num_tri = mesh.num_tri();
    mesh.face_normal.resize(num_tri, Vec3::new(0.0, 0.0, 1.0));

    // Struct for sorting triangles by area
    struct TriPriority {
        area2: f64,
        tri: usize,
    }

    // Compute face normals and priorities (sort largest faces first)
    let mut tri_priority: Vec<TriPriority> = (0..num_tri)
        .map(|tri| {
            // Mark coplanarID as unset
            if tri < mesh.mesh_relation.tri_ref.len() {
                mesh.mesh_relation.tri_ref[tri].coplanar_id = -1;
            }

            if mesh.halfedge[3 * tri].start_vert < 0 {
                return TriPriority { area2: 0.0, tri };
            }
            let v = mesh.vert_pos[mesh.halfedge[3 * tri].start_vert as usize];
            let n = cross(
                mesh.vert_pos[mesh.halfedge[3 * tri].end_vert as usize] - v,
                mesh.vert_pos[mesh.halfedge[3 * tri + 1].end_vert as usize] - v,
            );
            let normal = normalize(n);
            mesh.face_normal[tri] = if normal.x.is_nan() {
                Vec3::new(0.0, 0.0, 1.0)
            } else {
                normal
            };
            TriPriority { area2: length2(n), tri }
        })
        .collect();

    // Sort by area descending (largest triangles first → better coplanar seeds)
    tri_priority.sort_by(|a, b| b.area2.partial_cmp(&a.area2).unwrap_or(std::cmp::Ordering::Equal));

    // Flood-fill coplanar groups from each unassigned face
    let mut interior_halfedges: Vec<usize> = Vec::new();
    for tp in &tri_priority {
        let tri = tp.tri;
        if tri >= mesh.mesh_relation.tri_ref.len() {
            continue;
        }
        if mesh.mesh_relation.tri_ref[tri].coplanar_id >= 0 {
            continue;
        }

        mesh.mesh_relation.tri_ref[tri].coplanar_id = tri as i32;
        if mesh.halfedge[3 * tri].start_vert < 0 {
            continue;
        }

        let base = mesh.vert_pos[mesh.halfedge[3 * tri].start_vert as usize];
        let normal = mesh.face_normal[tri];

        interior_halfedges.clear();
        interior_halfedges.push(3 * tri);
        interior_halfedges.push(3 * tri + 1);
        interior_halfedges.push(3 * tri + 2);

        while let Some(h) = interior_halfedges.pop() {
            let paired = mesh.halfedge[h].paired_halfedge;
            if paired < 0 {
                continue;
            }
            let h2 = next_halfedge(paired) as usize;
            let h2_tri = h2 / 3;
            if h2_tri >= mesh.mesh_relation.tri_ref.len() {
                continue;
            }
            if mesh.mesh_relation.tri_ref[h2_tri].coplanar_id >= 0 {
                continue;
            }

            let v = mesh.vert_pos[mesh.halfedge[h2].end_vert as usize];
            if (dot(v - base, normal)).abs() < mesh.tolerance {
                mesh.mesh_relation.tri_ref[h2_tri].coplanar_id = tri as i32;
                mesh.face_normal[h2_tri] = normal;

                // Avoid re-pushing paired interior halfedges (cancel out)
                let last = interior_halfedges.last().copied();
                if last == Some(mesh.halfedge[h2].paired_halfedge as usize) {
                    interior_halfedges.pop();
                } else {
                    interior_halfedges.push(h2 as usize);
                }
                interior_halfedges.push(next_halfedge(h2 as i32) as usize);
            }
        }
    }

    calculate_vert_normals(mesh);
}

// -----------------------------------------------------------------------
// CalculateVertNormals
// -----------------------------------------------------------------------

/// Computes per-vertex normals as angle-weighted averages of incident face normals.
///
/// Mirrors `Manifold::Impl::CalculateVertNormals()` in `src/impl.cpp`.
pub fn calculate_vert_normals(mesh: &mut ManifoldImpl) {
    let num_vert = mesh.vert_pos.len();
    mesh.vert_normal.resize(num_vert, Vec3::new(0.0, 0.0, 0.0));

    // For each vertex, find the first halfedge that starts there
    let mut vert_first_edge = vec![i32::MAX; num_vert];
    for (i, edge) in mesh.halfedge.iter().enumerate() {
        let sv = edge.start_vert;
        if sv >= 0 && (sv as usize) < num_vert {
            let sv = sv as usize;
            if (i as i32) < vert_first_edge[sv] {
                vert_first_edge[sv] = i as i32;
            }
        }
    }

    for vert in 0..num_vert {
        let first_edge = vert_first_edge[vert];
        if first_edge == i32::MAX {
            mesh.vert_normal[vert] = Vec3::new(0.0, 0.0, 0.0);
            continue;
        }

        let mut normal = Vec3::new(0.0, 0.0, 0.0);
        let halfedge = &mesh.halfedge;
        let vert_pos = &mesh.vert_pos;
        let face_normal = &mesh.face_normal;

        // ForVert equivalent: walk CW around the vertex
        let mut current = first_edge as usize;
        loop {
            let h = &halfedge[current];
            let tri_verts = [
                h.start_vert as usize,
                h.end_vert as usize,
                halfedge[next_halfedge(current as i32) as usize].end_vert as usize,
            ];

            // Avoid degenerate triangles
            if tri_verts[0] < vert_pos.len()
                && tri_verts[1] < vert_pos.len()
                && tri_verts[2] < vert_pos.len()
            {
                let curr_edge_dir = vert_pos[tri_verts[1]] - vert_pos[tri_verts[0]];
                let prev_edge_dir = vert_pos[tri_verts[0]] - vert_pos[tri_verts[2]];
                let curr_len = length2(curr_edge_dir).sqrt();
                let prev_len = length2(prev_edge_dir).sqrt();

                if curr_len > 0.0 && prev_len > 0.0 {
                    let curr_norm = curr_edge_dir / curr_len;
                    let prev_norm = prev_edge_dir / prev_len;
                    if curr_norm.x.is_finite() && prev_norm.x.is_finite() {
                        let d = dot(prev_norm, curr_norm).clamp(-1.0, 1.0);
                        // Negate because prevEdge points into vert and currEdge points away
                        let phi = (-d).acos();
                        if phi.is_finite() && current / 3 < face_normal.len() {
                            normal = normal + face_normal[current / 3] * phi;
                        }
                    }
                }
            }

            // ForVert step: next_halfedge(halfedge[current].pairedHalfedge)
            let paired = halfedge[current].paired_halfedge;
            if paired < 0 {
                break;
            }
            current = next_halfedge(paired) as usize;
            if current == first_edge as usize {
                break;
            }
        }

        let len = length2(normal).sqrt();
        mesh.vert_normal[vert] = if len > 0.0 { normal / len } else { Vec3::new(0.0, 0.0, 0.0) };
    }
}

// -----------------------------------------------------------------------
// Tests
// -----------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linalg::Mat3x4;
    use crate::impl_mesh::ManifoldImpl;

    #[test]
    fn test_get_axis_aligned_projection_z() {
        // Normal primarily in Z: should project to XY plane
        let proj = get_axis_aligned_projection(Vec3::new(0.0, 0.0, 1.0));
        let v = Vec3::new(3.0, 4.0, 5.0);
        let p = proj.apply(v);
        assert!((p.x - 3.0).abs() < 1e-12);
        assert!((p.y - 4.0).abs() < 1e-12);
    }

    #[test]
    fn test_get_axis_aligned_projection_y() {
        // Normal primarily in Y: should project to ZX plane
        let proj = get_axis_aligned_projection(Vec3::new(0.0, 1.0, 0.0));
        let v = Vec3::new(3.0, 4.0, 5.0);
        let p = proj.apply(v);
        assert!((p.x - 5.0).abs() < 1e-12, "expected z=5, got {}", p.x);
        assert!((p.y - 3.0).abs() < 1e-12, "expected x=3, got {}", p.y);
    }

    #[test]
    fn test_get_axis_aligned_projection_x() {
        // Normal primarily in X: should project to YZ plane
        let proj = get_axis_aligned_projection(Vec3::new(1.0, 0.0, 0.0));
        let v = Vec3::new(3.0, 4.0, 5.0);
        let p = proj.apply(v);
        assert!((p.x - 4.0).abs() < 1e-12, "expected y=4, got {}", p.x);
        assert!((p.y - 5.0).abs() < 1e-12, "expected z=5, got {}", p.y);
    }

    #[test]
    fn test_get_axis_aligned_projection_negative_z() {
        // Normal primarily in -Z: row0 should be flipped
        let proj = get_axis_aligned_projection(Vec3::new(0.0, 0.0, -1.0));
        let v = Vec3::new(3.0, 4.0, 5.0);
        let p = proj.apply(v);
        // Flipped first row: x → -x
        assert!((p.x + 3.0).abs() < 1e-12, "expected -x=-3, got {}", p.x);
        assert!((p.y - 4.0).abs() < 1e-12);
    }

    #[test]
    fn test_set_normals_tetrahedron() {
        let mut m = ManifoldImpl::tetrahedron(&Mat3x4::identity());
        set_normals_and_coplanar(&mut m);
        // Every face normal should be unit length
        for n in &m.face_normal {
            let len = length2(*n).sqrt();
            assert!((len - 1.0).abs() < 1e-10, "face normal not unit: len={}", len);
        }
        // Every vert normal should be nonzero (tetrahedron has no degenerate verts)
        for n in &m.vert_normal {
            let len = length2(*n).sqrt();
            assert!(len > 0.0, "vert normal is zero");
        }
    }

    #[test]
    fn test_set_normals_cube() {
        let mut m = ManifoldImpl::cube(&Mat3x4::identity());
        set_normals_and_coplanar(&mut m);
        assert_eq!(m.face_normal.len(), m.num_tri());
        // Coplanar IDs should be assigned
        let any_coplanar_id_set = m.mesh_relation.tri_ref.iter()
            .any(|r| r.coplanar_id >= 0);
        assert!(any_coplanar_id_set, "no coplanar IDs were set");
    }
}
