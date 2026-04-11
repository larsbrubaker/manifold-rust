// face_op.rs — Phase 7a: Face normals, coplanarity, vertex normals
//
// Ports src/face_op.cpp, the face-normal and coplanarity portions of
// src/impl.cpp (SetNormalsAndCoplanar, CalculateVertNormals), and the
// GetAxisAlignedProjection utility from src/shared.h.

use std::collections::BTreeMap;

use crate::linalg::{Vec2, Vec3, IVec3, cross, dot, normalize, length2};
use crate::math;
use crate::polygon::{ccw, triangulate_idx};
use crate::types::{next_halfedge, Halfedge, PolyVert, PolygonsIdx, TriRef};
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
                        let phi = math::acos(-d);
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
// GetBarycentric — barycentric coordinates of point in triangle
// -----------------------------------------------------------------------

/// Compute barycentric coordinates of `v` with respect to triangle `tri_pos`.
/// Returns [u, v, w] where vertex i has weight uvw[i].
/// Returns exact 1.0 for vertices within `tolerance` of a triangle vertex,
/// and exact 0.0 for points within tolerance of an edge.
///
/// Mirrors `GetBarycentric` in `src/shared.h`.
pub fn get_barycentric(v: Vec3, tri_pos: [Vec3; 3], tolerance: f64) -> Vec3 {
    let edges = [
        tri_pos[2] - tri_pos[1],
        tri_pos[0] - tri_pos[2],
        tri_pos[1] - tri_pos[0],
    ];
    let d2 = [
        dot(edges[0], edges[0]),
        dot(edges[1], edges[1]),
        dot(edges[2], edges[2]),
    ];
    let long_side = if d2[0] > d2[1] && d2[0] > d2[2] {
        0
    } else if d2[1] > d2[2] {
        1
    } else {
        2
    };
    let cross_p = cross(edges[0], edges[1]);
    let area2 = dot(cross_p, cross_p);
    let tol2 = tolerance * tolerance;

    let mut uvw = Vec3::splat(0.0);
    for i in 0..3 {
        let dv = v - tri_pos[i];
        if dot(dv, dv) < tol2 {
            uvw = Vec3::splat(0.0);
            match i {
                0 => uvw.x = 1.0,
                1 => uvw.y = 1.0,
                _ => uvw.z = 1.0,
            }
            return uvw;
        }
    }

    if d2[long_side] < tol2 {
        // Degenerate point
        return Vec3::new(1.0, 0.0, 0.0);
    } else if area2 > d2[long_side] * tol2 {
        // Triangle case
        for i in 0..3 {
            let j = (i + 1) % 3;
            let cross_pv = cross(edges[i], v - tri_pos[j]);
            let area2v = dot(cross_pv, cross_pv);
            let val = if area2v < d2[i] * tol2 {
                0.0
            } else {
                dot(cross_pv, cross_p)
            };
            match i {
                0 => uvw.x = val,
                1 => uvw.y = val,
                _ => uvw.z = val,
            }
        }
        let sum = uvw.x + uvw.y + uvw.z;
        uvw = uvw / sum;
        return uvw;
    } else {
        // Line case
        let next_v = (long_side + 1) % 3;
        let alpha = dot(v - tri_pos[next_v], edges[long_side]) / d2[long_side];
        uvw = Vec3::splat(0.0);
        let last_v = (next_v + 1) % 3;
        match next_v {
            0 => uvw.x = 1.0 - alpha,
            1 => uvw.y = 1.0 - alpha,
            _ => uvw.z = 1.0 - alpha,
        }
        match last_v {
            0 => uvw.x = alpha,
            1 => uvw.y = alpha,
            _ => uvw.z = alpha,
        }
        return uvw;
    }
}

// -----------------------------------------------------------------------
// AssembleHalfedges — group halfedges into polygon loops
// -----------------------------------------------------------------------

/// Given a slice of halfedges (from a polygonal face), group them into polygon
/// loops by following start_vert → end_vert chains. Returns a vec of polygon
/// loops, where each loop is a vec of halfedge indices (offset by
/// `start_halfedge_idx`).
///
/// Mirrors `AssembleHalfedges` in `src/face_op.cpp`.
pub fn assemble_halfedges(
    halfedges: &[Halfedge],
    start_halfedge_idx: i32,
) -> Vec<Vec<i32>> {
    // Build multimap: start_vert → local edge index
    let mut vert_edge: BTreeMap<i32, Vec<usize>> = BTreeMap::new();
    for (i, he) in halfedges.iter().enumerate() {
        vert_edge.entry(he.start_vert).or_default().push(i);
    }

    let mut polys: Vec<Vec<i32>> = Vec::new();
    let mut start_edge: usize = 0;
    let mut this_edge: usize = start_edge;

    loop {
        if this_edge == start_edge {
            // Find next unvisited edge
            if vert_edge.is_empty() {
                break;
            }
            let (&_vert, edges) = vert_edge.iter().next().unwrap();
            start_edge = edges[0];
            this_edge = start_edge;
            polys.push(Vec::new());
        }
        polys.last_mut().unwrap().push(start_halfedge_idx + this_edge as i32);
        let end_vert = halfedges[this_edge].end_vert;
        let edges = vert_edge.get_mut(&end_vert).expect("non-manifold edge");
        // Remove the first occurrence
        let pos = 0; // take first available
        this_edge = edges.remove(pos);
        if edges.is_empty() {
            vert_edge.remove(&end_vert);
        }
    }
    polys
}

/// Project polygon loops into 2D using projection matrix and vertex positions.
///
/// Mirrors `ProjectPolygons` in `src/face_op.cpp`.
pub fn project_polygons(
    polys: &[Vec<i32>],
    halfedge: &[Halfedge],
    vert_pos: &[Vec3],
    projection: &Proj2x3,
) -> PolygonsIdx {
    let mut polygons: PolygonsIdx = Vec::new();
    for poly in polys {
        let mut simple_poly = Vec::new();
        for &edge in poly {
            let vert = halfedge[edge as usize].start_vert;
            simple_poly.push(PolyVert {
                pos: projection.apply(vert_pos[vert as usize]),
                idx: edge,
            });
        }
        polygons.push(simple_poly);
    }
    polygons
}

// -----------------------------------------------------------------------
// Face2Tri — triangulate polygonal faces into triangle halfedges
// -----------------------------------------------------------------------

/// Triangulates the faces represented by `face_edge` (polygon boundaries in
/// `self.halfedge`) and `halfedge_ref` (per-halfedge TriRef). After this call,
/// `self.halfedge` is reorganized into proper triangles and `self.mesh_relation.tri_ref`
/// is populated.
///
/// Mirrors `Manifold::Impl::Face2Tri` in `src/face_op.cpp`.
pub fn face2tri(mesh: &mut ManifoldImpl, face_edge: &[i32], halfedge_ref: &[TriRef]) {
    let mut tri_verts: Vec<IVec3> = Vec::new();
    let mut tri_normal: Vec<Vec3> = Vec::new();
    let mut tri_ref: Vec<TriRef> = Vec::new();

    let num_faces = face_edge.len() - 1;

    for face in 0..num_faces {
        let first_edge = face_edge[face] as usize;
        let last_edge = face_edge[face + 1] as usize;
        let num_edge = last_edge - first_edge;
        if num_edge == 0 {
            continue;
        }
        debug_assert!(num_edge >= 3, "face has less than three edges");
        let normal = mesh.face_normal[face];
        let ref_tri = halfedge_ref[first_edge];

        if num_edge == 3 {
            // Single triangle — just sort edges into correct winding order
            let mut tri_edge = [first_edge as i32, first_edge as i32 + 1, first_edge as i32 + 2];
            let mut tri = [
                mesh.halfedge[first_edge].start_vert,
                mesh.halfedge[first_edge + 1].start_vert,
                mesh.halfedge[first_edge + 2].start_vert,
            ];
            let ends = [
                mesh.halfedge[first_edge].end_vert,
                mesh.halfedge[first_edge + 1].end_vert,
                mesh.halfedge[first_edge + 2].end_vert,
            ];
            if ends[0] == tri[2] {
                tri_edge.swap(1, 2);
                tri.swap(1, 2);
            }
            tri_verts.push(IVec3::new(tri[0], tri[1], tri[2]));
            tri_normal.push(normal);
            tri_ref.push(ref_tri);
        } else if num_edge == 4 {
            // Quad — split into two triangles
            let projection = get_axis_aligned_projection(normal);
            let tri_ccw = |t: [i32; 3]| -> bool {
                ccw(
                    projection.apply(mesh.vert_pos[mesh.halfedge[t[0] as usize].start_vert as usize]),
                    projection.apply(mesh.vert_pos[mesh.halfedge[t[1] as usize].start_vert as usize]),
                    projection.apply(mesh.vert_pos[mesh.halfedge[t[2] as usize].start_vert as usize]),
                    mesh.epsilon,
                ) >= 0
            };

            let quad = assemble_halfedges(
                &mesh.halfedge[first_edge..last_edge],
                first_edge as i32,
            );
            let quad = &quad[0]; // Should be exactly one loop

            let tris0 = [
                [quad[0], quad[1], quad[2]],
                [quad[0], quad[2], quad[3]],
            ];
            let tris1 = [
                [quad[1], quad[2], quad[3]],
                [quad[0], quad[1], quad[3]],
            ];

            let choice = if !(tri_ccw(tris0[0]) && tri_ccw(tris0[1])) {
                1
            } else if tri_ccw(tris1[0]) && tri_ccw(tris1[1]) {
                let diag0 = mesh.vert_pos[mesh.halfedge[quad[0] as usize].start_vert as usize]
                    - mesh.vert_pos[mesh.halfedge[quad[2] as usize].start_vert as usize];
                let diag1 = mesh.vert_pos[mesh.halfedge[quad[1] as usize].start_vert as usize]
                    - mesh.vert_pos[mesh.halfedge[quad[3] as usize].start_vert as usize];
                if length2(diag0) > length2(diag1) { 1 } else { 0 }
            } else {
                0
            };

            let chosen = if choice == 0 { &tris0 } else { &tris1 };
            for tri in chosen {
                tri_verts.push(IVec3::new(
                    mesh.halfedge[tri[0] as usize].start_vert,
                    mesh.halfedge[tri[1] as usize].start_vert,
                    mesh.halfedge[tri[2] as usize].start_vert,
                ));
                tri_normal.push(normal);
                tri_ref.push(ref_tri);
            }
        } else {
            // General polygon — use full triangulator
            let projection = get_axis_aligned_projection(normal);
            let polys_loops = assemble_halfedges(
                &mesh.halfedge[first_edge..last_edge],
                first_edge as i32,
            );
            let polys = project_polygons(&polys_loops, &mesh.halfedge, &mesh.vert_pos, &projection);
            let tris = triangulate_idx(&polys, mesh.epsilon, true);
            for tri in &tris {
                tri_verts.push(IVec3::new(
                    mesh.halfedge[tri.x as usize].start_vert,
                    mesh.halfedge[tri.y as usize].start_vert,
                    mesh.halfedge[tri.z as usize].start_vert,
                ));
                tri_normal.push(normal);
                tri_ref.push(ref_tri);
            }
        }
    }

    // Rebuild halfedge array from triangulated faces
    let num_tri = tri_verts.len();
    let mut new_halfedge = vec![
        Halfedge {
            start_vert: -1,
            end_vert: -1,
            paired_halfedge: -1,
            prop_vert: -1,
        };
        num_tri * 3
    ];

    for (t, tri) in tri_verts.iter().enumerate() {
        let verts = [tri.x, tri.y, tri.z];
        for j in 0..3 {
            new_halfedge[3 * t + j].start_vert = verts[j];
            new_halfedge[3 * t + j].end_vert = verts[(j + 1) % 3];
        }
    }

    // Pair up halfedges
    // Build map: (start, end) → halfedge index
    let mut edge_map: std::collections::HashMap<(i32, i32), usize> = std::collections::HashMap::new();
    for i in 0..new_halfedge.len() {
        let he = &new_halfedge[i];
        if he.start_vert < 0 {
            continue;
        }
        let key = (he.end_vert, he.start_vert); // look for opposite
        if let Some(&pair_idx) = edge_map.get(&key) {
            new_halfedge[i].paired_halfedge = pair_idx as i32;
            new_halfedge[pair_idx].paired_halfedge = i as i32;
            edge_map.remove(&key);
        } else {
            edge_map.insert((he.start_vert, he.end_vert), i);
        }
    }

    mesh.halfedge = new_halfedge;
    mesh.face_normal = tri_normal;
    mesh.mesh_relation.tri_ref = tri_ref;
}

// -----------------------------------------------------------------------
// ReorderHalfedges — canonical ordering within each triangle
// -----------------------------------------------------------------------

/// Reorders halfedges within each face so the one with the smallest start_vert
/// is first, then fixes paired_halfedge references.
///
/// Mirrors `Manifold::Impl::ReorderHalfedges` in `src/sort.cpp`.
pub fn reorder_halfedges(mesh: &mut ManifoldImpl) {
    let num_tri = mesh.halfedge.len() / 3;

    // Step 1: rotate each triangle so smallest start_vert is first
    for tri in 0..num_tri {
        let base = tri * 3;
        let face = [
            mesh.halfedge[base],
            mesh.halfedge[base + 1],
            mesh.halfedge[base + 2],
        ];
        if face[0].start_vert < 0 {
            continue;
        }
        let mut index = 0;
        for i in 1..3 {
            if face[i].start_vert < face[index].start_vert {
                index = i;
            }
        }
        for i in 0..3 {
            mesh.halfedge[base + i] = face[(index + i) % 3];
        }
    }

    // Step 2: fix paired_halfedge references
    for tri in 0..num_tri {
        for i in 0..3 {
            let base = tri * 3 + i;
            let curr = mesh.halfedge[base];
            if curr.start_vert < 0 {
                break; // skip collapsed triangle
            }
            let opp_face = curr.paired_halfedge as usize / 3;
            let mut index = -1i32;
            for j in 0..3 {
                if curr.start_vert == mesh.halfedge[opp_face * 3 + j].end_vert {
                    index = j as i32;
                }
            }
            mesh.halfedge[base].paired_halfedge = opp_face as i32 * 3 + index;
        }
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
