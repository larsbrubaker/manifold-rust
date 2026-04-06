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

// Phase 11: Boolean Result — face assembly from intersection data
//
// C++ source: src/boolean_result.cpp (889 lines)
//
// Takes intersection data from Boolean3 and assembles the output mesh faces.
// The algorithm:
// 1. Convert winding numbers to inclusion values based on operation type
// 2. Compute vertex remapping via exclusive scan (with absolute sum)
// 3. Create output vertices (retained + new intersection verts)
// 4. Build edge maps for partial and new edges
// 5. Size output (count sides per face, allocate halfedges)
// 6. Assemble edges (partial, new, whole)
// 7. Triangulate polygonal faces (Face2Tri)
// 8. Create properties via barycentric interpolation
// 9. Finalize: simplify topology, sort geometry

use std::collections::BTreeMap;

use crate::boolean3::Boolean3;
use crate::edge_op::simplify_topology;
use crate::face_op::{face2tri, get_barycentric, reorder_halfedges};
use crate::impl_mesh::{reserve_ids, ManifoldImpl};
use crate::linalg::{Vec3, dot};
use crate::types::{Halfedge, OpType, TriRef};

// ---------------------------------------------------------------------------
// EdgePos — position of a vertex along an edge for pairing
// ---------------------------------------------------------------------------

#[derive(Clone)]
struct EdgePos {
    edge_pos: f64,
    vert: i32,
    collision_id: i32,
    is_start: bool,
}

impl EdgePos {
    fn sort_key(&self) -> (OrderedF64, i32) {
        (OrderedF64(self.edge_pos), self.collision_id)
    }
}

/// Wrapper for f64 that implements Ord (for sorting edge positions)
#[derive(Clone, Copy, PartialEq)]
struct OrderedF64(f64);

impl Eq for OrderedF64 {}

impl PartialOrd for OrderedF64 {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for OrderedF64 {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.partial_cmp(&other.0).unwrap_or(std::cmp::Ordering::Equal)
    }
}

// ---------------------------------------------------------------------------
// AbsSum — exclusive scan combiner
// ---------------------------------------------------------------------------

#[inline]
fn abs_sum(a: i32, b: i32) -> i32 {
    a.abs() + b.abs()
}

/// Exclusive scan with abs_sum combiner. Output[i] = init + sum(|input[0..i]|).
fn exclusive_scan_abs(input: &[i32], init: i32) -> Vec<i32> {
    let mut result = Vec::with_capacity(input.len());
    let mut acc = init;
    for &v in input {
        result.push(acc);
        acc = abs_sum(acc, v);
    }
    result
}

// ---------------------------------------------------------------------------
// SizeOutput — compute output face sizes and allocate halfedges
// ---------------------------------------------------------------------------

/// Counts sides per face and allocates output halfedges.
/// Returns (face_edge, face_pq2r).
fn size_output(
    out_r: &mut ManifoldImpl,
    in_p: &ManifoldImpl,
    in_q: &ManifoldImpl,
    i03: &[i32],
    i30: &[i32],
    i12: &[i32],
    i21: &[i32],
    p1q2: &[[i32; 2]],
    p2q1: &[[i32; 2]],
    invert_q: bool,
) -> (Vec<i32>, Vec<i32>) {
    let num_tri_p = in_p.num_tri();
    let num_tri_q = in_q.num_tri();
    let mut sides_per_face = vec![0i32; num_tri_p + num_tri_q];

    // Count retained vertex contributions
    for face in 0..num_tri_p {
        for j in 0..3 {
            let v = in_p.halfedge[3 * face + j].start_vert as usize;
            sides_per_face[face] += i03[v].abs();
        }
    }
    for face in 0..num_tri_q {
        for j in 0..3 {
            let v = in_q.halfedge[3 * face + j].start_vert as usize;
            sides_per_face[num_tri_p + face] += i30[v].abs();
        }
    }

    // Count new intersection vertex contributions
    for (idx, pair) in p1q2.iter().enumerate() {
        let edge_p = pair[0] as usize;
        let face_q = pair[1] as usize;
        let inclusion = i12[idx].abs();
        sides_per_face[num_tri_p + face_q] += inclusion;
        let half = in_p.halfedge[edge_p];
        sides_per_face[edge_p / 3] += inclusion;
        sides_per_face[half.paired_halfedge as usize / 3] += inclusion;
    }
    for (idx, pair) in p2q1.iter().enumerate() {
        let face_p = pair[0] as usize;
        let edge_q = pair[1] as usize;
        let inclusion = i21[idx].abs();
        sides_per_face[face_p] += inclusion;
        let half = in_q.halfedge[edge_q];
        sides_per_face[num_tri_p + edge_q / 3] += inclusion;
        sides_per_face[num_tri_p + half.paired_halfedge as usize / 3] += inclusion;
    }

    // Build face_pq2r: maps input face → output face index
    let mut face_pq2r = vec![0i32; num_tri_p + num_tri_q];
    let mut count = 0i32;
    for i in 0..sides_per_face.len() {
        face_pq2r[i] = count;
        if sides_per_face[i] > 0 {
            count += 1;
        }
    }
    let num_face_r = count as usize;

    // Build face normals for output
    out_r.face_normal.resize(num_face_r, Vec3::splat(0.0));
    let mut out_idx = 0;
    for i in 0..num_tri_p {
        if sides_per_face[i] > 0 {
            out_r.face_normal[out_idx] = in_p.face_normal[i];
            out_idx += 1;
        }
    }
    for i in 0..num_tri_q {
        if sides_per_face[num_tri_p + i] > 0 {
            let normal = in_q.face_normal[i];
            out_r.face_normal[out_idx] = if invert_q {
                Vec3::new(-normal.x, -normal.y, -normal.z)
            } else {
                normal
            };
            out_idx += 1;
        }
    }

    // Build face_edge: cumulative edge counts for each output face
    let active_sides: Vec<i32> = sides_per_face.iter().filter(|&&s| s > 0).copied().collect();
    let mut face_edge = vec![0i32; active_sides.len() + 1];
    for (i, &s) in active_sides.iter().enumerate() {
        face_edge[i + 1] = face_edge[i] + s;
    }
    let num_halfedge = *face_edge.last().unwrap() as usize;
    out_r.halfedge.resize(
        num_halfedge,
        Halfedge {
            start_vert: -1,
            end_vert: -1,
            paired_halfedge: -1,
            prop_vert: -1,
        },
    );

    (face_edge, face_pq2r)
}

// ---------------------------------------------------------------------------
// AddNewEdgeVerts — populate edge maps with intersection vertices
// ---------------------------------------------------------------------------

fn add_new_edge_verts(
    edges_p: &mut BTreeMap<i32, Vec<EdgePos>>,
    edges_new: &mut BTreeMap<(i32, i32), Vec<EdgePos>>,
    p1q2: &[[i32; 2]],
    i12: &[i32],
    v12r: &[i32],
    halfedge_p: &[Halfedge],
    forward: bool,
    offset: usize,
) {
    for i in 0..p1q2.len() {
        let edge_p = p1q2[i][if forward { 0 } else { 1 }];
        let face_q = p1q2[i][if forward { 1 } else { 0 }];
        let vert = v12r[i];
        let inclusion = i12[i];

        let halfedge = halfedge_p[edge_p as usize];
        let mut key_right = (halfedge.paired_halfedge / 3, face_q);
        if !forward {
            key_right = (key_right.1, key_right.0);
        }
        let mut key_left = (edge_p / 3, face_q);
        if !forward {
            key_left = (key_left.1, key_left.0);
        }

        let direction = inclusion < 0;
        let collision_id = (i + offset) as i32;

        // C++ captures all three is_start values at array creation time:
        // edgesP: direction
        // edgesNew[keyRight]: direction ^ !forward
        // edgesNew[keyLeft]: direction ^ forward
        let dir_p = direction;
        let dir_right = direction ^ !forward;
        let dir_left = direction ^ forward;

        // Add to edge P's map
        let ep = edges_p.entry(edge_p).or_default();
        for j in 0..inclusion.abs() {
            ep.push(EdgePos {
                edge_pos: 0.0,
                vert: vert + j,
                collision_id,
                is_start: dir_p,
            });
        }

        // Add to right new edge
        let er = edges_new.entry(key_right).or_default();
        for j in 0..inclusion.abs() {
            er.push(EdgePos {
                edge_pos: 0.0,
                vert: vert + j,
                collision_id,
                is_start: dir_right,
            });
        }

        // Add to left new edge
        let el = edges_new.entry(key_left).or_default();
        for j in 0..inclusion.abs() {
            el.push(EdgePos {
                edge_pos: 0.0,
                vert: vert + j,
                collision_id,
                is_start: dir_left,
            });
        }
    }
}

// ---------------------------------------------------------------------------
// PairUp — pair start/end vertices to form halfedges
// ---------------------------------------------------------------------------

fn pair_up<F: FnMut(Halfedge)>(edge_pos: &mut Vec<EdgePos>, mut f: F) {
    assert!(
        edge_pos.len() % 2 == 0,
        "Non-manifold edge! Not an even number of points. Got {} points: starts={}, ends={}",
        edge_pos.len(),
        edge_pos.iter().filter(|e| e.is_start).count(),
        edge_pos.iter().filter(|e| !e.is_start).count(),
    );
    let n_edges = edge_pos.len() / 2;

    // Partition: starts first, then ends
    let (starts, ends): (Vec<_>, Vec<_>) = edge_pos
        .drain(..)
        .partition(|e| e.is_start);

    // Stable sort both halves
    let mut starts: Vec<_> = starts;
    let mut ends: Vec<_> = ends;
    starts.sort_by_key(|e| e.sort_key());
    ends.sort_by_key(|e| e.sort_key());

    debug_assert_eq!(starts.len(), n_edges, "Non-manifold edge!");

    for i in 0..n_edges {
        f(Halfedge {
            start_vert: starts[i].vert,
            end_vert: ends[i].vert,
            paired_halfedge: -1,
            prop_vert: -1,
        });
    }
}

// ---------------------------------------------------------------------------
// AppendPartialEdges — edges of P/Q that have intersections
// ---------------------------------------------------------------------------

fn append_partial_edges(
    out_r: &mut ManifoldImpl,
    whole_halfedge_p: &mut Vec<bool>,
    face_ptr_r: &mut Vec<i32>,
    edges_p: &mut BTreeMap<i32, Vec<EdgePos>>,
    halfedge_ref: &mut Vec<TriRef>,
    in_p: &ManifoldImpl,
    i03: &[i32],
    vp2r: &[i32],
    face_p2r: &[i32],
    forward: bool,
) {
    for (&edge_p, edge_pos_p) in edges_p.iter_mut() {
        edge_pos_p.sort_by_key(|e| e.sort_key());

        let halfedge = in_p.halfedge[edge_p as usize];
        whole_halfedge_p[edge_p as usize] = false;
        whole_halfedge_p[halfedge.paired_halfedge as usize] = false;

        let v_start = halfedge.start_vert as usize;
        let v_end = halfedge.end_vert as usize;

        // Bounds check: vp2r must be valid for these vertices
        if vp2r[v_start] < 0 || vp2r[v_start] as usize >= out_r.vert_pos.len() {
            continue;
        }
        if vp2r[v_end] < 0 || vp2r[v_end] as usize >= out_r.vert_pos.len() {
            continue;
        }
        let edge_vec = out_r.vert_pos[vp2r[v_end] as usize] - out_r.vert_pos[vp2r[v_start] as usize];

        // Fill in edge positions of existing intersection verts
        for edge in edge_pos_p.iter_mut() {
            edge.edge_pos = dot(out_r.vert_pos[edge.vert as usize], edge_vec);
        }

        // Add start vertex
        let inclusion = i03[v_start];
        let start_pos = dot(out_r.vert_pos[vp2r[v_start] as usize], edge_vec);
        for j in 0..inclusion.abs() {
            edge_pos_p.push(EdgePos {
                edge_pos: start_pos,
                vert: vp2r[v_start] + j,
                collision_id: i32::MAX,
                is_start: inclusion > 0,
            });
        }

        // Add end vertex
        let inclusion = i03[v_end];
        let end_pos = dot(out_r.vert_pos[vp2r[v_end] as usize], edge_vec);
        for j in 0..inclusion.abs() {
            edge_pos_p.push(EdgePos {
                edge_pos: end_pos,
                vert: vp2r[v_end] + j,
                collision_id: i32::MAX,
                is_start: inclusion < 0,
            });
        }

        // Pair up and add halfedges to result
        let face_left_p = edge_p as usize / 3;
        let face_left = face_p2r[face_left_p];
        let face_right_p = halfedge.paired_halfedge as usize / 3;
        let face_right = face_p2r[face_right_p];

        let forward_ref = TriRef {
            mesh_id: if forward { 0 } else { 1 },
            original_id: -1,
            face_id: face_left_p as i32,
            coplanar_id: -1,
        };
        let backward_ref = TriRef {
            mesh_id: if forward { 0 } else { 1 },
            original_id: -1,
            face_id: face_right_p as i32,
            coplanar_id: -1,
        };

        pair_up(edge_pos_p, |mut e| {
            let forward_edge = face_ptr_r[face_left as usize];
            face_ptr_r[face_left as usize] += 1;
            let backward_edge = face_ptr_r[face_right as usize];
            face_ptr_r[face_right as usize] += 1;

            e.paired_halfedge = backward_edge;
            out_r.halfedge[forward_edge as usize] = e;
            halfedge_ref[forward_edge as usize] = forward_ref;

            let rev = Halfedge {
                start_vert: e.end_vert,
                end_vert: e.start_vert,
                paired_halfedge: forward_edge,
                prop_vert: -1,
            };
            out_r.halfedge[backward_edge as usize] = rev;
            halfedge_ref[backward_edge as usize] = backward_ref;
        });
    }
}

// ---------------------------------------------------------------------------
// AppendNewEdges — edges formed at face-face intersections
// ---------------------------------------------------------------------------

fn append_new_edges(
    out_r: &mut ManifoldImpl,
    face_ptr_r: &mut Vec<i32>,
    edges_new: &mut BTreeMap<(i32, i32), Vec<EdgePos>>,
    halfedge_ref: &mut Vec<TriRef>,
    face_pq2r: &[i32],
    num_face_p: usize,
) {
    for (&(face_p, face_q), edge_pos) in edges_new.iter_mut() {
        edge_pos.sort_by_key(|e| e.sort_key());

        // Compute bounding box to find longest dimension
        let mut min = Vec3::splat(f64::INFINITY);
        let mut max = Vec3::splat(f64::NEG_INFINITY);
        for edge in edge_pos.iter() {
            let p = out_r.vert_pos[edge.vert as usize];
            min.x = min.x.min(p.x);
            min.y = min.y.min(p.y);
            min.z = min.z.min(p.z);
            max.x = max.x.max(p.x);
            max.y = max.y.max(p.y);
            max.z = max.z.max(p.z);
        }
        let size = max - min;
        // Order points along longest dimension
        let dim = if size.x > size.y && size.x > size.z {
            0
        } else if size.y > size.z {
            1
        } else {
            2
        };

        for edge in edge_pos.iter_mut() {
            let p = out_r.vert_pos[edge.vert as usize];
            edge.edge_pos = match dim {
                0 => p.x,
                1 => p.y,
                _ => p.z,
            };
        }

        let face_left = face_pq2r[face_p as usize];
        let face_right = face_pq2r[num_face_p + face_q as usize];
        let forward_ref = TriRef {
            mesh_id: 0,
            original_id: -1,
            face_id: face_p,
            coplanar_id: -1,
        };
        let backward_ref = TriRef {
            mesh_id: 1,
            original_id: -1,
            face_id: face_q,
            coplanar_id: -1,
        };

        pair_up(edge_pos, |mut e| {
            let forward_edge = face_ptr_r[face_left as usize];
            face_ptr_r[face_left as usize] += 1;
            let backward_edge = face_ptr_r[face_right as usize];
            face_ptr_r[face_right as usize] += 1;

            e.paired_halfedge = backward_edge;
            out_r.halfedge[forward_edge as usize] = e;
            halfedge_ref[forward_edge as usize] = forward_ref;

            let rev = Halfedge {
                start_vert: e.end_vert,
                end_vert: e.start_vert,
                paired_halfedge: forward_edge,
                prop_vert: -1,
            };
            out_r.halfedge[backward_edge as usize] = rev;
            halfedge_ref[backward_edge as usize] = backward_ref;
        });
    }
}

// ---------------------------------------------------------------------------
// AppendWholeEdges — edges with no intersections (fully retained)
// ---------------------------------------------------------------------------

fn append_whole_edges(
    out_r: &mut ManifoldImpl,
    face_ptr_r: &mut Vec<i32>,
    halfedge_ref: &mut Vec<TriRef>,
    in_p: &ManifoldImpl,
    whole_halfedge_p: &[bool],
    i03: &[i32],
    vp2r: &[i32],
    face_p2r: &[i32],
    forward: bool,
) {
    for idx in 0..in_p.halfedge.len() {
        if !whole_halfedge_p[idx] {
            continue;
        }
        let mut halfedge = in_p.halfedge[idx];
        if !halfedge.is_forward() {
            continue;
        }

        let inclusion = i03[halfedge.start_vert as usize];
        if inclusion == 0 {
            continue;
        }
        if inclusion < 0 {
            // Reverse the halfedge
            std::mem::swap(&mut halfedge.start_vert, &mut halfedge.end_vert);
        }
        halfedge.start_vert = vp2r[halfedge.start_vert as usize];
        halfedge.end_vert = vp2r[halfedge.end_vert as usize];

        let face_left_p = idx / 3;
        let new_face = face_p2r[face_left_p];
        let face_right_p = halfedge.paired_halfedge as usize / 3;
        let face_right = face_p2r[face_right_p];

        let forward_ref = TriRef {
            mesh_id: if forward { 0 } else { 1 },
            original_id: -1,
            face_id: face_left_p as i32,
            coplanar_id: -1,
        };
        let backward_ref = TriRef {
            mesh_id: if forward { 0 } else { 1 },
            original_id: -1,
            face_id: face_right_p as i32,
            coplanar_id: -1,
        };

        for i in 0..inclusion.abs() {
            let forward_edge = face_ptr_r[new_face as usize];
            face_ptr_r[new_face as usize] += 1;
            let backward_edge = face_ptr_r[face_right as usize];
            face_ptr_r[face_right as usize] += 1;

            let he = Halfedge {
                start_vert: halfedge.start_vert + i,
                end_vert: halfedge.end_vert + i,
                paired_halfedge: backward_edge,
                prop_vert: -1,
            };
            out_r.halfedge[forward_edge as usize] = he;
            halfedge_ref[forward_edge as usize] = forward_ref;

            let rev = Halfedge {
                start_vert: halfedge.end_vert + i,
                end_vert: halfedge.start_vert + i,
                paired_halfedge: forward_edge,
                prop_vert: -1,
            };
            out_r.halfedge[backward_edge as usize] = rev;
            halfedge_ref[backward_edge as usize] = backward_ref;
        }
    }
}

// ---------------------------------------------------------------------------
// UpdateReference — map tri refs from input meshes to output
// ---------------------------------------------------------------------------

fn update_reference(out_r: &mut ManifoldImpl, in_p: &ManifoldImpl, in_q: &ManifoldImpl, invert_q: bool) {
    let offset_q = reserve_ids(in_q.mesh_relation.mesh_id_transform.len() as u32) as i32;

    for tri_ref in out_r.mesh_relation.tri_ref.iter_mut() {
        let tri = tri_ref.face_id;
        let pq = tri_ref.mesh_id == 0;
        if pq {
            if (tri as usize) < in_p.mesh_relation.tri_ref.len() {
                *tri_ref = in_p.mesh_relation.tri_ref[tri as usize];
            }
        } else {
            if (tri as usize) < in_q.mesh_relation.tri_ref.len() {
                *tri_ref = in_q.mesh_relation.tri_ref[tri as usize];
                tri_ref.mesh_id += offset_q;
            }
        }
    }

    for (&k, v) in &in_p.mesh_relation.mesh_id_transform {
        out_r.mesh_relation.mesh_id_transform.insert(k, v.clone());
    }
    for (&k, v) in &in_q.mesh_relation.mesh_id_transform {
        let mut rel = v.clone();
        rel.back_side ^= invert_q;
        out_r.mesh_relation.mesh_id_transform.insert(k + offset_q, rel);
    }
}

// ---------------------------------------------------------------------------
// CreateProperties — barycentric interpolation of properties
// ---------------------------------------------------------------------------

fn create_properties(out_r: &mut ManifoldImpl, in_p: &ManifoldImpl, in_q: &ManifoldImpl) {
    let num_prop_p = in_p.num_prop;
    let num_prop_q = in_q.num_prop;
    let num_prop = num_prop_p.max(num_prop_q);
    out_r.num_prop = num_prop;
    if num_prop == 0 {
        return;
    }

    let num_tri = out_r.num_tri();

    // Compute barycentric coordinates for each output halfedge
    let mut bary = vec![Vec3::splat(0.0); out_r.halfedge.len()];
    for tri in 0..num_tri {
        let ref_pq = out_r.mesh_relation.tri_ref[tri];
        if out_r.halfedge[3 * tri].start_vert < 0 {
            continue;
        }

        let tri_pq = ref_pq.face_id as usize;
        let pq = ref_pq.mesh_id == 0;
        let vert_pos = if pq { &in_p.vert_pos } else { &in_q.vert_pos };
        let halfedge = if pq { &in_p.halfedge } else { &in_q.halfedge };

        if 3 * tri_pq + 2 >= halfedge.len() {
            continue;
        }

        let tri_pos = [
            vert_pos[halfedge[3 * tri_pq].start_vert as usize],
            vert_pos[halfedge[3 * tri_pq + 1].start_vert as usize],
            vert_pos[halfedge[3 * tri_pq + 2].start_vert as usize],
        ];

        for i in 0..3 {
            let vert = out_r.halfedge[3 * tri + i].start_vert;
            if vert >= 0 && (vert as usize) < out_r.vert_pos.len() {
                bary[3 * tri + i] = get_barycentric(out_r.vert_pos[vert as usize], tri_pos, out_r.epsilon);
            }
        }
    }

    // Build properties
    out_r.properties.clear();
    out_r.properties.reserve(out_r.num_vert() * num_prop);
    let mut idx = 0i32;

    for tri in 0..num_tri {
        if out_r.halfedge[3 * tri].start_vert < 0 {
            continue;
        }
        let ref_pq = out_r.mesh_relation.tri_ref[tri];
        let pq = ref_pq.mesh_id == 0;
        let old_num_prop = if pq { num_prop_p } else { num_prop_q };
        let properties = if pq { &in_p.properties } else { &in_q.properties };
        let halfedge = if pq { &in_p.halfedge } else { &in_q.halfedge };

        for i in 0..3 {
            out_r.halfedge[3 * tri + i].prop_vert = idx;
            idx += 1;
            let uvw = bary[3 * tri + i];

            for p in 0..num_prop {
                if p < old_num_prop && 3 * ref_pq.face_id as usize + 2 < halfedge.len() {
                    let mut old_props = [0.0f64; 3];
                    for j in 0..3 {
                        let prop_vert = halfedge[3 * ref_pq.face_id as usize + j].prop_vert;
                        if prop_vert >= 0 {
                            let prop_idx = old_num_prop * prop_vert as usize + p;
                            if prop_idx < properties.len() {
                                old_props[j] = properties[prop_idx];
                            }
                        }
                    }
                    let val = uvw.x * old_props[0] + uvw.y * old_props[1] + uvw.z * old_props[2];
                    out_r.properties.push(val);
                } else {
                    out_r.properties.push(0.0);
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// boolean_result — the main entry point
// ---------------------------------------------------------------------------

/// Assemble the output mesh from Boolean3 intersection data.
///
/// This is the Rust port of `Boolean3::Result()` from `boolean_result.cpp`.
pub fn boolean_result(
    in_p: &ManifoldImpl,
    in_q: &ManifoldImpl,
    op: OpType,
    bool3: &Boolean3,
) -> ManifoldImpl {
    debug_assert!(
        bool3.expand_p == (op == OpType::Add),
        "Result op type not compatible with constructor op type."
    );

    let c1 = if op == OpType::Intersect { 0 } else { 1 };
    let c2 = if op == OpType::Add { 1 } else { 0 };
    let c3 = if op == OpType::Intersect { 1 } else { -1 };

    let invert_q = op == OpType::Subtract;

    // Convert winding numbers to inclusion values
    let i12: Vec<i32> = bool3.xv12.x12.iter().map(|&v| c3 * v).collect();
    let i21: Vec<i32> = bool3.xv21.x12.iter().map(|&v| c3 * v).collect();
    let i03: Vec<i32> = bool3.w03.iter().map(|&v| c1 + c3 * v).collect();
    let i30: Vec<i32> = bool3.w30.iter().map(|&v| c2 + c3 * v).collect();

    // Vertex remapping via exclusive scan with abs_sum
    let vp2r = exclusive_scan_abs(&i03, 0);
    let num_vert_r = if let Some(&last) = i03.last() {
        abs_sum(*vp2r.last().unwrap_or(&0), last)
    } else {
        0
    };
    let n_pv = num_vert_r;

    let vq2r = exclusive_scan_abs(&i30, num_vert_r);
    let num_vert_r = if let Some(&last) = i30.last() {
        abs_sum(*vq2r.last().unwrap_or(&num_vert_r), last)
    } else {
        num_vert_r
    };
    let n_qv = num_vert_r - n_pv;

    let v12r = if !bool3.xv12.v12.is_empty() {
        exclusive_scan_abs(&i12, num_vert_r)
    } else {
        Vec::new()
    };
    let num_vert_r = if !i12.is_empty() {
        abs_sum(*v12r.last().unwrap_or(&num_vert_r), *i12.last().unwrap())
    } else {
        num_vert_r
    };
    let n12 = num_vert_r - n_pv - n_qv;

    let v21r = if !bool3.xv21.v12.is_empty() {
        exclusive_scan_abs(&i21, num_vert_r)
    } else {
        Vec::new()
    };
    let num_vert_r = if !i21.is_empty() {
        abs_sum(*v21r.last().unwrap_or(&num_vert_r), *i21.last().unwrap())
    } else {
        num_vert_r
    };
    let _n21 = num_vert_r - n_pv - n_qv - n12;

    // Create the output Manifold
    let mut out_r = ManifoldImpl::new();
    if num_vert_r == 0 {
        return out_r;
    }

    out_r.epsilon = in_p.epsilon.max(in_q.epsilon);
    out_r.tolerance = in_p.tolerance.max(in_q.tolerance);

    // Allocate and populate output vertices
    out_r.vert_pos.resize(num_vert_r as usize, Vec3::splat(0.0));

    // DuplicateVerts: retained vertices from P
    for vert in 0..in_p.num_vert() {
        let n = i03[vert].abs();
        for i in 0..n {
            out_r.vert_pos[(vp2r[vert] + i) as usize] = in_p.vert_pos[vert];
        }
    }
    // Retained vertices from Q
    for vert in 0..in_q.num_vert() {
        let n = i30[vert].abs();
        for i in 0..n {
            out_r.vert_pos[(vq2r[vert] + i) as usize] = in_q.vert_pos[vert];
        }
    }
    // New vertices from P edges → Q faces
    for vert in 0..i12.len() {
        let n = i12[vert].abs();
        for i in 0..n {
            out_r.vert_pos[(v12r[vert] + i) as usize] = bool3.xv12.v12[vert];
        }
    }
    // New vertices from Q edges → P faces
    for vert in 0..i21.len() {
        let n = i21[vert].abs();
        for i in 0..n {
            out_r.vert_pos[(v21r[vert] + i) as usize] = bool3.xv21.v12[vert];
        }
    }

    // Build edge maps
    let mut edges_p: BTreeMap<i32, Vec<EdgePos>> = BTreeMap::new();
    let mut edges_q: BTreeMap<i32, Vec<EdgePos>> = BTreeMap::new();
    let mut edges_new: BTreeMap<(i32, i32), Vec<EdgePos>> = BTreeMap::new();

    add_new_edge_verts(
        &mut edges_p,
        &mut edges_new,
        &bool3.xv12.p1q2,
        &i12,
        &v12r,
        &in_p.halfedge,
        true,
        0,
    );
    add_new_edge_verts(
        &mut edges_q,
        &mut edges_new,
        &bool3.xv21.p1q2,
        &i21,
        &v21r,
        &in_q.halfedge,
        false,
        bool3.xv12.p1q2.len(),
    );

    // Size output
    let (face_edge, face_pq2r) = size_output(
        &mut out_r,
        in_p,
        in_q,
        &i03,
        &i30,
        &i12,
        &i21,
        &bool3.xv12.p1q2,
        &bool3.xv21.p1q2,
        invert_q,
    );

    // Assemble edges
    let mut face_ptr_r = face_edge.clone();
    let mut whole_halfedge_p = vec![true; in_p.halfedge.len()];
    let mut whole_halfedge_q = vec![true; in_q.halfedge.len()];
    let mut halfedge_ref = vec![
        TriRef {
            mesh_id: 0,
            original_id: -1,
            face_id: -1,
            coplanar_id: -1,
        };
        2 * out_r.num_edge()
    ];

    append_partial_edges(
        &mut out_r,
        &mut whole_halfedge_p,
        &mut face_ptr_r,
        &mut edges_p,
        &mut halfedge_ref,
        in_p,
        &i03,
        &vp2r,
        &face_pq2r[..in_p.num_tri()],
        true,
    );
    append_partial_edges(
        &mut out_r,
        &mut whole_halfedge_q,
        &mut face_ptr_r,
        &mut edges_q,
        &mut halfedge_ref,
        in_q,
        &i30,
        &vq2r,
        &face_pq2r[in_p.num_tri()..],
        false,
    );

    append_new_edges(
        &mut out_r,
        &mut face_ptr_r,
        &mut edges_new,
        &mut halfedge_ref,
        &face_pq2r,
        in_p.num_tri(),
    );

    append_whole_edges(
        &mut out_r,
        &mut face_ptr_r,
        &mut halfedge_ref,
        in_p,
        &whole_halfedge_p,
        &i03,
        &vp2r,
        &face_pq2r[..in_p.num_tri()],
        true,
    );
    append_whole_edges(
        &mut out_r,
        &mut face_ptr_r,
        &mut halfedge_ref,
        in_q,
        &whole_halfedge_q,
        &i30,
        &vq2r,
        &face_pq2r[in_p.num_tri()..],
        false,
    );

    // Triangulate polygonal faces
    face2tri(&mut out_r, &face_edge, &halfedge_ref);
    reorder_halfedges(&mut out_r);

    // Create properties via barycentric interpolation
    create_properties(&mut out_r, in_p, in_q);

    // Update references
    update_reference(&mut out_r, in_p, in_q, invert_q);

    // Simplify topology
    simplify_topology(&mut out_r, (n_pv + n_qv) as i32);
    out_r.remove_unreferenced_verts();

    // Finalize
    out_r.calculate_bbox();
    out_r.sort_geometry();
    // IncrementMeshIDs will be needed for CSG tree; skip for now

    out_r
}
