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

// Phase 11: Boolean Operations (Core)
//
// C++ sources: src/boolean3.cpp (531 lines), src/boolean_result.cpp (889 lines)
//
// This module implements the edge-face intersection detection algorithm from
// boolean3.cpp. The result is consumed by boolean_result.rs to assemble the
// output mesh.
//
// Key notation (from the C++ source):
// - P and Q are the two input manifolds, R is the output
// - Dimensions: vert=0, edge=1, face=2, solid=3
// - X = winding-number quantity, S = "shadow" subset of X
// - p1q2 = edges of P intersecting faces of Q
// - x12 = winding contribution at each intersection
// - v12 = 3D position of each intersection vertex

use std::collections::HashSet;

use crate::collider::Collider;
use crate::disjoint_sets::DisjointSets;
use crate::impl_mesh::ManifoldImpl;
use crate::linalg::{IVec3, Vec2, Vec3, Vec4};
use crate::sort::get_face_box_morton;
use crate::types::{Box as BBox, OpType};

// ---------------------------------------------------------------------------
// Intersections — sparse intersection data between two meshes
// ---------------------------------------------------------------------------

/// Stores the intersections of edges of one mesh with faces of the other.
/// In forward mode: edges of P with faces of Q.
/// In reverse mode: edges of Q with faces of P.
#[derive(Clone, Default)]
pub struct Intersections {
    /// Pairs [edge_idx, face_idx] — in forward mode [p1, q2], reverse [q1, p2]
    pub p1q2: Vec<[i32; 2]>,
    /// Winding number contribution at each intersection
    pub x12: Vec<i32>,
    /// 3D position of each intersection vertex
    pub v12: Vec<Vec3>,
}

// ---------------------------------------------------------------------------
// Boolean3 — the core intersection computation
// ---------------------------------------------------------------------------

/// Computes all edge-face intersections and winding numbers between two meshes.
pub struct Boolean3 {
    pub xv12: Intersections,
    pub xv21: Intersections,
    pub w03: Vec<i32>,
    pub w30: Vec<i32>,
    pub expand_p: bool,
    pub valid: bool,
}

// ---------------------------------------------------------------------------
// Geometric kernel functions
// ---------------------------------------------------------------------------
// These are the ONLY places where floating-point operations occur in the
// boolean algorithm. They are carefully designed to minimize rounding error
// and to eliminate it at edge cases. The branch structure must exactly match
// the C++ to produce identical results.

#[inline]
fn with_sign(pos: bool, v: f64) -> f64 {
    if pos { v } else { -v }
}

/// Interpolate along edge (aL, aR) at x-coordinate `x`.
/// Returns (y, z) at the interpolated point.
/// Uses the closer endpoint as the base to minimize rounding error.
fn interpolate(a_l: Vec3, a_r: Vec3, x: f64) -> Vec2 {
    let dx_l = x - a_l.x;
    let dx_r = x - a_r.x;
    debug_assert!(
        dx_l * dx_r <= 0.0,
        "Boolean manifold error: not in domain"
    );
    let use_l = dx_l.abs() < dx_r.abs();
    let d_lr = a_r - a_l;
    let lambda = (if use_l { dx_l } else { dx_r }) / d_lr.x;
    if !lambda.is_finite() || !d_lr.y.is_finite() || !d_lr.z.is_finite() {
        return Vec2::new(a_l.y, a_l.z);
    }
    Vec2::new(
        lambda * d_lr.y + if use_l { a_l.y } else { a_r.y },
        lambda * d_lr.z + if use_l { a_l.z } else { a_r.z },
    )
}

/// Find the intersection of two edges projected onto the yz-plane, parameterized
/// by their y-coordinates. Returns (x, y, z_a, z_b) at the intersection.
fn intersect_edges(a_l: Vec3, a_r: Vec3, b_l: Vec3, b_r: Vec3) -> Vec4 {
    let dy_l = b_l.y - a_l.y;
    let dy_r = b_r.y - a_r.y;
    debug_assert!(
        dy_l * dy_r <= 0.0,
        "Boolean manifold error: no intersection"
    );
    let use_l = dy_l.abs() < dy_r.abs();
    let dx = a_r.x - a_l.x;
    let mut lambda = (if use_l { dy_l } else { dy_r }) / (dy_l - dy_r);
    if !lambda.is_finite() {
        lambda = 0.0;
    }
    let x = lambda * dx + if use_l { a_l.x } else { a_r.x };
    let a_dy = a_r.y - a_l.y;
    let b_dy = b_r.y - b_l.y;
    let use_a = a_dy.abs() < b_dy.abs();
    let y = lambda * (if use_a { a_dy } else { b_dy })
        + if use_l {
            if use_a { a_l.y } else { b_l.y }
        } else if use_a {
            a_r.y
        } else {
            b_r.y
        };
    let z = lambda * (a_r.z - a_l.z) + if use_l { a_l.z } else { a_r.z };
    let w = lambda * (b_r.z - b_l.z) + if use_l { b_l.z } else { b_r.z };
    Vec4::new(x, y, z, w)
}

/// Symbolic perturbation shadow predicate.
/// When p == q, the tie is broken by the sign of dir.
#[inline]
fn shadows(p: f64, q: f64, dir: f64) -> bool {
    if p == q { dir < 0.0 } else { p < q }
}

// ---------------------------------------------------------------------------
// Shadow01 — vertex-edge shadow test
// ---------------------------------------------------------------------------
// Tests whether vertex a0 of mesh A shadows edge b1 of mesh B.
// Returns (winding contribution, (y,z) interpolated position).

fn shadow01(
    a0: usize,
    b1: usize,
    in_a: &ManifoldImpl,
    in_b: &ManifoldImpl,
    expand_p: bool,
    forward: bool,
) -> (i32, Vec2) {
    let b1s = in_b.halfedge[b1].start_vert as usize;
    let b1e = in_b.halfedge[b1].end_vert as usize;
    let a0x = in_a.vert_pos[a0].x;
    let b1sx = in_b.vert_pos[b1s].x;
    let b1ex = in_b.vert_pos[b1e].x;
    let a0xp = in_a.vert_normal[a0].x;
    let b1sxp = in_b.vert_normal[b1s].x;
    let b1exp = in_b.vert_normal[b1e].x;

    let mut s01 = if forward {
        shadows(a0x, b1ex, with_sign(expand_p, a0xp) - b1exp) as i32
            - shadows(a0x, b1sx, with_sign(expand_p, a0xp) - b1sxp) as i32
    } else {
        shadows(b1sx, a0x, with_sign(expand_p, b1sxp) - a0xp) as i32
            - shadows(b1ex, a0x, with_sign(expand_p, b1exp) - a0xp) as i32
    };

    let mut yz01 = Vec2::new(f64::NAN, f64::NAN);

    if s01 != 0 {
        yz01 = interpolate(in_b.vert_pos[b1s], in_b.vert_pos[b1e], in_a.vert_pos[a0].x);
        let b1pair = in_b.halfedge[b1].paired_halfedge as usize;
        let dir = in_b.face_normal[b1 / 3].y + in_b.face_normal[b1pair / 3].y;
        if forward {
            if !shadows(in_a.vert_pos[a0].y, yz01.x, -dir) {
                s01 = 0;
            }
        } else if !shadows(yz01.x, in_a.vert_pos[a0].y, with_sign(expand_p, dir)) {
            s01 = 0;
        }
    }
    (s01, yz01)
}

// ---------------------------------------------------------------------------
// Kernel11 — edge-edge intersection
// ---------------------------------------------------------------------------

fn kernel11(
    p1: usize,
    q1: usize,
    in_p: &ManifoldImpl,
    in_q: &ManifoldImpl,
    expand_p: bool,
) -> (i32, Vec4) {
    let mut xyzz11 = Vec4::splat(f64::NAN);
    let mut s11: i32 = 0;

    let mut k: usize = 0;
    let mut p_rl = [Vec3::splat(0.0); 2];
    let mut q_rl = [Vec3::splat(0.0); 2];
    let mut shadow_state = false;

    let p0 = [
        in_p.halfedge[p1].start_vert as usize,
        in_p.halfedge[p1].end_vert as usize,
    ];
    for i in 0..2 {
        let (s01, yz01) = shadow01(p0[i], q1, in_p, in_q, expand_p, true);
        if yz01.x.is_finite() {
            s11 += s01 * if i == 0 { -1 } else { 1 };
            if k < 2 && (k == 0 || (s01 != 0) != shadow_state) {
                shadow_state = s01 != 0;
                p_rl[k] = in_p.vert_pos[p0[i]];
                q_rl[k] = Vec3::new(p_rl[k].x, yz01.x, yz01.y);
                k += 1;
            }
        }
    }

    let q0 = [
        in_q.halfedge[q1].start_vert as usize,
        in_q.halfedge[q1].end_vert as usize,
    ];
    for i in 0..2 {
        let (s10, yz10) = shadow01(q0[i], p1, in_q, in_p, expand_p, false);
        if yz10.x.is_finite() {
            s11 += s10 * if i == 0 { -1 } else { 1 };
            if k < 2 && (k == 0 || (s10 != 0) != shadow_state) {
                shadow_state = s10 != 0;
                q_rl[k] = in_q.vert_pos[q0[i]];
                p_rl[k] = Vec3::new(q_rl[k].x, yz10.x, yz10.y);
                k += 1;
            }
        }
    }

    if s11 == 0 {
        xyzz11 = Vec4::splat(f64::NAN);
    } else {
        debug_assert_eq!(k, 2, "Boolean manifold error: s11");
        xyzz11 = intersect_edges(p_rl[0], p_rl[1], q_rl[0], q_rl[1]);

        let p1pair = in_p.halfedge[p1].paired_halfedge as usize;
        let dir_p = in_p.face_normal[p1 / 3].z + in_p.face_normal[p1pair / 3].z;
        let q1pair = in_q.halfedge[q1].paired_halfedge as usize;
        let dir_q = in_q.face_normal[q1 / 3].z + in_q.face_normal[q1pair / 3].z;
        if !shadows(xyzz11.z, xyzz11.w, with_sign(expand_p, dir_p) - dir_q) {
            s11 = 0;
        }
    }

    (s11, xyzz11)
}

// ---------------------------------------------------------------------------
// Kernel02 — vertex-face intersection
// ---------------------------------------------------------------------------

fn kernel02(
    a0: usize,
    b2: usize,
    in_a: &ManifoldImpl,
    in_b: &ManifoldImpl,
    expand_p: bool,
    forward: bool,
) -> (i32, f64) {
    let mut s02: i32 = 0;
    let mut z02: f64 = 0.0;

    let mut k: usize = 0;
    let mut yzz_rl = [Vec3::splat(0.0); 2];
    let mut shadow_state = false;

    for i in 0..3 {
        let b1 = 3 * b2 + i;
        let edge_b = in_b.halfedge[b1];
        let b1f = if edge_b.is_forward() {
            b1
        } else {
            edge_b.paired_halfedge as usize
        };

        let (s01, yz01) = shadow01(a0, b1f, in_a, in_b, expand_p, forward);
        if yz01.x.is_finite() {
            s02 += s01 * if forward == edge_b.is_forward() { -1 } else { 1 };
            if k < 2 && (k == 0 || (s01 != 0) != shadow_state) {
                shadow_state = s01 != 0;
                yzz_rl[k] = Vec3::new(yz01.x, yz01.y, yz01.y);
                k += 1;
            }
        }
    }

    if s02 == 0 {
        z02 = f64::NAN;
    } else {
        debug_assert_eq!(k, 2, "Boolean manifold error: s02");
        let vert_pos_a = in_a.vert_pos[a0];
        z02 = interpolate(yzz_rl[0], yzz_rl[1], vert_pos_a.y).y;
        if forward {
            if !shadows(vert_pos_a.z, z02, -in_b.face_normal[b2].z) {
                s02 = 0;
            }
        } else if !shadows(z02, vert_pos_a.z, with_sign(expand_p, in_b.face_normal[b2].z)) {
            s02 = 0;
        }
    }
    (s02, z02)
}

// ---------------------------------------------------------------------------
// Kernel12 — edge-face intersection
// ---------------------------------------------------------------------------

fn kernel12(
    a1: usize,
    b2: usize,
    in_a: &ManifoldImpl,
    in_b: &ManifoldImpl,
    in_p: &ManifoldImpl,
    in_q: &ManifoldImpl,
    expand_p: bool,
    forward: bool,
) -> (i32, Vec3) {
    let mut x12: i32 = 0;
    let mut v12 = Vec3::splat(f64::NAN);

    let mut k: usize = 0;
    let mut xzy_lr0 = [Vec3::splat(0.0); 2];
    let mut xzy_lr1 = [Vec3::splat(0.0); 2];
    let mut shadow_state = false;

    let edge_a = in_a.halfedge[a1];

    for vert_a in [edge_a.start_vert as usize, edge_a.end_vert as usize] {
        let (s, z) = kernel02(vert_a, b2, in_a, in_b, expand_p, forward);
        if z.is_finite() {
            x12 += s * if (vert_a == edge_a.start_vert as usize) == forward { 1 } else { -1 };
            if k < 2 && (k == 0 || (s != 0) != shadow_state) {
                shadow_state = s != 0;
                let pos = in_a.vert_pos[vert_a];
                xzy_lr0[k] = Vec3::new(pos.x, pos.z, pos.y);
                xzy_lr1[k] = xzy_lr0[k];
                xzy_lr1[k].y = z;
                k += 1;
            }
        }
    }

    for i in 0..3 {
        let b1 = 3 * b2 + i;
        let edge_b = in_b.halfedge[b1];
        let b1f = if edge_b.is_forward() {
            b1
        } else {
            edge_b.paired_halfedge as usize
        };
        let (s, xyzz) = if forward {
            kernel11(a1, b1f, in_p, in_q, expand_p)
        } else {
            kernel11(b1f, a1, in_p, in_q, expand_p)
        };
        if xyzz.x.is_finite() {
            x12 -= s * if edge_b.is_forward() { 1 } else { -1 };
            if k < 2 && (k == 0 || (s != 0) != shadow_state) {
                shadow_state = s != 0;
                xzy_lr0[k] = Vec3::new(xyzz.x, xyzz.z, xyzz.y);
                xzy_lr1[k] = xzy_lr0[k];
                xzy_lr1[k].y = xyzz.w;
                if !forward {
                    let tmp = xzy_lr0[k].y;
                    xzy_lr0[k].y = xzy_lr1[k].y;
                    xzy_lr1[k].y = tmp;
                }
                k += 1;
            }
        }
    }

    if x12 == 0 {
        v12 = Vec3::splat(f64::NAN);
    } else {
        debug_assert_eq!(k, 2, "Boolean manifold error: v12");
        let xzyy = intersect_edges(xzy_lr0[0], xzy_lr0[1], xzy_lr1[0], xzy_lr1[1]);
        v12.x = xzyy.x;
        v12.y = xzyy.z;
        v12.z = xzyy.y;
    }
    (x12, v12)
}

// ---------------------------------------------------------------------------
// Intersect12 — find all edge-face intersections using collider broadphase
// ---------------------------------------------------------------------------

fn intersect12(
    in_p: &ManifoldImpl,
    in_q: &ManifoldImpl,
    expand_p: bool,
    forward: bool,
) -> Intersections {
    // a: edge mesh, b: face mesh
    let (a, b) = if forward { (in_p, in_q) } else { (in_q, in_p) };

    let (face_box, face_morton) = get_face_box_morton(b);
    let collider = Collider::new(face_box, face_morton);

    let mut result = Intersections::default();

    // For each edge of a, generate its bounding box and test against faces of b
    let n = a.halfedge.len();
    for query_idx in 0..n {
        if !a.halfedge[query_idx].is_forward() {
            continue;
        }
        let edge_box = BBox::from_points(
            a.vert_pos[a.halfedge[query_idx].start_vert as usize],
            a.vert_pos[a.halfedge[query_idx].end_vert as usize],
        );
        if edge_box.is_empty() {
            continue;
        }
        for leaf_idx in 0..collider.leaf_count() {
            // Manual broadphase check since we need per-edge iteration
            let leaf_box = &collider.leaf_bbox()[leaf_idx];
            if !edge_box.does_overlap_box(leaf_box) {
                continue;
            }

            let (x, v) = kernel12(query_idx, leaf_idx, a, b, in_p, in_q, expand_p, forward);
            if v.x.is_finite() {
                if forward {
                    result.p1q2.push([query_idx as i32, leaf_idx as i32]);
                } else {
                    result.p1q2.push([leaf_idx as i32, query_idx as i32]);
                }
                result.x12.push(x);
                result.v12.push(v);
            }
        }
    }

    // Sort by edge index for deterministic results
    let mut indices: Vec<usize> = (0..result.p1q2.len()).collect();
    let sort_idx = if forward { 0 } else { 1 };
    indices.sort_by(|&a, &b| {
        let pa = result.p1q2[a];
        let pb = result.p1q2[b];
        pa[sort_idx]
            .cmp(&pb[sort_idx])
            .then(pa[1 - sort_idx].cmp(&pb[1 - sort_idx]))
    });

    let old_p1q2 = result.p1q2.clone();
    let old_x12 = result.x12.clone();
    let old_v12 = result.v12.clone();
    for (new_i, &old_i) in indices.iter().enumerate() {
        result.p1q2[new_i] = old_p1q2[old_i];
        result.x12[new_i] = old_x12[old_i];
        result.v12[new_i] = old_v12[old_i];
    }

    result
}

// ---------------------------------------------------------------------------
// Winding03 — compute winding numbers via flood-fill
// ---------------------------------------------------------------------------
// Groups vertices into connected components along unbroken edges (edges not
// cut by any intersection). For each component, picks a representative vertex
// and computes its winding number via kernel02 against all overlapping faces
// of the other mesh. Then flood-fills that winding number to all vertices in
// the component.

fn winding03(
    in_p: &ManifoldImpl,
    in_q: &ManifoldImpl,
    p1q2: &[[i32; 2]],
    expand_p: bool,
    forward: bool,
) -> Vec<i32> {
    let (a, b) = if forward { (in_p, in_q) } else { (in_q, in_p) };
    let sort_idx = if forward { 0 } else { 1 };

    // Build union-find: unite vertices along unbroken edges
    let u_a = DisjointSets::new(a.vert_pos.len() as u32);
    for edge in 0..a.halfedge.len() {
        let he = &a.halfedge[edge];
        if !he.is_forward() {
            continue;
        }
        // Check if this edge is broken (has an intersection)
        let is_broken = p1q2
            .binary_search_by(|pair| pair[sort_idx].cmp(&(edge as i32)))
            .is_ok();
        if !is_broken {
            u_a.unite(he.start_vert as u32, he.end_vert as u32);
        }
    }

    // Find unique component representatives
    let mut components = HashSet::new();
    for v in 0..a.vert_pos.len() {
        components.insert(u_a.find(v as u32));
    }
    let verts: Vec<usize> = components.into_iter().map(|v| v as usize).collect();

    // Build face collider for mesh b
    let (face_box, face_morton) = get_face_box_morton(b);
    let collider = Collider::new(face_box, face_morton);

    // For each representative vertex, compute winding number via kernel02
    let mut w03 = vec![0i32; a.vert_pos.len()];

    for &vi in &verts {
        let pt = a.vert_pos[vi];
        let mut total = 0i32;
        for leaf_idx in 0..collider.leaf_count() {
            let leaf_box = &collider.leaf_bbox()[leaf_idx];
            // Use XY-only overlap check: the winding number shoots a ray in Z
            // C++ Box::DoesOverlap(vec3) only checks x,y (projected in z)
            if !(pt.x <= leaf_box.max.x && pt.x >= leaf_box.min.x
                && pt.y <= leaf_box.max.y && pt.y >= leaf_box.min.y) {
                continue;
            }
            let (s02, z02) = kernel02(vi, leaf_idx, a, b, expand_p, forward);
            if z02.is_finite() {
                let contribution = s02 * if forward { 1 } else { -1 };
                total += contribution;
                w03[vi] += contribution;
            }
        }
    }

    // Flood fill: propagate representative's winding number to all component members
    for i in 0..w03.len() {
        let root = u_a.find(i as u32) as usize;
        if root != i {
            w03[i] = w03[root];
        }
    }

    w03
}

// ---------------------------------------------------------------------------
// Boolean3 constructor
// ---------------------------------------------------------------------------

impl Boolean3 {
    /// Compute all intersections between meshes inP and inQ for the given op.
    pub fn new(in_p: &ManifoldImpl, in_q: &ManifoldImpl, op: OpType) -> Self {
        let expand_p = op == OpType::Add;

        if in_p.is_empty() || in_q.is_empty() || !in_p.bbox.does_overlap_box(&in_q.bbox) {
            return Boolean3 {
                xv12: Intersections::default(),
                xv21: Intersections::default(),
                w03: vec![0; in_p.num_vert()],
                w30: vec![0; in_q.num_vert()],
                expand_p,
                valid: true,
            };
        }

        // Level 3: find all edge-face intersections in both directions
        let xv12 = intersect12(in_p, in_q, expand_p, true);
        let xv21 = intersect12(in_p, in_q, expand_p, false);

        if xv12.x12.len() > i32::MAX as usize || xv21.x12.len() > i32::MAX as usize {
            return Boolean3 {
                xv12: Intersections::default(),
                xv21: Intersections::default(),
                w03: Vec::new(),
                w30: Vec::new(),
                expand_p,
                valid: false,
            };
        }

        // Compute winding numbers via flood fill
        let w03 = winding03(in_p, in_q, &xv12.p1q2, expand_p, true);
        let w30 = winding03(in_p, in_q, &xv21.p1q2, expand_p, false);

        Boolean3 {
            xv12,
            xv21,
            w03,
            w30,
            expand_p,
            valid: true,
        }
    }
}

// ---------------------------------------------------------------------------
// compose_meshes — concatenate disjoint meshes (unchanged from before)
// ---------------------------------------------------------------------------

fn extract_tri_vert(mesh: &ManifoldImpl) -> Vec<IVec3> {
    (0..mesh.num_tri())
        .map(|tri| {
            IVec3::new(
                mesh.halfedge[3 * tri].start_vert,
                mesh.halfedge[3 * tri + 1].start_vert,
                mesh.halfedge[3 * tri + 2].start_vert,
            )
        })
        .collect()
}

fn extract_tri_prop(mesh: &ManifoldImpl) -> Vec<IVec3> {
    (0..mesh.num_tri())
        .map(|tri| {
            IVec3::new(
                mesh.halfedge[3 * tri].prop_vert,
                mesh.halfedge[3 * tri + 1].prop_vert,
                mesh.halfedge[3 * tri + 2].prop_vert,
            )
        })
        .collect()
}

fn property_row(mesh: &ManifoldImpl, row: usize, width: usize) -> Vec<f64> {
    if mesh.num_prop == 0 {
        vec![0.0; width]
    } else {
        let mut out = vec![0.0; width];
        let src = &mesh.properties[row * mesh.num_prop..(row + 1) * mesh.num_prop];
        out[..src.len()].copy_from_slice(src);
        out
    }
}

/// Concatenate multiple disjoint meshes into one. This is a genuine utility
/// used by both boolean operations and CSG compose. It does NOT perform any
/// boolean intersection — the meshes must be non-overlapping for correct results.
pub fn compose_meshes(meshes: &[ManifoldImpl]) -> ManifoldImpl {
    if meshes.is_empty() {
        return ManifoldImpl::new();
    }
    if meshes.len() == 1 {
        return meshes[0].clone();
    }

    let num_prop = meshes.iter().map(|m| m.num_prop).max().unwrap_or(0);
    let mut vert_pos = Vec::new();
    let mut properties = Vec::new();
    let mut tri_vert = Vec::new();
    let mut tri_prop = Vec::new();
    let mut vert_offset = 0i32;
    let mut prop_offset = 0i32;

    for mesh in meshes {
        vert_pos.extend_from_slice(&mesh.vert_pos);

        let old_tri_vert = extract_tri_vert(mesh);
        let old_tri_prop = extract_tri_prop(mesh);
        tri_vert.extend(old_tri_vert.into_iter().map(|t| {
            IVec3::new(t.x + vert_offset, t.y + vert_offset, t.z + vert_offset)
        }));
        tri_prop.extend(old_tri_prop.into_iter().map(|t| {
            IVec3::new(t.x + prop_offset, t.y + prop_offset, t.z + prop_offset)
        }));

        if num_prop > 0 {
            let prop_rows = mesh.num_prop_vert();
            for row in 0..prop_rows {
                properties.extend(property_row(mesh, row, num_prop));
            }
            prop_offset += prop_rows as i32;
        } else {
            prop_offset += mesh.num_prop_vert() as i32;
        }
        vert_offset += mesh.num_vert() as i32;
    }

    let mut out = ManifoldImpl::new();
    out.vert_pos = vert_pos;
    out.num_prop = num_prop;
    out.properties = properties;
    out.create_halfedges(&tri_prop, &tri_vert);
    out.initialize_original();
    out.calculate_bbox();
    out.set_epsilon(-1.0, false);
    // required to remove parts that are smaller than the tolerance (matches C++)
    crate::edge_op::remove_degenerates(&mut out, 0);
    out.sort_geometry();
    out.increment_mesh_ids();
    out.set_normals_and_coplanar();
    out
}

// ---------------------------------------------------------------------------
// boolean — public entry point
// ---------------------------------------------------------------------------

/// Perform a 3D boolean operation on two manifold meshes.
///
/// For overlapping meshes, uses the full Boolean3 intersection algorithm.
/// For disjoint meshes, uses fast-path shortcuts.
pub fn boolean(mesh_a: &ManifoldImpl, mesh_b: &ManifoldImpl, op: OpType) -> ManifoldImpl {
    if mesh_a.is_empty() {
        return match op {
            OpType::Add => mesh_b.clone(),
            OpType::Intersect => ManifoldImpl::new(),
            OpType::Subtract => ManifoldImpl::new(),
        };
    }
    if mesh_b.is_empty() {
        return match op {
            OpType::Add | OpType::Subtract => mesh_a.clone(),
            OpType::Intersect => ManifoldImpl::new(),
        };
    }

    if !mesh_a.bbox.does_overlap_box(&mesh_b.bbox) {
        // Non-overlapping fast paths
        return match op {
            OpType::Add => compose_meshes(&[mesh_a.clone(), mesh_b.clone()]),
            OpType::Intersect => ManifoldImpl::new(),
            OpType::Subtract => mesh_a.clone(),
        };
    }

    // Full boolean — compute intersections
    let bool3 = Boolean3::new(mesh_a, mesh_b, op);
    if !bool3.valid {
        return ManifoldImpl::new();
    }

    let result = crate::boolean_result::boolean_result(mesh_a, mesh_b, op, &bool3);
    result
}

#[cfg(test)]
#[path = "boolean3_tests.rs"]
mod tests;
