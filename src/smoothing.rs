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

// Phase 9: smoothing/tangent generation — ported from cpp-reference/manifold/src/smoothing.cpp

use std::collections::{BTreeMap, HashMap};

use crate::impl_mesh::ManifoldImpl;
use crate::linalg::{cross, dot, length, normalize, Vec3, Vec4};
use crate::types::{degrees, radians, K_PI, K_PRECISION, K_TWO_PI, Smoothness, TriRef};

#[inline]
fn vec3_from_vec4(v: Vec4) -> Vec3 {
    Vec3::new(v.x, v.y, v.z)
}

#[inline]
fn safe_normalize(v: Vec3) -> Vec3 {
    let len2 = dot(v, v);
    if len2 <= 0.0 || !len2.is_finite() {
        Vec3::new(0.0, 0.0, 0.0)
    } else {
        v / len2.sqrt()
    }
}

fn wrap(radians: f64) -> f64 {
    if radians < -K_PI {
        radians + K_TWO_PI
    } else if radians > K_PI {
        radians - K_TWO_PI
    } else {
        radians
    }
}

fn angle_between(a: Vec3, b: Vec3) -> f64 {
    let d = dot(a, b);
    if d >= 1.0 {
        0.0
    } else if d <= -1.0 {
        K_PI
    } else {
        d.acos()
    }
}

pub fn circular_tangent(tangent: Vec3, edge_vec: Vec3) -> Vec4 {
    let dir = safe_normalize(tangent);
    let weight = 0.5f64.max(dot(dir, safe_normalize(edge_vec)));
    let bz2 = Vec4::new(
        dir.x * 0.5 * length(edge_vec),
        dir.y * 0.5 * length(edge_vec),
        dir.z * 0.5 * length(edge_vec),
        weight,
    );
    let bz3 = Vec4::new(
        bz2.x * (2.0 / 3.0),
        bz2.y * (2.0 / 3.0),
        bz2.z * (2.0 / 3.0),
        1.0 + (bz2.w - 1.0) * (2.0 / 3.0),
    );
    Vec4::new(bz3.x / bz3.w, bz3.y / bz3.w, bz3.z / bz3.w, bz3.w)
}

fn prev3(i: usize) -> usize {
    if i == 0 { 2 } else { i - 1 }
}

fn collect_vertex_cycle(mesh: &ManifoldImpl, start: usize) -> Vec<usize> {
    let mut cycle = Vec::new();
    let mut current = start;
    loop {
        cycle.push(current);
        let paired = mesh.halfedge[current].paired_halfedge;
        if paired < 0 {
            break;
        }
        current = crate::impl_mesh::next_halfedge(paired) as usize;
        if current == start {
            break;
        }
    }
    cycle
}

impl ManifoldImpl {
    pub fn get_normal(&self, halfedge: usize, normal_idx: usize) -> Vec3 {
        let prop = self.halfedge[halfedge].prop_vert as usize;
        let base = prop * self.num_prop + normal_idx;
        Vec3::new(
            self.properties[base],
            self.properties[base + 1],
            self.properties[base + 2],
        )
    }

    pub fn tangent_from_normal(&self, normal: Vec3, halfedge: usize) -> Vec4 {
        let edge = self.halfedge[halfedge];
        let edge_vec = self.vert_pos[edge.end_vert as usize] - self.vert_pos[edge.start_vert as usize];
        let edge_normal =
            self.face_normal[halfedge / 3] + self.face_normal[edge.paired_halfedge as usize / 3];
        let dir = cross(cross(edge_normal, edge_vec), normal);
        circular_tangent(dir, edge_vec)
    }

    pub fn is_inside_quad(&self, halfedge: usize) -> bool {
        if !self.halfedge_tangent.is_empty() {
            return self.halfedge_tangent[halfedge].w < 0.0;
        }
        let tri = halfedge / 3;
        let ref_tri = self.mesh_relation.tri_ref[tri];
        let pair = self.halfedge[halfedge].paired_halfedge as usize;
        let pair_tri = pair / 3;
        let pair_ref = self.mesh_relation.tri_ref[pair_tri];
        if !ref_tri.same_face(&pair_ref) {
            return false;
        }

        let same_face = |edge_idx: usize, reference: TriRef| -> bool {
            let pair = self.halfedge[edge_idx].paired_halfedge as usize / 3;
            reference.same_face(&self.mesh_relation.tri_ref[pair])
        };

        let mut neighbor = crate::impl_mesh::next_halfedge(halfedge as i32) as usize;
        if same_face(neighbor, ref_tri) {
            return false;
        }
        neighbor = crate::impl_mesh::next_halfedge(neighbor as i32) as usize;
        if same_face(neighbor, ref_tri) {
            return false;
        }
        neighbor = crate::impl_mesh::next_halfedge(pair as i32) as usize;
        if same_face(neighbor, pair_ref) {
            return false;
        }
        neighbor = crate::impl_mesh::next_halfedge(neighbor as i32) as usize;
        if same_face(neighbor, pair_ref) {
            return false;
        }
        true
    }

    pub fn is_marked_inside_quad(&self, halfedge: usize) -> bool {
        !self.halfedge_tangent.is_empty() && self.halfedge_tangent[halfedge].w < 0.0
    }

    pub fn update_sharpened_edges(&self, sharpened_edges: &[Smoothness]) -> Vec<Smoothness> {
        let mut old_halfedge_to_new = HashMap::new();
        for tri in 0..self.num_tri() {
            let old_tri = self.mesh_relation.tri_ref[tri].face_id;
            for i in 0..3 {
                old_halfedge_to_new.insert((3 * old_tri + i as i32) as usize, 3 * tri + i);
            }
        }
        let mut out = sharpened_edges.to_vec();
        for edge in &mut out {
            if let Some(&new_edge) = old_halfedge_to_new.get(&edge.halfedge) {
                edge.halfedge = new_edge;
            }
        }
        out
    }

    pub fn flat_faces(&self) -> Vec<bool> {
        let num_tri = self.num_tri();
        let mut tri_is_flat_face = vec![false; num_tri];
        for tri in 0..num_tri {
            let reference = self.mesh_relation.tri_ref[tri];
            let mut face_neighbors = 0;
            let mut face_tris = [-1; 3];
            for j in 0..3 {
                let neighbor_tri = self.halfedge[3 * tri + j].paired_halfedge as usize / 3;
                let j_ref = self.mesh_relation.tri_ref[neighbor_tri];
                if j_ref.same_face(&reference) {
                    face_neighbors += 1;
                    face_tris[j] = neighbor_tri as i32;
                }
            }
            if face_neighbors > 1 {
                tri_is_flat_face[tri] = true;
                for j in 0..3 {
                    if face_tris[j] >= 0 {
                        tri_is_flat_face[face_tris[j] as usize] = true;
                    }
                }
            }
        }
        tri_is_flat_face
    }

    pub fn vert_flat_face(&self, flat_faces: &[bool]) -> Vec<i32> {
        let mut vert_flat_face = vec![-1; self.num_vert()];
        let mut vert_ref = vec![TriRef::default(); self.num_vert()];
        for tri in 0..self.num_tri() {
            if flat_faces[tri] {
                for j in 0..3 {
                    let vert = self.halfedge[3 * tri + j].start_vert as usize;
                    if vert_ref[vert].same_face(&self.mesh_relation.tri_ref[tri]) {
                        continue;
                    }
                    vert_ref[vert] = self.mesh_relation.tri_ref[tri];
                    vert_flat_face[vert] = if vert_flat_face[vert] == -1 { tri as i32 } else { -2 };
                }
            }
        }
        vert_flat_face
    }

    pub fn vert_halfedge(&self) -> Vec<i32> {
        let mut vert_halfedge = vec![-1; self.num_vert()];
        for (idx, edge) in self.halfedge.iter().enumerate() {
            let start = edge.start_vert as usize;
            if vert_halfedge[start] < 0 {
                vert_halfedge[start] = idx as i32;
            }
        }
        vert_halfedge
    }

    pub fn sharpen_edges(&self, min_sharp_angle: f64, min_smoothness: f64) -> Vec<Smoothness> {
        let mut out = Vec::new();
        let min_radians = radians(min_sharp_angle);
        for e in 0..self.halfedge.len() {
            if !self.halfedge[e].is_forward() {
                continue;
            }
            let pair = self.halfedge[e].paired_halfedge as usize;
            let dihedral = dot(self.face_normal[e / 3], self.face_normal[pair / 3]).acos();
            if dihedral > min_radians {
                out.push(Smoothness { halfedge: e, smoothness: min_smoothness });
                out.push(Smoothness { halfedge: pair, smoothness: min_smoothness });
            }
        }
        out
    }

    pub fn sharpen_tangent(&mut self, halfedge: usize, smoothness: f64) {
        self.halfedge_tangent[halfedge] = Vec4::new(
            smoothness * self.halfedge_tangent[halfedge].x,
            smoothness * self.halfedge_tangent[halfedge].y,
            smoothness * self.halfedge_tangent[halfedge].z,
            if smoothness == 0.0 { 0.0 } else { self.halfedge_tangent[halfedge].w },
        );
    }

    pub fn linearize_flat_tangents(&mut self) {
        for halfedge in 0..self.halfedge_tangent.len() {
            if !self.halfedge[halfedge].is_forward() {
                continue;
            }
            let pair = self.halfedge[halfedge].paired_halfedge as usize;
            let tangent = self.halfedge_tangent[halfedge];
            let other = self.halfedge_tangent[pair];
            let flat = [tangent.w == 0.0, other.w == 0.0];
            if !flat[0] && !flat[1] {
                continue;
            }
            let edge_vec = self.vert_pos[self.halfedge[halfedge].end_vert as usize]
                - self.vert_pos[self.halfedge[halfedge].start_vert as usize];

            if flat[0] && flat[1] {
                self.halfedge_tangent[halfedge] = Vec4::new(edge_vec.x / 3.0, edge_vec.y / 3.0, edge_vec.z / 3.0, 1.0);
                self.halfedge_tangent[pair] = Vec4::new(-edge_vec.x / 3.0, -edge_vec.y / 3.0, -edge_vec.z / 3.0, 1.0);
            } else if flat[0] {
                let other_v = vec3_from_vec4(other);
                let v = (edge_vec + other_v) / 2.0;
                self.halfedge_tangent[halfedge] = Vec4::new(v.x, v.y, v.z, 1.0);
            } else {
                let tan_v = vec3_from_vec4(tangent);
                let v = (-edge_vec + tan_v) / 2.0;
                self.halfedge_tangent[pair] = Vec4::new(v.x, v.y, v.z, 1.0);
            }
        }
    }

    pub fn distribute_tangents(&mut self, fixed_halfedges: &[bool]) {
        for halfedge in 0..fixed_halfedges.len() {
            if !fixed_halfedges[halfedge] {
                continue;
            }

            let mut start = halfedge;
            if self.is_marked_inside_quad(start) {
                start = crate::impl_mesh::next_halfedge(self.halfedge[start].paired_halfedge) as usize;
            }

            let mut normal = Vec3::new(0.0, 0.0, 0.0);
            let mut current_angle = Vec::new();
            let mut desired_angle = Vec::new();

            let approx_normal = self.vert_normal[self.halfedge[start].start_vert as usize];
            let center = self.vert_pos[self.halfedge[start].start_vert as usize];
            let mut last_edge_vec =
                safe_normalize(self.vert_pos[self.halfedge[start].end_vert as usize] - center);
            let first_tangent = safe_normalize(vec3_from_vec4(self.halfedge_tangent[start]));
            let mut last_tangent = first_tangent;
            let mut current = start;
            let mut guard = 0usize;

            loop {
                guard += 1;
                if guard > self.halfedge.len() + 1 {
                    break;
                }
                current = crate::impl_mesh::next_halfedge(self.halfedge[current].paired_halfedge) as usize;
                if self.is_marked_inside_quad(current) {
                    if current == start {
                        break;
                    }
                    continue;
                }
                let this_edge_vec =
                    safe_normalize(self.vert_pos[self.halfedge[current].end_vert as usize] - center);
                let this_tangent = safe_normalize(vec3_from_vec4(self.halfedge_tangent[current]));
                normal = normal + cross(this_tangent, last_tangent);

                let cumulative = angle_between(this_edge_vec, last_edge_vec)
                    + desired_angle.last().copied().unwrap_or(0.0);
                desired_angle.push(cumulative);

                if current == start {
                    current_angle.push(K_TWO_PI);
                } else {
                    let mut angle = angle_between(this_tangent, first_tangent);
                    if dot(approx_normal, cross(this_tangent, first_tangent)) < 0.0 {
                        angle = K_TWO_PI - angle;
                    }
                    current_angle.push(angle);
                }

                last_edge_vec = this_edge_vec;
                last_tangent = this_tangent;
                if fixed_halfedges[current] {
                    break;
                }
            }

            if current_angle.len() == 1 || dot(normal, normal) == 0.0 {
                continue;
            }

            let scale = current_angle.last().copied().unwrap_or(K_TWO_PI)
                / desired_angle.last().copied().unwrap_or(K_TWO_PI);
            let mut offset = 0.0;
            if current == start {
                for i in 0..current_angle.len() {
                    offset += wrap(current_angle[i] - scale * desired_angle[i]);
                }
                offset /= current_angle.len() as f64;
            }

            current = start;
            let axis = safe_normalize(normal);
            let mut i = 0usize;
            let mut guard = 0usize;
            loop {
                guard += 1;
                if guard > self.halfedge.len() + 1 {
                    break;
                }
                current = crate::impl_mesh::next_halfedge(self.halfedge[current].paired_halfedge) as usize;
                if self.is_marked_inside_quad(current) {
                    if current == start {
                        break;
                    }
                    continue;
                }
                desired_angle[i] *= scale;
                let last_angle = if i > 0 { desired_angle[i - 1] } else { 0.0 };
                if desired_angle[i] - last_angle > K_PI {
                    desired_angle[i] = last_angle + K_PI;
                } else if i + 1 < desired_angle.len() && scale * desired_angle[i + 1] - desired_angle[i] > K_PI {
                    desired_angle[i] = scale * desired_angle[i + 1] - K_PI;
                }

                let angle = current_angle[i] - desired_angle[i] - offset;
                let tangent = vec3_from_vec4(self.halfedge_tangent[current]);
                let q = crate::linalg::rotation_quat_axis_angle(axis, angle);
                let rotated = crate::linalg::qrot(q, tangent);
                self.halfedge_tangent[current] = Vec4::new(
                    rotated.x,
                    rotated.y,
                    rotated.z,
                    self.halfedge_tangent[current].w,
                );
                i += 1;
                if fixed_halfedges[current] {
                    break;
                }
            }
        }
    }

    pub fn create_tangents_from_normals(&mut self, normal_idx: usize) {
        if self.is_empty() {
            return;
        }
        let num_vert = self.num_vert();
        let num_halfedge = self.halfedge.len();
        let mut tangent = vec![Vec4::new(0.0, 0.0, 0.0, 0.0); num_halfedge];
        let mut fixed_halfedge = vec![false; num_halfedge];
        let vert_halfedge = self.vert_halfedge();

        for &e in vert_halfedge.iter().take(num_vert) {
            if e < 0 {
                continue;
            }
            let e = e as usize;
            let cycle = collect_vertex_cycle(self, e);
            let mut face_edges = [-1isize, -1isize];

            for idx in 0..cycle.len() {
                let halfedge = cycle[idx];
                let next = cycle[(idx + 1) % cycle.len()];
                let here_normal = self.get_normal(halfedge, normal_idx);
                let next_normal = self.get_normal(next, normal_idx);
                let here_diff = self.face_normal[halfedge / 3] - here_normal;
                let next_diff = self.face_normal[next / 3] - next_normal;
                let here_is_flat = dot(here_diff, here_diff) < K_PRECISION * K_PRECISION;
                let next_is_flat = dot(next_diff, next_diff) < K_PRECISION * K_PRECISION;

                if self.is_inside_quad(halfedge) {
                    tangent[halfedge] = Vec4::new(0.0, 0.0, 0.0, -1.0);
                    continue;
                }

                let diff = next_normal - here_normal;
                let different_normals = dot(diff, diff) > K_PRECISION * K_PRECISION;
                if different_normals || here_is_flat != next_is_flat {
                    fixed_halfedge[halfedge] = true;
                    if face_edges[0] == -1 {
                        face_edges[0] = halfedge as isize;
                    } else if face_edges[1] == -1 {
                        face_edges[1] = halfedge as isize;
                    } else {
                        face_edges[0] = -2;
                    }
                }

                tangent[halfedge] = if different_normals {
                    let edge_vec = self.vert_pos[self.halfedge[halfedge].end_vert as usize]
                        - self.vert_pos[self.halfedge[halfedge].start_vert as usize];
                    let dir = cross(here_normal, next_normal);
                    let signed_dir = if dot(dir, edge_vec) < 0.0 { -dir } else { dir };
                    circular_tangent(signed_dir, edge_vec)
                } else {
                    self.tangent_from_normal(here_normal, halfedge)
                };
            }

            if face_edges[0] >= 0 && face_edges[1] >= 0 {
                let f0 = face_edges[0] as usize;
                let f1 = face_edges[1] as usize;
                let edge0 =
                    self.vert_pos[self.halfedge[f0].end_vert as usize] - self.vert_pos[self.halfedge[f0].start_vert as usize];
                let edge1 =
                    self.vert_pos[self.halfedge[f1].end_vert as usize] - self.vert_pos[self.halfedge[f1].start_vert as usize];
                let new_tangent = normalize(edge0) - normalize(edge1);
                tangent[f0] = circular_tangent(new_tangent, edge0);
                tangent[f1] = circular_tangent(-new_tangent, edge1);
            } else if face_edges[0] == -1 && face_edges[1] == -1 {
                fixed_halfedge[e] = true;
            }
        }

        self.halfedge_tangent = tangent;
        self.distribute_tangents(&fixed_halfedge);
    }

    pub fn create_tangents(&mut self, mut sharpened_edges: Vec<Smoothness>) {
        if self.is_empty() {
            return;
        }
        let num_halfedge = self.halfedge.len();
        let vert_halfedge = self.vert_halfedge();
        let tri_is_flat_face = self.flat_faces();
        let vert_flat_face = self.vert_flat_face(&tri_is_flat_face);
        let mut vert_normal = self.vert_normal.clone();
        for v in 0..self.num_vert() {
            if vert_flat_face[v] >= 0 {
                vert_normal[v] = self.face_normal[vert_flat_face[v] as usize];
            }
        }

        let mut tangent = vec![Vec4::new(0.0, 0.0, 0.0, 0.0); num_halfedge];
        let mut fixed_halfedge = vec![false; num_halfedge];
        for (edge_idx, tan) in tangent.iter_mut().enumerate() {
            *tan = if self.is_inside_quad(edge_idx) {
                Vec4::new(0.0, 0.0, 0.0, -1.0)
            } else {
                self.tangent_from_normal(vert_normal[self.halfedge[edge_idx].start_vert as usize], edge_idx)
            };
        }
        self.halfedge_tangent = tangent;

        for tri in 0..self.num_tri() {
            if !tri_is_flat_face[tri] {
                continue;
            }
            for j in 0..3 {
                let tri2 = self.halfedge[3 * tri + j].paired_halfedge as usize / 3;
                if !tri_is_flat_face[tri2]
                    || !self.mesh_relation.tri_ref[tri].same_face(&self.mesh_relation.tri_ref[tri2])
                {
                    sharpened_edges.push(Smoothness { halfedge: 3 * tri + j, smoothness: 0.0 });
                }
            }
        }

        type Pair = (Smoothness, Smoothness);
        let mut edges: BTreeMap<usize, Pair> = BTreeMap::new();
        for edge in sharpened_edges {
            if edge.smoothness >= 1.0 {
                continue;
            }
            let forward = self.halfedge[edge.halfedge].is_forward();
            let pair = self.halfedge[edge.halfedge].paired_halfedge as usize;
            let idx = if forward { edge.halfedge } else { pair };
            edges.entry(idx)
                .and_modify(|existing| {
                    let e = if forward { &mut existing.0 } else { &mut existing.1 };
                    e.smoothness = e.smoothness.min(edge.smoothness);
                })
                .or_insert_with(|| {
                    let mut pair_entry = (edge, Smoothness { halfedge: pair, smoothness: 1.0 });
                    if !forward {
                        pair_entry = (pair_entry.1, pair_entry.0);
                    }
                    pair_entry
                });
        }

        let mut vert_tangents: BTreeMap<usize, Vec<Pair>> = BTreeMap::new();
        for edge in edges.values() {
            vert_tangents
                .entry(self.halfedge[edge.0.halfedge].start_vert as usize)
                .or_default()
                .push(*edge);
            vert_tangents
                .entry(self.halfedge[edge.1.halfedge].start_vert as usize)
                .or_default()
                .push((edge.1, edge.0));
        }

        for v in 0..self.num_vert() {
            let Some(vert) = vert_tangents.get(&v) else {
                if vert_halfedge[v] >= 0 {
                    fixed_halfedge[vert_halfedge[v] as usize] = true;
                }
                continue;
            };

            if vert.len() == 1 {
                continue;
            }
            if vert.len() == 2 {
                let first = vert[0].0.halfedge;
                let second = vert[1].0.halfedge;
                fixed_halfedge[first] = true;
                fixed_halfedge[second] = true;
                let new_tangent =
                    normalize(vec3_from_vec4(self.halfedge_tangent[first]) - vec3_from_vec4(self.halfedge_tangent[second]));
                let pos = self.vert_pos[self.halfedge[first].start_vert as usize];
                self.halfedge_tangent[first] =
                    circular_tangent(new_tangent, self.vert_pos[self.halfedge[first].end_vert as usize] - pos);
                self.halfedge_tangent[second] =
                    circular_tangent(-new_tangent, self.vert_pos[self.halfedge[second].end_vert as usize] - pos);

                let mut smoothness = (vert[0].1.smoothness + vert[1].0.smoothness) / 2.0;
                for current in collect_vertex_cycle(self, first) {
                    if current == second {
                        smoothness = (vert[1].1.smoothness + vert[0].0.smoothness) / 2.0;
                    } else if current != first && !self.is_marked_inside_quad(current) {
                        self.sharpen_tangent(current, smoothness);
                    }
                }
            } else {
                let mut smoothness = 0.0;
                let mut denom = 0.0;
                for pair in vert {
                    smoothness += pair.0.smoothness + pair.1.smoothness;
                    denom += if pair.0.smoothness == 0.0 { 0.0 } else { 1.0 };
                    denom += if pair.1.smoothness == 0.0 { 0.0 } else { 1.0 };
                }
                if denom > 0.0 {
                    smoothness /= denom;
                }

                for current in collect_vertex_cycle(self, vert[0].0.halfedge) {
                    if !self.is_marked_inside_quad(current) {
                        let pair = self.halfedge[current].paired_halfedge as usize;
                        let s = if tri_is_flat_face[current / 3] || tri_is_flat_face[pair / 3] {
                            0.0
                        } else {
                            smoothness
                        };
                        self.sharpen_tangent(current, s);
                    }
                }
            }
        }

        self.linearize_flat_tangents();
        self.distribute_tangents(&fixed_halfedge);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linalg::Mat3x4;

    #[test]
    fn test_circular_tangent_is_finite() {
        let t = circular_tangent(Vec3::new(0.0, 1.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
        assert!(t.x.is_finite());
        assert!(t.y.is_finite());
        assert!(t.z.is_finite());
        assert!(t.w.is_finite());
        assert!(t.w > 0.0);
    }

    #[test]
    fn test_sharpen_edges_cube() {
        let m = ManifoldImpl::cube(&Mat3x4::identity());
        let edges = m.sharpen_edges(45.0, 0.0);
        assert_eq!(edges.len(), 24);
    }

    #[test]
    fn test_create_tangents_cube() {
        let mut m = ManifoldImpl::cube(&Mat3x4::identity());
        let sharp = m.sharpen_edges(45.0, 0.0);
        m.create_tangents(sharp);
        assert_eq!(m.halfedge_tangent.len(), m.num_halfedge());
        assert!(m.halfedge_tangent.iter().all(|t| t.x.is_finite() && t.y.is_finite() && t.z.is_finite() && t.w.is_finite()));
        assert!(m.halfedge_tangent.iter().any(|t| t.w < 0.0));
    }

    #[test]
    fn test_create_tangents_from_normals_cube() {
        let mut m = ManifoldImpl::cube(&Mat3x4::identity());
        m.num_prop = 3;
        m.properties = m
            .vert_normal
            .iter()
            .flat_map(|n| [n.x, n.y, n.z])
            .collect();
        m.create_tangents_from_normals(0);
        assert_eq!(m.halfedge_tangent.len(), m.num_halfedge());
        assert!(m.halfedge_tangent.iter().all(|t| t.w.is_finite()));
    }
}
