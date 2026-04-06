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

// Phase 10: collider — API-compatible foundation for overlap queries and
// geometric distance checks. The current implementation uses a deterministic
// brute-force broadphase while preserving the same entry points the later
// boolean phases will rely on.

use crate::impl_mesh::ManifoldImpl;
use crate::linalg::{cross, distance2, dot, Mat3x4, Vec3};
use crate::sort::{get_face_box_morton, morton_code};
use crate::types::Box as BBox;

#[derive(Clone, Debug, Default)]
pub struct Collider {
    leaf_bbox: Vec<BBox>,
    leaf_morton: Vec<u32>,
}

impl Collider {
    pub fn new(leaf_bbox: Vec<BBox>, leaf_morton: Vec<u32>) -> Self {
        debug_assert_eq!(leaf_bbox.len(), leaf_morton.len());
        Self { leaf_bbox, leaf_morton }
    }

    pub fn update_boxes(&mut self, leaf_bbox: Vec<BBox>) {
        debug_assert_eq!(leaf_bbox.len(), self.leaf_bbox.len());
        self.leaf_bbox = leaf_bbox;
    }

    pub fn transform(&mut self, transform: &Mat3x4) {
        debug_assert!(Self::is_axis_aligned(transform));
        for bbox in &mut self.leaf_bbox {
            *bbox = bbox.transform(transform);
        }
    }

    pub fn collisions_with_boxes<F: FnMut(usize, usize)>(
        &self,
        queries: &[BBox],
        self_collision: bool,
        mut record: F,
    ) {
        for (query_idx, query) in queries.iter().enumerate() {
            for (leaf_idx, leaf) in self.leaf_bbox.iter().enumerate() {
                if self_collision && query_idx == leaf_idx {
                    continue;
                }
                if query.does_overlap_box(leaf) {
                    record(query_idx, leaf_idx);
                }
            }
        }
    }

    /// Function-based collision query. For each query index 0..n, calls
    /// `query_box_fn(i)` to generate a bounding box on the fly, then tests it
    /// against all leaf boxes. Non-empty overlaps are reported via `record`.
    /// This is the API used by Boolean3: edge bounding boxes are generated
    /// from halfedge indices without pre-materializing an array.
    pub fn collisions_fn<F, R>(
        &self,
        query_box_fn: F,
        n: usize,
        mut record: R,
    ) where
        F: Fn(usize) -> BBox,
        R: FnMut(usize, usize),
    {
        for query_idx in 0..n {
            let query = query_box_fn(query_idx);
            if query.is_empty() {
                continue;
            }
            for (leaf_idx, leaf) in self.leaf_bbox.iter().enumerate() {
                if query.does_overlap_box(leaf) {
                    record(query_idx, leaf_idx);
                }
            }
        }
    }

    /// Point-based collision query. For each query index in 0..n, calls
    /// `point_fn(i)` to get a point, wraps it in a zero-volume box, and tests
    /// against all leaf boxes. Used by Winding03 for vertex queries.
    pub fn collisions_point<F, R>(
        &self,
        point_fn: F,
        n: usize,
        mut record: R,
    ) where
        F: Fn(usize) -> Vec3,
        R: FnMut(usize, usize),
    {
        for query_idx in 0..n {
            let pt = point_fn(query_idx);
            let query = BBox::from_point(pt);
            for (leaf_idx, leaf) in self.leaf_bbox.iter().enumerate() {
                if query.does_overlap_box(leaf) {
                    record(query_idx, leaf_idx);
                }
            }
        }
    }

    pub fn leaf_count(&self) -> usize {
        self.leaf_bbox.len()
    }

    pub fn leaf_bbox(&self) -> &[BBox] {
        &self.leaf_bbox
    }

    pub fn morton_code(position: Vec3, bbox: &BBox) -> u32 {
        morton_code(position, bbox)
    }

    pub fn is_axis_aligned(transform: &Mat3x4) -> bool {
        for row in 0..3 {
            let mut zero_count = 0;
            for col in 0..3 {
                if transform[col][row] == 0.0 {
                    zero_count += 1;
                }
            }
            if zero_count != 2 {
                return false;
            }
        }
        true
    }

    pub fn leaf_morton(&self) -> &[u32] {
        &self.leaf_morton
    }
}

pub fn edge_edge_dist(p: Vec3, a: Vec3, q: Vec3, b: Vec3) -> (Vec3, Vec3) {
    let t_vec = q - p;
    let a_dot_a = dot(a, a);
    let b_dot_b = dot(b, b);
    let a_dot_b = dot(a, b);
    let a_dot_t = dot(a, t_vec);
    let b_dot_t = dot(b, t_vec);

    let denom = a_dot_a * b_dot_b - a_dot_b * a_dot_b;
    let mut t = if denom != 0.0 {
        ((a_dot_t * b_dot_b - b_dot_t * a_dot_b) / denom).clamp(0.0, 1.0)
    } else {
        0.0
    };

    let u = if b_dot_b != 0.0 {
        let mut u = (t * a_dot_b - b_dot_t) / b_dot_b;
        if u < 0.0 {
            t = if a_dot_a != 0.0 { (a_dot_t / a_dot_a).clamp(0.0, 1.0) } else { 0.0 };
            0.0
        } else if u > 1.0 {
            t = if a_dot_a != 0.0 {
                ((a_dot_b + a_dot_t) / a_dot_a).clamp(0.0, 1.0)
            } else {
                0.0
            };
            1.0
        } else {
            u
        }
    } else {
        t = if a_dot_a != 0.0 { (a_dot_t / a_dot_a).clamp(0.0, 1.0) } else { 0.0 };
        0.0
    };

    (p + a * t, q + b * u)
}

pub fn distance_triangle_triangle_squared(p: [Vec3; 3], q: [Vec3; 3]) -> f64 {
    let sv = [p[1] - p[0], p[2] - p[1], p[0] - p[2]];
    let tv = [q[1] - q[0], q[2] - q[1], q[0] - q[2]];

    let mut shown_disjoint = false;
    let mut mindd = f64::MAX;

    for i in 0..3 {
        for j in 0..3 {
            let (cp, cq) = edge_edge_dist(p[i], sv[i], q[j], tv[j]);
            let v = cq - cp;
            let dd = dot(v, v);

            if dd <= mindd {
                mindd = dd;

                let mut id = i + 2;
                if id >= 3 {
                    id -= 3;
                }
                let mut z = p[id] - cp;
                let mut a = dot(z, v);

                id = j + 2;
                if id >= 3 {
                    id -= 3;
                }
                z = q[id] - cq;
                let mut b = dot(z, v);

                if a <= 0.0 && b >= 0.0 {
                    return dot(v, v);
                }

                if a <= 0.0 {
                    a = 0.0;
                } else if b > 0.0 {
                    b = 0.0;
                }

                if mindd - a + b > 0.0 {
                    shown_disjoint = true;
                }
            }
        }
    }

    let sn = cross(sv[0], sv[1]);
    let snl = dot(sn, sn);
    if snl > 1e-15 {
        let tp = Vec3::new(dot(p[0] - q[0], sn), dot(p[0] - q[1], sn), dot(p[0] - q[2], sn));
        let mut index = None;
        if tp.x > 0.0 && tp.y > 0.0 && tp.z > 0.0 {
            let mut idx = if tp.x < tp.y { 0 } else { 1 };
            if tp.z < tp[idx] {
                idx = 2;
            }
            index = Some(idx);
        } else if tp.x < 0.0 && tp.y < 0.0 && tp.z < 0.0 {
            let mut idx = if tp.x > tp.y { 0 } else { 1 };
            if tp.z > tp[idx] {
                idx = 2;
            }
            index = Some(idx);
        }

        if let Some(index) = index {
            shown_disjoint = true;
            let q_index = q[index];
            let mut v = q_index - p[0];
            let mut z = cross(sn, sv[0]);
            if dot(v, z) > 0.0 {
                v = q_index - p[1];
                z = cross(sn, sv[1]);
                if dot(v, z) > 0.0 {
                    v = q_index - p[2];
                    z = cross(sn, sv[2]);
                    if dot(v, z) > 0.0 {
                        let cp = q_index + sn * (tp[index] / snl);
                        let cq = q_index;
                        return dot(cp - cq, cp - cq);
                    }
                }
            }
        }
    }

    let tn = cross(tv[0], tv[1]);
    let tnl = dot(tn, tn);
    if tnl > 1e-15 {
        let sp = Vec3::new(dot(q[0] - p[0], tn), dot(q[0] - p[1], tn), dot(q[0] - p[2], tn));
        let mut index = None;
        if sp.x > 0.0 && sp.y > 0.0 && sp.z > 0.0 {
            let mut idx = if sp.x < sp.y { 0 } else { 1 };
            if sp.z < sp[idx] {
                idx = 2;
            }
            index = Some(idx);
        } else if sp.x < 0.0 && sp.y < 0.0 && sp.z < 0.0 {
            let mut idx = if sp.x > sp.y { 0 } else { 1 };
            if sp.z > sp[idx] {
                idx = 2;
            }
            index = Some(idx);
        }

        if let Some(index) = index {
            shown_disjoint = true;
            let p_index = p[index];
            let mut v = p_index - q[0];
            let mut z = cross(tn, tv[0]);
            if dot(v, z) > 0.0 {
                v = p_index - q[1];
                z = cross(tn, tv[1]);
                if dot(v, z) > 0.0 {
                    v = p_index - q[2];
                    z = cross(tn, tv[2]);
                    if dot(v, z) > 0.0 {
                        let cp = p_index;
                        let cq = p_index + tn * (sp[index] / tnl);
                        return dot(cp - cq, cp - cq);
                    }
                }
            }
        }
    }

    if shown_disjoint { mindd } else { 0.0 }
}

pub fn ray_triangle_intersection(
    origin: Vec3,
    direction: Vec3,
    tri: [Vec3; 3],
) -> Option<f64> {
    let eps = 1e-9;
    let edge1 = tri[1] - tri[0];
    let edge2 = tri[2] - tri[0];
    let h = cross(direction, edge2);
    let a = dot(edge1, h);
    if a.abs() < eps {
        return None;
    }
    let f = 1.0 / a;
    let s = origin - tri[0];
    let u = f * dot(s, h);
    if !(0.0..=1.0).contains(&u) {
        return None;
    }
    let q = cross(s, edge1);
    let v = f * dot(direction, q);
    if v < 0.0 || u + v > 1.0 {
        return None;
    }
    let t = f * dot(edge2, q);
    if t > eps { Some(t) } else { None }
}

impl ManifoldImpl {
    pub fn is_self_intersecting(&self) -> bool {
        let ep = 2.0 * self.epsilon;
        let epsilon_sq = ep * ep;
        let (face_box, face_morton) = get_face_box_morton(self);
        let collider = Collider::new(face_box.clone(), face_morton);
        let mut intersecting = false;

        collider.collisions_with_boxes(&face_box, true, |tri0, tri1| {
            if intersecting {
                return;
            }
            let tri_verts0 = self.face_triangle_vertices(tri0);
            let tri_verts1 = self.face_triangle_vertices(tri1);

            for a in &tri_verts0 {
                for b in &tri_verts1 {
                    if distance2(*a, *b) <= epsilon_sq {
                        return;
                    }
                }
            }

            if distance_triangle_triangle_squared(tri_verts0, tri_verts1) == 0.0 {
                let mut tmp0 = tri_verts0;
                let mut tmp1 = tri_verts1;
                for i in 0..3 {
                    tmp0[i] = tri_verts0[i] + self.face_normal[tri1] * ep;
                }
                if distance_triangle_triangle_squared(tmp0, tri_verts1) > 0.0 {
                    return;
                }
                for i in 0..3 {
                    tmp0[i] = tri_verts0[i] - self.face_normal[tri1] * ep;
                }
                if distance_triangle_triangle_squared(tmp0, tri_verts1) > 0.0 {
                    return;
                }
                for i in 0..3 {
                    tmp1[i] = tri_verts1[i] + self.face_normal[tri0] * ep;
                }
                if distance_triangle_triangle_squared(tri_verts0, tmp1) > 0.0 {
                    return;
                }
                for i in 0..3 {
                    tmp1[i] = tri_verts1[i] - self.face_normal[tri0] * ep;
                }
                if distance_triangle_triangle_squared(tri_verts0, tmp1) > 0.0 {
                    return;
                }
                intersecting = true;
            }
        });

        intersecting
    }

    pub fn min_gap(&self, other: &ManifoldImpl, search_length: f64) -> f64 {
        let (self_box, self_morton) = get_face_box_morton(self);
        let (mut other_box, _) = get_face_box_morton(other);
        for bbox in &mut other_box {
            bbox.min = bbox.min - Vec3::splat(search_length);
            bbox.max = bbox.max + Vec3::splat(search_length);
        }

        let collider = Collider::new(self_box, self_morton);
        let mut min_distance = f64::INFINITY;
        collider.collisions_with_boxes(&other_box, false, |tri_other, tri| {
            let p = self.face_triangle_vertices(tri);
            let q = other.face_triangle_vertices(tri_other);
            min_distance = min_distance.min(distance_triangle_triangle_squared(p, q));
        });

        min_distance.min(search_length * search_length).sqrt()
    }

    fn face_triangle_vertices(&self, tri: usize) -> [Vec3; 3] {
        [
            self.vert_pos[self.halfedge[3 * tri].start_vert as usize],
            self.vert_pos[self.halfedge[3 * tri + 1].start_vert as usize],
            self.vert_pos[self.halfedge[3 * tri + 2].start_vert as usize],
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linalg::{mat4_to_mat3x4, translation_matrix};

    #[test]
    fn test_collider_box_overlap() {
        let boxes = vec![
            BBox::from_points(Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 1.0, 1.0)),
            BBox::from_points(Vec3::new(2.0, 2.0, 2.0), Vec3::new(3.0, 3.0, 3.0)),
        ];
        let collider = Collider::new(boxes.clone(), vec![0, 1]);
        let mut hits = Vec::new();
        collider.collisions_with_boxes(&boxes, true, |a, b| hits.push((a, b)));
        assert!(hits.is_empty());

        let queries = vec![BBox::from_points(Vec3::new(0.5, 0.5, 0.5), Vec3::new(2.5, 2.5, 2.5))];
        collider.collisions_with_boxes(&queries, false, |a, b| hits.push((a, b)));
        assert_eq!(hits, vec![(0, 0), (0, 1)]);
    }

    #[test]
    fn test_ray_triangle_intersection() {
        let tri = [
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ];
        let hit = ray_triangle_intersection(Vec3::new(0.25, 0.25, -1.0), Vec3::new(0.0, 0.0, 1.0), tri);
        assert!(hit.is_some());
    }

    #[test]
    fn test_triangle_triangle_distance_zero_for_intersection() {
        let a = [
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ];
        let b = [
            Vec3::new(0.25, 0.25, -1.0),
            Vec3::new(0.25, 0.25, 1.0),
            Vec3::new(0.75, 0.25, 0.0),
        ];
        assert_eq!(distance_triangle_triangle_squared(a, b), 0.0);
    }

    #[test]
    fn test_cube_not_self_intersecting() {
        let m = ManifoldImpl::cube(&Mat3x4::identity());
        assert!(!m.is_self_intersecting());
    }

    #[test]
    fn test_min_gap_between_cubes() {
        let a = ManifoldImpl::cube(&Mat3x4::identity());
        let b = ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(Vec3::new(2.0, 0.0, 0.0))));
        let gap = a.min_gap(&b, 5.0);
        assert!((gap - 1.0).abs() < 1e-8, "gap = {}", gap);
    }
}
