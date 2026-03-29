// Phase 3: Polygon Triangulation — ported from src/polygon.cpp, src/tree2d.h/cpp, src/utils.h
//
// Key algorithm: ear-clipping with 2D KD-tree acceleration, convex fast-path,
// hole key-holing. Must produce identical triangulations to the C++ version.

use std::collections::HashMap;
use crate::linalg::{Vec2, IVec3};
use crate::types::{PolyVert, PolygonsIdx, SimplePolygonIdx, Polygons, Rect, K_PRECISION};

type IVec3Out = IVec3;

const INVALID: usize = usize::MAX;
const K_BEST: f64 = f64::NEG_INFINITY;

// ---------------------------------------------------------------------------
// CCW (from utils.h)
// ---------------------------------------------------------------------------

/// Determines if p0, p1, p2 are wound CCW, CW, or colinear within tolerance.
/// Returns 1 (CCW), -1 (CW), or 0 (colinear within tol).
pub fn ccw(p0: Vec2, p1: Vec2, p2: Vec2, tol: f64) -> i32 {
    let v1 = p1 - p0;
    let v2 = p2 - p0;
    let area = v1.x * v2.y - v1.y * v2.x;
    let base2 = (v1.x * v1.x + v1.y * v1.y).max(v2.x * v2.x + v2.y * v2.y);
    if area * area * 4.0 <= base2 * tol * tol {
        0
    } else if area > 0.0 {
        1
    } else {
        -1
    }
}

#[inline]
fn determinant2x2(a: Vec2, b: Vec2) -> f64 {
    a.x * b.y - a.y * b.x
}

#[inline]
fn safe_normalize_2d(v: Vec2) -> Vec2 {
    let len = (v.x * v.x + v.y * v.y).sqrt();
    if len == 0.0 || !len.is_finite() {
        Vec2::new(0.0, 0.0)
    } else {
        Vec2::new(v.x / len, v.y / len)
    }
}

#[inline]
fn dot2d(a: Vec2, b: Vec2) -> f64 {
    a.x * b.x + a.y * b.y
}

// ---------------------------------------------------------------------------
// 2D KD-tree (from tree2d.h/cpp)
// ---------------------------------------------------------------------------

/// Recursive in-place kd-tree construction on a PolyVert slice.
/// Alternates between sorting by x and y.
fn build_two_d_tree_impl(points: &mut [PolyVert], sort_x: bool) {
    if sort_x {
        points.sort_by(|a, b| a.pos.x.partial_cmp(&b.pos.x).unwrap_or(std::cmp::Ordering::Equal));
    } else {
        points.sort_by(|a, b| a.pos.y.partial_cmp(&b.pos.y).unwrap_or(std::cmp::Ordering::Equal));
    }
    if points.len() < 2 {
        return;
    }
    let mid = points.len() / 2;
    build_two_d_tree_impl(&mut points[..mid], !sort_x);
    build_two_d_tree_impl(&mut points[mid + 1..], !sort_x);
}

fn build_two_d_tree(points: &mut Vec<PolyVert>) {
    if points.len() <= 8 {
        return;
    }
    build_two_d_tree_impl(points.as_mut_slice(), true);
}

/// Query the 2D kd-tree, calling `f` for every point inside rect `r`.
fn query_two_d_tree<F: FnMut(PolyVert)>(points: &[PolyVert], r: Rect, mut f: F) {
    if points.len() <= 8 {
        for p in points {
            if r.contains_point(p.pos) {
                f(*p);
            }
        }
        return;
    }

    // Stack-based traversal
    let mut stack: Vec<(Rect, usize, usize, i32)> = Vec::with_capacity(64); // (current_rect, start, len, level)

    // Initial rect: infinite
    let inf_rect = Rect {
        min: Vec2::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
        max: Vec2::new(f64::INFINITY, f64::INFINITY),
    };
    stack.push((inf_rect, 0, points.len(), 0));

    while let Some((current, start, len, level)) = stack.pop() {
        if len <= 8 {
            for p in &points[start..start + len] {
                if r.contains_point(p.pos) {
                    f(*p);
                }
            }
            continue;
        }

        let mid = len / 2;
        let middle = points[start + mid];

        let mut left = current;
        let mut right = current;
        if level % 2 == 0 {
            left.max.x = middle.pos.x;
            right.min.x = middle.pos.x;
        } else {
            left.max.y = middle.pos.y;
            right.min.y = middle.pos.y;
        }

        if r.contains_point(middle.pos) {
            f(middle);
        }

        if left.does_overlap(&r) {
            stack.push((left, start, mid, level + 1));
        }
        if right.does_overlap(&r) {
            stack.push((right, start + mid + 1, len - mid - 1, level + 1));
        }
    }
}

// ---------------------------------------------------------------------------
// IsConvex fast-path check
// ---------------------------------------------------------------------------

fn is_convex(polys: &PolygonsIdx, epsilon: f64) -> bool {
    for poly in polys {
        if poly.is_empty() {
            continue;
        }
        let first_edge = poly[0].pos - poly[poly.len() - 1].pos;
        let mut last_edge = safe_normalize_2d(first_edge);
        for v in 0..poly.len() {
            let edge = if v + 1 < poly.len() {
                poly[v + 1].pos - poly[v].pos
            } else {
                first_edge
            };
            let det = determinant2x2(last_edge, edge);
            if det <= 0.0 || (det.abs() < epsilon && dot2d(last_edge, edge) < 0.0) {
                return false;
            }
            last_edge = safe_normalize_2d(edge);
        }
    }
    true
}

// ---------------------------------------------------------------------------
// TriangulateConvex — fast alternating triangulation for convex polygons
// ---------------------------------------------------------------------------

fn triangulate_convex(polys: &PolygonsIdx) -> Vec<IVec3Out> {
    let num_tri: usize = polys.iter().map(|p| p.len().saturating_sub(2)).sum();
    let mut triangles = Vec::with_capacity(num_tri);
    for poly in polys {
        let mut i = 0usize;
        let mut k = poly.len().saturating_sub(1);
        let mut right = true;
        while i + 1 < k {
            let j = if right { i + 1 } else { k - 1 };
            triangles.push(IVec3Out::new(poly[i].idx, poly[j].idx, poly[k].idx));
            if right {
                i = j;
            } else {
                k = j;
            }
            right = !right;
        }
    }
    triangles
}

// ---------------------------------------------------------------------------
// EarClip — main triangulator
// ---------------------------------------------------------------------------

struct Vert {
    mesh_idx: i32,
    cost: f64,
    /// Incremented each time the vert is removed from the ear queue (lazy deletion).
    ear_version: u32,
    pos: Vec2,
    right_dir: Vec2,
    left: usize,
    right: usize,
}

impl Clone for Vert {
    fn clone(&self) -> Self {
        Vert {
            mesh_idx: self.mesh_idx,
            cost: self.cost,
            ear_version: self.ear_version,
            pos: self.pos,
            right_dir: self.right_dir,
            left: self.left,
            right: self.right,
        }
    }
}

/// Entry in the min-heap ear queue.
#[derive(Clone, Copy, PartialEq)]
struct EarEntry {
    cost: f64,
    idx: usize,
    version: u32,
}

impl Eq for EarEntry {}

impl PartialOrd for EarEntry {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for EarEntry {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // Min-heap: lower cost = higher priority (appears first in pop)
        // BinaryHeap is max-heap, so "greater" = higher priority = lower cost
        other
            .cost
            .partial_cmp(&self.cost)
            .unwrap_or(std::cmp::Ordering::Equal)
    }
}

struct IdxCollider {
    points: Vec<PolyVert>,
    itr: Vec<usize>, // itr[i] = polygon vert index for points[i].idx
}

struct EarClip {
    polygon: Vec<Vert>,
    holes: Vec<usize>,             // start indices of hole polygons (negative area), sorted by pos.x desc
    outers: Vec<usize>,            // start indices of outer polygons
    simples: Vec<usize>,           // start indices of all simple polygons (after keyholing)
    hole2bbox: HashMap<usize, Rect>,
    ears_queue: std::collections::BinaryHeap<EarEntry>,
    triangles: Vec<IVec3Out>,
    bbox: Rect,
    epsilon: f64,
}

impl EarClip {
    fn new(polys: &PolygonsIdx, epsilon: f64) -> Self {
        let num_vert: usize = polys.iter().map(|p| p.len()).sum();
        let mut ec = EarClip {
            polygon: Vec::with_capacity(num_vert + 2 * polys.len()),
            holes: Vec::new(),
            outers: Vec::new(),
            simples: Vec::new(),
            hole2bbox: HashMap::new(),
            ears_queue: std::collections::BinaryHeap::new(),
            triangles: Vec::with_capacity(num_vert + 2 * polys.len()),
            bbox: Rect::new(),
            epsilon,
        };

        let starts = ec.initialize(polys);

        // Clip degenerate verts
        for i in 0..ec.polygon.len() {
            ec.clip_if_degenerate(i);
        }

        // Find start/classify each polygon
        for start in starts {
            ec.find_start(start);
        }

        ec
    }

    fn triangulate(mut self) -> (Vec<IVec3Out>, f64) {
        // Key-hole all holes into outer polygons
        let holes: Vec<usize> = self.holes.clone();
        for start in holes {
            self.cut_keyhole(start);
        }

        // Ear-clip each simple polygon
        let simples: Vec<usize> = self.simples.clone();
        for start in simples {
            self.triangulate_poly(start);
        }

        let eps = self.epsilon;
        (self.triangles, eps)
    }

    // -----------------------------------------------------------------------
    // Linked-list helpers
    // -----------------------------------------------------------------------

    fn clipped(&self, v: usize) -> bool {
        self.polygon[self.polygon[v].right].left != v
    }

    fn link(&mut self, left: usize, right: usize) {
        self.polygon[left].right = right;
        self.polygon[right].left = left;
        let dir = self.polygon[right].pos - self.polygon[left].pos;
        self.polygon[left].right_dir = safe_normalize_2d(dir);
    }

    /// Apply `f` to each unclipped vert in the polygon ring starting from `first`.
    /// Returns `Some(last_v)` on success (last_v == first), or `None` if degenerate.
    fn loop_verts(&self, first: usize) -> Option<Vec<usize>> {
        let mut result = Vec::new();
        let mut v = first;
        let mut cur_first = first;
        loop {
            if self.clipped(v) {
                cur_first = self.polygon[self.polygon[v].right].left;
                if !self.clipped(cur_first) {
                    v = cur_first;
                    if self.polygon[v].right == self.polygon[v].left {
                        return None;
                    }
                    result.push(v);
                }
            } else {
                if self.polygon[v].right == self.polygon[v].left {
                    return None;
                }
                result.push(v);
            }
            v = self.polygon[v].right;
            if v == cur_first {
                break;
            }
        }
        Some(result)
    }

    // -----------------------------------------------------------------------
    // Vert predicate methods (take vert index `v`)
    // -----------------------------------------------------------------------

    fn vert_is_short(&self, v: usize) -> bool {
        let edge = self.polygon[self.polygon[v].right].pos - self.polygon[v].pos;
        dot2d(edge, edge) * 4.0 < self.epsilon * self.epsilon
    }

    /// Returns 1 if p is inside the angle at v, -1 outside, 0 within epsilon.
    fn vert_interior(&self, v: usize, p: Vec2) -> i32 {
        let left = self.polygon[v].left;
        let right = self.polygon[v].right;
        let diff = p - self.polygon[v].pos;
        if dot2d(diff, diff) < self.epsilon * self.epsilon {
            return 0;
        }
        ccw(self.polygon[v].pos, self.polygon[left].pos, self.polygon[right].pos, self.epsilon)
            + ccw(self.polygon[v].pos, self.polygon[right].pos, p, self.epsilon)
            + ccw(self.polygon[v].pos, p, self.polygon[left].pos, self.epsilon)
    }

    /// Returns true if vert `v` is on the inside of the edge tail → tail.right.
    /// to_left: walk v's polygon edges to the left (vs right).
    fn vert_inside_edge(&self, v: usize, tail: usize, to_left: bool) -> bool {
        let p2 = self.epsilon * self.epsilon;
        let mut next_l = self.polygon[self.polygon[v].left].right;
        let mut next_r = self.polygon[tail].right;
        let mut center = tail;
        let mut last = center;

        let v_stop = if to_left { self.polygon[v].right } else { self.polygon[v].left };

        loop {
            if next_l == next_r || tail == next_r || next_l == v_stop {
                break;
            }

            let edge_l = self.polygon[next_l].pos - self.polygon[center].pos;
            let l2 = dot2d(edge_l, edge_l);
            if l2 <= p2 {
                next_l = if to_left { self.polygon[next_l].left } else { self.polygon[next_l].right };
                continue;
            }

            let edge_r = self.polygon[next_r].pos - self.polygon[center].pos;
            let r2 = dot2d(edge_r, edge_r);
            if r2 <= p2 {
                next_r = self.polygon[next_r].right;
                continue;
            }

            let vec_lr = self.polygon[next_r].pos - self.polygon[next_l].pos;
            let lr2 = dot2d(vec_lr, vec_lr);
            if lr2 <= p2 {
                last = center;
                center = next_l;
                next_l = if to_left { self.polygon[next_l].left } else { self.polygon[next_l].right };
                if next_l == next_r {
                    break;
                }
                next_r = self.polygon[next_r].right;
                continue;
            }

            let mut convexity = ccw(
                self.polygon[next_l].pos,
                self.polygon[center].pos,
                self.polygon[next_r].pos,
                self.epsilon,
            );
            if center != last {
                convexity += ccw(self.polygon[last].pos, self.polygon[center].pos, self.polygon[next_l].pos, self.epsilon)
                    + ccw(self.polygon[next_r].pos, self.polygon[center].pos, self.polygon[last].pos, self.epsilon);
            }
            if convexity != 0 {
                return convexity > 0;
            }

            if l2 < r2 {
                center = next_l;
                next_l = if to_left { self.polygon[next_l].left } else { self.polygon[next_l].right };
            } else {
                center = next_r;
                next_r = self.polygon[next_r].right;
            }
            last = center;
        }
        true
    }

    fn vert_is_convex(&self, v: usize, epsilon: f64) -> bool {
        let left = self.polygon[v].left;
        let right = self.polygon[v].right;
        ccw(self.polygon[left].pos, self.polygon[v].pos, self.polygon[right].pos, epsilon) >= 0
    }

    fn vert_is_reflex(&self, v: usize) -> bool {
        let left = self.polygon[v].left;
        !self.vert_inside_edge(left, self.polygon[left].right, true)
    }

    /// x-value on this vert's right edge at start.y, or NaN if no crossing.
    fn vert_interp_y2x(&self, v: usize, start: Vec2, on_top: i32) -> f64 {
        let pos = self.polygon[v].pos;
        let rpos = self.polygon[self.polygon[v].right].pos;
        let eps = self.epsilon;
        if (pos.y - start.y).abs() <= eps {
            if rpos.y <= start.y + eps || on_top == 1 {
                f64::NAN
            } else {
                pos.x
            }
        } else if pos.y < start.y - eps {
            if rpos.y > start.y + eps {
                pos.x + (start.y - pos.y) * (rpos.x - pos.x) / (rpos.y - pos.y)
            } else if rpos.y < start.y - eps || on_top == -1 {
                f64::NAN
            } else {
                rpos.x
            }
        } else {
            f64::NAN
        }
    }

    /// Signed distance of vert `other` relative to the edge at `v` in direction `unit`.
    fn vert_signed_dist(&self, v: usize, other: usize, unit: Vec2) -> f64 {
        let eps = self.epsilon;
        let d = determinant2x2(unit, self.polygon[other].pos - self.polygon[v].pos);
        if d.abs() < eps {
            let d_r = determinant2x2(unit, self.polygon[self.polygon[other].right].pos - self.polygon[v].pos);
            if d_r.abs() > eps {
                return d_r;
            }
            let d_l = determinant2x2(unit, self.polygon[self.polygon[other].left].pos - self.polygon[v].pos);
            if d_l.abs() > eps {
                return d_l;
            }
        }
        d
    }

    /// Cost of vert `other` within ear `v`, where `open_side` is the unit vector left→right.
    fn vert_cost(&self, v: usize, other: usize, open_side: Vec2) -> f64 {
        let left = self.polygon[v].left;
        let right = self.polygon[v].right;
        let cost = self.vert_signed_dist(v, other, self.polygon[v].right_dir)
            .min(self.vert_signed_dist(left, other, self.polygon[left].right_dir));
        let open_cost = determinant2x2(open_side, self.polygon[other].pos - self.polygon[right].pos);
        cost.min(open_cost)
    }

    fn delaunay_cost(diff: Vec2, scale: f64, epsilon: f64) -> f64 {
        -epsilon - scale * dot2d(diff, diff)
    }

    fn vert_ear_cost(&self, v: usize, collider: &IdxCollider) -> f64 {
        let left = self.polygon[v].left;
        let right = self.polygon[v].right;
        let open_side_vec = self.polygon[left].pos - self.polygon[right].pos;
        let center = (self.polygon[left].pos + self.polygon[right].pos) * 0.5;
        let denom = dot2d(open_side_vec, open_side_vec);
        let scale = if denom > 0.0 { 4.0 / denom } else { 0.0 };
        let radius = denom.sqrt() * 0.5;
        let open_side = safe_normalize_2d(open_side_vec);

        let mut total_cost = dot2d(self.polygon[left].right_dir, self.polygon[v].right_dir) - 1.0 - self.epsilon;

        // Folded ears: clip first
        if ccw(self.polygon[v].pos, self.polygon[left].pos, self.polygon[right].pos, self.epsilon) == 0 {
            return total_cost;
        }

        // Build ear bounding box, expanded to include pos, then by epsilon
        let cx = center.x;
        let cy = center.y;
        let mut ear_box = Rect::from_points(
            Vec2::new(cx - radius, cy - radius),
            Vec2::new(cx + radius, cy + radius),
        );
        ear_box.union_point(self.polygon[v].pos);
        ear_box.min = ear_box.min - Vec2::splat(self.epsilon);
        ear_box.max = ear_box.max + Vec2::splat(self.epsilon);

        let lid = self.polygon[left].mesh_idx;
        let rid = self.polygon[right].mesh_idx;
        let mid_id = self.polygon[v].mesh_idx;

        let mut tc = total_cost;
        query_two_d_tree(&collider.points, ear_box, |point| {
            let test = collider.itr[point.idx as usize];
            if !self.clipped(test)
                && self.polygon[test].mesh_idx != mid_id
                && self.polygon[test].mesh_idx != lid
                && self.polygon[test].mesh_idx != rid
            {
                let mut cost = self.vert_cost(v, test, open_side);
                if cost < -self.epsilon {
                    cost = Self::delaunay_cost(self.polygon[test].pos - center, scale, self.epsilon);
                }
                if cost > tc {
                    tc = cost;
                }
            }
        });
        tc
    }

    // -----------------------------------------------------------------------
    // Core algorithm methods
    // -----------------------------------------------------------------------

    /// Remove ear vert from polygon and emit a triangle.
    fn clip_ear(&mut self, ear: usize) {
        let left = self.polygon[ear].left;
        let right = self.polygon[ear].right;
        self.link(left, right);
        if self.polygon[left].mesh_idx != self.polygon[ear].mesh_idx
            && self.polygon[ear].mesh_idx != self.polygon[right].mesh_idx
            && self.polygon[right].mesh_idx != self.polygon[left].mesh_idx
        {
            self.triangles.push(IVec3Out::new(
                self.polygon[left].mesh_idx,
                self.polygon[ear].mesh_idx,
                self.polygon[right].mesh_idx,
            ));
        }
    }

    /// Clip degenerate ears (zero-area or very short edges).
    fn clip_if_degenerate(&mut self, ear: usize) {
        if self.clipped(ear) {
            return;
        }
        if self.polygon[ear].left == self.polygon[ear].right {
            return;
        }
        let right = self.polygon[ear].right;
        let left = self.polygon[ear].left;
        let is_short = self.vert_is_short(ear);
        let is_colinear_fold = {
            let lp = self.polygon[left].pos;
            let ep = self.polygon[ear].pos;
            let rp = self.polygon[right].pos;
            let diff_l = lp - ep;
            let diff_r = rp - ep;
            ccw(lp, ep, rp, self.epsilon) == 0 && dot2d(diff_l, diff_r) > 0.0
        };
        if is_short || is_colinear_fold {
            self.clip_ear(ear);
            let l = self.polygon[ear].left;
            let r = self.polygon[ear].right;
            self.clip_if_degenerate(l);
            self.clip_if_degenerate(r);
        }
    }

    /// Build the circular polygon list from input polygons. Returns start indices.
    fn initialize(&mut self, polys: &PolygonsIdx) -> Vec<usize> {
        let mut starts = Vec::new();
        for poly in polys {
            if poly.is_empty() {
                continue;
            }
            let first_idx = self.polygon.len();
            let vert = &poly[0];
            self.polygon.push(Vert {
                mesh_idx: vert.idx,
                cost: 0.0,
                ear_version: 0,
                pos: vert.pos,
                right_dir: Vec2::new(0.0, 0.0),
                left: INVALID,
                right: INVALID,
            });
            self.bbox.union_point(vert.pos);
            let mut last = first_idx;
            starts.push(first_idx);

            for vert in &poly[1..] {
                let next_idx = self.polygon.len();
                self.polygon.push(Vert {
                    mesh_idx: vert.idx,
                    cost: 0.0,
                    ear_version: 0,
                    pos: vert.pos,
                    right_dir: Vec2::new(0.0, 0.0),
                    left: INVALID,
                    right: INVALID,
                });
                self.bbox.union_point(vert.pos);
                self.link(last, next_idx);
                last = next_idx;
            }
            self.link(last, first_idx);
        }

        if self.epsilon < 0.0 {
            self.epsilon = self.bbox.scale() * K_PRECISION;
        }

        starts
    }

    /// Classify a polygon as hole or outer. For holes, find the rightmost reflex vert.
    fn find_start(&mut self, first: usize) {
        let origin = self.polygon[first].pos;
        let mut start = first;
        let mut max_x = f64::NEG_INFINITY;
        let mut bbox = Rect::new();
        // Kahan summation
        let mut area: f64 = 0.0;
        let mut area_comp: f64 = 0.0;

        let verts = match self.loop_verts(first) {
            None => return,
            Some(v) => v,
        };

        for &v in &verts {
            bbox.union_point(self.polygon[v].pos);
            let a1 = determinant2x2(
                self.polygon[v].pos - origin,
                self.polygon[self.polygon[v].right].pos - origin,
            );
            let t1 = area + a1;
            area_comp += (area - t1) + a1;
            area = t1;

            if self.polygon[v].pos.x > max_x && self.vert_is_reflex(v) {
                max_x = self.polygon[v].pos.x;
                start = v;
            }
        }

        area += area_comp;
        let size = bbox.size();
        let min_area = self.epsilon * size.x.max(size.y);

        if max_x.is_finite() && area < -min_area {
            // Hole (negative area)
            // Insert into holes sorted by max_x descending
            self.holes.push(start);
            self.hole2bbox.insert(start, bbox);
        } else {
            self.simples.push(start);
            if area > min_area {
                self.outers.push(start);
            }
        }
    }

    /// Attach a hole to an outer polygon via a keyhole.
    fn cut_keyhole(&mut self, start: usize) {
        let bbox = *self.hole2bbox.get(&start).unwrap();
        let start_pos = self.polygon[start].pos;
        let on_top: i32 = if start_pos.y >= bbox.max.y - self.epsilon {
            1
        } else if start_pos.y <= bbox.min.y + self.epsilon {
            -1
        } else {
            0
        };
        let mut connector: usize = INVALID;

        let outers: Vec<usize> = self.outers.clone();
        for outer_start in &outers {
            let verts = match self.loop_verts(*outer_start) {
                None => continue,
                Some(v) => v,
            };
            for &edge in &verts {
                let x = self.vert_interp_y2x(edge, start_pos, on_top);
                if x.is_finite() && self.vert_inside_edge(start, self.polygon[start].left, true) {
                    // Simplified connector selection: first valid edge
                    if connector == INVALID {
                        connector = edge;
                    } else {
                        // Pick the edge that results in a more leftward connection
                        let edge_x = x;
                        let conn_x = self.vert_interp_y2x(connector, start_pos, on_top);
                        let conn_right = self.polygon[connector].right;
                        let inside_edge_check = ccw(
                            Vec2::new(edge_x, start_pos.y),
                            self.polygon[connector].pos,
                            self.polygon[conn_right].pos,
                            self.epsilon,
                        );
                        if inside_edge_check == 1 {
                            connector = edge;
                        } else if inside_edge_check == 0 {
                            // Use vertical ordering as tiebreak
                            if self.polygon[connector].pos.y < self.polygon[edge].pos.y {
                                if self.vert_inside_edge(edge, connector, false) {
                                    connector = edge;
                                }
                            } else if !self.vert_inside_edge(connector, edge, false) {
                                connector = edge;
                            }
                        }
                    }
                }
            }
        }

        if connector == INVALID {
            self.simples.push(start);
            return;
        }

        connector = self.find_closer_bridge(start, connector);
        self.join_polygons(start, connector);
    }

    /// Refine keyhole connector: find any reflex vert closer to start.
    fn find_closer_bridge(&self, start: usize, edge: usize) -> usize {
        let start_pos = self.polygon[start].pos;
        let edge_right = self.polygon[edge].right;
        let mut connector = if self.polygon[edge].pos.x < start_pos.x {
            edge_right
        } else if self.polygon[edge_right].pos.x < start_pos.x {
            edge
        } else if self.polygon[edge_right].pos.y - start_pos.y > start_pos.y - self.polygon[edge].pos.y {
            edge
        } else {
            edge_right
        };

        if (self.polygon[connector].pos.y - start_pos.y).abs() <= self.epsilon {
            return connector;
        }
        let above: f64 = if self.polygon[connector].pos.y > start_pos.y {
            1.0
        } else {
            -1.0
        };

        let outers: Vec<usize> = self.outers.clone();
        for outer_start in &outers {
            let verts = match self.loop_verts(*outer_start) {
                None => continue,
                Some(v) => v,
            };
            for &vert in &verts {
                let inside = above
                    * ccw(start_pos, self.polygon[vert].pos, self.polygon[connector].pos, self.epsilon) as f64;
                let vp = self.polygon[vert].pos;
                let cp = self.polygon[connector].pos;
                if vp.x > start_pos.x - self.epsilon
                    && vp.y * above > start_pos.y * above - self.epsilon
                    && (inside > 0.0
                        || (inside == 0.0 && vp.x < cp.x && vp.y * above < cp.y * above))
                    && self.vert_inside_edge(vert, edge, true)
                    && self.vert_is_reflex(vert)
                {
                    connector = vert;
                }
            }
        }

        connector
    }

    /// Create a keyhole between hole `start` and outer polygon `connector`.
    fn join_polygons(&mut self, start: usize, connector: usize) {
        let new_start = self.polygon.len();
        self.polygon.push(self.polygon[start].clone());
        let new_connector = self.polygon.len();
        self.polygon.push(self.polygon[connector].clone());

        let start_right = self.polygon[start].right;
        self.polygon[start_right].left = new_start;
        let connector_left = self.polygon[connector].left;
        self.polygon[connector_left].right = new_connector;

        self.link(start, connector);
        self.link(new_connector, new_start);

        self.clip_if_degenerate(start);
        self.clip_if_degenerate(new_start);
        self.clip_if_degenerate(connector);
        self.clip_if_degenerate(new_connector);
    }

    /// Update ear queue entry for vert v.
    fn process_ear(&mut self, v: usize, collider: &IdxCollider) {
        // Lazy-delete existing queue entry
        self.polygon[v].ear_version = self.polygon[v].ear_version.wrapping_add(1);

        if self.vert_is_short(v) {
            self.polygon[v].cost = K_BEST;
            let version = self.polygon[v].ear_version;
            self.ears_queue.push(EarEntry { cost: K_BEST, idx: v, version });
        } else if self.vert_is_convex(v, 2.0 * self.epsilon) {
            let cost = self.vert_ear_cost(v, collider);
            self.polygon[v].cost = cost;
            let version = self.polygon[v].ear_version;
            self.ears_queue.push(EarEntry { cost, idx: v, version });
        } else {
            self.polygon[v].cost = 1.0; // reflex, not an ear
        }
    }

    /// Build a 2D KD-tree collider of all polygon verts for ear cost queries.
    fn vert_collider(&self, start: usize) -> IdxCollider {
        let verts = match self.loop_verts(start) {
            None => return IdxCollider { points: Vec::new(), itr: Vec::new() },
            Some(v) => v,
        };

        let mut itr = Vec::with_capacity(verts.len());
        let mut points = Vec::with_capacity(verts.len());
        for (k, &v) in verts.iter().enumerate() {
            points.push(PolyVert { pos: self.polygon[v].pos, idx: k as i32 });
            itr.push(v);
        }

        build_two_d_tree(&mut points);
        IdxCollider { points, itr }
    }

    /// Main ear-clipping loop for a simple polygon.
    fn triangulate_poly(&mut self, start: usize) {
        let collider = self.vert_collider(start);
        if collider.itr.is_empty() {
            return;
        }

        self.ears_queue.clear();
        let mut num_tri: i32 = -2;

        // Queue all verts and count
        let verts = match self.loop_verts(start) {
            None => return,
            Some(v) => v,
        };

        for &v in &verts {
            self.process_ear(v, &collider);
            num_tri += 1;
        }

        let mut last_v = *verts.last().unwrap();

        while num_tri > 0 {
            // Pop the cheapest valid ear
            let ear = loop {
                match self.ears_queue.pop() {
                    None => {
                        // Fallback: use last_v
                        break last_v;
                    }
                    Some(entry) => {
                        if entry.version == self.polygon[entry.idx].ear_version {
                            break entry.idx;
                        }
                        // Stale entry, try next
                    }
                }
            };

            self.clip_ear(ear);
            num_tri -= 1;

            let ear_left = self.polygon[ear].left;
            let ear_right = self.polygon[ear].right;
            self.process_ear(ear_left, &collider);
            self.process_ear(ear_right, &collider);
            last_v = ear_right;
        }
    }
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Triangulates indexed polygons. Returns triangle indices referencing original vertex indices.
pub fn triangulate_idx(polys: &PolygonsIdx, epsilon: f64, allow_convex: bool) -> Vec<IVec3Out> {
    if polys.is_empty() || polys.iter().all(|p| p.is_empty()) {
        return Vec::new();
    }
    if allow_convex && is_convex(polys, epsilon) {
        return triangulate_convex(polys);
    }
    let (triangles, _updated_eps) = EarClip::new(polys, epsilon).triangulate();
    triangles
}

/// Triangulates unindexed polygons. Vertices are indexed sequentially across all contours.
pub fn triangulate(polygons: &Polygons, epsilon: f64, allow_convex: bool) -> Vec<IVec3Out> {
    let mut idx = 0i32;
    let mut polygons_indexed: PolygonsIdx = Vec::new();
    for poly in polygons {
        let simple: SimplePolygonIdx = poly
            .iter()
            .map(|&pos| {
                let v = PolyVert { pos, idx };
                idx += 1;
                v
            })
            .collect();
        polygons_indexed.push(simple);
    }
    triangulate_idx(&polygons_indexed, epsilon, allow_convex)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linalg::Vec2;
    use crate::types::PolyVert;

    fn make_poly(coords: &[(f64, f64)]) -> SimplePolygonIdx {
        coords
            .iter()
            .enumerate()
            .map(|(i, &(x, y))| PolyVert { pos: Vec2::new(x, y), idx: i as i32 })
            .collect()
    }

    #[test]
    fn test_ccw_basic() {
        // CCW triangle
        assert_eq!(ccw(Vec2::new(0.0, 0.0), Vec2::new(1.0, 0.0), Vec2::new(0.5, 1.0), 1e-10), 1);
        // CW triangle
        assert_eq!(ccw(Vec2::new(0.0, 0.0), Vec2::new(0.5, 1.0), Vec2::new(1.0, 0.0), 1e-10), -1);
        // Colinear
        assert_eq!(ccw(Vec2::new(0.0, 0.0), Vec2::new(1.0, 0.0), Vec2::new(2.0, 0.0), 1e-10), 0);
    }

    #[test]
    fn test_triangle_ccw() {
        let poly = vec![make_poly(&[(0.0, 0.0), (1.0, 0.0), (0.5, 1.0)])];
        let tris = triangulate_idx(&poly, 1e-10, true);
        assert_eq!(tris.len(), 1);
        assert_eq!(tris[0], IVec3Out::new(0, 1, 2));
    }

    #[test]
    fn test_square_ccw() {
        // Unit square, CCW
        let poly = vec![make_poly(&[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)])];
        let tris = triangulate_idx(&poly, 1e-10, true);
        assert_eq!(tris.len(), 2);
        // Total triangles should cover the whole square
        for tri in &tris {
            assert!(tri.x >= 0 && tri.x < 4);
            assert!(tri.y >= 0 && tri.y < 4);
            assert!(tri.z >= 0 && tri.z < 4);
        }
    }

    #[test]
    fn test_pentagon_ccw() {
        use std::f64::consts::PI;
        let n = 5usize;
        let coords: Vec<(f64, f64)> = (0..n)
            .map(|i| {
                let a = 2.0 * PI * i as f64 / n as f64;
                (a.cos(), a.sin())
            })
            .collect();
        let poly = vec![make_poly(&coords)];
        let tris = triangulate_idx(&poly, 1e-10, false);
        assert_eq!(tris.len(), n - 2);
    }

    #[test]
    fn test_convex_fast_path() {
        use std::f64::consts::PI;
        let n = 8usize;
        let coords: Vec<(f64, f64)> = (0..n)
            .map(|i| {
                let a = 2.0 * PI * i as f64 / n as f64;
                (a.cos(), a.sin())
            })
            .collect();
        let poly = vec![make_poly(&coords)];
        // is_convex should return true → triangulate_convex used
        assert!(is_convex(&poly, 1e-10));
        let tris = triangulate_idx(&poly, 1e-10, true);
        assert_eq!(tris.len(), n - 2);
    }

    #[test]
    fn test_triangulate_unindexed() {
        let poly = vec![vec![
            Vec2::new(0.0, 0.0),
            Vec2::new(1.0, 0.0),
            Vec2::new(1.0, 1.0),
            Vec2::new(0.0, 1.0),
        ]];
        let tris = triangulate(&poly, 1e-10, true);
        assert_eq!(tris.len(), 2);
    }

    #[test]
    fn test_kd_tree_query() {
        let mut points: Vec<PolyVert> = (0..20)
            .map(|i| PolyVert {
                pos: Vec2::new(i as f64, (i % 5) as f64),
                idx: i,
            })
            .collect();
        build_two_d_tree(&mut points);

        let query_rect = Rect::from_points(Vec2::new(3.0, 0.0), Vec2::new(8.0, 3.0));
        let mut found = Vec::new();
        query_two_d_tree(&points, query_rect, |p| found.push(p.idx));
        // Points with x in [3,8] and y in [0,3]
        found.sort();
        // i=3: x=3, y=3 ✓; i=4: x=4, y=4 ✗; i=5: x=5, y=0 ✓; i=6: x=6, y=1 ✓;
        // i=7: x=7, y=2 ✓; i=8: x=8, y=3 ✓
        // (y = i%5: 3→3, 5→0, 6→1, 7→2, 8→3)
        assert!(found.contains(&3));
        assert!(found.contains(&5));
        assert!(found.contains(&6));
        assert!(found.contains(&7));
        assert!(found.contains(&8));
        assert!(!found.contains(&4)); // y=4 outside
    }
}
