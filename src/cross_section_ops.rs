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

// cross_section_ops.rs — the Clipper2-backed operations on CrossSection:
// two-operand booleans, batch booleans and compose, decompose, simplify,
// offset, Minkowski sum, warp and convex hull.
//
// Ports the corresponding members of src/cross_section/cross_section.cpp.
// A child module of cross_section.rs (which owns the struct, constructors,
// transforms and queries) so these methods keep access to the private
// `polygons` field and callers keep the same `CrossSection::...` paths.

use clipper2_rust::{
    difference_d, inflate_paths_d, intersect_d, minkowski_sum_d, simplify_paths, union_d, EndType,
    FillRule, JoinType, PathsD,
};

use super::{contour_area, from_paths, path_area, to_paths, CrossSection, PRECISION};
use crate::linalg::Vec2;
use crate::math;
use crate::types::{OpType, Quality, Rect};

impl CrossSection {
    pub fn union(&self, other: &Self) -> Self {
        Self::new(from_paths(&union_d(
            &to_paths(&self.polygons),
            &to_paths(&other.polygons),
            FillRule::Positive,
            PRECISION,
        )))
    }

    pub fn intersection(&self, other: &Self) -> Self {
        Self::new(from_paths(&intersect_d(
            &to_paths(&self.polygons),
            &to_paths(&other.polygons),
            FillRule::Positive,
            PRECISION,
        )))
    }

    pub fn difference(&self, other: &Self) -> Self {
        Self::new(from_paths(&difference_d(
            &to_paths(&self.polygons),
            &to_paths(&other.polygons),
            FillRule::Positive,
            PRECISION,
        )))
    }
    /// Decompose into connected components. Each component maintains its
    /// contours (outer boundary + holes).
    pub fn decompose(&self) -> Vec<Self> {
        // Simple decomposition: use clipper union to normalize, then separate
        // non-overlapping groups by bounding box.
        let normalized = self.union(&Self::default());
        let polys = &normalized.polygons;
        if polys.is_empty() {
            return vec![];
        }

        // Group polygons: outer polygons are CCW (positive area), holes are CW.
        // Each outer polygon starts a new component, holes are assigned to the
        // outer polygon whose bbox contains them.
        let mut outers: Vec<(usize, Rect)> = Vec::new();
        let mut holes: Vec<(usize, Vec2)> = Vec::new();

        for (i, poly) in polys.iter().enumerate() {
            if poly.len() < 3 {
                continue;
            }
            let sa = contour_area(poly);
            if sa >= 0.0 {
                // Outer (CCW in our convention)
                let mut r = Rect::new();
                for &p in poly {
                    r.union_point(p);
                }
                outers.push((i, r));
            } else {
                // Hole — use first point as representative
                holes.push((i, poly[0]));
            }
        }

        let mut components: Vec<Vec<usize>> = outers.iter().map(|(i, _)| vec![*i]).collect();

        for (hole_idx, pt) in &holes {
            // Find smallest outer bbox that contains this hole's representative point
            let mut best = None;
            let mut best_area = f64::MAX;
            for (ci, (_, rect)) in outers.iter().enumerate() {
                if rect.contains_point(*pt) {
                    let a = (rect.max.x - rect.min.x) * (rect.max.y - rect.min.y);
                    if a < best_area {
                        best_area = a;
                        best = Some(ci);
                    }
                }
            }
            if let Some(ci) = best {
                components[ci].push(*hole_idx);
            }
        }

        components
            .into_iter()
            .map(|indices| {
                let component_polys = indices.into_iter().map(|i| polys[i].clone()).collect();
                Self::new(component_polys)
            })
            .collect()
    }

    /// Simplify contours by removing near-collinear vertices.
    /// Mirrors C++ CrossSection::Simplify(epsilon=1e-6): normalizes via union,
    /// filters tiny polygons, then applies SimplifyPaths with epsilon.
    pub fn simplify(&self, epsilon: f64) -> Self {
        if self.polygons.is_empty() {
            return Self::default();
        }
        // Normalize via union (removes overlaps/inversions).
        let paths = to_paths(&self.polygons);
        let unified = union_d(&paths, &PathsD::new(), FillRule::Positive, PRECISION);
        // Filter out contours smaller than epsilon (area vs bounding box).
        let filtered: PathsD = unified
            .into_iter()
            .filter(|poly| {
                let a = path_area(poly).abs();
                // Compute bounding box max extent
                let (mut min_x, mut min_y) = (f64::MAX, f64::MAX);
                let (mut max_x, mut max_y) = (f64::MIN, f64::MIN);
                for p in poly {
                    if p.x < min_x {
                        min_x = p.x;
                    }
                    if p.x > max_x {
                        max_x = p.x;
                    }
                    if p.y < min_y {
                        min_y = p.y;
                    }
                    if p.y > max_y {
                        max_y = p.y;
                    }
                }
                let max_size = (max_x - min_x).max(max_y - min_y);
                a > max_size * epsilon
            })
            .collect();
        let simplified = simplify_paths(&filtered, epsilon, true);
        Self::new(from_paths(&simplified))
    }

    /// Offset with the C++ `CrossSection::Offset` defaults: Round joins,
    /// miter_limit 2.0, circularSegments 0 (segments from Quality).
    pub fn offset(&self, delta: f64) -> Self {
        self.offset_with_params(delta, 1, 2.0, 0)
    }

    /// Offset with explicit join type and segment count.
    /// join_type: 0=Square, 1=Round, 2=Miter, 3=Bevel (the C++
    /// `CrossSection::JoinType` enumerator order). Other codes fall through to
    /// Square, the value C++ `jt()` starts from before its switch.
    pub fn offset_with_params(
        &self,
        delta: f64,
        join_type: i32,
        miter_limit: f64,
        circular_segments: i32,
    ) -> Self {
        let jt = match join_type {
            1 => JoinType::Round,
            2 => JoinType::Miter,
            3 => JoinType::Bevel,
            _ => JoinType::Square,
        };
        // For round joins, compute arc_tolerance from circular_segments (or,
        // when it is <= 2, Quality's count for radius delta) to get the exact
        // segment count. Matches C++ CrossSection::Offset:
        //   arc_tol = (math::cos(π/n) - 1) * -|delta|
        let arc_tol = if jt == JoinType::Round {
            let n = if circular_segments > 2 {
                circular_segments
            } else {
                Quality::get_circular_segments(delta)
            };
            let abs_delta = delta.abs();
            (math::cos(std::f64::consts::PI / n as f64) - 1.0) * -abs_delta
        } else {
            0.0
        };
        Self::new(from_paths(&inflate_paths_d(
            &to_paths(&self.polygons),
            delta,
            jt,
            EndType::Polygon,
            miter_limit,
            PRECISION,
            arc_tol,
        )))
    }

    pub fn minkowski_sum(&self, other: &Self) -> Self {
        let mut result = Vec::new();
        for a in to_paths(&self.polygons) {
            for b in to_paths(&other.polygons) {
                result.extend(minkowski_sum_d(&a, &b, true, PRECISION));
            }
        }
        Self::new(from_paths(&result))
    }

    /// Apply a function to every vertex in-place.
    pub fn warp<F: FnMut(&mut Vec2)>(&self, mut f: F) -> Self {
        let polys = self
            .polygons
            .iter()
            .map(|poly| {
                poly.iter()
                    .map(|&v| {
                        let mut v2 = v;
                        f(&mut v2);
                        v2
                    })
                    .collect()
            })
            .collect();
        Self { polygons: polys }
    }

    /// Batch boolean operation on a slice of CrossSections.
    /// OpType::Add = union, Subtract = difference, Intersect = intersection.
    pub fn batch_boolean(sections: &[Self], op: OpType) -> Self {
        if sections.is_empty() {
            return Self::default();
        }
        match op {
            OpType::Add => {
                let mut paths = PathsD::new();
                for s in sections {
                    for p in to_paths(&s.polygons) {
                        paths.push(p);
                    }
                }
                let empty = PathsD::new();
                Self {
                    polygons: from_paths(&union_d(&paths, &empty, FillRule::Positive, PRECISION)),
                }
            }
            OpType::Subtract => {
                let mut result = sections[0].clone();
                for s in &sections[1..] {
                    result = result.difference(s);
                }
                result
            }
            OpType::Intersect => {
                let mut result = sections[0].clone();
                for s in &sections[1..] {
                    result = result.intersection(s);
                }
                result
            }
        }
    }

    /// Compute convex hull of all vertices in a slice of CrossSections.
    pub fn hull_cross_sections(sections: &[Self]) -> Self {
        let points: Vec<Vec2> = sections
            .iter()
            .flat_map(|s| s.polygons.iter().flat_map(|p| p.iter().cloned()))
            .collect();
        Self::hull_points(&points)
    }

    /// Compute convex hull of a set of 2D points (Andrew's monotone chain).
    pub fn hull_points(points: &[Vec2]) -> Self {
        if points.len() < 3 {
            return Self::default();
        }
        let mut pts: Vec<Vec2> = points.to_vec();
        pts.sort_by(|a, b| {
            a.x.partial_cmp(&b.x)
                .unwrap()
                .then(a.y.partial_cmp(&b.y).unwrap())
        });
        pts.dedup_by(|a, b| (a.x - b.x).abs() < 1e-10 && (a.y - b.y).abs() < 1e-10);

        let cross = |o: Vec2, a: Vec2, b: Vec2| -> f64 {
            (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x)
        };

        let n = pts.len();
        if n < 3 {
            return Self::default();
        }
        let mut hull: Vec<Vec2> = Vec::with_capacity(2 * n);
        // Lower hull
        for &p in &pts {
            while hull.len() >= 2 && cross(hull[hull.len() - 2], hull[hull.len() - 1], p) <= 0.0 {
                hull.pop();
            }
            hull.push(p);
        }
        // Upper hull
        let lower_len = hull.len();
        for &p in pts.iter().rev() {
            while hull.len() > lower_len
                && cross(hull[hull.len() - 2], hull[hull.len() - 1], p) <= 0.0
            {
                hull.pop();
            }
            hull.push(p);
        }
        hull.pop(); // last point == first
        if hull.len() < 3 {
            return Self::default();
        }
        Self::new(vec![hull])
    }

    /// Compose (merge) multiple CrossSections by combining all their contours.
    /// Matches C++ CrossSection::Compose(vector<CrossSection>) which unions all polygons.
    pub fn compose(sections: &[Self]) -> Self {
        let all: Vec<Vec<Vec2>> = sections
            .iter()
            .flat_map(|s| s.polygons.iter().cloned())
            .collect();
        if all.is_empty() {
            return Self::default();
        }
        Self::from_polygons_fill(all)
    }
}
