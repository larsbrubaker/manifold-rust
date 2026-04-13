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

use clipper2_rust::{difference_d, inflate_paths_d, intersect_d, minkowski_sum_d, union_d, EndType, FillRule, JoinType, PathD, PathsD, Point};

use crate::linalg::Vec2;
use crate::math;
use crate::types::{Polygons, Rect};

#[derive(Clone, Debug, Default)]
pub struct CrossSection {
    polygons: Polygons,
}

fn to_paths(polygons: &Polygons) -> PathsD {
    polygons
        .iter()
        .map(|poly| poly.iter().map(|p| Point::new(p.x, p.y)).collect::<PathD>())
        .collect()
}

fn from_paths(paths: &PathsD) -> Polygons {
    paths.iter()
        .map(|path| path.iter().map(|p| Vec2::new(p.x, p.y)).collect())
        .collect()
}

fn signed_area(poly: &[Vec2]) -> f64 {
    let mut area = 0.0;
    for i in 0..poly.len() {
        let j = (i + 1) % poly.len();
        area += poly[i].x * poly[j].y - poly[j].x * poly[i].y;
    }
    area * 0.5
}

impl CrossSection {
    pub fn new(polygons: Polygons) -> Self {
        Self { polygons }
    }

    /// Create a CrossSection from polygons, normalizing via Clipper2 Union.
    /// Mirrors C++ CrossSection(Polygons, FillRule) constructor which runs
    /// the polygons through C2::Union to merge overlapping regions.
    pub fn from_polygons_fill(polygons: Polygons) -> Self {
        if polygons.is_empty() {
            return Self::default();
        }
        let paths = to_paths(&polygons);
        let empty = PathsD::new();
        let result = union_d(&paths, &empty, FillRule::NonZero, 6);
        Self { polygons: from_paths(&result) }
    }

    pub fn square(size: f64) -> Self {
        if size <= 0.0 {
            return Self { polygons: vec![] };
        }
        Self::new(vec![vec![
            Vec2::new(0.0, 0.0),
            Vec2::new(size, 0.0),
            Vec2::new(size, size),
            Vec2::new(0.0, size),
        ]])
    }

    pub fn circle(radius: f64, segments: i32) -> Self {
        if radius <= 0.0 {
            return Self { polygons: vec![] };
        }
        let segments = segments.max(3) as usize;
        let poly = (0..segments)
            .map(|i| {
                let a = (i as f64 / segments as f64) * std::f64::consts::TAU;
                Vec2::new(radius * math::cos(a), radius * math::sin(a))
            })
            .collect();
        Self::new(vec![poly])
    }

    pub fn to_polygons(&self) -> Polygons {
        self.polygons.clone()
    }

    pub fn translate(&self, v: Vec2) -> Self {
        Self::new(
            self.polygons
                .iter()
                .map(|poly| poly.iter().map(|p| *p + v).collect())
                .collect(),
        )
    }

    pub fn area(&self) -> f64 {
        self.polygons.iter().map(|p| signed_area(p).abs()).sum()
    }

    pub fn bounds(&self) -> Rect {
        let mut rect = Rect::new();
        for poly in &self.polygons {
            for &p in poly {
                rect.union_point(p);
            }
        }
        rect
    }

    pub fn union(&self, other: &Self) -> Self {
        Self::new(from_paths(&union_d(&to_paths(&self.polygons), &to_paths(&other.polygons), FillRule::NonZero, 6)))
    }

    pub fn intersection(&self, other: &Self) -> Self {
        Self::new(from_paths(&intersect_d(&to_paths(&self.polygons), &to_paths(&other.polygons), FillRule::NonZero, 6)))
    }

    pub fn difference(&self, other: &Self) -> Self {
        Self::new(from_paths(&difference_d(&to_paths(&self.polygons), &to_paths(&other.polygons), FillRule::NonZero, 6)))
    }

    pub fn scale(&self, v: Vec2) -> Self {
        Self::new(
            self.polygons
                .iter()
                .map(|poly| poly.iter().map(|p| Vec2::new(p.x * v.x, p.y * v.y)).collect())
                .collect(),
        )
    }

    pub fn rotate(&self, degrees: f64) -> Self {
        let rad = degrees.to_radians();
        let c = math::cos(rad);
        let s = math::sin(rad);
        Self::new(
            self.polygons
                .iter()
                .map(|poly| {
                    poly.iter()
                        .map(|p| Vec2::new(p.x * c - p.y * s, p.x * s + p.y * c))
                        .collect()
                })
                .collect(),
        )
    }

    /// Mirror through a line perpendicular to the given axis vector.
    /// Matches C++ `CrossSection::Mirror(ax)` which uses `I - 2*n*n^T`.
    pub fn mirror(&self, axis: Vec2) -> Self {
        let len_sq = axis.x * axis.x + axis.y * axis.y;
        if len_sq < 1e-20 {
            return Self::default();
        }
        // Reflection matrix: R = I - 2*n*n^T where n = normalize(axis)
        let nx = axis.x / len_sq.sqrt();
        let ny = axis.y / len_sq.sqrt();
        let r00 = 1.0 - 2.0 * nx * nx;
        let r01 = -2.0 * nx * ny;
        let r10 = -2.0 * nx * ny;
        let r11 = 1.0 - 2.0 * ny * ny;
        Self::new(
            self.polygons
                .iter()
                .map(|poly| {
                    // Mirror reverses winding, so reverse the polygon
                    poly.iter()
                        .rev()
                        .map(|p| Vec2::new(r00 * p.x + r01 * p.y, r10 * p.x + r11 * p.y))
                        .collect()
                })
                .collect(),
        )
    }

    pub fn is_empty(&self) -> bool {
        self.polygons.is_empty() || self.polygons.iter().all(|p| p.len() < 3)
    }

    pub fn num_vert(&self) -> usize {
        self.polygons.iter().map(|p| p.len()).sum()
    }

    pub fn num_contour(&self) -> usize {
        self.polygons.iter().filter(|p| p.len() >= 3).count()
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
            let sa = signed_area(poly);
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

    pub fn offset(&self, delta: f64) -> Self {
        Self::new(from_paths(&inflate_paths_d(&to_paths(&self.polygons), delta, JoinType::Round, EndType::Polygon, 2.0, 6, 0.0)))
    }

    /// Offset with explicit join type and segment count.
    /// join_type: 0=Square, 1=Round, 2=Miter
    pub fn offset_with_params(
        &self,
        delta: f64,
        join_type: i32,
        miter_limit: f64,
        circular_segments: i32,
    ) -> Self {
        let jt = match join_type {
            0 => JoinType::Square,
            2 => JoinType::Miter,
            _ => JoinType::Round,
        };
        Self::new(from_paths(&inflate_paths_d(
            &to_paths(&self.polygons),
            delta,
            jt,
            EndType::Polygon,
            miter_limit,
            6,
            0.0,
        )))
    }

    pub fn minkowski_sum(&self, other: &Self) -> Self {
        let mut result = Vec::new();
        for a in to_paths(&self.polygons) {
            for b in to_paths(&other.polygons) {
                result.extend(minkowski_sum_d(&a, &b, true, 6));
            }
        }
        Self::new(from_paths(&result))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cross_section_area_bounds() {
        let cs = CrossSection::square(2.0);
        assert!((cs.area() - 4.0).abs() < 1e-10);
        let b = cs.bounds();
        assert!((b.max.x - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_cross_section_boolean() {
        let a = CrossSection::square(2.0);
        let b = CrossSection::square(2.0).translate(Vec2::new(1.0, 0.0));
        assert!(a.intersection(&b).area() > 0.9);
        assert!(a.union(&b).area() > a.area());
        assert!(a.difference(&b).area() < a.area());
    }

    #[test]
    fn test_cross_section_offset() {
        let a = CrossSection::square(1.0);
        let b = a.offset(0.25);
        assert!(b.area() > a.area());
    }

    /// C++ TEST(CrossSection, Square) — cube from extrusion matches cube
    #[test]
    fn test_cpp_cross_section_square() {
        let cs = CrossSection::square(5.0);
        let a = crate::manifold::Manifold::cube(
            crate::linalg::Vec3::new(5.0, 5.0, 5.0),
            false,
        );
        let b = crate::manifold::Manifold::extrude(
            &cs.to_polygons(),
            5.0,
            0,
            0.0,
            crate::linalg::Vec2::new(1.0, 1.0),
        );
        let diff = a.difference(&b);
        assert!(
            diff.volume().abs() < 1e-6,
            "CrossSection square extrusion should match cube, diff volume: {}",
            diff.volume()
        );
    }

    /// C++ TEST(CrossSection, Empty) — empty cross section from empty polygons
    #[test]
    fn test_cpp_cross_section_empty() {
        let polys: crate::types::Polygons = vec![vec![], vec![]];
        let cs = CrossSection::new(polys);
        assert!(cs.area().abs() < 1e-10, "CrossSection from empty polygons should have zero area");
    }
}
