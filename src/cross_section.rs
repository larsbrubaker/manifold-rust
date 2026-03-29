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

    pub fn square(size: f64) -> Self {
        Self::new(vec![vec![
            Vec2::new(0.0, 0.0),
            Vec2::new(size, 0.0),
            Vec2::new(size, size),
            Vec2::new(0.0, size),
        ]])
    }

    pub fn circle(radius: f64, segments: i32) -> Self {
        let segments = segments.max(3) as usize;
        let poly = (0..segments)
            .map(|i| {
                let a = (i as f64 / segments as f64) * std::f64::consts::TAU;
                Vec2::new(radius * a.cos(), radius * a.sin())
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

    pub fn offset(&self, delta: f64) -> Self {
        Self::new(from_paths(&inflate_paths_d(&to_paths(&self.polygons), delta, JoinType::Round, EndType::Polygon, 2.0, 6, 0.0)))
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
}
