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

// cross_section.rs — the 2D CrossSection type: a set of contours (outer
// boundaries CCW, holes CW) stored as Clipper2-ready polygons. Ports
// src/cross_section/cross_section.cpp of the C++ reference.
//
// This file owns the struct, the Clipper2 path conversions and area helpers,
// the constructors, the affine transforms and the read-only queries. The
// Clipper2-backed operations (booleans, decompose, simplify, offset, warp,
// hull, batch booleans) live in cross_section_ops.rs, a child module so they
// keep access to the private `polygons` field; tests are in
// cross_section_tests.rs. Manifold::slice / project (manifold.rs) and the
// extrude/revolve constructors (constructors.rs) consume CrossSections.

use clipper2_rust::{union_d, FillRule, PathD, PathsD, Point};

use crate::linalg::Vec2;
use crate::math;
use crate::types::{Polygons, Rect};

/// Decimal places Clipper2 keeps when scaling to integer coordinates; mirrors
/// `precision_` in C++ cross_section.cpp, passed to every Clipper2 call.
const PRECISION: i32 = 8;

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
    paths
        .iter()
        .map(|path| path.iter().map(|p| Vec2::new(p.x, p.y)).collect())
        .collect()
}

/// Exact port of Clipper2's `Area(const Path<T>&)` (clipper.core.h at commit
/// 46f6391, the version C++ Manifold pins). Clipper2 walks the trapezoid form
/// over edges (n-1,0), (0,1), ..., (n-2,n-1), accumulating
/// `(prev.y + cur.y) * (prev.x - cur.x)` in that order; its two-edges-per-step
/// unrolling does not change the order. This differs in the last bits from a
/// shoelace sum (and from clipper2_rust's `area`), so every place C++ calls
/// `C2::Area` uses this instead.
fn clipper2_area_by<F: Fn(usize) -> (f64, f64)>(cnt: usize, pt: F) -> f64 {
    if cnt < 3 {
        return 0.0;
    }
    let mut a = 0.0;
    let mut prev = cnt - 1;
    for cur in 0..cnt {
        let (px, py) = pt(prev);
        let (cx, cy) = pt(cur);
        a += (py + cy) * (px - cx);
        prev = cur;
    }
    a * 0.5
}

fn contour_area(poly: &[Vec2]) -> f64 {
    clipper2_area_by(poly.len(), |i| (poly[i].x, poly[i].y))
}

fn path_area(path: &PathD) -> f64 {
    clipper2_area_by(path.len(), |i| (path[i].x, path[i].y))
}

impl CrossSection {
    pub fn new(polygons: Polygons) -> Self {
        Self { polygons }
    }

    /// Create a CrossSection from polygons, normalizing via Clipper2 Union.
    /// Mirrors C++ CrossSection(Polygons, FillRule) constructor with its
    /// default FillRule::Positive, which runs the polygons through C2::Union
    /// to merge overlapping regions.
    pub fn from_polygons_fill(polygons: Polygons) -> Self {
        if polygons.is_empty() {
            return Self::default();
        }
        let paths = to_paths(&polygons);
        let empty = PathsD::new();
        let result = union_d(&paths, &empty, FillRule::Positive, PRECISION);
        Self {
            polygons: from_paths(&result),
        }
    }

    /// Create a CrossSection from a Rect (axis-aligned rectangle).
    /// Matches C++ CrossSection(Rect) constructor.
    pub fn from_rect(rect: &Rect) -> Self {
        if rect.is_empty() {
            return Self::default();
        }
        Self::new(vec![vec![
            Vec2::new(rect.min.x, rect.min.y),
            Vec2::new(rect.max.x, rect.min.y),
            Vec2::new(rect.max.x, rect.max.y),
            Vec2::new(rect.min.x, rect.max.y),
        ]])
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

    /// Create a rectangle of size (w, h), optionally centered at origin.
    /// Matches C++ CrossSection::Square(vec2, center).
    pub fn square_vec2(size: Vec2, center: bool) -> Self {
        let (w, h) = (size.x, size.y);
        if w <= 0.0 || h <= 0.0 {
            return Self { polygons: vec![] };
        }
        let (x0, y0, x1, y1) = if center {
            (-w / 2.0, -h / 2.0, w / 2.0, h / 2.0)
        } else {
            (0.0, 0.0, w, h)
        };
        Self::new(vec![vec![
            Vec2::new(x0, y0),
            Vec2::new(x1, y0),
            Vec2::new(x1, y1),
            Vec2::new(x0, y1),
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

    /// Net enclosed area: the sum of signed contour areas, so CCW outers add
    /// and CW holes subtract. Mirrors C++ `CrossSection::Area`, i.e.
    /// Clipper2's `Area(Paths)`: an explicit fold from +0.0 in contour order,
    /// so an empty section yields +0.0 rather than `.sum()`'s -0.0.
    pub fn area(&self) -> f64 {
        self.polygons.iter().fold(0.0, |a, p| a + contour_area(p))
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

    pub fn scale(&self, v: Vec2) -> Self {
        Self::new(
            self.polygons
                .iter()
                .map(|poly| {
                    poly.iter()
                        .map(|p| Vec2::new(p.x * v.x, p.y * v.y))
                        .collect()
                })
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

    /// Create CrossSection from a simple polygon with a specified fill rule.
    /// fill_rule: 0=EvenOdd, 1=NonZero, 2=Positive, 3=Negative (the C++
    /// `CrossSection::FillRule` enumerator order). Other codes fall through to
    /// EvenOdd, the value C++ `fr()` starts from before its switch.
    pub fn from_polygon_with_fill_rule(polygon: Vec<Vec2>, fill_rule: i32) -> Self {
        let fr = match fill_rule {
            1 => FillRule::NonZero,
            2 => FillRule::Positive,
            3 => FillRule::Negative,
            _ => FillRule::EvenOdd,
        };
        let path: PathD = polygon.iter().map(|v| Point::new(v.x, v.y)).collect();
        let paths = PathsD::from(vec![path]);
        let empty = PathsD::new();
        let result = union_d(&paths, &empty, fr, PRECISION);
        Self {
            polygons: from_paths(&result),
        }
    }
}

#[path = "cross_section_ops.rs"]
mod ops;

#[cfg(test)]
#[path = "cross_section_tests.rs"]
mod tests;
