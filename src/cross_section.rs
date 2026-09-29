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

use crate::types::OpType;
use clipper2_rust::{
    difference_d, inflate_paths_d, intersect_d, minkowski_sum_d, simplify_paths, union_d, EndType,
    FillRule, JoinType, PathD, PathsD, Point,
};

use crate::linalg::Vec2;
use crate::math;
use crate::types::{Polygons, Quality, Rect};

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

    /// A 10x10 square minus an inner 4x4 square yields an outer contour and a
    /// hole; the hole's area is subtracted, giving 100 - 16 = 84.
    #[test]
    fn test_cross_section_area_subtracts_holes() {
        let outer = CrossSection::square(10.0);
        let hole = CrossSection::square(4.0).translate(Vec2::new(3.0, 3.0));
        let ring = outer.difference(&hole);
        assert_eq!(
            ring.num_contour(),
            2,
            "difference should yield outer + hole"
        );
        assert_eq!(ring.area(), 84.0);
    }

    /// C++ TEST(CrossSection, Square) — cube from extrusion matches cube
    #[test]
    fn test_cpp_cross_section_square() {
        let cs = CrossSection::square(5.0);
        let a = crate::manifold::Manifold::cube(crate::linalg::Vec3::new(5.0, 5.0, 5.0), false);
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
        assert!(
            cs.area().abs() < 1e-10,
            "CrossSection from empty polygons should have zero area"
        );
    }

    /// C++ `CrossSection::Area` is `C2::Area(paths)`, which starts from
    /// `a = 0.0` and adds each contour, so a section with no contours reports
    /// +0.0 (an iterator `.sum()` of no f64s yields -0.0).
    #[test]
    fn test_cross_section_area_empty_is_positive_zero() {
        assert_eq!(CrossSection::default().area().to_bits(), 0.0f64.to_bits());
    }

    /// Off-origin polygons separate Clipper2's trapezoid Area from a shoelace
    /// sum in the last bits. Expected bits come from compiling Clipper2 commit
    /// 46f6391's `clipper.core.h` `Area` (MSVC /O2) on these exact
    /// coordinates: odd count (33) and even count (first 32 points).
    #[test]
    fn test_cross_section_area_matches_clipper2_bits() {
        let cs = CrossSection::circle(1.0, 33).translate(Vec2::new(100.0, -50.0));
        assert_eq!(cs.area().to_bits(), 0x4008fb2d94a5b1f1);
        let mut even = cs.to_polygons();
        even[0].pop();
        assert_eq!(CrossSection::new(even).area().to_bits(), 0x4008f42c81cc8074);
    }

    /// C++ runs every Clipper2 op at `precision_ = 8` decimal places. ClipperD
    /// scales by the power of two above 10^precision (2^27 at 8, 2^20 at 6),
    /// so x = 1.00000012 snaps to 1 + 16 * 2^-27 = 1 + 2^-23, where precision
    /// 6 would round the 1.2e-7 feature away to 1.0.
    #[test]
    fn test_cross_section_union_keeps_eighth_decimal() {
        let x = 1.000_000_12;
        let snapped = 1.0 + 2f64.powi(-23);
        let a = CrossSection::new(vec![vec![
            Vec2::new(0.0, 0.0),
            Vec2::new(x, 0.0),
            Vec2::new(x, 1.0),
            Vec2::new(0.0, 1.0),
        ]]);
        let u = a.union(&CrossSection::default());
        assert_eq!(u.bounds().max.x, snapped);
        let f = CrossSection::from_polygons_fill(a.to_polygons());
        assert_eq!(f.bounds().max.x, snapped);
    }

    /// C++ booleans use FillRule::Positive and the Polygons constructor
    /// defaults to Positive, so a clockwise (negative) contour fills nothing.
    #[test]
    fn test_cross_section_booleans_use_positive_fill() {
        let cw = CrossSection::new(vec![vec![
            Vec2::new(0.0, 0.0),
            Vec2::new(0.0, 1.0),
            Vec2::new(1.0, 1.0),
            Vec2::new(1.0, 0.0),
        ]]);
        assert!(cw.union(&CrossSection::default()).is_empty());
        assert!(CrossSection::from_polygons_fill(cw.to_polygons()).is_empty());
        let batch = CrossSection::batch_boolean(&[cw.clone(), cw.clone()], OpType::Add);
        assert!(batch.is_empty());
        let sq = CrossSection::square(1.0);
        assert!(sq.intersection(&cw).is_empty());
        assert_eq!(sq.difference(&cw).area(), 1.0);
    }

    /// C++ `fr()` starts from EvenOdd and only overrides it for the three
    /// other enumerators, and `jt()` likewise starts from Square; unknown
    /// integer codes fall through to those initial values.
    #[test]
    fn test_cross_section_unknown_codes_match_cpp_defaults() {
        let star: Vec<Vec2> = (0..5)
            .map(|i| {
                let a = (i as f64) * 4.0 * std::f64::consts::PI / 5.0;
                Vec2::new(10.0 * math::cos(a), 10.0 * math::sin(a))
            })
            .collect();
        let even_odd = CrossSection::from_polygon_with_fill_rule(star.clone(), 0);
        let positive = CrossSection::from_polygon_with_fill_rule(star.clone(), 2);
        let unknown = CrossSection::from_polygon_with_fill_rule(star, 99);
        assert!(even_odd.area() < positive.area());
        assert_eq!(unknown.to_polygons(), even_odd.to_polygons());
        let sq = CrossSection::square(1.0);
        assert_eq!(
            sq.offset_with_params(0.5, 99, 2.0, 0).to_polygons(),
            sq.offset_with_params(0.5, 0, 2.0, 0).to_polygons()
        );
    }

    /// C++ `Offset` defaults to Round joins with `circularSegments = 0`, which
    /// derives the arc tolerance from `Quality::GetCircularSegments(delta)`.
    #[test]
    fn test_cross_section_offset_default_segments_match_quality() {
        let sq = CrossSection::square(1.0);
        let n = crate::types::Quality::get_circular_segments(3.0);
        let expected = sq.offset_with_params(3.0, 1, 2.0, n).to_polygons();
        assert_eq!(sq.offset(3.0).to_polygons(), expected);
        assert_eq!(
            sq.offset_with_params(3.0, 1, 2.0, 0).to_polygons(),
            expected
        );
    }
}
