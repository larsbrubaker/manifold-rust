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

// cross_section_tests.rs — unit tests for CrossSection (cross_section.rs and
// its child module cross_section_ops.rs). Expected values marked as
// C++-derived come from the C++ reference and Clipper2 at the pinned commit.

use super::*;
use crate::types::OpType;
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
