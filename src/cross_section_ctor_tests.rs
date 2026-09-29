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

// cross_section_ctor_tests.rs — tests pinning CrossSection's constructors
// (cross_section.rs) and Manifold::slice / project, which wrap their
// polygons in one, to the C++ reference. Expected bit patterns come from the
// C++ reference compiled with MSVC against Clipper2 46f6391.

use super::*;
use crate::manifold::Manifold;

type Bits = &'static [&'static [(u64, u64)]];

fn bits(p: &Polygons) -> Vec<Vec<(u64, u64)>> {
    p.iter()
        .map(|c| c.iter().map(|v| (v.x.to_bits(), v.y.to_bits())).collect())
        .collect()
}

fn want(b: Bits) -> Vec<Vec<(u64, u64)>> {
    b.iter().map(|c| c.to_vec()).collect()
}

/// C++ `Circle` steps `360.0 / n` degrees through `cosd` / `sind`, which
/// land exactly on the axes (with the signed zeros below), and falls back
/// to `Quality::GetCircularSegments(radius)` when `circularSegments <= 2`.
#[test]
fn test_circle_matches_cpp_bits() {
    const CIRCLE_1_8: Bits = &[&[
        (0x3ff0000000000000, 0x0000000000000000),
        (0x3fe6a09e667f3bcc, 0x3fe6a09e667f3bcc),
        (0x8000000000000000, 0x3ff0000000000000),
        (0xbfe6a09e667f3bcc, 0x3fe6a09e667f3bcc),
        (0xbff0000000000000, 0x8000000000000000),
        (0xbfe6a09e667f3bcc, 0xbfe6a09e667f3bcc),
        (0x0000000000000000, 0xbff0000000000000),
        (0x3fe6a09e667f3bcc, 0xbfe6a09e667f3bcc),
    ]];
    const CIRCLE_2_0: Bits = &[&[
        (0x4000000000000000, 0x0000000000000000),
        (0x3ffbb67ae8584cab, 0x3fefffffffffffff),
        (0x3fefffffffffffff, 0x3ffbb67ae8584cab),
        (0x8000000000000000, 0x4000000000000000),
        (0xbfefffffffffffff, 0x3ffbb67ae8584cab),
        (0xbffbb67ae8584cab, 0x3fefffffffffffff),
        (0xc000000000000000, 0x8000000000000000),
        (0xbffbb67ae8584cab, 0xbfefffffffffffff),
        (0xbfefffffffffffff, 0xbffbb67ae8584cab),
        (0x0000000000000000, 0xc000000000000000),
        (0x3fefffffffffffff, 0xbffbb67ae8584cab),
        (0x3ffbb67ae8584cab, 0xbfefffffffffffff),
    ]];
    assert_eq!(
        bits(&CrossSection::circle(1.0, 8).to_polygons()),
        want(CIRCLE_1_8)
    );
    assert_eq!(
        bits(&CrossSection::circle(2.0, 0).to_polygons()),
        want(CIRCLE_2_0)
    );
    assert_eq!(
        bits(&CrossSection::circle(2.0, 2).to_polygons()),
        want(CIRCLE_2_0)
    );
    assert!(CrossSection::circle(0.0, 8).is_empty());
}

/// C++ `Square` returns empty only for a negative dimension or a zero-length
/// size vector, so a zero-height rectangle is one degenerate contour; the
/// centered square starts at (+w/2, +h/2) and runs counter-clockwise.
#[test]
fn test_square_matches_cpp() {
    let p =
        |c: &[(f64, f64)]| -> Polygons { vec![c.iter().map(|&(x, y)| Vec2::new(x, y)).collect()] };
    assert_eq!(
        CrossSection::square_vec2(Vec2::new(2.0, 4.0), true).to_polygons(),
        p(&[(1.0, 2.0), (-1.0, 2.0), (-1.0, -2.0), (1.0, -2.0)])
    );
    assert_eq!(
        CrossSection::square_vec2(Vec2::new(5.0, 0.0), false).to_polygons(),
        p(&[(0.0, 0.0), (5.0, 0.0), (5.0, 0.0), (0.0, 0.0)])
    );
    assert_eq!(
        CrossSection::square_vec2(Vec2::new(0.0, 0.0), false).num_vert(),
        0
    );
    assert_eq!(
        CrossSection::square_vec2(Vec2::new(-1.0, 1.0), false).num_vert(),
        0
    );
    assert_eq!(CrossSection::square(0.0).num_vert(), 0);
}

/// C++ `CrossSection(const Polygons&, FillRule = Positive)` always runs
/// `C2::Union`, so overlapping contours merge and coordinates snap to
/// Clipper2's grid at `precision_`.
#[test]
fn test_new_unions_like_cpp_polygons_ctor() {
    const TRI: Bits = &[&[
        (0x0000000000000000, 0x3ff0000000000000),
        (0x0000000000000000, 0x0000000000000000),
        (0x3ff0000020000000, 0x0000000000000000),
    ]];
    let tri = vec![vec![
        Vec2::new(0.0, 0.0),
        Vec2::new(1.000_000_12, 0.0),
        Vec2::new(0.0, 1.0),
    ]];
    assert_eq!(
        bits(&CrossSection::new(tri.clone()).to_polygons()),
        want(TRI)
    );
    assert_eq!(
        bits(&CrossSection::from_polygons_fill(tri).to_polygons()),
        want(TRI)
    );
    let sq = |x: f64, y: f64| {
        vec![
            Vec2::new(x, y),
            Vec2::new(x + 2.0, y),
            Vec2::new(x + 2.0, y + 2.0),
            Vec2::new(x, y + 2.0),
        ]
    };
    let merged: Polygons = vec![[
        (2.0, 1.0),
        (3.0, 1.0),
        (3.0, 3.0),
        (1.0, 3.0),
        (1.0, 2.0),
        (0.0, 2.0),
        (0.0, 0.0),
        (2.0, 0.0),
    ]
    .iter()
    .map(|&(x, y)| Vec2::new(x, y))
    .collect()];
    assert_eq!(
        CrossSection::new(vec![sq(0.0, 0.0), sq(1.0, 1.0)]).to_polygons(),
        merged
    );
}

/// C++ `CrossSection(const Rect&)` wraps the four corners with no union and
/// no emptiness check, so a default (inverted, infinite) Rect gives one
/// contour of infinities.
#[test]
fn test_from_rect_matches_cpp() {
    let r = Rect {
        min: Vec2::new(0.0, 0.0),
        max: Vec2::new(2.0, 1.0),
    };
    let p =
        |c: &[(f64, f64)]| -> Polygons { vec![c.iter().map(|&(x, y)| Vec2::new(x, y)).collect()] };
    assert_eq!(
        CrossSection::from_rect(&r).to_polygons(),
        p(&[(0.0, 0.0), (2.0, 0.0), (2.0, 1.0), (0.0, 1.0)])
    );
    let inf = f64::INFINITY;
    assert_eq!(
        CrossSection::from_rect(&Rect::new()).to_polygons(),
        p(&[(inf, inf), (-inf, inf), (-inf, -inf), (inf, -inf)])
    );
}

/// C++ `Slice` / `Project` return raw `Polygons`, and a CrossSection is
/// made from them only through the Positive-union Polygons constructor
/// (`CrossSection bottom = cube.Slice();` in manifold_test.cpp), so the
/// wrapped sections must equal C++ `CrossSection(m.Slice())` /
/// `CrossSection(m.Project())`. The raw projection already agrees bit for
/// bit; the raw slice does not (C++ starts from `*unordered_set::begin()`
/// and a few interpolated coordinates differ by one ULP), but the union's
/// snap to Clipper2's grid and its normalized start vertex erase both.
#[test]
fn test_slice_and_project_wrap_like_cpp() {
    const SLICE_CS: Bits = &[&[
        (0x3fda0e0998000000, 0xbfe6a09e68000000),
        (0x3fe6a09e68000000, 0xbfda0e0998000000),
        (0x3fec06075c000000, 0x0000000000000000),
        (0x3fe6a09e68000000, 0x3fda0e0998000000),
        (0x3fda0e0998000000, 0x3fe6a09e68000000),
        (0x0000000000000000, 0x3fec06075c000000),
        (0xbfda0e0998000000, 0x3fe6a09e68000000),
        (0xbfe6a09e68000000, 0x3fda0e0998000000),
        (0xbfec06075c000000, 0x0000000000000000),
        (0xbfe6a09e68000000, 0xbfda0e0998000000),
        (0xbfda0e0998000000, 0xbfe6a09e68000000),
        (0x0000000000000000, 0xbfec06075c000000),
    ]];
    const PROJ_RAW: Bits = &[&[
        (0xbfe6a09e667f3bcd, 0xbfe6a09e667f3bcd),
        (0x3c91a62633145c07, 0xbff0000000000000),
        (0x3fe6a09e667f3bcd, 0xbfe6a09e667f3bcc),
        (0x3ff0000000000000, 0x3c91a62633145c07),
        (0x3fe6a09e667f3bcd, 0x3fe6a09e667f3bcd),
        (0x3c91a62633145c07, 0x3ff0000000000000),
        (0xbfe6a09e667f3bcc, 0x3fe6a09e667f3bcd),
        (0xbff0000000000000, 0x3c91a62633145c07),
    ]];
    const PROJ_CS: Bits = &[&[
        (0x3fe6a09e68000000, 0xbfe6a09e68000000),
        (0x3ff0000000000000, 0x0000000000000000),
        (0x3fe6a09e68000000, 0x3fe6a09e68000000),
        (0x0000000000000000, 0x3ff0000000000000),
        (0xbfe6a09e68000000, 0x3fe6a09e68000000),
        (0xbff0000000000000, 0x0000000000000000),
        (0xbfe6a09e68000000, 0xbfe6a09e68000000),
        (0x0000000000000000, 0xbff0000000000000),
    ]];
    let s = Manifold::sphere(1.0, 8);
    assert_eq!(bits(&s.as_impl().project()), want(PROJ_RAW));
    assert_eq!(bits(&s.slice(0.3).to_polygons()), want(SLICE_CS));
    assert_eq!(bits(&s.project().to_polygons()), want(PROJ_CS));
}
