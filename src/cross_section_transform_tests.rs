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

// cross_section_transform_tests.rs — tests pinning CrossSection's lazy affine
// transforms (translate / rotate / scale / mirror in cross_section.rs) to the
// C++ reference: the composed mat2x3, its application as `m * vec3(x, y, 1)`,
// the identity shortcut, winding reversal on a negative determinant, and the
// point at which a read materializes the transform. Expected bit patterns
// come from the C++ reference compiled with MSVC against Clipper2 46f6391.

use super::*;

type Bits = &'static [&'static [(u64, u64)]];

fn bits(p: &Polygons) -> Vec<Vec<(u64, u64)>> {
    p.iter()
        .map(|c| c.iter().map(|v| (v.x.to_bits(), v.y.to_bits())).collect())
        .collect()
}

fn want(b: Bits) -> Vec<Vec<(u64, u64)>> {
    b.iter().map(|c| c.to_vec()).collect()
}

const CHAIN: Bits = &[&[
    (0x3ff120ab8f1b81cd, 0xbfe77420ce509266),
    (0x3fdb35858b1ea3f2, 0x3fc10e6311d12e36),
    (0xbfdf32889fab9bc0, 0x3ff07e32c3f80733),
    (0xbff21fb18c247aeb, 0x3ff6f18ced3e5df7),
    (0xbff21f87fb86eff7, 0x3ff1b483e64af552),
    (0xbfdf30f73ccc5c9a, 0x3fcec538e74432c6),
    (0x3fdb3716edfde318, 0xbfe5077e89aab626),
    (0x3ff120d51fb90cc1, 0xbff0f7196e1bb1d7),
]];

const CHAIN_MAT: Bits = &[&[
    (0x3ff120ab8f1b81cd, 0xbfe77420ce509267),
    (0x3fdb35858b1ea3f3, 0x3fc10e6311d12e36),
    (0xbfdf32889fab9bc1, 0x3ff07e32c3f80732),
    (0xbff21fb18c247aec, 0x3ff6f18ced3e5df7),
    (0xbff21f87fb86eff7, 0x3ff1b483e64af552),
    (0xbfdf30f73ccc5c9b, 0x3fcec538e74432c2),
    (0x3fdb3716edfde319, 0xbfe5077e89aab627),
    (0x3ff120d51fb90cc1, 0xbff0f7196e1bb1d7),
]];

const THERE_BACK: Bits = &[&[
    (0x3ff0000000000000, 0x0000000000000000),
    (0x3fe6a09e667f3bcc, 0x3fe6a09e667f3bcc),
    (0x8000000000000000, 0x3ff0000000000000),
    (0xbfe6a09e667f3bcc, 0x3fe6a09e667f3bcc),
    (0xbff0000000000000, 0x8000000000000000),
    (0xbfe6a09e667f3bcc, 0xbfe6a09e667f3bcc),
    (0x0000000000000000, 0xbff0000000000000),
    (0x3fe6a09e667f3bcc, 0xbfe6a09e667f3bcc),
]];

const ROT90: Bits = &[&[
    (0x0000000000000000, 0x0000000000000000),
    (0x0000000000000000, 0x4000000000000000),
    (0xbff0000000000000, 0x4000000000000000),
    (0xbff0000000000000, 0x0000000000000000),
]];

const ROT180: Bits = &[&[
    (0x0000000000000000, 0x0000000000000000),
    (0xc000000000000000, 0x0000000000000000),
    (0xc000000000000000, 0xbff0000000000000),
    (0x0000000000000000, 0xbff0000000000000),
]];

const ROTM90: Bits = &[&[
    (0x0000000000000000, 0x0000000000000000),
    (0x0000000000000000, 0xc000000000000000),
    (0x3ff0000000000000, 0xc000000000000000),
    (0x3ff0000000000000, 0x0000000000000000),
]];

const ROT45: Bits = &[&[
    (0x0000000000000000, 0x0000000000000000),
    (0x3ff6a09e667f3bcc, 0x3ff6a09e667f3bcc),
    (0x3fe6a09e667f3bcc, 0x4000f876ccdf6cd9),
    (0xbfe6a09e667f3bcc, 0x3fe6a09e667f3bcc),
]];

const ROT30_C8: Bits = &[&[
    (0x3febb67ae8584cab, 0x3fdfffffffffffff),
    (0x3fd0907dc1930691, 0x3feee8dd4748bf14),
    (0xbfdfffffffffffff, 0x3febb67ae8584cab),
    (0xbfeee8dd4748bf14, 0x3fd0907dc1930691),
    (0xbfebb67ae8584cab, 0xbfdfffffffffffff),
    (0xbfd0907dc1930691, 0xbfeee8dd4748bf14),
    (0x3fdfffffffffffff, 0xbfebb67ae8584cab),
    (0x3feee8dd4748bf14, 0xbfd0907dc1930691),
]];

const MIRROR_1E_11: Bits = &[&[
    (0x0000000000000000, 0x3ff0000000000000),
    (0xc000000000000000, 0x3ff0000000000000),
    (0xc000000000000000, 0x0000000000000000),
    (0x0000000000000000, 0x0000000000000000),
]];

const MIRROR_1E_160: Bits = &[&[
    (0x0000000000000000, 0x3ff0000000000000),
    (0xc000001758f3cba8, 0x3ff0000000000000),
    (0xc000001758f3cba8, 0x0000000000000000),
    (0x0000000000000000, 0x0000000000000000),
]];

const MIRROR_Y: Bits = &[&[
    (0x0000000000000000, 0xbff0000000000000),
    (0x4000000000000000, 0xbff0000000000000),
    (0x4000000000000000, 0x0000000000000000),
    (0x0000000000000000, 0x0000000000000000),
]];

const SCALE_NEG: Bits = &[&[
    (0x0000000000000000, 0x3ff0000000000000),
    (0xc000000000000000, 0x3ff0000000000000),
    (0xc000000000000000, 0x0000000000000000),
    (0x0000000000000000, 0x0000000000000000),
]];

const SCALE_ZERO: Bits = &[&[
    (0x0000000000000000, 0x0000000000000000),
    (0x0000000000000000, 0x0000000000000000),
    (0x0000000000000000, 0x3ff0000000000000),
    (0x0000000000000000, 0x3ff0000000000000),
]];

const UNION_TRANSFORMED: Bits = &[&[
    (0x3fe16daed8000000, 0xbfead663a8000000),
    (0x3fef4cfc34000000, 0xbfca9cd9b0000000),
    (0x3fec95bd30000000, 0x3fd0000000000000),
    (0x4004000000000000, 0x3fd0000000000000),
    (0x4004000000000000, 0x3ff4000000000000),
    (0x3fe0000000000000, 0x3ff4000000000000),
    (0x3fe0000000000000, 0x3fe8e077c8000000),
    (0x3fca9cd9b0000000, 0x3fef4cfc34000000),
    (0xbfe16daed8000000, 0x3fead663a8000000),
    (0xbfef4cfc34000000, 0x3fca9cd9b0000000),
    (0xbfead663a8000000, 0xbfe16daed8000000),
    (0xbfca9cd9b0000000, 0xbfef4cfc34000000),
]];

/// C++ composes every transform into one mat2x3 (`m * Mat3(transform_)`) and
/// applies it once on read, so a chain rounds differently from applying each
/// step to the vertices.
#[test]
fn test_chained_transforms_compose_lazily_like_cpp() {
    let cs = CrossSection::circle(1.0, 8)
        .translate(Vec2::new(0.1, 0.2))
        .rotate(33.0)
        .scale(Vec2::new(1.7, -0.3))
        .mirror(Vec2::new(0.3, 0.7))
        .translate(Vec2::new(-0.05, 0.11));
    assert_eq!(bits(&cs.to_polygons()), want(CHAIN));
    assert_eq!(cs.area().to_bits(), 0x3ff714789bbf37de);
}

/// A read (C++ `GetPaths`) bakes the pending transform into the paths and
/// resets it to identity, so later transforms compose from the baked paths.
#[test]
fn test_read_materializes_transform_like_cpp() {
    let x = CrossSection::circle(1.0, 8)
        .translate(Vec2::new(0.1, 0.2))
        .rotate(33.0);
    x.area();
    let cs = x
        .scale(Vec2::new(1.7, -0.3))
        .mirror(Vec2::new(0.3, 0.7))
        .translate(Vec2::new(-0.05, 0.11));
    assert_eq!(bits(&cs.to_polygons()), want(CHAIN_MAT));
    assert_eq!(cs.area().to_bits(), 0x3ff714789bbf37df);
}

/// Translating there and back composes to exactly the identity, which C++
/// `GetPaths` skips, so the vertices (including the circle's -0.0s) are
/// returned untouched.
#[test]
fn test_identity_composite_leaves_paths_untouched() {
    let cs = CrossSection::circle(1.0, 8)
        .translate(Vec2::new(0.1, 0.3))
        .translate(Vec2::new(-0.1, -0.3));
    assert_eq!(bits(&cs.to_polygons()), want(THERE_BACK));
}

/// C++ `Rotate` takes `sind` / `cosd`, which are exact at multiples of 90
/// degrees.
#[test]
fn test_rotate_uses_sind_cosd_like_cpp() {
    let r = CrossSection::square_vec2(Vec2::new(2.0, 1.0), false);
    assert_eq!(bits(&r.rotate(90.0).to_polygons()), want(ROT90));
    assert_eq!(bits(&r.rotate(180.0).to_polygons()), want(ROT180));
    assert_eq!(bits(&r.rotate(-90.0).to_polygons()), want(ROTM90));
    assert_eq!(bits(&r.rotate(45.0).to_polygons()), want(ROT45));
    assert_eq!(r.rotate(45.0).area().to_bits(), 0x3ffffffffffffffe);
    let c = CrossSection::circle(1.0, 8).rotate(30.0);
    assert_eq!(bits(&c.to_polygons()), want(ROT30_C8));
}

/// C++ `Mirror` returns an empty section only when `la::length(ax) == 0`
/// (which includes lengths that underflow); any other axis is normalized and
/// reflected, with the winding reversed by the negative determinant.
#[test]
fn test_mirror_guard_and_matrix_like_cpp() {
    let r = CrossSection::square_vec2(Vec2::new(2.0, 1.0), false);
    assert_eq!(r.mirror(Vec2::new(0.0, 0.0)).num_vert(), 0);
    assert_eq!(
        bits(&r.mirror(Vec2::new(1e-11, 0.0)).to_polygons()),
        want(MIRROR_1E_11)
    );
    assert_eq!(
        bits(&r.mirror(Vec2::new(1e-160, 0.0)).to_polygons()),
        want(MIRROR_1E_160)
    );
    assert_eq!(r.mirror(Vec2::new(1e-170, 0.0)).num_vert(), 0);
    assert_eq!(
        bits(&r.mirror(Vec2::new(0.0, 1.0)).to_polygons()),
        want(MIRROR_Y)
    );
}

/// Scale goes through the same matrix: a negative determinant reverses each
/// contour, a zero one does not.
#[test]
fn test_scale_matrix_and_winding_like_cpp() {
    let r = CrossSection::square_vec2(Vec2::new(2.0, 1.0), false);
    assert_eq!(
        bits(&r.scale(Vec2::new(-1.0, 1.0)).to_polygons()),
        want(SCALE_NEG)
    );
    assert_eq!(
        bits(&r.scale(Vec2::new(0.0, 1.0)).to_polygons()),
        want(SCALE_ZERO)
    );
}

/// `m * vec3(x, y, 1)` multiplies every column, so translating the infinite
/// corners of `CrossSection::from_rect(&Rect::new())` gives `0 * inf` = NaN in
/// every coordinate, as in C++. (NaN sign and payload are not specified by
/// Rust, so only NaN-ness is compared.)
#[test]
fn test_translate_infinite_rect_gives_nan_like_cpp() {
    let cs = CrossSection::from_rect(&Rect::new()).translate(Vec2::new(1.0, 2.0));
    let p = cs.to_polygons();
    assert_eq!(p.len(), 1);
    assert_eq!(p[0].len(), 4);
    assert!(p[0].iter().all(|v| v.x.is_nan() && v.y.is_nan()));
    assert!(cs.area().is_nan());
}

/// Bounds and booleans read through the pending transform too.
#[test]
fn test_readers_apply_pending_transform_like_cpp() {
    let b = CrossSection::circle(1.0, 8)
        .translate(Vec2::new(0.1, 0.2))
        .rotate(33.0)
        .bounds();
    assert_eq!(
        (
            b.min.x.to_bits(),
            b.min.y.to_bits(),
            b.max.x.to_bits(),
            b.max.y.to_bits()
        ),
        (
            0xbff00d243325f02b,
            0xbfe830bd2e65849d,
            0x3fee7faffea820a9,
            0x3ff3349d9b473e31
        )
    );
    let r = CrossSection::square_vec2(Vec2::new(2.0, 1.0), false);
    let u = CrossSection::circle(1.0, 8)
        .rotate(33.0)
        .union(&r.translate(Vec2::new(0.5, 0.25)));
    assert_eq!(bits(&u.to_polygons()), want(UNION_TRANSFORMED));
    assert_eq!(u.area().to_bits(), 0x4012b987c182097b);
}
