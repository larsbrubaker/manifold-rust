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
    assert_eq!(
        CrossSection::from_raw(even).area().to_bits(),
        0x4008f42c81cc8074
    );
}

/// C++ runs every Clipper2 op at `precision_ = 8` decimal places. ClipperD
/// scales by the power of two above 10^precision (2^27 at 8, 2^20 at 6),
/// so x = 1.00000012 snaps to 1 + 16 * 2^-27 = 1 + 2^-23, where precision
/// 6 would round the 1.2e-7 feature away to 1.0.
#[test]
fn test_cross_section_union_keeps_eighth_decimal() {
    let x = 1.000_000_12;
    let snapped = 1.0 + 2f64.powi(-23);
    let a = CrossSection::from_raw(vec![vec![
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
    let cw = CrossSection::from_raw(vec![vec![
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
            Vec2::new(10.0 * crate::math::cos(a), 10.0 * crate::math::sin(a))
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

/// Build Polygons from coordinate-pair literals.
fn polys(contours: &[&[(f64, f64)]]) -> Polygons {
    contours
        .iter()
        .map(|c| c.iter().map(|&(x, y)| Vec2::new(x, y)).collect())
        .collect()
}

/// The bar [0,10]x[0,2] with hole [8,9]x[0.5,1.5], unioned with a U whose
/// bbox [7,11]x[-0.5,2.5] covers the hole's vertices while its opening
/// embraces the bar's right end.
fn bar_and_u() -> CrossSection {
    let bar = CrossSection::square_vec2(Vec2::new(10.0, 2.0), false).difference(
        &CrossSection::square_vec2(Vec2::new(1.0, 1.0), false).translate(Vec2::new(8.0, 0.5)),
    );
    let u = CrossSection::square_vec2(Vec2::new(4.0, 3.0), false)
        .translate(Vec2::new(7.0, -0.5))
        .difference(
            &CrossSection::square_vec2(Vec2::new(3.5, 2.5), false).translate(Vec2::new(7.0, -0.25)),
        );
    bar.union(&u)
}

const U_OUTLINE: &[(f64, f64)] = &[
    (11.0, 2.5),
    (7.0, 2.5),
    (7.0, 2.25),
    (10.5, 2.25),
    (10.5, -0.25),
    (7.0, -0.25),
    (7.0, -0.5),
    (11.0, -0.5),
];

/// C++ `Decompose` groups holes by Clipper2's PolyTree containment, so the
/// bar keeps its hole even though the U's bounding box also covers it.
/// Expected contours and order from the C++ reference compiled against
/// Clipper2 46f6391.
#[test]
fn test_decompose_keeps_hole_with_its_outline() {
    let cs = bar_and_u();
    let bar = &[(10.0, 2.0), (0.0, 2.0), (0.0, 0.0), (10.0, 0.0)][..];
    let hole = &[(8.0, 1.5), (9.0, 1.5), (9.0, 0.5), (8.0, 0.5)][..];
    assert_eq!(cs.to_polygons(), polys(&[U_OUTLINE, bar, hole]));
    let comps: Vec<Polygons> = cs.decompose().iter().map(|c| c.to_polygons()).collect();
    assert_eq!(comps, vec![polys(&[bar, hole]), polys(&[U_OUTLINE])]);
}

/// C++ emits the reversed stack of its outline/hole recursion: an island
/// inside a hole is pushed before its enclosing outline, later siblings
/// after. Expected order from the compiled C++ reference.
#[test]
fn test_decompose_order_matches_cpp() {
    let ring = |outer: f64, inner: f64| {
        CrossSection::square_vec2(Vec2::new(outer, outer), true)
            .difference(&CrossSection::square_vec2(Vec2::new(inner, inner), true))
    };
    let nest = ring(10.0, 8.0)
        .union(&ring(4.0, 2.0))
        .union(&CrossSection::square(1.0).translate(Vec2::new(20.0, 0.0)));
    let comps: Vec<Polygons> = nest.decompose().iter().map(|c| c.to_polygons()).collect();
    let sq = |h: f64| [(h, h), (-h, h), (-h, -h), (h, -h)];
    let hole = |h: f64| [(-h, h), (h, h), (h, -h), (-h, -h)];
    assert_eq!(
        comps,
        vec![
            polys(&[&[(21.0, 1.0), (20.0, 1.0), (20.0, 0.0), (21.0, 0.0)]]),
            polys(&[&sq(5.0), &hole(4.0)]),
            polys(&[&sq(2.0), &hole(1.0)]),
        ]
    );
}

/// C++ returns `*this` unchanged when `NumContour() < 2`: an empty section
/// decomposes to one empty section, and a single contour is not pushed
/// through Clipper2 (which would snap it to the 2^-27 grid).
#[test]
fn test_decompose_short_circuits_below_two_contours() {
    let empty = CrossSection::default().decompose();
    assert_eq!(empty.len(), 1);
    assert!(empty[0].is_empty());
    let circ = CrossSection::circle(1.0, 8).translate(Vec2::new(0.1, 0.2));
    let comps = circ.decompose();
    assert_eq!(comps.len(), 1);
    let bits = |p: &Polygons| -> Vec<(u64, u64)> {
        p.iter()
            .flatten()
            .map(|v| (v.x.to_bits(), v.y.to_bits()))
            .collect()
    };
    assert_eq!(bits(&comps[0].to_polygons()), bits(&circ.to_polygons()));
}

/// C++ `Simplify` unions into a PolyTree and `flatten`s it, pushing each
/// node's descendants before the node itself, so holes precede their
/// outline. Expected contours from the compiled C++ reference.
#[test]
fn test_simplify_flattens_polytree_like_cpp() {
    let ring = CrossSection::square(10.0)
        .difference(&CrossSection::square(4.0).translate(Vec2::new(3.0, 3.0)));
    assert_eq!(
        ring.simplify(1e-6).to_polygons(),
        polys(&[
            &[(3.0, 7.0), (7.0, 7.0), (7.0, 3.0), (3.0, 3.0)],
            &[(10.0, 10.0), (0.0, 10.0), (0.0, 0.0), (10.0, 0.0)],
        ])
    );
    let ring = |outer: f64, inner: f64| {
        CrossSection::square_vec2(Vec2::new(outer, outer), true)
            .difference(&CrossSection::square_vec2(Vec2::new(inner, inner), true))
    };
    let nest = ring(10.0, 8.0)
        .union(&ring(4.0, 2.0))
        .union(&CrossSection::square(1.0).translate(Vec2::new(20.0, 0.0)));
    let sq = |h: f64| [(h, h), (-h, h), (-h, -h), (h, -h)];
    let hole = |h: f64| [(-h, h), (h, h), (h, -h), (-h, -h)];
    assert_eq!(
        nest.simplify(1e-6).to_polygons(),
        polys(&[
            &hole(1.0),
            &sq(2.0),
            &hole(4.0),
            &sq(5.0),
            &[(21.0, 1.0), (20.0, 1.0), (20.0, 0.0), (21.0, 0.0)],
        ])
    );
}

fn three_squares() -> Vec<CrossSection> {
    vec![
        CrossSection::square(2.0),
        CrossSection::square(2.0).translate(Vec2::new(1.0, 1.0)),
        CrossSection::square(2.0).translate(Vec2::new(-1.0, 1.5)),
    ]
}

/// C++ `BatchBoolean` Add/Subtract run one `BooleanOp` with the first
/// section as subject and the rest as clips; Intersect folds pairwise, and
/// `Compose` is BatchBoolean Add. Expected contours from the compiled C++.
#[test]
fn test_batch_boolean_matches_cpp() {
    let secs = three_squares();
    let add = polys(&[&[
        (2.0, 1.0),
        (3.0, 1.0),
        (3.0, 3.0),
        (1.0, 3.0),
        (1.0, 3.5),
        (-1.0, 3.5),
        (-1.0, 1.5),
        (0.0, 1.5),
        (0.0, 0.0),
        (2.0, 0.0),
    ]]);
    assert_eq!(
        CrossSection::batch_boolean(&secs, OpType::Add).to_polygons(),
        add
    );
    assert_eq!(CrossSection::compose(&secs).to_polygons(), add);
    assert_eq!(
        CrossSection::batch_boolean(&secs, OpType::Subtract).to_polygons(),
        polys(&[&[
            (2.0, 1.0),
            (1.0, 1.0),
            (1.0, 1.5),
            (0.0, 1.5),
            (0.0, 0.0),
            (2.0, 0.0)
        ]])
    );
    assert!(CrossSection::batch_boolean(&secs, OpType::Intersect).is_empty());
}

/// C++ `BatchBoolean` returns `crossSections[0]` itself for a single input,
/// so neither it nor `Compose` snaps the contours through Clipper2.
#[test]
fn test_batch_boolean_single_section_is_unchanged() {
    let circ = CrossSection::circle(1.0, 8).translate(Vec2::new(0.1, 0.2));
    let one = [circ.clone()];
    let bits = |p: &Polygons| -> Vec<(u64, u64)> {
        p.iter()
            .flatten()
            .map(|v| (v.x.to_bits(), v.y.to_bits()))
            .collect()
    };
    let want = bits(&circ.to_polygons());
    for op in [OpType::Add, OpType::Subtract, OpType::Intersect] {
        assert_eq!(
            bits(&CrossSection::batch_boolean(&one, op).to_polygons()),
            want
        );
    }
    assert_eq!(bits(&CrossSection::compose(&one).to_polygons()), want);
    assert!(CrossSection::compose(&[]).is_empty());
}

/// Subtract runs one `BooleanOp` with every tail contour as a clip; a
/// pairwise fold reaches the same region with its contours in another
/// order. Clip triangles as C++ `Hull` emits them; expected contours from
/// the compiled C++ reference (its pairwise fold gives c1/c2 swapped).
#[test]
fn test_batch_subtract_is_one_boolean_op() {
    let tri = |p: [(f64, f64); 3]| CrossSection::from_raw(polys(&[&p]));
    let secs = [
        CrossSection::square(8.0).translate(Vec2::new(1.0, 1.0)),
        tri([
            (5.505859375, 9.8291015625),
            (6.05078125, 0.2421875),
            (9.4619140625, 1.42578125),
        ]),
        tri([
            (2.4287109375, 5.0869140625),
            (5.408203125, 0.8125),
            (4.029296875, 9.94140625),
        ]),
    ];
    assert_eq!(
        CrossSection::batch_boolean(&secs, OpType::Subtract).to_polygons(),
        polys(&[
            &[
                (2.4287109375, 5.0869140625),
                (3.718903623521328, 9.0),
                (1.0, 9.0),
                (1.0, 1.0),
                (5.2775057330727577, 1.0),
            ],
            &[
                (5.5529856532812119, 9.0),
                (4.1714947372674942, 9.0),
                (5.3798815608024597, 1.0),
                (6.0077070519328117, 1.0),
            ],
            &[
                (9.0, 9.0),
                (5.8961778432130814, 9.0),
                (9.0, 2.4069637954235077)
            ],
            &[
                (9.0, 1.2655064538121223),
                (8.2348068803548813, 1.0),
                (9.0, 1.0)
            ],
        ])
    );
}

/// C++ `Warp` goes through `WarpBatch`, which re-unions the moved contours
/// with FillRule::Positive at `precision_`: a warp that twists a square into
/// a bowtie keeps only the positively wound lobe (expected contour from the
/// compiled C++ reference), and moved vertices land on Clipper2's grid.
#[test]
fn test_warp_unions_like_cpp() {
    let bowtie = CrossSection::square(2.0).warp(|v| {
        if v.y > 1.0 {
            v.x = 2.0 - v.x;
        }
    });
    assert_eq!(
        bowtie.to_polygons(),
        polys(&[&[(1.0, 1.0), (0.0, 0.0), (2.0, 0.0)]])
    );
    assert_eq!(bowtie.area(), 1.0);
    let stretched = CrossSection::square(1.0).warp(|v| v.x *= 1.000_000_12);
    assert_eq!(stretched.bounds().max.x, 1.0 + 2f64.powi(-23));
}

/// C++ `IsEmpty` is `paths_.empty()` and `NumContour` is `paths_.size()`:
/// contours with fewer than three vertices still count. C++ `Hull` produces
/// exactly such sections — one empty contour for fewer than three points
/// (`h_2pts`: contours=1 nvert=0 empty=0) and a two-vertex contour for
/// collinear points (`h_collinear`: contours=1 nvert=2 empty=0) — through the
/// private no-union constructor that `from_raw` mirrors.
#[test]
fn test_is_empty_and_num_contour_count_every_path_like_cpp() {
    let one_empty = CrossSection::from_raw(vec![vec![]]);
    assert!(!one_empty.is_empty());
    assert_eq!(one_empty.num_contour(), 1);
    assert_eq!(one_empty.num_vert(), 0);
    let degenerate = CrossSection::from_raw(vec![
        vec![],
        vec![Vec2::new(0.0, 0.0), Vec2::new(3.0, 0.0)],
    ]);
    assert!(!degenerate.is_empty());
    assert_eq!(degenerate.num_contour(), 2);
    assert_eq!(degenerate.num_vert(), 2);
    let none = CrossSection::default();
    assert!(none.is_empty());
    assert_eq!(none.num_contour(), 0);
    assert_eq!(none.num_vert(), 0);
}
