use super::*;

/// C++ TEST(CrossSection, BevelOffset) — bevel join offset of a square
#[test]
fn test_cpp_cross_section_bevel_offset() {
    let a = CrossSection::square(20.0).translate(Vec2::new(-10.0, -10.0));
    let segments = 20;
    // join_type=3 is Bevel
    let rounded = a.offset_with_params(5.0, 3, 2.0, segments);
    let result = Manifold::extrude(&rounded.to_polygons(), 5.0, 0, 0.0, Vec2::new(1.0, 1.0));

    assert_eq!(result.genus(), 0, "BevelOffset genus={}", result.genus());
    // Volume = height * (outer_area - corner_cuts); with bevel, corners cut off triangles of size 5*5
    let expected_vol = 5.0 * ((20.0 + 2.0 * 5.0) * (20.0 + 2.0 * 5.0) - 2.0 * 5.0 * 5.0);
    assert!((result.volume() - expected_vol).abs() < 1.0,
        "BevelOffset vol={} expected~{}", result.volume(), expected_vol);
    assert_eq!(rounded.num_vert(), 4 + 4, "BevelOffset NumVert={}", rounded.num_vert());
}

/// C++ TEST(CrossSection, Warp) — warp function applied to square vertices
#[test]
fn test_cpp_cross_section_warp() {
    let sq = CrossSection::square(10.0);
    let _a = sq.scale(Vec2::new(2.0, 3.0)).translate(Vec2::new(4.0, 5.0));
    let _b = sq.warp(|v| {
        v.x = v.x * 2.0 + 4.0;
        v.y = v.y * 3.0 + 5.0;
    });
    assert_eq!(sq.num_vert(), 4, "sq NumVert={}", sq.num_vert());
    assert_eq!(sq.num_contour(), 1, "sq NumContour={}", sq.num_contour());
}

/// C++ TEST(CrossSection, FillRule) — fill rule affects area of self-intersecting polygon
#[test]
fn test_cpp_cross_section_fill_rule() {
    let polygon = vec![
        Vec2::new(-7.0, 13.0),
        Vec2::new(-7.0, 12.0),
        Vec2::new(-5.0, 9.0),
        Vec2::new(-5.0, 8.1),
        Vec2::new(-4.8, 8.0),
    ];

    // fill_rule: 0=EvenOdd, 1=NonZero, 2=Positive, 3=Negative
    let positive = CrossSection::from_polygon_with_fill_rule(polygon.clone(), 2);
    assert!((positive.area() - 0.683).abs() < 0.001,
        "Positive area={}", positive.area());

    let negative = CrossSection::from_polygon_with_fill_rule(polygon.clone(), 3);
    assert!((negative.area() - 0.193).abs() < 0.001,
        "Negative area={}", negative.area());

    let even_odd = CrossSection::from_polygon_with_fill_rule(polygon.clone(), 0);
    assert!((even_odd.area() - 0.875).abs() < 0.001,
        "EvenOdd area={}", even_odd.area());

    let non_zero = CrossSection::from_polygon_with_fill_rule(polygon.clone(), 1);
    assert!((non_zero.area() - 0.875).abs() < 0.001,
        "NonZero area={}", non_zero.area());
}

/// C++ TEST(CrossSection, HullError) — rounded rectangle via hull of circles
#[test]
fn test_cpp_cross_section_hull_error() {
    let rounded_rectangle = |x: f64, y: f64, radius: f64, segments: i32| {
        let circ = CrossSection::circle(radius, segments);
        let vl = vec![
            circ.translate(Vec2::new(radius, radius)),
            circ.translate(Vec2::new(x - radius, radius)),
            circ.translate(Vec2::new(x - radius, y - radius)),
            circ.translate(Vec2::new(radius, y - radius)),
        ];
        CrossSection::hull_cross_sections(&vl)
    };
    let rr = rounded_rectangle(51.0, 36.0, 9.0, 36);
    assert!((rr.area() - 1765.179_037_5).abs() < 1.0,
        "HullError area={}", rr.area());
    assert_eq!(rr.num_vert(), 40, "HullError NumVert={}", rr.num_vert());
}

/// C++ TEST(CrossSection, BatchBoolean) — batch boolean ops on cross sections
#[test]
fn test_cpp_cross_section_batch_boolean() {
    let square = CrossSection::square(100.0);
    let circle1 = CrossSection::circle(30.0, 30).translate(Vec2::new(-10.0, 30.0));
    let circle2 = CrossSection::circle(20.0, 30).translate(Vec2::new(110.0, 20.0));
    let circle3 = CrossSection::circle(40.0, 30).translate(Vec2::new(50.0, 110.0));

    let sections = vec![square.clone(), circle1.clone(), circle2.clone(), circle3.clone()];

    let intersect = CrossSection::batch_boolean(&sections, OpType::Intersect);
    assert!((intersect.area()).abs() < 1e-4,
        "BatchBoolean intersect area={}", intersect.area());
    assert_eq!(intersect.num_vert(), 0,
        "BatchBoolean intersect NumVert={}", intersect.num_vert());

    let add = CrossSection::batch_boolean(&sections, OpType::Add);
    assert!((add.area() - 16278.637_002).abs() < 1.0,
        "BatchBoolean add area={}", add.area());
    assert_eq!(add.num_vert(), 66,
        "BatchBoolean add NumVert={}", add.num_vert());

    let subtract = CrossSection::batch_boolean(&sections, OpType::Subtract);
    assert!((subtract.area() - 7234.478_452).abs() < 1.0,
        "BatchBoolean subtract area={}", subtract.area());
    assert_eq!(subtract.num_vert(), 42,
        "BatchBoolean subtract NumVert={}", subtract.num_vert());
}
