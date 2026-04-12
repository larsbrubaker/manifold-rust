use super::*;

/// C++ TEST(SDF, Resize) — Layers SDF produces correct bounds and genus
#[test]
fn test_cpp_sdf_resize() {
    let size = 20.0;
    let layers = Manifold::level_set(
        |p: Vec3| {
            let a = ((2.0 * p.z).round() % 4.0) as i32;
            if a == 0 { 1.0 } else if a == 2 { -1.0 } else { 0.0 }
        },
        crate::types::Box {
            min: Vec3::new(0.0, 0.0, 0.0),
            max: Vec3::new(size, size, size),
        },
        1.0,
    );

    assert_eq!(layers.status(), crate::types::Error::NoError);
    assert_eq!(layers.genus(), -8, "Resize: genus should be -8, got {}", layers.genus());
    let epsilon = layers.get_tolerance();
    let bounds = layers.bounding_box();
    assert!((bounds.min.x - 0.0).abs() < epsilon, "min.x");
    assert!((bounds.min.y - 0.0).abs() < epsilon, "min.y");
    assert!((bounds.min.z - 1.5).abs() < epsilon, "min.z={}", bounds.min.z);
    assert!((bounds.max.x - size).abs() < epsilon, "max.x");
    assert!((bounds.max.y - size).abs() < epsilon, "max.y");
    assert!((bounds.max.z - (size - 1.5)).abs() < epsilon, "max.z={}", bounds.max.z);
}

/// C++ TEST(SDF, SineSurface) — SDF sine surface with smooth + refine
#[test]
#[ignore = "Requires SmoothOut implementation"]
fn test_cpp_sdf_sine_surface() {
    let pi = std::f64::consts::PI;
    let surface = Manifold::level_set(
        move |p: Vec3| {
            let mid = p.x.sin() + p.y.sin();
            if p.z > mid - 0.5 && p.z < mid + 0.5 { 1.0 } else { -1.0 }
        },
        crate::types::Box {
            min: Vec3::new(-1.75 * pi, -1.75 * pi, -1.75 * pi),
            max: Vec3::new(1.75 * pi, 1.75 * pi, 1.75 * pi),
        },
        1.0,
    ).simplify(0.0);
    // smoothed = surface.SmoothOut(180).RefineToLength(0.05)
    // EXPECT_EQ(smoothed.Genus(), 38);
    // EXPECT_NEAR(smoothed.Volume(), 107.4, 0.1);
    assert_eq!(surface.status(), crate::types::Error::NoError);
}

/// C++ TEST(SDF, Blobs) — metaball SDF using smoothstep
#[test]
#[ignore = "Slow: 263s in debug mode due to fine edge_length=0.05"]
fn test_cpp_sdf_blobs() {
    let blend = 1.0f64;
    let balls: Vec<[f64; 4]> = vec![
        [0.0, 0.0, 0.0, 2.0],
        [1.0, 2.0, 3.0, 2.0],
        [-2.0, 2.0, -2.0, 1.0],
        [-2.0, -3.0, -2.0, 2.0],
        [-3.0, -1.0, -3.0, 1.0],
        [2.0, -3.0, -2.0, 2.0],
        [-2.0, 3.0, 2.0, 2.0],
        [-2.0, -3.0, 2.0, 2.0],
        [1.0, -1.0, 1.0, -2.0],
        [-4.0, -3.0, -2.0, 1.0],
    ];
    let blobs = Manifold::level_set_with_level(
        move |p: Vec3| {
            let mut d = 0.0;
            for ball in &balls {
                let center = Vec3::new(ball[0], ball[1], ball[2]);
                let w = ball[3];
                let sign = if w > 0.0 { 1.0 } else { -1.0 };
                let diff = p - center;
                let dist = (diff.x * diff.x + diff.y * diff.y + diff.z * diff.z).sqrt();
                d += sign * crate::types::smoothstep(-blend, blend, w.abs() - dist);
            }
            d
        },
        crate::types::Box {
            min: Vec3::new(-5.0, -5.0, -5.0),
            max: Vec3::new(5.0, 5.0, 5.0),
        },
        0.05,
        0.5,
    );

    assert_eq!(blobs.status(), crate::types::Error::NoError);
    assert!(!blobs.is_empty(), "Blobs should not be empty");
    // C++ computes genus = 1 - chi/2 where chi = NumVert - NumTri/2
    let chi = blobs.num_vert() as i32 - blobs.num_tri() as i32 / 2;
    let genus = 1 - chi / 2;
    assert_eq!(genus, 0, "Blobs genus should be 0, got {}", genus);
}
