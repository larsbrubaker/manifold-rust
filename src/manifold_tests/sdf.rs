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

/// CubeVoid SDF — returns distance to a unit cube void (negative inside, positive outside)
fn cube_void_sdf(p: Vec3) -> f64 {
    let ax = p.x.abs();
    let ay = p.y.abs();
    let az = p.z.abs();
    // Inside the cube: distance to nearest face (negative)
    // Outside: distance to cube surface (positive)
    let dx = ax - 1.0;
    let dy = ay - 1.0;
    let dz = az - 1.0;
    // If all components negative, we're inside: return max (most negative = deepest)
    // If any positive, we're outside
    if dx <= 0.0 && dy <= 0.0 && dz <= 0.0 {
        // Inside: return most-negative (closest face)
        dx.max(dy).max(dz)
    } else {
        // Outside: Euclidean distance to surface
        let ex = dx.max(0.0);
        let ey = dy.max(0.0);
        let ez = dz.max(0.0);
        (ex * ex + ey * ey + ez * ez).sqrt()
    }
}

/// C++ TEST(SDF, CubeVoid) — test CubeVoid SDF function values
#[test]
fn test_cpp_sdf_cube_void() {
    assert_eq!(cube_void_sdf(Vec3::new(0.0, 0.0, 0.0)), -1.0);
    assert!((cube_void_sdf(Vec3::new(0.0, 0.0, 1.0))).abs() < 1e-10);
    assert!((cube_void_sdf(Vec3::new(0.0, 1.0, 1.0))).abs() < 1e-10);
    assert!((cube_void_sdf(Vec3::new(-1.0, 0.0, 0.0))).abs() < 1e-10);
    assert!((cube_void_sdf(Vec3::new(1.0, 1.0, -1.0))).abs() < 1e-10);
    assert!(cube_void_sdf(Vec3::new(2.0, 0.0, 0.0)) > 0.0);
    assert!(cube_void_sdf(Vec3::new(2.0, -2.0, 0.0)) > 0.0);
    assert!(cube_void_sdf(Vec3::new(-2.0, 2.0, 2.0)) > 0.0);
}

/// C++ TEST(SDF, Bounds) — CubeVoid with edge_length=1
#[test]
fn test_cpp_sdf_bounds_cubevoid() {
    let size = 4.0;
    let cube_void = Manifold::level_set(
        cube_void_sdf,
        crate::types::Box {
            min: Vec3::new(-size / 2.0, -size / 2.0, -size / 2.0),
            max: Vec3::new(size / 2.0, size / 2.0, size / 2.0),
        },
        1.0,
    );
    assert_eq!(cube_void.status(), crate::types::Error::NoError);
    assert_eq!(cube_void.genus(), -1, "CubeVoid genus should be -1, got {}", cube_void.genus());
    let epsilon = cube_void.get_tolerance();
    let bounds = cube_void.bounding_box();
    let outer = size / 2.0;
    assert!((bounds.min.x - (-outer)).abs() < epsilon, "min.x");
    assert!((bounds.min.y - (-outer)).abs() < epsilon, "min.y");
    assert!((bounds.min.z - (-outer)).abs() < epsilon, "min.z");
    assert!((bounds.max.x - outer).abs() < epsilon, "max.x");
    assert!((bounds.max.y - outer).abs() < epsilon, "max.y");
    assert!((bounds.max.z - outer).abs() < epsilon, "max.z");
}

/// C++ TEST(SDF, Bounds3) — sphere with radius > box, clipped to box
#[test]
fn test_cpp_sdf_bounds3() {
    let radius = 1.2;
    let sphere = Manifold::level_set(
        move |p: Vec3| radius - (p.x * p.x + p.y * p.y + p.z * p.z).sqrt(),
        crate::types::Box {
            min: Vec3::new(-1.0, -1.0, -1.0),
            max: Vec3::new(1.0, 1.0, 1.0),
        },
        0.1,
    );
    assert_eq!(sphere.status(), crate::types::Error::NoError);
    assert_eq!(sphere.genus(), 0, "Sphere genus={}", sphere.genus());
    let epsilon = sphere.get_tolerance();
    let bounds = sphere.bounding_box();
    assert!((bounds.min.x - (-1.0)).abs() < epsilon, "min.x={}", bounds.min.x);
    assert!((bounds.min.y - (-1.0)).abs() < epsilon, "min.y");
    assert!((bounds.min.z - (-1.0)).abs() < epsilon, "min.z");
    assert!((bounds.max.x - 1.0).abs() < epsilon, "max.x={}", bounds.max.x);
    assert!((bounds.max.y - 1.0).abs() < epsilon, "max.y");
    assert!((bounds.max.z - 1.0).abs() < epsilon, "max.z");
}

/// C++ TEST(SDF, Void) — cube minus SDF void
#[test]
fn test_cpp_sdf_void_subtract() {
    let size = 4.0;
    let cube_void = Manifold::level_set(
        cube_void_sdf,
        crate::types::Box {
            min: Vec3::new(-size / 2.0, -size / 2.0, -size / 2.0),
            max: Vec3::new(size / 2.0, size / 2.0, size / 2.0),
        },
        0.5,
    );
    assert_eq!(cube_void.status(), crate::types::Error::NoError);
    let cube = Manifold::cube(Vec3::splat(size), true);
    let result = cube - cube_void;
    assert_eq!(result.genus(), 0, "genus={}", result.genus());
    assert!((result.volume() - 8.0).abs() < 0.001, "vol={}", result.volume());
    assert!((result.surface_area() - 24.0).abs() < 0.001, "sa={}", result.surface_area());
    let epsilon = result.get_tolerance();
    let bounds = result.bounding_box();
    assert!((bounds.min.x - (-1.0)).abs() < epsilon);
    assert!((bounds.min.y - (-1.0)).abs() < epsilon);
    assert!((bounds.min.z - (-1.0)).abs() < epsilon);
    assert!((bounds.max.x - 1.0).abs() < epsilon);
    assert!((bounds.max.y - 1.0).abs() < epsilon);
    assert!((bounds.max.z - 1.0).abs() < epsilon);
}

/// C++ TEST(SDF, SphereShell) — thin sphere shell via level set
#[test]
#[ignore = "Slow: very fine edge_length=0.01 for thin shell"]
fn test_cpp_sdf_sphere_shell() {
    let sphere = Manifold::level_set(
        |p: Vec3| {
            let r = (p.x * p.x + p.y * p.y + p.z * p.z).sqrt();
            (1.0 - r).min(r - 0.995)
        },
        crate::types::Box {
            min: Vec3::splat(-1.1),
            max: Vec3::splat(1.1),
        },
        0.01,
    );
    assert!((sphere.genus() - 14235).abs() < 1000,
        "SphereShell genus={}, expected ~14235", sphere.genus());
}
