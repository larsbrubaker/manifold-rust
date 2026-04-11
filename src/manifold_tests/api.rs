use super::*;

/// C++ TEST(Properties, MinGapCubeSphereOverlapping) — overlapping returns 0
#[test]
fn test_cpp_min_gap_cube_sphere_overlapping() {
    let a = Manifold::cube(Vec3::new(1.0, 1.0, 1.0), false);
    let b = Manifold::sphere(1.0, 0);
    let distance = a.min_gap(&b, 0.1);
    assert_eq!(distance, 0.0, "MinGapCubeSphereOverlapping: {} expected 0", distance);
}

/// C++ TEST(Properties, MinGapSphereSphereOutOfBounds) — returns search_length
#[test]
fn test_cpp_min_gap_sphere_sphere_out_of_bounds() {
    let a = Manifold::sphere(1.0, 0);
    let b = Manifold::sphere(1.0, 0).translate(Vec3::new(2.0, 2.0, 0.0));
    let distance = a.min_gap(&b, 0.8);
    assert_eq!(distance, 0.8,
        "MinGapSphereSphereOutOfBounds: {} expected 0.8 (search_length)", distance);
}

/// C++ TEST(Properties, MingapAfterTransformations) — rotated/scaled spheres
#[test]
#[ignore = "Slow in debug: 512-segment spheres"]
fn test_cpp_min_gap_after_transformations() {
    let a = Manifold::sphere(1.0, 512).rotate(30.0, 30.0, 30.0);
    let b = Manifold::sphere(1.0, 512)
        .scale(Vec3::new(3.0, 1.0, 1.0))
        .rotate(0.0, 90.0, 45.0)
        .translate(Vec3::new(3.0, 0.0, 0.0));
    let distance = a.min_gap(&b, 1.1);
    assert!((distance - 1.0).abs() < 0.001,
        "MingapAfterTransformations: {} expected ~1.0", distance);
}

/// C++ TEST(Manifold, ValidInputOneRunIndex) — empty mesh with runIndex={0}
#[test]
fn test_cpp_valid_input_one_run_index() {
    let mut empty_mesh = MeshGL::default();
    empty_mesh.run_index = vec![0];
    let empty = Manifold::from_mesh_gl(&empty_mesh);
    assert!(empty.is_empty(), "ValidInputOneRunIndex: should be empty");
}

/// C++ TEST(Manifold, Empty) — default manifold is empty
#[test]
fn test_cpp_manifold_empty() {
    let empty = Manifold::empty();
    assert!(empty.is_empty());
    assert_eq!(empty.num_vert(), 0);
    assert_eq!(empty.num_tri(), 0);
    assert_eq!(empty.volume(), 0.0);
}

/// C++ TEST(Manifold, Simplify) from manifold_test.cpp — simplify cube
#[test]
fn test_cpp_manifold_simplify() {
    // Cube from MeshGL should survive simplify
    let cube = Manifold::cube(Vec3::new(1.0, 1.0, 1.0), false);
    let simplified = cube.as_original();
    assert!(!simplified.is_empty(), "Simplified cube should not be empty");
    assert_eq!(simplified.num_vert(), 8);
    assert_eq!(simplified.num_tri(), 12);
}
