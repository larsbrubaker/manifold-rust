use super::*;
use crate::linalg::Vec3;

#[test]
fn test_convex_hull_tetrahedron() {
    let pts = vec![
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    ];
    let hull = convex_hull(&pts);
    assert_eq!(hull.num_vert(), 4);
    assert_eq!(hull.num_tri(), 4);
}

#[test]
fn test_convex_hull_cube_points() {
    let pts = vec![
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(1.0, 1.0, 0.0),
        Vec3::new(1.0, 0.0, 1.0),
        Vec3::new(0.0, 1.0, 1.0),
        Vec3::new(1.0, 1.0, 1.0),
    ];
    let hull = convex_hull(&pts);
    assert_eq!(hull.num_vert(), 8);
    assert_eq!(hull.num_tri(), 12);
}

#[test]
fn test_convex_hull_empty() {
    let hull = convex_hull(&[]);
    assert!(hull.is_empty());
}

#[test]
fn test_convex_hull_single_point() {
    let hull = convex_hull(&[Vec3::new(1.0, 2.0, 3.0)]);
    // Degenerate: should produce something (possibly degenerate mesh)
    // Just check it doesn't panic
    let _ = hull.num_tri();
}

#[test]
fn test_convex_hull_two_points() {
    let hull = convex_hull(&[Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0)]);
    let _ = hull.num_tri();
}

#[test]
fn test_convex_hull_coplanar_points() {
    let pts = vec![
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(1.0, 1.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.5, 0.5, 0.0),
    ];
    let hull = convex_hull(&pts);
    // Planar case -- should still produce a valid mesh
    assert!(hull.num_tri() > 0);
}

#[test]
fn test_convex_hull_sphere_points() {
    // Generate points on a sphere
    let mut pts = Vec::new();
    let n = 20;
    for i in 0..n {
        let phi = std::f64::consts::PI * (i as f64) / (n as f64 - 1.0);
        for j in 0..n {
            let theta = 2.0 * std::f64::consts::PI * (j as f64) / n as f64;
            pts.push(Vec3::new(
                phi.sin() * theta.cos(),
                phi.sin() * theta.sin(),
                phi.cos(),
            ));
        }
    }
    let hull = convex_hull(&pts);
    assert!(hull.num_tri() > 0);
    // All vertices should be at distance ~1 from origin
    for v in &hull.vert_pos {
        let r = (v.x * v.x + v.y * v.y + v.z * v.z).sqrt();
        assert!((r - 1.0).abs() < 0.01, "vertex not on unit sphere: r={}", r);
    }
}

#[test]
fn test_convex_hull_interior_points_excluded() {
    // Cube corners + interior point
    let pts = vec![
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(1.0, 1.0, 0.0),
        Vec3::new(1.0, 0.0, 1.0),
        Vec3::new(0.0, 1.0, 1.0),
        Vec3::new(1.0, 1.0, 1.0),
        Vec3::new(0.5, 0.5, 0.5), // interior point
    ];
    let hull = convex_hull(&pts);
    assert_eq!(hull.num_vert(), 8); // interior point should be excluded
    assert_eq!(hull.num_tri(), 12);
}

#[test]
fn test_convex_hull_is_convex() {
    let pts = vec![
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(1.0, 1.0, 0.0),
        Vec3::new(1.0, 0.0, 1.0),
        Vec3::new(0.0, 1.0, 1.0),
        Vec3::new(1.0, 1.0, 1.0),
    ];
    let hull = convex_hull(&pts);
    assert!(hull.is_convex());
}

// A defect found by the convex-dilation union tree (manifold-sharp) on Thingi10K
// 63451 dilated by Sphere(0.3, 8): the hull of one flat triangle swept by the sphere
// was not convex. At iteration 29 the apex lay on the line through two hull points,
// so in the plane of both faces on that edge; the float plane distance put one face
// at 5.6e-17 (visible) and the other at 0 (hidden), so the edge became a horizon edge
// and the new face had zero area and a noise normal. Divergence ledger entry 11
// decides visibility with the exact orientation instead (`quickhull::is_above`).

/// The worst distance any input point lies outside any face of its hull, after
/// sweeping `triangle` by `tool` the way the Minkowski sum does (every corner plus
/// every tool vertex, in the same order).
fn worst_outside_of_swept_hull(triangle: &[Vec3; 3], tool: &crate::manifold::Manifold) -> f64 {
    let tool = tool.as_impl();
    let mut points = Vec::with_capacity(3 * tool.vert_pos.len());
    for &corner in triangle {
        for &tool_vert in &tool.vert_pos {
            points.push(corner + tool_vert);
        }
    }

    let hull = convex_hull(&points);
    assert_eq!(hull.status, crate::types::Error::NoError);

    let mut worst_outside = 0.0f64;
    for tri in 0..hull.num_tri() {
        let p0 = hull.vert_pos[hull.halfedge[3 * tri].start_vert as usize];
        let p1 = hull.vert_pos[hull.halfedge[3 * tri + 1].start_vert as usize];
        let p2 = hull.vert_pos[hull.halfedge[3 * tri + 2].start_vert as usize];
        let normal = crate::linalg::normalize(crate::linalg::cross(p1 - p0, p2 - p0));
        for &point in &points {
            worst_outside = worst_outside.max(crate::linalg::dot(normal, point - p0));
        }
    }
    worst_outside
}

#[test]
fn test_hull_of_a_flat_triangle_swept_by_sphere_is_convex() {
    // Triangle 163 of Thingi10K 63451 after the demo import, bit-exact.
    let triangle = [
        Vec3::new(-0.373046875, 0.33203125, 0.0234375),
        Vec3::new(-0.5078125, 0.466796875, 0.0234375),
        Vec3::new(-0.5078125, 0.197265625, 0.0234375),
    ];
    let worst = worst_outside_of_swept_hull(&triangle, &crate::manifold::Manifold::sphere(0.3, 8));
    assert!(worst < 1e-12, "a hull vertex lies {worst} outside a face");
}

#[test]
fn test_thingi641145_triangle109_swept_hull_is_convex() {
    // Triangle 109 of Thingi10K 641145 after the demo import, bit-exact; the sweep
    // radius is 0.02 of the part's diagonal. Before the exact visibility test, input
    // points lay 0.234 outside the hull.
    let triangle = [
        Vec3::new(-0.6162518858909607, 0.6155887842178345, 0.3000994920730591),
        Vec3::new(-0.7850430011749268, 0.6208153367042542, 0.3000994920730591),
        Vec3::new(-0.6162518858909607, -0.6155887246131897, 0.3000994920730591),
    ];
    let worst = worst_outside_of_swept_hull(
        &triangle,
        &crate::manifold::Manifold::sphere(0.05781898171099809, 12),
    );
    assert!(worst < 1e-12, "a hull vertex lies {worst} outside a face");
}
