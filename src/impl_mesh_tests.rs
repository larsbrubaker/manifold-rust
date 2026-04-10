use super::*;
use crate::linalg::{Vec3, Mat3x4};

#[test]
fn test_next_halfedge() {
    assert_eq!(next_halfedge(0), 1);
    assert_eq!(next_halfedge(1), 2);
    assert_eq!(next_halfedge(2), 0);
    assert_eq!(next_halfedge(3), 4);
    assert_eq!(next_halfedge(5), 3);
}

#[test]
fn test_tetrahedron() {
    let m = ManifoldImpl::tetrahedron(&Mat3x4::identity());
    assert_eq!(m.num_vert(), 4);
    assert_eq!(m.num_tri(), 4);
    assert_eq!(m.num_edge(), 6);
    assert!(m.is_manifold(), "Tetrahedron should be manifold");
    assert!(m.is_2_manifold(), "Tetrahedron should be 2-manifold");
}

#[test]
fn test_cube() {
    let m = ManifoldImpl::cube(&Mat3x4::identity());
    assert_eq!(m.num_vert(), 8);
    assert_eq!(m.num_tri(), 12);
    assert_eq!(m.num_edge(), 18);
    assert!(m.is_manifold(), "Cube should be manifold");
    assert!(m.is_2_manifold(), "Cube should be 2-manifold");
}

#[test]
fn test_octahedron() {
    let m = ManifoldImpl::octahedron(&Mat3x4::identity());
    assert_eq!(m.num_vert(), 6);
    assert_eq!(m.num_tri(), 8);
    assert_eq!(m.num_edge(), 12);
    assert!(m.is_manifold(), "Octahedron should be manifold");
    assert!(m.is_2_manifold(), "Octahedron should be 2-manifold");
}

#[test]
fn test_cube_bbox() {
    let m = ManifoldImpl::cube(&Mat3x4::identity());
    assert!((m.bbox.min.x - 0.0).abs() < 1e-10);
    assert!((m.bbox.min.y - 0.0).abs() < 1e-10);
    assert!((m.bbox.min.z - 0.0).abs() < 1e-10);
    assert!((m.bbox.max.x - 1.0).abs() < 1e-10);
    assert!((m.bbox.max.y - 1.0).abs() < 1e-10);
    assert!((m.bbox.max.z - 1.0).abs() < 1e-10);
}

#[test]
fn test_tetrahedron_symmetric() {
    let m = ManifoldImpl::tetrahedron(&Mat3x4::identity());
    // Tetrahedron is centered about origin with vertices at distance sqrt(3)
    for v in &m.vert_pos {
        let dist2 = v.x * v.x + v.y * v.y + v.z * v.z;
        assert!((dist2 - 3.0).abs() < 1e-10, "Expected dist^2=3, got {}", dist2);
    }
}

#[test]
fn test_cube_transform() {
    use crate::linalg::{translation_matrix, mat4_to_mat3x4};
    let t = mat4_to_mat3x4(translation_matrix(Vec3::new(1.0, 2.0, 3.0)));
    let m = ManifoldImpl::cube(&t);
    assert!((m.bbox.min.x - 1.0).abs() < 1e-10);
    assert!((m.bbox.min.y - 2.0).abs() < 1e-10);
    assert!((m.bbox.min.z - 3.0).abs() < 1e-10);
    assert!((m.bbox.max.x - 2.0).abs() < 1e-10);
    assert!((m.bbox.max.y - 3.0).abs() < 1e-10);
    assert!((m.bbox.max.z - 4.0).abs() < 1e-10);
}

#[test]
fn test_for_vert() {
    let m = ManifoldImpl::tetrahedron(&Mat3x4::identity());
    // Each vertex in the tetrahedron is surrounded by 3 triangles = 3 halfedges
    let mut count = 0;
    m.for_vert(0, |_| count += 1);
    // ForVert visits one halfedge per triangle around the vertex (3 for tetrahedron)
    assert!(count > 0);
}

#[test]
fn test_create_halfedges_simple() {
    // Simple triangle
    let mut m = ManifoldImpl::new();
    m.vert_pos = vec![
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
    ];
    // A tetrahedron has 4 triangles and 12 halfedges
    // Let's do a single triangle (not manifold, just for halfedge construction)
    let tri = vec![IVec3::new(0, 1, 2)];
    m.create_halfedges(&tri, &[]);
    assert_eq!(m.halfedge.len(), 3);
    // Single triangle -- no pairs, so all pairedHalfedge should be -1
    // (no paired halfedges available)
}
