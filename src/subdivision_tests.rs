use super::*;
use crate::linalg::{Mat3x4, Vec3};
use crate::manifold::Manifold;
use crate::types::Error;

#[test]
fn test_subdivide_cube_once() {
    let cube = ManifoldImpl::cube(&Mat3x4::identity());
    let sub = subdivide_impl(&cube, 1);
    assert_eq!(sub.num_tri(), cube.num_tri() * 4);
    assert!(sub.num_vert() > cube.num_vert());
}

#[test]
fn test_partition_single_triangle() {
    // Simplest case: 1 division per edge = no subdivision
    let part = Partition::get_partition(IVec4::new(1, 1, 1, 0));
    assert_eq!(part.tri_vert.len(), 1);
    assert_eq!(part.vert_bary.len(), 3);
}

#[test]
fn test_partition_two_divisions() {
    // 2 divisions on each edge of a triangle
    let part = Partition::get_partition(IVec4::new(2, 2, 2, 0));
    assert_eq!(part.tri_vert.len(), 4); // 4 sub-triangles
}

#[test]
fn test_subdivide_edge_divisions() {
    // Test that the full Subdivide method works with edge_divisions callback
    let mut cube = ManifoldImpl::cube(&Mat3x4::identity());
    let _vert_bary = cube.subdivide(&|_vec, _t0, _t1| 1, false);
    assert_eq!(cube.num_tri(), 12 * 4); // each of 12 tris becomes 4
}

#[test]
fn test_create_tmp_edges() {
    let cube = ManifoldImpl::cube(&Mat3x4::identity());
    let edges = create_tmp_edges(&cube.halfedge);
    // A cube has 12 triangles = 36 halfedges = 18 edges
    assert_eq!(edges.len(), cube.halfedge.len() / 2);
}

#[test]
fn test_partition_asymmetric() {
    // Asymmetric divisions: 3,2,1
    let part = Partition::get_partition(IVec4::new(3, 2, 1, 0));
    assert!(part.tri_vert.len() > 1);
    // All barycentric coordinates should sum to 1
    for bary in &part.vert_bary {
        let sum = bary.x + bary.y + bary.z + bary.w;
        assert!(
            (sum - 1.0).abs() < 1e-10,
            "Barycentric coords should sum to 1, got {}",
            sum
        );
    }
}

// ---------------------------------------------------------------------------
// Regression: subdivide_impl must run the whole of refine's finishing tail.
//
// It used to run two of its six steps, calculate_bbox and set_epsilon, leaving
// behind a collider built over the pre-subdivision faces, a vert_normal list
// shorter than vert_pos, and any tangents the input carried still sized to the
// old halfedges. subdivide appends vertices and faces, so all three are stale by
// construction, and the boolean reads vert_normal per vertex.
// ---------------------------------------------------------------------------

/// The symptom: `subdivide_impl` must hand back a mesh a boolean can consume.
fn subdivided_cube_is_usable_in_a_boolean(levels: usize, expected_tri: usize) {
    let cube = ManifoldImpl::cube(&Mat3x4::identity());
    let subdivided = Manifold::from_impl(subdivide_impl(&cube, levels));
    let bore = Manifold::cube(Vec3::new(0.5, 0.5, 4.0), true).translate(Vec3::new(0.5, 0.5, 0.5));

    let drilled = subdivided.difference(&bore);

    assert_eq!(drilled.status(), Error::NoError);
    assert_eq!(drilled.num_tri(), expected_tri);
    // The subdivided cube is still the unit cube, so volume and area are the
    // subdivision-independent answer: 1 minus the 0.25 bore, and the six faces
    // minus the two 0.25 holes plus the four 0.5x1 bore walls.
    assert_eq!(drilled.volume(), 0.75);
    assert_eq!(drilled.surface_area(), 7.5);
}

#[test]
fn test_subdivided_cube_is_usable_in_a_boolean_one_level() {
    subdivided_cube_is_usable_in_a_boolean(1, 64);
}

#[test]
fn test_subdivided_cube_is_usable_in_a_boolean_two_levels() {
    subdivided_cube_is_usable_in_a_boolean(2, 176);
}

/// Root cause 1: subdivision appends faces, so the collider the pre-subdivision
/// mesh carried describes triangles that no longer exist. Before the repair the
/// cube's 12 leaves stood against 48 faces and 45 of them found nothing.
#[test]
fn test_subdivided_cube_collider_matches_its_own_faces() {
    for levels in 1..=2 {
        let m = subdivide_impl(&ManifoldImpl::cube(&Mat3x4::identity()), levels);
        assert_eq!(
            crate::sort::collider_self_misses(&m),
            0,
            "every face must overlap its own leaf box in the cached collider (levels={})",
            levels
        );
    }
}

/// Root cause 2: the boolean reads `vert_normal` per vertex as its symbolic
/// perturbation direction (`boolean3_kernels::shadow01`), so a vertex list that
/// outgrew its normals indexes off the end. This is the half a bare
/// `sort_geometry` does not repair — `sort_verts` permutes `vert_normal` only
/// when its length already equals the vertex count, so a stale list stays stale.
#[test]
fn test_subdivided_cube_has_a_normal_for_every_vertex() {
    for levels in 1..=2 {
        let m = subdivide_impl(&ManifoldImpl::cube(&Mat3x4::identity()), levels);
        assert_eq!(
            m.vert_normal.len(),
            m.num_vert(),
            "shadow01 indexes vert_normal by vertex (levels={})",
            levels
        );
    }
}

/// Root cause 3: subdivision multiplies the halfedges, so tangents the input
/// carried no longer describe them. Left in place, `gather_faces` copies one
/// tangent per *new* halfedge out of the old array and indexes past its end —
/// which is what would turn the collider repair into a panic on a smoothed mesh.
#[test]
fn test_subdivide_drops_tangents_it_can_no_longer_describe() {
    let smooth = Manifold::smooth(&Manifold::tetrahedron().get_mesh_gl(-1), &[]);
    let smooth = smooth.as_impl();
    assert_eq!(
        smooth.halfedge_tangent.len(),
        smooth.halfedge.len(),
        "the fixture must actually carry tangents for this to test anything"
    );

    let m = subdivide_impl(smooth, 1);

    assert_eq!(
        m.halfedge_tangent.len(),
        0,
        "tangents that no longer match the halfedge count must be dropped, not kept"
    );
}

#[test]
fn test_subdivide_preserves_manifold() {
    let cube = ManifoldImpl::cube(&Mat3x4::identity());
    let sub = subdivide_impl(&cube, 1);
    // Every halfedge should have a valid pair
    for (i, he) in sub.halfedge.iter().enumerate() {
        assert!(he.paired_halfedge >= 0, "Halfedge {} has no pair", i);
        let pair = he.paired_halfedge as usize;
        assert!(
            pair < sub.halfedge.len(),
            "Halfedge {} pair {} out of range",
            i,
            pair
        );
    }
}
