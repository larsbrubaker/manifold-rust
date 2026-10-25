// boolean3_compose_tests.rs — Unit tests for compose_meshes' private
// helpers (boolean3.rs), kept apart from boolean3_tests.rs, which is at the
// file-length limit.

use super::*;
use crate::linalg::{mat4_to_mat3x4, translation_matrix, Vec3};

/// compose_meshes' face-normal concatenation: an input with a full set is
/// copied verbatim; one without falls back to geometric normals, where a
/// removed triangle (start_vert < 0) and a degenerate one both get (0, 0, 1),
/// as set_normals_and_coplanar assigns them, instead of an out-of-bounds
/// index or a NaN.
#[test]
fn test_concat_face_normals_fallback_handles_removed_and_degenerate_tris() {
    let full = ManifoldImpl::cube(&mat4_to_mat3x4(translation_matrix(Vec3::new(
        0.0, 0.0, 0.0,
    ))));
    let mut bare = full.clone();
    bare.face_normal.clear();
    for h in 0..3 {
        bare.halfedge[h].start_vert = -1;
        bare.halfedge[h].end_vert = -1;
        bare.halfedge[3 + h].start_vert = 0;
    }
    let normals = concat_face_normals(&[full.clone(), bare]);
    let n = full.num_tri();
    assert_eq!(normals.len(), 2 * n);
    assert_eq!(&normals[..n], &full.face_normal[..]);
    let up = Vec3::new(0.0, 0.0, 1.0);
    assert_eq!(normals[n], up, "removed triangle");
    assert_eq!(normals[n + 1], up, "degenerate triangle");
    for (i, v) in normals[n + 2..].iter().enumerate() {
        assert_eq!(*v, full.face_normal[i + 2], "triangle {}", i + 2);
    }
}
