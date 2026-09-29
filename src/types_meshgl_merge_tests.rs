// types_meshgl_merge_tests.rs — regression tests for MeshGLP::merge
// (types_meshgl.rs), the port of C++ MergeMeshGLP in src/sort.cpp.
//
// The cases here pin the open-edge bookkeeping that used to diverge from the
// C++ (docs/CPP_DIVERGENCES.md, retired entry 6): C++ keeps open halfedges in
// a std::multiset, so a halfedge listed twice survives one reverse match, and
// every remaining open halfedge contributes one open-vertex entry.

use super::{MeshGL, MeshGL64, MeshGLP, MeshIndex, MeshPrecision};

/// Unit tetrahedron (verts 0..3) with face (0,2,1) listed twice first, plus a
/// separate open triangle (4,5,6) whose vert 4 coincides with vert 0.
fn doubled_face_probe<P: MeshPrecision, I: MeshIndex>() -> MeshGLP<P, I> {
    let pos: [[f64; 3]; 7] = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, -1.0, 0.0],
    ];
    let tris: [usize; 18] = [
        0, 2, 1, // listed twice before its reverse halfedges arrive
        0, 2, 1, //
        0, 1, 3, //
        1, 2, 3, //
        0, 3, 2, //
        4, 5, 6, // open triangle
    ];
    MeshGLP {
        num_prop: I::from_usize(3),
        vert_properties: pos.iter().flatten().map(|&v| P::from_f64(v)).collect(),
        tri_verts: tris.iter().map(|&v| I::from_usize(v)).collect(),
        ..Default::default()
    }
}

/// C++ MergeMeshGLP: the doubled face leaves one copy of each of its three
/// halfedges open, so verts 0, 1, 2 are open alongside 4, 5, 6 and the
/// coincident vert 4 merges into 0. The old BTreeSet port deduplicated the
/// doubled halfedges, left only 4, 5, 6 open and found no merge.
#[test]
fn merge_keeps_duplicate_open_halfedges_like_cpp_multiset() {
    let mut mesh: MeshGL = doubled_face_probe();
    assert!(mesh.merge());
    assert_eq!(mesh.merge_from_vert, vec![4u32]);
    assert_eq!(mesh.merge_to_vert, vec![0u32]);
}

/// C++ MergeMeshGLP is a template instantiated for MeshGL64 too; the same
/// probe must merge identically at double precision / 64-bit indices.
#[test]
fn merge_meshgl64_matches_meshgl_on_doubled_face_probe() {
    let mut mesh: MeshGL64 = doubled_face_probe();
    assert!(mesh.merge());
    assert_eq!(mesh.merge_from_vert, vec![4u64]);
    assert_eq!(mesh.merge_to_vert, vec![0u64]);
    // Second pass: the doubled face and the open triangle still leave open
    // edges (true), but no new coincident pairs, so the result is unchanged.
    assert!(mesh.merge());
    assert_eq!(mesh.merge_from_vert, vec![4u64]);
    assert_eq!(mesh.merge_to_vert, vec![0u64]);
}

/// MeshGL64 on a closed mesh: no open edges, so merge reports false and
/// leaves the merge vectors untouched, as C++ returns before clearing them.
#[test]
fn merge_meshgl64_closed_mesh_returns_false() {
    let mut mesh: MeshGL64 = doubled_face_probe();
    // Drop the duplicate face and the open triangle: a clean tetrahedron.
    mesh.tri_verts = mesh.tri_verts[3..15].to_vec();
    assert!(!mesh.merge());
    assert!(mesh.merge_from_vert.is_empty());
    assert!(mesh.merge_to_vert.is_empty());
}

/// Pinched boundary, no duplicate halfedges: verts 2 and 4 each start two
/// open halfedges, so C++ lists them twice in openVerts. The partition is the
/// same either way ({0, 4, 5} merge; 2 and 3 are just over tolerance apart),
/// but the duplicate leaves change the BVH shape and hence the order in which
/// union-by-rank sees the pairs, so the representative differs: C++ picks 0,
/// the old deduplicating port picked 4 (from=[0 5], to=[4 4]). Expected values
/// come from compiling C++ v3.5.2 MergeMeshGLP (src/sort.cpp) against the
/// reference headers and running it on this exact input (found by a
/// randomized C++-vs-Rust sweep of 3000 open meshes, all of which now agree).
#[test]
fn merge_pinched_boundary_picks_cpp_representative() {
    let pos: [f32; 18] = [
        1.02016, 1.03536, 0.0, //
        0.0, -0.06496, 0.0, //
        -0.05568, 0.01168, 0.0, //
        0.07056, 0.01072, 0.0, //
        1.03808, 0.99952, 0.0, //
        0.95856, 0.98976, 0.0, //
    ];
    let mut mesh = MeshGL {
        num_prop: 3,
        vert_properties: pos.to_vec(),
        tri_verts: vec![0, 4, 2, 5, 2, 3, 5, 3, 4],
        tolerance: 0.1,
        ..Default::default()
    };
    assert!(mesh.merge());
    assert_eq!(mesh.merge_from_vert, vec![4u32, 5]);
    assert_eq!(mesh.merge_to_vert, vec![0u32, 0]);
}
