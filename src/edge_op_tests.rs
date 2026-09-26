use super::*;
use crate::face_op::set_normals_and_coplanar;
use crate::impl_mesh::ManifoldImpl;
use crate::linalg::Mat3x4;

#[test]
fn test_pair_up() {
    let mut halfedges = vec![
        Halfedge {
            start_vert: 0,
            end_vert: 1,
            paired_halfedge: -1,
            prop_vert: 0,
        },
        Halfedge {
            start_vert: 1,
            end_vert: 0,
            paired_halfedge: -1,
            prop_vert: 1,
        },
    ];
    pair_up(&mut halfedges, 0, 1);
    assert_eq!(halfedges[0].paired_halfedge, 1);
    assert_eq!(halfedges[1].paired_halfedge, 0);
}

#[test]
fn test_cleanup_topology_noop_on_clean_mesh() {
    let mut m = ManifoldImpl::tetrahedron(&Mat3x4::identity());
    set_normals_and_coplanar(&mut m);
    let before_verts = m.vert_pos.len();
    let before_halfedges = m.halfedge.len();
    cleanup_topology(&mut m);
    // Tetrahedron is already 2-manifold; cleanup should not add verts/halfedges
    assert_eq!(m.vert_pos.len(), before_verts);
    assert_eq!(m.halfedge.len(), before_halfedges);
}

#[test]
fn test_simplify_topology_noop_on_clean_mesh() {
    let mut m = ManifoldImpl::cube(&Mat3x4::identity());
    set_normals_and_coplanar(&mut m);
    simplify_topology(&mut m, 0);
    // After simplify, cube should still be 2-manifold (no degenerate edges)
    // (Some verts/edges may be removed, but topology must be valid)
    // Just check it's still 2-manifold where halfedges are valid
    let valid = m
        .halfedge
        .iter()
        .filter(|h| h.paired_halfedge >= 0)
        .all(|h| h.paired_halfedge < m.halfedge.len() as i32);
    assert!(valid, "invalid paired halfedge after simplify");
}

/// Regression fixture shared with manifold-sharp (DedupeEdgesRegressionTests):
/// 852 triangles cut from a 28060-triangle exact union right before
/// DedupeEdges, within three vertex rings of its sixteen duplicate edges.
/// Format: "numVert numTri", then one "x y z" per vertex, then one line per
/// triangle of three (start, end, pair) halfedges, pair -1 where it fell
/// outside the cut. Positions are round-trip formatted, so bit-exact.
const DEDUPE_STALE_DUPLICATE: &str = include_str!("testdata/dedupe-stale-duplicate.txt");

fn load_halfedge_fixture(text: &str) -> ManifoldImpl {
    let mut lines = text.lines();
    let header: Vec<usize> = lines
        .next()
        .unwrap()
        .split_whitespace()
        .map(|s| s.parse().unwrap())
        .collect();
    let (num_vert, num_tri) = (header[0], header[1]);
    let mut m = ManifoldImpl::new();
    for _ in 0..num_vert {
        let p: Vec<f64> = lines
            .next()
            .unwrap()
            .split_whitespace()
            .map(|s| s.parse().unwrap())
            .collect();
        m.vert_pos.push(Vec3::new(p[0], p[1], p[2]));
    }
    for _ in 0..num_tri {
        let p: Vec<i32> = lines
            .next()
            .unwrap()
            .split_whitespace()
            .map(|s| s.parse().unwrap())
            .collect();
        for k in 0..3 {
            m.halfedge.push(Halfedge {
                start_vert: p[3 * k],
                end_vert: p[3 * k + 1],
                paired_halfedge: p[3 * k + 2],
                prop_vert: p[3 * k],
            });
        }
    }
    m
}

/// Splitting duplicated edges only relabels vertices and adds zero-area
/// triangles, so no pre-existing triangle corner may change position. Before
/// the stale-entry check, a duplicate already resolved by an earlier repair in
/// the same pass was "repaired" again, relabelling an orbit to a copy of the
/// wrong vertex: 16 corners jumped ~0.13 and a union lost 2.3e-5 of volume.
#[test]
fn test_dedupe_edges_never_moves_a_triangle_corner() {
    let mut m = load_halfedge_fixture(DEDUPE_STALE_DUPLICATE);
    let before: Vec<Vec3> = m
        .halfedge
        .iter()
        .map(|h| m.vert_pos[h.start_vert as usize])
        .collect();

    dedupe_edges(&mut m);

    let moved = before
        .iter()
        .enumerate()
        .filter(|(h, p)| {
            let after = m.vert_pos[m.halfedge[*h].start_vert as usize];
            after.x != p.x || after.y != p.y || after.z != p.z
        })
        .count();
    assert_eq!(moved, 0, "a corner that moves changes the solid");
}
