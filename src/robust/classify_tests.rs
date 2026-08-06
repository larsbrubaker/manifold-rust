// robust/classify_tests.rs — Integration tests for the classification
// pipeline: intersection_graph → classify → propagate → ray_shoot on
// hand-built solids. The strongest oracle is exact and independent: a piece
// of P tagged Union must have its centroid outside Q (winding 0), a piece
// tagged Inter inside Q — checked per piece wherever the piece does not lie
// on the other mesh's surface.

use num_rational::BigRational;
use num_traits::{Signed, Zero};

use crate::linalg::Vec3;

use super::super::intersection_graph::{build_graph, IntersectionGraph};
use super::super::propagate::propagate;
use super::super::ray_shoot::{piece_centroid, winding_number};
use super::{classify_rings, Classification, Tag};

fn v(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z)
}

/// Axis-aligned cube [lo,hi]³ as 12 outward-wound triangles.
fn cube(lo: f64, hi: f64) -> Vec<[Vec3; 3]> {
    cube_at([lo, lo, lo], [hi, hi, hi])
}

fn cube_at(lo: [f64; 3], hi: [f64; 3]) -> Vec<[Vec3; 3]> {
    let p = |x: f64, y: f64, z: f64| v(x, y, z);
    let (x0, y0, z0) = (lo[0], lo[1], lo[2]);
    let (x1, y1, z1) = (hi[0], hi[1], hi[2]);
    let quads = [
        // -z, +z, -y, +y, -x, +x (outward normals)
        [p(x0, y0, z0), p(x0, y1, z0), p(x1, y1, z0), p(x1, y0, z0)],
        [p(x0, y0, z1), p(x1, y0, z1), p(x1, y1, z1), p(x0, y1, z1)],
        [p(x0, y0, z0), p(x1, y0, z0), p(x1, y0, z1), p(x0, y0, z1)],
        [p(x0, y1, z0), p(x0, y1, z1), p(x1, y1, z1), p(x1, y1, z0)],
        [p(x0, y0, z0), p(x0, y0, z1), p(x0, y1, z1), p(x0, y1, z0)],
        [p(x1, y0, z0), p(x1, y1, z0), p(x1, y1, z1), p(x1, y0, z1)],
    ];
    let mut out = Vec::new();
    for q in quads {
        out.push([q[0], q[1], q[2]]);
        out.push([q[0], q[2], q[3]]);
    }
    out
}

/// Run the full classification pipeline; returns (graph, per-piece final
/// tags, discarded flags).
fn run(
    p: &[[Vec3; 3]],
    q: &[[Vec3; 3]],
    q_complement: bool,
) -> (IntersectionGraph, Vec<Option<Tag>>, Vec<bool>) {
    let graph = build_graph(p, q);
    let cls: Classification = classify_rings(&graph);
    let prop = propagate(&graph, &cls);
    let mut tags = prop.tags.clone();
    // Ray-shoot untagged components against the *other* mesh.
    for &(root, rep) in &prop.untagged {
        let piece = &graph.pieces[rep];
        let other: &[[Vec3; 3]] = if piece.mesh == 0 { q } else { p };
        let other_complement = if piece.mesh == 0 { q_complement } else { false };
        let centroid = piece_centroid(graph.piece_verts(rep));
        let w = winding_number(&centroid, other);
        let inside = if other_complement { w == 0 } else { w != 0 };
        let tag = if inside { Tag::Inter } else { Tag::Union };
        for pi in 0..graph.pieces.len() {
            if !cls.discarded[pi] && prop.component[pi] == root {
                tags[pi] = Some(tag);
            }
        }
    }
    (graph, tags, cls.discarded)
}

/// Independent oracle: every tagged piece not lying on the other surface
/// must agree with an exact winding query of its centroid.
fn check_against_winding_oracle(
    graph: &IntersectionGraph,
    tags: &[Option<Tag>],
    discarded: &[bool],
    p: &[[Vec3; 3]],
    q: &[[Vec3; 3]],
) {
    for (pi, piece) in graph.pieces.iter().enumerate() {
        if discarded[pi] {
            continue;
        }
        let other: &[[Vec3; 3]] = if piece.mesh == 0 { q } else { p };
        let centroid = piece_centroid(graph.piece_verts(pi));
        // Skip pieces on the other mesh's surface (winding undefined there):
        // detect via any zero orient3d against a containing plane — cheap
        // proxy: coplanar-overlap pieces are exactly the ones both meshes
        // share, and their centroids sit on the other surface. We detect by
        // testing the winding parity stability under a tiny exact offset
        // along the piece normal instead — simpler: skip when the centroid
        // lies on any triangle plane of `other` within that triangle.
        if centroid_on_surface(&centroid, other) {
            continue;
        }
        let tag = tags[pi].expect("piece left untagged");
        let inside = winding_number(&centroid, other) != 0;
        let expect = if inside { Tag::Inter } else { Tag::Union };
        assert_eq!(
            tag, expect,
            "piece {pi} (mesh {}, tri {}) misclassified",
            piece.mesh, piece.tri
        );
    }
}

fn centroid_on_surface(c: &super::super::exact::rational::R3, tris: &[[Vec3; 3]]) -> bool {
    use super::super::exact::predicates::{orient3d_r, point_in_tri_2d, tri_normal_r, TriLoc};
    use super::super::exact::rational::R3;
    use super::super::exact::Sign;
    use super::super::tri_tri::dominant_axis;
    for t in tris {
        let a = R3::from_vec3(t[0]);
        let b = R3::from_vec3(t[1]);
        let cc = R3::from_vec3(t[2]);
        if orient3d_r(&a, &b, &cc, c) != Sign::Zero {
            continue;
        }
        let n = tri_normal_r(&a, &b, &cc);
        let axis = dominant_axis(&n);
        let loc = point_in_tri_2d(
            &c.project_drop(axis),
            &a.project_drop(axis),
            &b.project_drop(axis),
            &cc.project_drop(axis),
        );
        if loc != TriLoc::Outside {
            return true;
        }
    }
    false
}

/// Exact double-area of a piece, valid for axis-plane-aligned geometry used
/// in these tests (cross product has one nonzero component).
fn piece_area2(v: [&super::super::exact::rational::R3; 3]) -> BigRational {
    let n = v[1].sub(v[0]).cross(&v[2].sub(v[0]));
    let comps = [n.x.abs(), n.y.abs(), n.z.abs()];
    let mut nonzero: Vec<&BigRational> = comps.iter().filter(|c| !c.is_zero()).collect();
    assert!(nonzero.len() <= 1, "test helper requires axis-aligned pieces");
    nonzero.pop().cloned().unwrap_or_else(BigRational::zero)
}

#[test]
fn overlapping_cubes_classify_correctly() {
    let p = cube(0.0, 2.0);
    let q_tris = cube_at([1.0, 1.0, 1.0], [3.0, 3.0, 3.0]);
    let (graph, tags, discarded) = run(&p, &q_tris, false);
    assert!(graph.any_intersections);
    assert!(discarded.iter().all(|d| !d));
    assert!(tags.iter().zip(&discarded).all(|(t, d)| *d || t.is_some()));
    check_against_winding_oracle(&graph, &tags, &discarded, &p, &q_tris);
}

#[test]
fn skewered_cube_classifies_correctly() {
    // A long thin box punching straight through a larger cube — every face
    // of the small box is cut by the big cube's faces.
    let p = cube(0.0, 4.0);
    let q_tris = cube_at([1.5, 1.5, -2.0], [2.5, 2.5, 6.0]);
    let (graph, tags, discarded) = run(&p, &q_tris, false);
    assert!(graph.any_intersections);
    check_against_winding_oracle(&graph, &tags, &discarded, &p, &q_tris);
}

#[test]
fn nested_cube_resolved_by_ray_shooting() {
    // Small cube strictly inside the big one — no intersections at all;
    // classification must come entirely from winding numbers.
    let p = cube(0.0, 6.0);
    let q_tris = cube(2.0, 4.0);
    let (graph, tags, discarded) = run(&p, &q_tris, false);
    assert!(!graph.any_intersections);
    assert!(discarded.iter().all(|d| !d));
    for (pi, piece) in graph.pieces.iter().enumerate() {
        let expect = if piece.mesh == 0 { Tag::Union } else { Tag::Inter };
        assert_eq!(tags[pi], Some(expect), "piece {pi}");
    }
}

#[test]
fn face_touching_cubes_discard_shared_face() {
    // Stacked cubes sharing the z=2 plane: P's top face and Q's bottom face
    // coincide with opposite orientations → regularization discards both.
    let p = cube_at([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
    let q_tris = cube_at([0.0, 0.0, 2.0], [2.0, 2.0, 4.0]);
    let (graph, tags, discarded) = run(&p, &q_tris, false);
    assert!(graph.any_intersections);
    // Discarded area on each side equals the full shared face (2×2 = area
    // 4, double-area 8).
    let mut discarded_area2 = [BigRational::zero(), BigRational::zero()];
    for (pi, piece) in graph.pieces.iter().enumerate() {
        if discarded[pi] {
            discarded_area2[piece.mesh as usize] =
                &discarded_area2[piece.mesh as usize] + piece_area2(graph.piece_verts(pi));
        }
    }
    let eight = BigRational::from_integer(8.into());
    assert_eq!(discarded_area2[0], eight, "P-side discarded area");
    assert_eq!(discarded_area2[1], eight, "Q-side discarded area");
    // Everything not discarded is outside the other solid → Union.
    for (pi, piece) in graph.pieces.iter().enumerate() {
        if discarded[pi] {
            continue;
        }
        assert_eq!(tags[pi], Some(Tag::Union), "piece {pi} of mesh {}", piece.mesh);
    }
}

#[test]
fn identical_cubes_share_every_face_once_per_output() {
    // Equivalent operands: every face coplanar with the twin's, same
    // orientation → each region gets one Union and one Inter copy; nothing
    // discarded, nothing untagged.
    let p = cube(0.0, 2.0);
    let q_tris = cube(0.0, 2.0);
    let (graph, tags, discarded) = run(&p, &q_tris, false);
    assert!(graph.any_intersections);
    assert!(discarded.iter().all(|d| !d));
    let full_surface2 = BigRational::from_integer(48.into()); // 6 faces × area 4 × 2
    let mut union_area2 = BigRational::zero();
    let mut inter_area2 = BigRational::zero();
    for (pi, piece) in graph.pieces.iter().enumerate() {
        match tags[pi].expect("all pieces must be tagged for identical cubes") {
            Tag::Union => union_area2 = &union_area2 + piece_area2(graph.piece_verts(pi)),
            Tag::Inter => inter_area2 = &inter_area2 + piece_area2(graph.piece_verts(pi)),
        }
    }
    assert_eq!(union_area2, full_surface2, "union output must cover the cube once");
    assert_eq!(inter_area2, full_surface2, "intersection output must cover the cube once");
}

#[test]
fn edge_touching_cubes_have_no_discards_and_all_union() {
    // Two cubes sharing exactly one edge (classic non-manifold union input):
    // the contact is a segment, not an area — nothing discarded, everything
    // outside the other solid.
    let p = cube_at([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
    let q_tris = cube_at([2.0, 2.0, 0.0], [4.0, 4.0, 2.0]);
    let (graph, tags, discarded) = run(&p, &q_tris, false);
    assert!(graph.any_intersections);
    assert!(discarded.iter().all(|d| !d));
    for (pi, _) in graph.pieces.iter().enumerate() {
        assert_eq!(tags[pi], Some(Tag::Union), "piece {pi}");
    }
    check_against_winding_oracle(&graph, &tags, &discarded, &p, &q_tris);
}

#[test]
fn subtract_style_flipped_operand() {
    // P minus Q via flipped Q: overlapping cubes; ray-shooting not needed
    // here, but ring classification must handle the reversed orientations.
    let p = cube(0.0, 2.0);
    let q_flipped: Vec<[Vec3; 3]> = cube_at([1.0, 1.0, 1.0], [3.0, 3.0, 3.0])
        .iter()
        .map(|t| [t[0], t[2], t[1]])
        .collect();
    let (graph, tags, _discarded) = run(&p, &q_flipped, true);
    assert!(graph.any_intersections);
    // Difference = Inter pieces of (P, Q^c): P pieces outside Q plus Q
    // pieces inside P (flipped). Sanity: at least some of each tag on P.
    let p_union = tags
        .iter()
        .zip(&graph.pieces)
        .filter(|(t, pc)| pc.mesh == 0 && **t == Some(Tag::Union))
        .count();
    let p_inter = tags
        .iter()
        .zip(&graph.pieces)
        .filter(|(t, pc)| pc.mesh == 0 && **t == Some(Tag::Inter))
        .count();
    assert!(p_union > 0 && p_inter > 0);
}
