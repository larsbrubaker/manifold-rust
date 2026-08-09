// robust/arrangement_tests.rs — Unit tests for the per-triangle 2D
// arrangement (robust/arrangement.rs): segment splitting at crossings and
// at points, provenance tracking, constraint preservation through the CDT,
// and exact area conservation.

use std::collections::BTreeSet;

use super::super::exact::backend::{rat_abs, rat_new, rat_zero, Rational};

use crate::linalg::Vec3;

use super::super::exact::rational::{R2, R3};
use super::super::exact::Sign;
use super::{build, Arrangement, ArrangementInput};

fn v(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z)
}

fn r3(x: f64, y: f64, z: f64) -> R3 {
    R3::from_vec3(v(x, y, z))
}

const TRI: [Vec3; 3] = [
    Vec3 { x: 0.0, y: 0.0, z: 0.0 },
    Vec3 { x: 8.0, y: 0.0, z: 0.0 },
    Vec3 { x: 0.0, y: 8.0, z: 0.0 },
];

fn area2(a: &R2, b: &R2, c: &R2) -> Rational {
    b.sub(a).cross(&c.sub(a))
}

/// Shared checks: CCW sub-triangles, exact area conservation, every
/// constraint edge realized in the triangulation.
fn validate(arr: &Arrangement) {
    let mut total = rat_zero();
    for t in &arr.tris {
        let a2 = area2(
            &arr.points2[t[0]],
            &arr.points2[t[1]],
            &arr.points2[t[2]],
        );
        assert_eq!(Sign::of_rat(&a2), Sign::Pos, "sub-triangle not CCW");
        total = total + a2;
    }
    let base = rat_abs(&area2(&arr.points2[0], &arr.points2[1], &arr.points2[2]));
    assert_eq!(total, base, "area not conserved");

    let edges: BTreeSet<(usize, usize)> = arr
        .tris
        .iter()
        .flat_map(|t| (0..3).map(|e| {
            let (u, w) = (t[e], t[(e + 1) % 3]);
            (u.min(w), u.max(w))
        }))
        .collect();
    for edge in arr.constraints.keys() {
        assert!(edges.contains(edge), "constraint edge {edge:?} not realized");
    }
}

#[test]
fn empty_input_yields_single_triangle() {
    let arr = build(TRI, &ArrangementInput::default(), None).unwrap();
    assert_eq!(arr.tris.len(), 1);
    assert_eq!(arr.points3.len(), 3);
    assert!(arr.constraints.is_empty());
    validate(&arr);
}

#[test]
fn single_segment_becomes_one_constraint() {
    let input = ArrangementInput {
        points: vec![],
        segments: vec![(r3(1.0, 1.0, 0.0), r3(3.0, 2.0, 0.0), 7)],
    };
    let arr = build(TRI, &input, None).unwrap();
    assert_eq!(arr.constraints.len(), 1);
    let provs = arr.constraints.values().next().unwrap();
    assert_eq!(provs.as_slice(), &[7]);
    validate(&arr);
}

#[test]
fn crossing_segments_split_at_intersection() {
    // Two segments crossing at (2,2); each splits in half → 4 constraints,
    // and the crossing point must exist as an arrangement vertex.
    let input = ArrangementInput {
        points: vec![],
        segments: vec![
            (r3(1.0, 1.0, 0.0), r3(3.0, 3.0, 0.0), 0),
            (r3(1.0, 3.0, 0.0), r3(3.0, 1.0, 0.0), 1),
        ],
    };
    let arr = build(TRI, &input, None).unwrap();
    assert_eq!(arr.constraints.len(), 4);
    assert!(
        arr.points3.contains(&r3(2.0, 2.0, 0.0)),
        "crossing point missing"
    );
    // Each sub-constraint carries exactly its parent's provenance.
    for (edge, provs) in &arr.constraints {
        assert_eq!(provs.len(), 1, "edge {edge:?} has provenance {provs:?}");
    }
    validate(&arr);
}

#[test]
fn point_primitive_splits_segment() {
    // An isolated point primitive lying on a segment splits it.
    let input = ArrangementInput {
        points: vec![(r3(2.0, 1.0, 0.0), 5)],
        segments: vec![(r3(1.0, 1.0, 0.0), r3(3.0, 1.0, 0.0), 9)],
    };
    let arr = build(TRI, &input, None).unwrap();
    assert_eq!(arr.constraints.len(), 2, "point on segment must split it");
    for provs in arr.constraints.values() {
        assert_eq!(provs.as_slice(), &[9]);
    }
    validate(&arr);
}

#[test]
fn collinear_overlapping_segments_merge_provenance() {
    // Two overlapping collinear segments: [1,4] and [2,6] on y=1. Points
    // 1,2,4,6 → sub-edges [1,2](prov 0), [2,4](both), [4,6](prov 1).
    let input = ArrangementInput {
        points: vec![],
        segments: vec![
            (r3(1.0, 1.0, 0.0), r3(4.0, 1.0, 0.0), 0),
            (r3(2.0, 1.0, 0.0), r3(6.0, 1.0, 0.0), 1),
        ],
    };
    let arr = build(TRI, &input, None).unwrap();
    assert_eq!(arr.constraints.len(), 3);
    let mut prov_sets: Vec<Vec<usize>> = arr.constraints.values().cloned().collect();
    for p in &mut prov_sets {
        p.sort();
    }
    prov_sets.sort();
    assert_eq!(prov_sets, vec![vec![0], vec![0, 1], vec![1]]);
    validate(&arr);
}

#[test]
fn segment_touching_triangle_edge_and_corner() {
    // Segment from a corner to the middle of the opposite edge — endpoints
    // dedup against the corner registry and split the hull edge.
    let input = ArrangementInput {
        points: vec![],
        segments: vec![(r3(0.0, 0.0, 0.0), r3(4.0, 4.0, 0.0), 3)],
    };
    let arr = build(TRI, &input, None).unwrap();
    assert_eq!(arr.points3.len(), 4, "only the midpoint is new");
    assert_eq!(arr.constraints.len(), 1);
    validate(&arr);
}

#[test]
fn polygon_boundary_ring() {
    // A coplanar-overlap hexagon boundary as six segments with one
    // provenance — mimics how the pipeline feeds Coplanar results.
    let hex = [
        r3(2.0, 1.0, 0.0),
        r3(3.0, 1.0, 0.0),
        r3(4.0, 2.0, 0.0),
        r3(3.0, 3.0, 0.0),
        r3(2.0, 3.0, 0.0),
        r3(1.0, 2.0, 0.0),
    ];
    let segments = (0..6)
        .map(|i| (hex[i].clone(), hex[(i + 1) % 6].clone(), 11))
        .collect();
    let arr = build(TRI, &ArrangementInput {
        points: vec![],
        segments,
    }, None)
    .unwrap();
    assert_eq!(arr.constraints.len(), 6);
    validate(&arr);
}

#[test]
fn skew_plane_arrangement_lifts_exactly() {
    // Non-axis-aligned triangle: verify all 3D points still lie exactly on
    // its plane after crossing construction + lifting.
    let tri = [v(0.1, 0.2, 0.3), v(6.7, -0.4, 0.9), v(-0.6, 6.1, 2.2)];
    // Two crossing segments built from midpoints (guaranteed on the plane).
    let a = R3::from_vec3(tri[0]);
    let b = R3::from_vec3(tri[1]);
    let c = R3::from_vec3(tri[2]);
    let half = rat_new(1.into(), 2.into());
    let mab = a.add(&b).scale(&half);
    let mbc = b.add(&c).scale(&half);
    let mca = c.add(&a).scale(&half);
    let input = ArrangementInput {
        points: vec![],
        segments: vec![(mab.clone(), mbc.clone(), 0), (mbc.clone(), mca.clone(), 1), (mca, mab, 2)],
    };
    let arr = build(tri, &input, None).unwrap();
    use crate::robust::exact::predicates::orient3d_r;
    for p in &arr.points3 {
        assert_eq!(orient3d_r(&a, &b, &c, p), Sign::Zero, "point off plane");
    }
    validate(&arr);
    // The midpoint triangle splits the face into 4 regions ≥ 4 sub-tris.
    assert!(arr.tris.len() >= 4);
}

#[test]
fn constraint_edges_reference_valid_points() {
    let input = ArrangementInput {
        points: vec![],
        segments: vec![
            (r3(1.0, 1.0, 0.0), r3(5.0, 1.0, 0.0), 0),
            (r3(2.0, 0.0, 0.0), r3(2.0, 4.0, 0.0), 1),
        ],
    };
    let arr = build(TRI, &input, None).unwrap();
    for (u, w) in arr.constraints.keys() {
        assert!(*u < arr.points3.len() && *w < arr.points3.len());
        assert_ne!(arr.points2[*u], arr.points2[*w]);
    }
    validate(&arr);
}
