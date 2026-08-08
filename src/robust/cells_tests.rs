// robust/cells_tests.rs — Cell complex and winding propagation tests.
//
// The invariants under test are the ones that make the arrangement
// formulation robust: the number of cells matches the geometry, every cell
// reachable from the outside gets a winding, and the winding in each region
// equals how many operands actually contain it.

use super::*;
use crate::linalg::Vec3;
use crate::robust::intersection_graph::build_graph;
use crate::types::{Error, OpType};

/// An axis-aligned box as 12 outward-oriented triangles.
fn cube(lo: Vec3, hi: Vec3) -> Vec<[Vec3; 3]> {
    let v = [
        Vec3::new(lo.x, lo.y, lo.z),
        Vec3::new(hi.x, lo.y, lo.z),
        Vec3::new(hi.x, hi.y, lo.z),
        Vec3::new(lo.x, hi.y, lo.z),
        Vec3::new(lo.x, lo.y, hi.z),
        Vec3::new(hi.x, lo.y, hi.z),
        Vec3::new(hi.x, hi.y, hi.z),
        Vec3::new(lo.x, hi.y, hi.z),
    ];
    let idx = [
        [0, 3, 2],
        [0, 2, 1], // −z
        [4, 5, 6],
        [4, 6, 7], // +z
        [0, 1, 5],
        [0, 5, 4], // −y
        [3, 7, 6],
        [3, 6, 2], // +y
        [0, 4, 7],
        [0, 7, 3], // −x
        [1, 2, 6],
        [1, 6, 5], // +x
    ];
    idx.iter().map(|t| [v[t[0]], v[t[1]], v[t[2]]]).collect()
}

/// Winding numbers for every cell, each component anchored by an exact
/// query.
fn all_windings(
    graph: &IntersectionGraph,
    complex: &CellComplex,
    p: &[[Vec3; 3]],
    q: &[[Vec3; 3]],
) -> Windings {
    windings(graph, complex, [p, q])
}

/// Winding of the region just inside `mesh`'s surface, by looking at the
/// anti side of any piece belonging to it.
fn inside_winding(
    graph: &IntersectionGraph,
    complex: &CellComplex,
    wind: &Windings,
    mesh: u8,
) -> Vec<[i32; 2]> {
    let mut seen = Vec::new();
    for (pi, piece) in graph.pieces.iter().enumerate() {
        if piece.mesh != mesh {
            continue;
        }
        let c = complex.cell(pi, ANTI);
        if wind.known[c] && !seen.contains(&wind.w[c]) {
            seen.push(wind.w[c]);
        }
    }
    seen
}

#[test]
fn disjoint_cubes_have_three_cells() {
    let p = cube(Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 1.0, 1.0));
    let q = cube(Vec3::new(2.0, 0.0, 0.0), Vec3::new(3.0, 1.0, 1.0));
    let graph = build_graph(&p, &q);
    let complex = build_cells(&graph);

    // Each shell bounds its own inside/outside pair. The two outer cells are
    // one region in space, but nothing connects them combinatorially — that
    // is exactly what the per-component seeding resolves.
    assert_eq!(complex.num_cells, 4, "two disjoint shells, two cells each");

    let wind = all_windings(&graph, &complex, &p, &q);
    assert!(wind.known.iter().all(|&k| k), "seeding reaches every cell");
    assert!(
        wind.w.iter().any(|&w| w == [0, 0]),
        "the region outside both must be present"
    );
    assert_eq!(inside_winding(&graph, &complex, &wind, 0), vec![[1, 0]]);
    assert_eq!(inside_winding(&graph, &complex, &wind, 1), vec![[0, 1]]);
}

#[test]
fn overlapping_cubes_wind_to_two_in_the_lens() {
    let p = cube(Vec3::new(0.0, 0.0, 0.0), Vec3::new(2.0, 2.0, 2.0));
    let q = cube(Vec3::new(1.0, 1.0, 1.0), Vec3::new(3.0, 3.0, 3.0));
    let graph = build_graph(&p, &q);
    let complex = build_cells(&graph);

    let wind = all_windings(&graph, &complex, &p, &q);

    // The three material regions: P only, Q only, and the shared lens.
    let mut regions: Vec<[i32; 2]> = Vec::new();
    for c in 0..complex.num_cells {
        if wind.known[c] && !regions.contains(&wind.w[c]) {
            regions.push(wind.w[c]);
        }
    }
    regions.sort();
    assert_eq!(
        regions,
        vec![[0, 0], [0, 1], [1, 0], [1, 1]],
        "outside, Q only, P only, and the overlap where both wind to one"
    );
}

/// Run one operation end to end through the cell complex and assemble it.
fn boolean_via_cells(p: &[[Vec3; 3]], q: &[[Vec3; 3]], op: OpType) -> crate::manifold::Manifold {
    let graph = build_graph(p, q);
    let complex = build_cells(&graph);
    let wind = all_windings(&graph, &complex, p, q);
    let pieces = extract(&graph, &complex, &wind, op);
    crate::robust::assemble::assemble(&pieces, &graph.verts, &graph.verts_f64, |_| true, None)
}

/// Two 2³ cubes overlapping in a 1³ corner: union 15, intersection 1,
/// difference 7. Each result must also be a closed manifold — that is the
/// property derived orientation is supposed to guarantee.
#[test]
fn extracted_booleans_have_correct_volume() {
    let p = cube(Vec3::new(0.0, 0.0, 0.0), Vec3::new(2.0, 2.0, 2.0));
    let q = cube(Vec3::new(1.0, 1.0, 1.0), Vec3::new(3.0, 3.0, 3.0));

    for (op, want) in [
        (OpType::Add, 15.0),
        (OpType::Intersect, 1.0),
        (OpType::Subtract, 7.0),
    ] {
        let m = boolean_via_cells(&p, &q, op);
        assert_eq!(m.status(), Error::NoError, "{op:?} status");
        assert!(!m.as_impl().is_soup, "{op:?} must close into a manifold");
        assert!(
            (m.volume() - want).abs() < 1e-9,
            "{op:?} volume {}, want {want}",
            m.volume()
        );
    }
}

/// Inverting one operand's winding must not change the result: the cell
/// labels decide orientation, so a reversed input shell still yields the
/// same solid. This is the property that fixes the inverted-orientation
/// class of NotClosed failures.
#[test]
fn inverted_operand_yields_the_same_union() {
    let p = cube(Vec3::new(0.0, 0.0, 0.0), Vec3::new(2.0, 2.0, 2.0));
    let q = cube(Vec3::new(1.0, 1.0, 1.0), Vec3::new(3.0, 3.0, 3.0));
    let flipped: Vec<[Vec3; 3]> = q.iter().map(|t| [t[0], t[2], t[1]]).collect();

    let m = boolean_via_cells(&p, &flipped, OpType::Add);
    assert_eq!(m.status(), Error::NoError, "inverted-operand union status");
    // The reversed shell bounds {w <= -1}, so it contributes no material;
    // the union is P alone, and critically the result still closes.
    assert!(!m.as_impl().is_soup, "result must still be a manifold");
    assert!((m.volume() - 8.0).abs() < 1e-9, "volume {}", m.volume());
}

/// A doubled shell must step the winding by two, not one — this is the
/// multiplicity behaviour that lets self-overlapping scans classify
/// correctly without an explicit regularization pass.
#[test]
fn doubled_shell_steps_winding_by_two() {
    let mut p = cube(Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 1.0, 1.0));
    let dup = p.clone();
    p.extend(dup); // every facet present twice, same orientation
    let q = cube(Vec3::new(5.0, 0.0, 0.0), Vec3::new(6.0, 1.0, 1.0));
    let graph = build_graph(&p, &q);
    let complex = build_cells(&graph);

    let wind = all_windings(&graph, &complex, &p, &q);
    assert_eq!(
        inside_winding(&graph, &complex, &wind, 0),
        vec![[2, 0]],
        "a double cover winds to two inside"
    );
}
