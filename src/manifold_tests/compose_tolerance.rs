// Tests that composing disjoint meshes keeps the inputs' tolerance and
// epsilon — boolean3::compose_meshes, behind Manifold::compose, the exact
// engine's disjoint-union fast path (and so the pairwise batch_boolean fold)
// and csg_tree's batch_union.
//
// C++ CsgLeafNode::Compose (csg_tree.cpp:245-268) gives the combined impl
// the largest of its nodes' epsilons and tolerances, as boolean_result
// (boolean_result.cpp:832) does for an overlapping boolean. The port used to
// reset both to the automatic kPrecision * bbox.scale(), so a user's
// set_tolerance or imported MeshGL.tolerance vanished whenever a union
// happened to be disjoint. The geometry, counts and run tables are unchanged
// by the carry; every test here also pins that against the same inputs at
// default tolerance.

use super::*;
use crate::csg_tree::{CsgLeafNode, CsgNode};
use crate::linalg::Mat3x4;
use crate::types::MeshGL64;

/// Everything a composed MeshGL64 carries except its tolerance and the
/// process-wide original IDs (which differ between freshly imported
/// fixtures): the geometry, the counts and the run table.
fn shape(gl: &MeshGL64) -> (Vec<u64>, Vec<u64>, Vec<u64>, Vec<u64>, Vec<u64>, usize) {
    (
        gl.vert_properties.iter().map(|x| x.to_bits()).collect(),
        gl.tri_verts.clone(),
        gl.run_index.clone(),
        gl.run_transform.iter().map(|x| x.to_bits()).collect(),
        gl.face_id.clone(),
        gl.run_original_id.len(),
    )
}

fn unit_cube() -> Manifold {
    Manifold::cube(Vec3::new(1.0, 1.0, 1.0), false)
}

/// A unit cube imported with `MeshGL.tolerance` — which raises epsilon to
/// the tolerance as well (C++ `SetEpsilon`).
fn imported_cube(tolerance: f64) -> Manifold {
    let mut gl = unit_cube().get_mesh_gl64(-1);
    gl.tolerance = tolerance;
    Manifold::from_mesh_gl64(&gl)
}

/// A unit cube moved out to x = 1000 and back: its epsilon scaled with the
/// far bounding box (kPrecision * 1001) and a transform never lowers it, so it
/// carries an epsilon well above its own automatic 1e-12.
fn travelled_cube() -> Manifold {
    unit_cube()
        .translate(Vec3::new(1000.0, 0.0, 0.0))
        .translate(Vec3::new(-1000.0, 0.0, 0.0))
}

/// The rule the overlapping exact boolean applies, which the disjoint paths
/// must agree with: max of the operands' epsilon and of their tolerance.
fn expected(parts: &[&Manifold]) -> (f64, f64) {
    let eps = parts.iter().map(|m| m.get_epsilon()).fold(-1.0, f64::max);
    let tol = parts.iter().map(|m| m.get_tolerance()).fold(-1.0, f64::max);
    (eps, tol)
}

fn assert_carried(m: &Manifold, parts: &[&Manifold], what: &str) {
    let (eps, tol) = expected(parts);
    assert_eq!(m.get_epsilon(), eps, "{what}: epsilon");
    assert_eq!(m.get_tolerance(), tol, "{what}: tolerance");
    assert_eq!(
        m.get_mesh_gl64(-1).tolerance,
        tol,
        "{what}: MeshGL64 tolerance"
    );
}

/// The overlapping exact boolean really does keep `max` — the reference the
/// disjoint paths are held to.
#[test]
fn overlapping_exact_union_keeps_max_tolerance_and_epsilon() {
    let a = imported_cube(1e-3);
    let b = travelled_cube().translate(Vec3::new(0.5, 0.5, 0.5));
    let u = a.boolean_with_engine(&b, OpType::Add, crate::types::BooleanEngine::Exact);
    assert_carried(&u, &[&a, &b], "overlapping");
}

/// Every disjoint-union entry point, for one input pair: the result keeps
/// the max epsilon and tolerance, and its shape equals that of the same
/// union of `plain` inputs (default tolerance) — the carry changes nothing
/// else. The run table orders runs by original ID, so each pair must be
/// created in the same order as its plain counterpart (`a` before `b`).
fn check_entry_points(
    a: &Manifold,
    b: &Manifold,
    plain_a: &Manifold,
    plain_b: &Manifold,
    what: &str,
) {
    let b = b.translate(Vec3::new(5.0, 0.0, 0.0));
    let plain_b = plain_b.translate(Vec3::new(5.0, 0.0, 0.0));
    let parts = [a, &b];
    type Entry = fn(&Manifold, &Manifold) -> Manifold;
    let entries: [(&str, Entry); 4] = [
        ("compose", |a, b| Manifold::compose(&[a.clone(), b.clone()])),
        ("exact union", |a, b| {
            a.boolean_with_engine(b, OpType::Add, crate::types::BooleanEngine::Exact)
        }),
        ("batch_boolean", |a, b| {
            Manifold::batch_boolean(&[a.clone(), b.clone()], OpType::Add)
        }),
        ("csg batch_union", |a, b| {
            let leaf = |m: &Manifold| CsgNode::leaf(m.as_impl().clone());
            Manifold::from_impl(CsgNode::op_n(OpType::Add, vec![leaf(a), leaf(b)]).evaluate())
        }),
    ];
    for (name, entry) in entries {
        let r = entry(a, &b);
        assert_eq!(r.status(), crate::types::Error::NoError, "{what} {name}");
        assert_eq!(r.volume(), 2.0, "{what} {name}: volume");
        assert_carried(&r, &parts, &format!("{what} {name}"));
        let p = entry(plain_a, &plain_b);
        assert_eq!(
            shape(&r.get_mesh_gl64(-1)),
            shape(&p.get_mesh_gl64(-1)),
            "{what} {name}: geometry and runs"
        );
    }
}

#[test]
fn disjoint_union_keeps_set_tolerance() {
    let a = unit_cube().set_tolerance(1e-3);
    check_entry_points(&a, &unit_cube(), &unit_cube(), &unit_cube(), "a 1e-3");
}

#[test]
fn disjoint_union_keeps_the_larger_tolerance_either_side() {
    for (ta, tb) in [(0.0, 1e-3), (1e-4, 1e-3), (1e-3, 1e-4)] {
        let with = |t: f64| {
            let c = unit_cube();
            if t > 0.0 {
                c.set_tolerance(t)
            } else {
                c
            }
        };
        let (a, b) = (with(ta), with(tb));
        let (plain_a, plain_b) = (unit_cube(), unit_cube());
        check_entry_points(&a, &b, &plain_a, &plain_b, &format!("a {ta}, b {tb}"));
    }
}

#[test]
fn disjoint_union_keeps_imported_tolerance_and_epsilon() {
    let a = imported_cube(1e-3);
    assert_eq!(a.get_epsilon(), 1e-3);
    let b = unit_cube();
    let (plain_a, plain_b) = (imported_cube(0.0), unit_cube());
    check_entry_points(&a, &b, &plain_a, &plain_b, "a imported 1e-3");
    let a = unit_cube();
    let b = imported_cube(1e-3);
    let (plain_a, plain_b) = (unit_cube(), imported_cube(0.0));
    check_entry_points(&a, &b, &plain_a, &plain_b, "b imported 1e-3");
}

#[test]
fn disjoint_union_keeps_a_transformed_inputs_epsilon() {
    let a = travelled_cube();
    assert!(a.get_epsilon() > 1e-9, "{}", a.get_epsilon());
    let b = unit_cube().set_tolerance(1e-3);
    check_entry_points(&a, &b, &a, &unit_cube(), "travelled + 1e-3");
}

/// The CSG tree's leaves carry lazy transforms; batch_union composes their
/// materialized impls, whose epsilon is what the overlapping boolean of the
/// same leaves would see. A leaf scaled by 1000 carries 1000x its epsilon.
#[test]
fn csg_batch_union_keeps_scaled_leaf_epsilon() {
    let cube = unit_cube().set_tolerance(1e-3);
    let mut scale = Mat3x4::identity();
    scale.x.x = 1000.0;
    scale.y.y = 1000.0;
    scale.z.z = 1000.0;
    let scaled = CsgLeafNode::with_transform(cube.as_impl().clone(), scale);
    let mut far = Mat3x4::identity();
    far.w.x = 5000.0;
    let moved = CsgLeafNode::with_transform(unit_cube().as_impl().clone(), far);
    let (sa, sb) = (
        Manifold::from_impl(scaled.get_impl()),
        Manifold::from_impl(moved.get_impl()),
    );
    let r = Manifold::from_impl(
        CsgNode::op_n(
            OpType::Add,
            vec![CsgNode::leaf_node(scaled), CsgNode::leaf_node(moved)],
        )
        .evaluate(),
    );
    assert_eq!(sa.get_epsilon(), 1e-9);
    assert_carried(&r, &[&sa, &sb], "scaled leaf");
    assert!((r.volume() - (1e9 + 1.0)).abs() < 1e-6, "{}", r.volume());
}

/// The four disjoint-union entry points, by name.
fn union_entries() -> [(&'static str, fn(&Manifold, &Manifold) -> Manifold); 4] {
    [
        ("compose", |a, b| Manifold::compose(&[a.clone(), b.clone()])),
        ("exact union", |a, b| {
            a.boolean_with_engine(b, OpType::Add, crate::types::BooleanEngine::Exact)
        }),
        ("batch_boolean", |a, b| {
            Manifold::batch_boolean(&[a.clone(), b.clone()], OpType::Add)
        }),
        ("csg batch_union", |a, b| {
            let leaf = |m: &Manifold| CsgNode::leaf(m.as_impl().clone());
            Manifold::from_impl(CsgNode::op_n(OpType::Add, vec![leaf(a), leaf(b)]).evaluate())
        }),
    ]
}

/// A part's larger epsilon and tolerance are metadata for the composed
/// result, never a licence to coarsen the other parts: a fine sphere (edges
/// about 5e-4, under the cube's imported 1e-3) composed with that cube keeps
/// every vertex and triangle it has when composed with a default cube.
/// Collapsing it at the carried 1e-3 instead indexed the composed impl's
/// still-empty face normals and panicked.
#[test]
fn compose_does_not_coarsen_a_fine_part_at_another_parts_tolerance() {
    let coarse = imported_cube(1e-3);
    let plain = imported_cube(0.0);
    let sphere = Manifold::sphere(0.01, 128).translate(Vec3::new(5.0, 0.0, 0.0));
    for (name, entry) in union_entries() {
        let r = entry(&coarse, &sphere);
        let p = entry(&plain, &sphere);
        assert_eq!(r.status(), crate::types::Error::NoError, "{name}");
        assert_eq!(r.num_tri(), 12 + sphere.num_tri(), "{name}: triangles");
        assert_eq!(r.num_vert(), 8 + sphere.num_vert(), "{name}: vertices");
        assert_eq!(
            shape(&r.get_mesh_gl64(-1)),
            shape(&p.get_mesh_gl64(-1)),
            "{name}: geometry and runs"
        );
        assert_carried(&r, &[&coarse, &sphere], name);
    }
}

/// A part finer than the *combined* automatic epsilon — a tiny sphere at the
/// origin composed with a cube a million units away, whose box makes
/// `kPrecision * scale` 1e-6 against the sphere's ~5e-7 edges. This panicked
/// before the compose carried face normals (any default-tolerance input
/// could reach it; no tolerance is set here). The clean-up still collapses
/// the sphere at that epsilon, as C++ Compose's RemoveDegenerates does —
/// what is pinned here is only that it no longer panics and that the result
/// is a valid manifold keeping the cube.
#[test]
fn compose_with_a_part_finer_than_the_combined_epsilon_does_not_panic() {
    let tiny = Manifold::sphere(1e-5, 128);
    let far = unit_cube().translate(Vec3::new(1e6, 0.0, 0.0));
    for (name, entry) in union_entries() {
        let r = entry(&tiny, &far);
        assert_eq!(r.status(), crate::types::Error::NoError, "{name}");
        assert!(r.as_impl().is_manifold(), "{name}: manifold");
        assert!(r.num_tri() >= 12, "{name}: the cube survives");
        assert!(
            r.volume() >= 1.0 - 1e-9,
            "{name}: cube volume {}",
            r.volume()
        );
    }
}
