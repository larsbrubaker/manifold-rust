// robust/engine_tests.rs — End-to-end tests for the robust boolean engine
// and its selection plumbing: primitive booleans on the Robust engine
// checked against exact-engine volumes/areas, BooleanConfig default
// handling, Auto dispatch, fast paths, cancellation, and chained robust
// operations.

use crate::cancel::CancelToken;
use crate::linalg::Vec3;
use crate::manifold::Manifold;
use crate::types::{BooleanConfig, BooleanEngine, Error, OpType};

// STL fixture loading, shared with any other robust test module that wants
// it; declared here so it needs no entry in `robust/mod.rs`.
#[path = "stl_fixtures.rs"]
mod stl_fixtures;
use stl_fixtures::import_stl_like_demo;

fn v(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z)
}

fn assert_close(a: f64, b: f64, tol: f64, what: &str) {
    assert!(
        (a - b).abs() <= tol * b.abs().max(1.0),
        "{what}: {a} vs {b}"
    );
}

/// Both engines on the same manifold inputs must agree on volume and area
/// to near-f64 precision, and the robust result must be a true manifold.
fn check_ops_match(a: &Manifold, b: &Manifold, what: &str) {
    for op in [OpType::Add, OpType::Subtract, OpType::Intersect] {
        let exact = a.boolean_with_engine(b, op, BooleanEngine::Exact);
        let robust = a.boolean_with_engine(b, op, BooleanEngine::Robust);
        assert_eq!(robust.status(), Error::NoError, "{what} {op:?}");
        assert!(!robust.as_impl().is_soup, "{what} {op:?}: output must be manifold");
        assert_close(robust.volume(), exact.volume(), 1e-9, &format!("{what} {op:?} volume"));
        assert_close(
            robust.surface_area(),
            exact.surface_area(),
            1e-9,
            &format!("{what} {op:?} area"),
        );
        assert_eq!(robust.genus(), exact.genus(), "{what} {op:?} genus");
    }
}

#[test]
fn cube_cube_overlap_matches_exact() {
    let a = Manifold::cube(v(2.0, 2.0, 2.0), false);
    let b = a.translate(v(1.0, 1.0, 1.0));
    check_ops_match(&a, &b, "cube/cube");
}

#[test]
fn cube_sphere_matches_exact() {
    // Modest tessellation: debug-build BigRational time grows quickly with
    // intersection count; battery_extended covers denser meshes in release.
    let a = Manifold::cube(v(2.0, 2.0, 2.0), true);
    let b = Manifold::sphere(1.3, 12);
    check_ops_match(&a, &b, "cube/sphere");
}

#[test]
fn sphere_cylinder_matches_exact() {
    let a = Manifold::sphere(1.0, 12);
    let b = Manifold::cylinder(3.0, 0.6, 0.6, 8).translate(v(0.0, 0.0, -1.5));
    check_ops_match(&a, &b, "sphere/cylinder");
}

#[test]
fn tetra_tetra_matches_exact() {
    let a = Manifold::tetrahedron();
    let b = Manifold::tetrahedron().rotate(0.0, 0.0, 45.0).translate(v(0.2, 0.1, 0.3));
    check_ops_match(&a, &b, "tetra/tetra");
}

#[test]
fn disjoint_and_empty_fast_paths() {
    let a = Manifold::cube(v(1.0, 1.0, 1.0), false);
    let b = a.translate(v(5.0, 0.0, 0.0));
    let u = a.union_with_engine(&b, BooleanEngine::Robust);
    assert_close(u.volume(), 2.0, 1e-12, "disjoint union volume");
    assert!(a
        .intersection_with_engine(&b, BooleanEngine::Robust)
        .is_empty());
    assert_close(
        a.difference_with_engine(&b, BooleanEngine::Robust).volume(),
        1.0,
        1e-12,
        "disjoint difference",
    );
    let empty = Manifold::empty();
    assert_close(
        a.union_with_engine(&empty, BooleanEngine::Robust).volume(),
        1.0,
        1e-12,
        "union with empty",
    );
    assert!(a
        .intersection_with_engine(&empty, BooleanEngine::Robust)
        .is_empty());
}

#[test]
fn global_default_engine_config() {
    assert_eq!(BooleanConfig::default_engine(), BooleanEngine::Exact);
    BooleanConfig::set_default_engine(BooleanEngine::Auto);
    assert_eq!(BooleanConfig::default_engine(), BooleanEngine::Auto);
    // Plain boolean on manifold inputs under Auto → exact path → identical
    // to the explicit Exact call.
    let a = Manifold::cube(v(2.0, 2.0, 2.0), false);
    let b = a.translate(v(1.0, 0.0, 0.0));
    let via_auto = a.union(&b);
    let via_exact = a.union_with_engine(&b, BooleanEngine::Exact);
    assert_eq!(via_auto.num_vert(), via_exact.num_vert());
    assert_eq!(via_auto.num_tri(), via_exact.num_tri());
    assert_eq!(via_auto.volume(), via_exact.volume());
    BooleanConfig::reset_to_defaults();
    assert_eq!(BooleanConfig::default_engine(), BooleanEngine::Exact);
}

#[test]
fn auto_dispatches_soup_operands_to_robust() {
    // Non-manifold soup (edge-sharing cubes) + a manifold cutter.
    let mut tris = Vec::new();
    let cube = |lo: [f64; 3], hi: [f64; 3]| {
        Manifold::cube(v(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]), false)
            .translate(v(lo[0], lo[1], lo[2]))
    };
    let c1 = cube([0.0; 3], [2.0; 3]);
    let c2 = cube([2.0, 2.0, 0.0], [4.0, 4.0, 2.0]);
    for m in [&c1, &c2] {
        let gl = m.get_mesh_gl64(-1);
        for t in 0..gl.num_tri() {
            let tv = gl.get_tri_verts(t);
            let p = |i: u64| {
                let o = i as usize * 3;
                v(
                    gl.vert_properties[o],
                    gl.vert_properties[o + 1],
                    gl.vert_properties[o + 2],
                )
            };
            tris.push([p(tv[0]), p(tv[1]), p(tv[2])]);
        }
    }
    let mut mesh = crate::types::MeshGL64::default();
    mesh.num_prop = 3;
    for t in &tris {
        for p in t {
            mesh.vert_properties.extend([p.x, p.y, p.z]);
        }
    }
    mesh.tri_verts = (0..3 * tris.len() as u64).collect();
    let soup = Manifold::from_mesh_gl64_robust(&mesh);
    assert!(soup.as_impl().is_soup);
    assert_close(soup.volume(), 16.0, 1e-12, "soup volume");

    let cutter = Manifold::cube(v(1.0, 1.0, 1.0), false).translate(v(0.5, 0.5, 0.5));
    // Exact refuses.
    assert_eq!(
        soup.boolean_with_engine(&cutter, OpType::Subtract, BooleanEngine::Exact)
            .status(),
        Error::NotManifold
    );
    // Auto falls through to robust and produces the right volume.
    let diff = soup.boolean_with_engine(&cutter, OpType::Subtract, BooleanEngine::Auto);
    assert_eq!(diff.status(), Error::NoError);
    assert_close(diff.volume(), 15.0, 1e-9, "auto soup difference volume");
    // The shared edge persists in the output; whether it re-imports as a
    // soup or as an index-paired halfedge mesh depends on how the strict
    // pairing resolves the 4-halfedge fan — both are valid. What matters is
    // that chained booleans keep working either way (below).

    // Chained robust op on the soup result.
    let cutter2 = Manifold::cube(v(1.0, 1.0, 1.0), false).translate(v(2.5, 2.5, 0.5));
    let diff2 = diff.boolean_with_engine(&cutter2, OpType::Subtract, BooleanEngine::Auto);
    assert_eq!(diff2.status(), Error::NoError);
    assert_close(diff2.volume(), 14.0, 1e-9, "chained difference volume");
}

/// Thingi10K #92068 welds into a topologically manifold mesh, but its shells
/// are triple-wound duplicates of one another, so the exact engine
/// mis-integrates the union: it reports 1.7699 where the Monte-Carlo ground
/// truth is 1.5333 +/- 0.0111 and the robust engine gives 1.5259. `Auto`
/// must therefore pick Robust on the strength of the self-intersection test
/// alone — neither operand is a soup.
#[test]
fn auto_dispatches_self_intersecting_operands_to_robust() {
    let a = import_stl_like_demo(include_bytes!("testdata/92068.stl"));
    let b = import_stl_like_demo(include_bytes!("testdata/39926.stl"))
        .translate(v(0.3, 0.0, 0.0));
    assert_eq!(a.status(), Error::NoError, "operand A import");
    assert_eq!(b.status(), Error::NoError, "operand B import");
    assert!(!a.as_impl().is_soup, "operand A welds to a manifold");
    assert!(!b.as_impl().is_soup, "operand B welds to a manifold");
    assert!(
        a.has_self_intersections(),
        "92068's shells are coincident duplicates"
    );

    let auto = a.union_with_engine(&b, BooleanEngine::Auto);
    let robust = a.union_with_engine(&b, BooleanEngine::Robust);
    assert_eq!(auto.status(), Error::NoError, "auto union status");
    assert_eq!(
        auto.volume(),
        robust.volume(),
        "Auto must resolve to the robust engine"
    );
}

#[test]
fn robust_cancellation() {
    let a = Manifold::sphere(1.0, 16);
    let b = Manifold::sphere(1.0, 16).translate(v(0.5, 0.0, 0.0));
    let token = CancelToken::new();
    token.cancel();
    let r = a.boolean_with_engine_and_token(&b, OpType::Add, BooleanEngine::Robust, Some(&token));
    assert_eq!(r.status(), Error::Cancelled);
}

#[test]
fn coincident_cubes_union_and_subtract() {
    // Identical operands: union == either, difference == empty — classic
    // failure cases for tolerance-based engines, exact here.
    let a = Manifold::cube(v(2.0, 2.0, 2.0), false);
    let u = a.union_with_engine(&a.clone(), BooleanEngine::Robust);
    assert_close(u.volume(), 8.0, 1e-12, "self union volume");
    assert_close(u.surface_area(), 24.0, 1e-12, "self union area");
    let d = a.difference_with_engine(&a.clone(), BooleanEngine::Robust);
    assert!(d.is_empty() || d.volume().abs() < 1e-12, "self difference must vanish");
    let i = a.intersection_with_engine(&a.clone(), BooleanEngine::Robust);
    assert_close(i.volume(), 8.0, 1e-12, "self intersection volume");
}

#[test]
fn face_touching_union_merges_cleanly() {
    // Stacked cubes sharing a full face: the union must dissolve the shared
    // wall (volume 16, area of a 2x2x4 box).
    let a = Manifold::cube(v(2.0, 2.0, 2.0), false);
    let b = a.translate(v(0.0, 0.0, 2.0));
    let u = a.union_with_engine(&b, BooleanEngine::Robust);
    assert_eq!(u.status(), Error::NoError);
    assert_close(u.volume(), 16.0, 1e-12, "stacked union volume");
    assert_close(u.surface_area(), 40.0, 1e-12, "stacked union area");
}
