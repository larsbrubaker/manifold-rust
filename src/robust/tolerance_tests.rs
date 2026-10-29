// robust/tolerance_tests.rs — The robust engine keeps its operands'
// tolerance and epsilon.
//
// The exact engine (boolean_result_assemble.rs, C++ boolean_result.cpp:832)
// gives a boolean result `max` of its operands' epsilon and tolerance. The
// robust engine re-imports its assembled result through a MeshGL64 that
// carries no tolerance, so without carrying them across explicitly a user's
// `set_tolerance` (or an imported `MeshGL.tolerance`) silently resets to the
// automatic `kPrecision * bbox.scale()` after any robust boolean — and
// everything downstream (`simplify`, further booleans, export) loses it.
// These tests pin the robust result to the exact engine's values on the same
// inputs, through every robust entry: the pipeline, `Auto` routed to robust,
// soup operands, the disjoint fast paths, compose, and chained operations.
// The last two pin the other half of the contract: the carried values are
// metadata, so the robust geometry is bit-identical whatever the operands'
// tolerance (the header of robust/assemble.rs says why, with the measurements).

use crate::linalg::Vec3;
use crate::manifold::Manifold;
use crate::types::{BooleanEngine, Error, MeshGL, OpType, WindingRule};

const OPS: [OpType; 3] = [OpType::Add, OpType::Subtract, OpType::Intersect];

fn v(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z)
}

/// Overlapping operand pair: `a` with tolerance `tol_a`, `b` (a translated
/// cube) with `tol_b`; a tolerance of 0 leaves that operand at its default.
fn pair(tol_a: f64, tol_b: f64) -> (Manifold, Manifold) {
    let a = Manifold::cube(v(2.0, 2.0, 2.0), false);
    let b = Manifold::cube(v(2.0, 2.0, 2.0), false).translate(v(1.0, 1.0, 1.0));
    let a = if tol_a > 0.0 {
        a.set_tolerance(tol_a)
    } else {
        a
    };
    let b = if tol_b > 0.0 {
        b.set_tolerance(tol_b)
    } else {
        b
    };
    (a, b)
}

/// Robust result's tolerance and epsilon must equal the exact engine's on
/// the same overlapping inputs, and the f64 export must report it.
fn assert_matches_exact(a: &Manifold, b: &Manifold, what: &str) {
    for op in OPS {
        let exact = a.boolean_with_engine(b, op, BooleanEngine::Exact);
        let robust = a.boolean_with_engine(b, op, BooleanEngine::Robust);
        assert_eq!(robust.status(), Error::NoError, "{what} {op:?}");
        assert!(
            !robust.is_empty(),
            "{what} {op:?}: result must not be empty"
        );
        assert_eq!(
            robust.get_tolerance(),
            exact.get_tolerance(),
            "{what} {op:?}: tolerance"
        );
        assert_eq!(
            robust.get_epsilon(),
            exact.get_epsilon(),
            "{what} {op:?}: epsilon"
        );
        assert_eq!(
            robust.get_mesh_gl64(-1).tolerance,
            exact.get_mesh_gl64(-1).tolerance,
            "{what} {op:?}: exported MeshGL64 tolerance"
        );
    }
}

#[test]
fn robust_boolean_keeps_set_tolerance() {
    let (a, b) = pair(1e-3, 1e-3);
    assert_matches_exact(&a, &b, "both 1e-3");
    let r = a.boolean_with_engine(&b, OpType::Add, BooleanEngine::Robust);
    assert_eq!(r.get_tolerance(), 1e-3);
}

#[test]
fn robust_boolean_keeps_the_larger_operand_tolerance() {
    let (a, b) = pair(1e-3, 0.0);
    assert_matches_exact(&a, &b, "a 1e-3, b default");
    let (a, b) = pair(0.0, 1e-3);
    assert_matches_exact(&a, &b, "a default, b 1e-3");
    let r = a.boolean_with_engine(&b, OpType::Subtract, BooleanEngine::Robust);
    assert_eq!(r.get_tolerance(), 1e-3);
}

/// An imported `MeshGL.tolerance` raises epsilon too (C++ `SetEpsilon`), so
/// this pins the epsilon carry as well as the tolerance one.
#[test]
fn robust_boolean_keeps_imported_tolerance_and_epsilon() {
    let mut gl = Manifold::cube(v(2.0, 2.0, 2.0), false).get_mesh_gl64(-1);
    gl.tolerance = 1e-3;
    let a = Manifold::from_mesh_gl64(&gl);
    assert_eq!(a.get_epsilon(), 1e-3);
    let b = Manifold::cube(v(2.0, 2.0, 2.0), false).translate(v(1.0, 1.0, 1.0));
    assert_matches_exact(&a, &b, "imported 1e-3");
}

#[test]
fn auto_routed_to_robust_keeps_tolerance() {
    let (a, b) = pair(1e-3, 0.0);
    for op in OPS {
        // `Nonzero` forces `Auto` onto the robust engine.
        let r = a.boolean_with_engine_and_rule(&b, op, BooleanEngine::Auto, WindingRule::Nonzero);
        let exact = a.boolean_with_engine(&b, op, BooleanEngine::Exact);
        assert_eq!(r.get_tolerance(), exact.get_tolerance(), "{op:?}");
        assert_eq!(r.get_epsilon(), exact.get_epsilon(), "{op:?}");
    }
}

#[test]
fn chained_robust_operations_keep_tolerance() {
    let (a, b) = pair(1e-3, 0.0);
    let r = a.boolean_with_engine(&b, OpType::Add, BooleanEngine::Robust);
    // `simplify(0)` means "at the current tolerance" — it must not drop it.
    assert_eq!(r.simplify(0.0).get_tolerance(), 1e-3);
    let c = Manifold::cube(v(1.0, 1.0, 1.0), false).translate(v(2.5, 2.5, 2.5));
    let r2 = r.boolean_with_engine(&c, OpType::Subtract, BooleanEngine::Robust);
    assert_eq!(r2.status(), Error::NoError);
    assert_eq!(r2.get_tolerance(), 1e-3);
    let gl = r2.get_mesh_gl(-1);
    assert_eq!(gl.tolerance, 1e-3);
    // And through a float round trip.
    assert_eq!(Manifold::from_mesh_gl(&gl).get_tolerance(), 1e-3f32 as f64);
}

/// The bbox-disjoint fast paths keep the operands' larger tolerance too:
/// the difference (which returns A) as the exact engine does, the union
/// (a concatenation) as C++ `CsgLeafNode::Compose` does.
#[test]
fn robust_disjoint_fast_paths_keep_tolerance() {
    let a = Manifold::cube(v(1.0, 1.0, 1.0), false);
    let b = Manifold::cube(v(1.0, 1.0, 1.0), false)
        .translate(v(5.0, 0.0, 0.0))
        .set_tolerance(1e-3);
    let sub = a.boolean_with_engine(&b, OpType::Subtract, BooleanEngine::Robust);
    let exact_sub = a.boolean_with_engine(&b, OpType::Subtract, BooleanEngine::Exact);
    assert_eq!(sub.get_tolerance(), exact_sub.get_tolerance());
    assert_eq!(sub.get_tolerance(), 1e-3);
    let add = a.boolean_with_engine(&b, OpType::Add, BooleanEngine::Robust);
    assert_eq!(add.volume(), 2.0);
    assert_eq!(add.get_tolerance(), 1e-3);
    let add = b.boolean_with_engine(&a, OpType::Add, BooleanEngine::Robust);
    assert_eq!(add.get_tolerance(), 1e-3);
}

/// Cube [lo,hi]³ as a fully disconnected 12-triangle soup MeshGL.
fn soup_cube_gl(lo: f64, hi: f64, tolerance: f32) -> MeshGL {
    let c = |x: bool, y: bool, z: bool| {
        let p = |b: bool| if b { hi as f32 } else { lo as f32 };
        [p(x), p(y), p(z)]
    };
    let quads = [
        [
            c(false, false, false),
            c(false, true, false),
            c(true, true, false),
            c(true, false, false),
        ],
        [
            c(false, false, true),
            c(true, false, true),
            c(true, true, true),
            c(false, true, true),
        ],
        [
            c(false, false, false),
            c(true, false, false),
            c(true, false, true),
            c(false, false, true),
        ],
        [
            c(false, true, false),
            c(false, true, true),
            c(true, true, true),
            c(true, true, false),
        ],
        [
            c(false, false, false),
            c(false, false, true),
            c(false, true, true),
            c(false, true, false),
        ],
        [
            c(true, false, false),
            c(true, true, false),
            c(true, true, true),
            c(true, false, true),
        ],
    ];
    let mut mesh = MeshGL::default();
    mesh.num_prop = 3;
    for q in quads {
        for p in [q[0], q[1], q[2], q[0], q[2], q[3]] {
            mesh.vert_properties.extend(p);
        }
    }
    mesh.tri_verts = (0..36).collect();
    mesh.tolerance = tolerance;
    mesh
}

#[test]
fn soup_operand_tolerance_survives_robust_boolean() {
    let soup = Manifold::from_mesh_gl_robust(&soup_cube_gl(0.0, 2.0, 1e-3));
    assert!(soup.as_impl().is_soup);
    let soup_tol = soup.get_tolerance();
    assert_eq!(soup_tol, 1e-3f32 as f64);
    let b = Manifold::cube(v(2.0, 2.0, 2.0), false).translate(v(1.0, 1.0, 1.0));
    for op in OPS {
        for engine in [BooleanEngine::Robust, BooleanEngine::Auto] {
            let r = soup.boolean_with_engine(&b, op, engine);
            assert_eq!(r.status(), Error::NoError, "{op:?} {engine:?}");
            assert_eq!(r.get_tolerance(), soup_tol, "{op:?} {engine:?}");
            assert_eq!(r.get_epsilon(), soup.get_epsilon(), "{op:?} {engine:?}");
        }
    }
}

/// `Manifold::compose` with a soup part runs the robust concatenation
/// (`robust::assemble_all`); it keeps the parts' largest tolerance, as C++
/// `CsgLeafNode::Compose` does.
#[test]
fn soup_compose_keeps_tolerance() {
    let soup = Manifold::from_mesh_gl_robust(&soup_cube_gl(0.0, 1.0, 1e-3));
    let cube = Manifold::cube(v(1.0, 1.0, 1.0), false).translate(v(3.0, 0.0, 0.0));
    let composed = Manifold::compose(&[cube, soup.clone()]);
    assert_eq!(composed.status(), Error::NoError);
    assert_eq!(composed.get_tolerance(), soup.get_tolerance());
    assert_eq!(composed.get_epsilon(), soup.get_epsilon());
}

/// Geometry of a robust result, bit for bit.
fn geometry(m: &Manifold) -> (Vec<u64>, Vec<u64>) {
    let gl = m.get_mesh_gl64(-1);
    (
        gl.vert_properties.iter().map(|x| x.to_bits()).collect(),
        gl.tri_verts,
    )
}

/// The carried tolerance is metadata: the robust result's geometry is the
/// same as at the default tolerance, even where the exact engine's
/// tolerance collapse would coarsen it. Here `b`'s top corners sit 1e-4
/// above `a`'s top face, so at tolerance 1e-3 the exact engine collapses the
/// created vertices there (32 triangles for the union, and the difference's
/// volume moves from 6.0001 to 5.9999667); the robust result stays exact.
#[test]
fn carried_tolerance_does_not_coarsen_robust_geometry() {
    let b = Manifold::cube(v(2.0, 2.0, 2.0), false).translate(v(1.0, 1.0, 1e-4));
    let plain = Manifold::cube(v(2.0, 2.0, 2.0), false);
    let coarse = plain.set_tolerance(1e-3);
    for op in OPS {
        let r_plain = plain.boolean_with_engine(&b, op, BooleanEngine::Robust);
        let r_coarse = coarse.boolean_with_engine(&b, op, BooleanEngine::Robust);
        assert_eq!(geometry(&r_coarse), geometry(&r_plain), "{op:?}");
        assert_eq!(r_coarse.get_tolerance(), 1e-3, "{op:?}");
    }
    let sub = coarse.boolean_with_engine(&b, OpType::Subtract, BooleanEngine::Robust);
    assert_eq!(sub.volume(), 8.0 - (2.0 - 1e-4));
}

/// An imported tolerance raises epsilon as well; the robust result keeps
/// both but still simplifies at the automatic epsilon. Collapsing at the
/// carried 1e-3 instead broke this thin-shell difference (a sphere minus a
/// rotated copy shifted 5e-4): it came out self-intersecting with genus 14
/// instead of the clean genus-20 shell the automatic epsilon gives.
#[test]
fn imported_tolerance_does_not_coarsen_thin_shell_difference() {
    let a = Manifold::sphere(1.0, 48);
    let b = Manifold::sphere(1.0, 48)
        .rotate(13.0, 7.0, 3.0)
        .translate(v(0.0005, 0.0, 0.0));
    let mut gl = a.get_mesh_gl64(-1);
    gl.tolerance = 1e-3;
    let coarse = Manifold::from_mesh_gl64(&gl);
    assert_eq!(coarse.get_epsilon(), 1e-3);
    let plain = Manifold::from_mesh_gl64(&a.get_mesh_gl64(-1));
    let r_plain = plain.boolean_with_engine(&b, OpType::Subtract, BooleanEngine::Robust);
    let r_coarse = coarse.boolean_with_engine(&b, OpType::Subtract, BooleanEngine::Robust);
    assert_eq!(geometry(&r_coarse), geometry(&r_plain));
    assert!(!r_coarse.has_self_intersections());
    assert_eq!(r_coarse.genus(), 20);
    assert_eq!(r_coarse.get_epsilon(), 1e-3);
    assert_eq!(r_coarse.get_tolerance(), 1e-3);
}
