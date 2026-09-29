// Tests ported from C++ TEST(Manifold, Normals*) in manifold_test.cpp (#1718).
// Verify that CalculateNormals records world-frame normals on the Manifold and
// that GetMeshGL() auto-substitutes / round-trips them, including across
// transforms and Boolean subtraction (cavity inversion).

use super::*;
use crate::linalg::{dot, length, normalize, Vec3};
use crate::types::Error;

// C++ CalculateNormals() / CalculateNormals(0) default minSharpAngle is 52.5.
const DEFAULT_ANGLE: f64 = 52.5;

/// Count verts on the sphere surface (|pos| ~ radius) whose stored normal at
/// channel 3..5 aligns with `expected(pos)` (dot > 0.9). Returns (good, bad).
fn count_sphere_normal_alignment(
    gl: &crate::types::MeshGL,
    radius: f64,
    expected: impl Fn(Vec3) -> Vec3,
) -> (i32, i32) {
    let np = gl.num_prop as usize;
    let (mut good, mut bad) = (0, 0);
    let num_vert = gl.vert_properties.len() / np;
    for v in 0..num_vert {
        let pos = Vec3::new(
            gl.vert_properties[v * np] as f64,
            gl.vert_properties[v * np + 1] as f64,
            gl.vert_properties[v * np + 2] as f64,
        );
        if (length(pos) - radius).abs() > 0.1 {
            continue;
        }
        let n = Vec3::new(
            gl.vert_properties[v * np + 3] as f64,
            gl.vert_properties[v * np + 4] as f64,
            gl.vert_properties[v * np + 5] as f64,
        );
        if dot(n, expected(pos)) > 0.9 {
            good += 1;
        } else {
            bad += 1;
        }
    }
    (good, bad)
}

/// C++ TEST(Manifold, NormalsCavity) — inner-sphere normals from a Boolean diff
/// should point toward the origin (outward from the surrounding solid).
#[test]
fn test_cpp_normals_cavity() {
    let mesh = Manifold::sphere(10.0, 32)
        .difference(&Manifold::sphere(3.0, 32))
        .calculate_normals(0, DEFAULT_ANGLE)
        .get_mesh_gl(-1);
    assert!(mesh.num_prop >= 6, "numProp={}", mesh.num_prop);
    let (good, bad) = count_sphere_normal_alignment(&mesh, 3.0, |pos| normalize(pos * -1.0));
    assert!(good > 0, "NormalsCavity: good={}", good);
    assert_eq!(bad, 0, "NormalsCavity: bad={}", bad);
}

/// C++ TEST(Manifold, NormalsRotateBeforeCalc) — rotation before CalculateNormals:
/// SetNormals computes from already-rotated face normals and stores world-frame.
#[test]
fn test_cpp_normals_rotate_before_calc() {
    let mesh = Manifold::sphere(10.0, 32)
        .rotate(45.0, 0.0, 0.0)
        .calculate_normals(0, DEFAULT_ANGLE)
        .get_mesh_gl(-1);
    let (_, bad) = count_sphere_normal_alignment(&mesh, 10.0, normalize);
    assert_eq!(bad, 0, "NormalsRotateBeforeCalc: bad={}", bad);
}

/// C++ TEST(Manifold, NormalsRotateAfterCalc) — rotation *after* CalculateNormals:
/// Transform eager-transforms the stored slot 0..2 so it tracks the new orientation.
#[test]
fn test_cpp_normals_rotate_after_calc() {
    let mesh = Manifold::sphere(10.0, 32)
        .calculate_normals(0, DEFAULT_ANGLE)
        .rotate(45.0, 0.0, 0.0)
        .get_mesh_gl(-1);
    let (_, bad) = count_sphere_normal_alignment(&mesh, 10.0, normalize);
    assert_eq!(bad, 0, "NormalsRotateAfterCalc: bad={}", bad);
}

/// C++ TEST(Manifold, NormalsAutoSubstitute) — no-arg invocation defaults to
/// slot 0 and sets the per-run hasNormals bit on every output run.
#[test]
fn test_cpp_normals_auto_substitute() {
    let mesh = Manifold::sphere(10.0, 32)
        .calculate_normals(0, DEFAULT_ANGLE)
        .get_mesh_gl(-1);
    assert!(mesh.num_prop >= 6, "numProp={}", mesh.num_prop);
    assert!(!mesh.run_original_id.is_empty(), "no runs");
    assert!(mesh.has_normals(0), "run 0 should have hasNormals bit");
}

/// C++ TEST(Manifold, NormalsRoundTrip) — getMesh -> ofMesh -> getMesh preserves
/// the per-run flag, so the second getMesh still emits world-frame normals.
#[test]
fn test_cpp_normals_round_trip() {
    let round = Manifold::sphere(10.0, 32)
        .difference(&Manifold::sphere(3.0, 32))
        .calculate_normals(0, DEFAULT_ANGLE);
    let out1 = round.get_mesh_gl(-1);
    assert!(out1.has_normals(0), "out1 run 0 hasNormals");
    let out2 = Manifold::from_mesh_gl(&out1).get_mesh_gl(-1);
    assert!(out2.has_normals(0), "out2 run 0 hasNormals");
    let (good, bad) = count_sphere_normal_alignment(&out2, 3.0, |pos| normalize(pos * -1.0));
    assert!(good > 0, "NormalsRoundTrip: good={}", good);
    assert_eq!(bad, 0, "NormalsRoundTrip: bad={}", bad);
}

/// C++ TEST(Manifold, NormalsRefinePreserved) — Refine keeps the recording
/// (linearly-interpolated normals at new verts, but the flag survives).
#[test]
fn test_cpp_normals_refine_preserved() {
    let mesh = Manifold::sphere(10.0, 32)
        .calculate_normals(0, DEFAULT_ANGLE)
        .refine(2)
        .get_mesh_gl(-1);
    assert!(!mesh.run_original_id.is_empty(), "no runs");
    assert!(mesh.has_normals(0), "Refine should preserve hasNormals");
}

/// C++ TEST(Manifold, NormalsSmoothByNormalsNoArg) — no-arg SmoothByNormals reads
/// the recorded slot 0 and produces a valid manifold.
#[test]
fn test_cpp_normals_smooth_by_normals_no_arg() {
    let smoothed = Manifold::sphere(10.0, 32)
        .calculate_normals(0, DEFAULT_ANGLE)
        .smooth_by_normals(0);
    assert_eq!(smoothed.status(), Error::NoError);
}

/// C++ TEST(Manifold, NormalsNonStandardSlotNotRecorded) — CalculateNormals(3)
/// does NOT set the recording, since a non-standard slot can't be safely
/// auto-substituted on GetMeshGL(-1).
#[test]
fn test_cpp_normals_non_standard_slot_not_recorded() {
    let mesh = Manifold::sphere(10.0, 32)
        .calculate_normals(3, DEFAULT_ANGLE)
        .get_mesh_gl(-1);
    assert!(!mesh.run_original_id.is_empty(), "no runs");
    assert!(
        !mesh.has_normals(0),
        "non-standard slot must not record hasNormals"
    );
}

/// C++ TEST(Manifold, MirroredNormals) — upstream 422ab6fc (issue #1781). A
/// mirror flips triangle winding; the per-corner prop indices must follow their
/// start verts, or normals land on the wrong corners and GetMeshGL duplicates
/// vertices. Vert count must be unchanged and normals must stay outward.
#[test]
fn test_cpp_mirrored_normals() {
    let s = Manifold::sphere(1.0, 32).calculate_normals(0, 180.0);
    let mesh = s.get_mesh_gl(-1);
    let mirrored = s.mirror(Vec3::new(1.0, 0.0, 0.0)).get_mesh_gl(-1);
    assert_eq!(mesh.num_vert(), mirrored.num_vert());
    assert!(mirrored.has_normals(0));
    let (good, bad) = count_sphere_normal_alignment(&mirrored, 1.0, normalize);
    assert_eq!(bad, 0, "MirroredNormals: bad={}", bad);
    assert_eq!(good, mirrored.num_vert() as i32, "MirroredNormals: good={}", good);
}

// Property-adding ops after a mirror: with num_prop == 0 the flip must still
// keep prop_vert == start_vert, because set_properties, calculate_curvature and
// calculate_normals index existing props by prop_vert. Upstream 422ab6fc;
// see docs/CPP_DIVERGENCES.md entry 4.

#[test]
fn test_mirrored_cone_set_properties_shares_verts() {
    // radius_low == 0 builds the apex-bottom cone by mirroring internally.
    let cone = Manifold::cylinder(2.0, 0.0, 1.0, 16);
    let gl = cone.set_properties(1, |p, pos, _| p[0] = pos.z).get_mesh_gl(-1);
    assert_eq!(gl.num_vert(), 17);
}

#[test]
fn test_mirrored_cube_set_properties_shares_verts() {
    let cube = Manifold::cube(Vec3::new(1.0, 2.0, 3.0), true).mirror(Vec3::new(1.0, 0.0, 0.0));
    let gl = cube.set_properties(1, |p, pos, _| p[0] = pos.x).get_mesh_gl(-1);
    assert_eq!(gl.num_vert(), 8);
}

#[test]
fn test_mirrored_cone_calculate_normals_vert_count() {
    let cone = Manifold::cylinder(2.0, 0.0, 1.0, 16);
    let gl = cone.calculate_normals(0, 60.0).get_mesh_gl(-1);
    assert_eq!(gl.num_vert(), 33);
}
