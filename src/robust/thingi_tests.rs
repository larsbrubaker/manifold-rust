// robust/thingi_tests.rs — Regression tests against real-world Thingi10K
// meshes (checked into src/robust/testdata/). These reproduce failures found
// through the WASM demo's Boolean Gallery: both operands import cleanly via
// `Manifold::from_mesh_gl_robust`, yet a robust boolean of the pair returned
// `Error::NotClosed`. The fixtures follow the demo's exact pipeline: binary
// STL -> f32 triangle soup -> normalize (center at bbox center, max side =
// 2.0, f32 storage with f64 arithmetic, as in JS) -> MeshGL::merge ->
// robust import.

use crate::linalg::Vec3;
use crate::manifold::Manifold;
use crate::types::{BooleanEngine, Error, MeshGL};

/// Parse a binary STL (80-byte header, u32 face count, 50-byte records of
/// 12 f32 = normal + 3 vertices, then u16 attribute count) into a raw f32
/// position soup, one vertex per corner.
fn parse_binary_stl(data: &[u8]) -> Vec<f32> {
    assert!(data.len() >= 84, "STL too short");
    let n_faces = u32::from_le_bytes(data[80..84].try_into().unwrap()) as usize;
    assert!(
        84 + n_faces * 50 <= data.len(),
        "STL truncated: header declares {n_faces} faces"
    );
    let mut positions = Vec::with_capacity(n_faces * 9);
    for f in 0..n_faces {
        let rec = 84 + f * 50;
        // Skip the 12-byte facet normal; read the 3 corner vertices.
        for i in 0..9 {
            let off = rec + 12 + i * 4;
            positions.push(f32::from_le_bytes(data[off..off + 4].try_into().unwrap()));
        }
    }
    positions
}

/// Normalize exactly as the demo's `fetchMesh` does: bbox center to origin,
/// uniform scale so the longest bbox side is 2.0. JS stores positions in a
/// Float32Array but computes in doubles, so: read f32 -> widen to f64 ->
/// (x - c) * s in f64 -> truncate to f32 on store.
fn normalize_positions(positions: &mut [f32]) {
    let n_verts = positions.len() / 3;
    let (mut min, mut max) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for i in 0..n_verts {
        for k in 0..3 {
            let x = positions[i * 3 + k] as f64;
            if x < min[k] {
                min[k] = x;
            }
            if x > max[k] {
                max[k] = x;
            }
        }
    }
    let center = [
        (min[0] + max[0]) / 2.0,
        (min[1] + max[1]) / 2.0,
        (min[2] + max[2]) / 2.0,
    ];
    let max_side = (max[0] - min[0]).max(max[1] - min[1]).max(max[2] - min[2]);
    let scale = if max_side > 0.0 { 2.0 / max_side } else { 1.0 };
    for i in 0..n_verts {
        for k in 0..3 {
            positions[i * 3 + k] = ((positions[i * 3 + k] as f64 - center[k]) * scale) as f32;
        }
    }
}

/// The demo's import path: soup MeshGL (num_prop = 3, sequential indices),
/// weld with `MeshGL::merge`, import through the robust entry point.
fn import_stl_like_demo(stl: &[u8]) -> Manifold {
    let mut positions = parse_binary_stl(stl);
    normalize_positions(&mut positions);
    let n_verts = (positions.len() / 3) as u32;
    let mut mesh = MeshGL::default();
    mesh.num_prop = 3;
    mesh.vert_properties = positions;
    mesh.tri_verts = (0..n_verts).collect();
    mesh.merge();
    Manifold::from_mesh_gl_robust(&mesh)
}

const TENTACLE_939888: &[u8] = include_bytes!("testdata/939888.stl");
const PICKAXE_93557: &[u8] = include_bytes!("testdata/93557.stl");

/// Thingi10K #939888 union #93557 (demo repro): both operands import as
/// closed manifolds, so the robust union must produce a valid result — the
/// engine returned NotClosed.
#[test]
fn thingi_939888_union_93557_is_closed() {
    let a = import_stl_like_demo(TENTACLE_939888);
    assert_eq!(a.status(), Error::NoError, "operand A import");
    assert!(!a.as_impl().is_soup, "operand A should weld to a manifold");

    let b = import_stl_like_demo(PICKAXE_93557)
        .rotate(356.0, 140.0, 322.0)
        .translate(Vec3::new(0.3, 0.0, 0.0));
    assert_eq!(b.status(), Error::NoError, "operand B import");
    assert!(!b.as_impl().is_soup, "operand B should weld to a manifold");

    let result = a.union_with_engine(&b, BooleanEngine::Robust);
    assert_eq!(result.status(), Error::NoError, "robust union status");
    assert!(!result.is_empty(), "robust union should be non-empty");
    // The union of two closed solids contains at least each operand's
    // volume-maximum and must exceed either operand alone minus overlap;
    // just sanity-check positivity here — status is the regression.
    assert!(result.volume() > 0.0, "union volume must be positive");
}
