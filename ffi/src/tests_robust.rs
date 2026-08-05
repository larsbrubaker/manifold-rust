// Copyright 2026 Lars Brubaker
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// C-ABI tests for the robust (non-manifold) import and the boolean-engine
// selector: status codes 0/2/15, the engine set/get round trip, and an
// Auto-engine boolean on soup geometry through manifold_rs_batch_boolean.

use super::tests::cube_mesh;
use super::*;

/// The cube mesh re-expressed as fully disconnected triangle soup (three
/// duplicated verts per triangle) — geometrically identical, guaranteed to
/// fail strict pairing.
fn cube_soup(origin: [f32; 3], size: f32) -> (Vec<f32>, Vec<u32>) {
    let (verts, tris) = cube_mesh(origin, size);
    let mut soup_verts = Vec::with_capacity(9 * tris.len() / 3);
    for &i in &tris {
        let o = 3 * i as usize;
        soup_verts.extend_from_slice(&verts[o..o + 3]);
    }
    let soup_tris: Vec<u32> = (0..tris.len() as u32).collect();
    (soup_verts, soup_tris)
}

#[test]
fn robust_import_accepts_soup_strict_rejects() {
    let (verts, tris) = cube_soup([0.0; 3], 2.0);
    let strict = unsafe {
        manifold_rs_from_mesh(verts.as_ptr(), verts.len(), tris.as_ptr(), tris.len(), 3)
    };
    assert!(!strict.is_null());
    assert_eq!(unsafe { manifold_rs_status(strict) }, 2, "strict: NotManifold");
    unsafe { manifold_rs_destroy(strict) };

    let robust = unsafe {
        manifold_rs_from_mesh_robust(verts.as_ptr(), verts.len(), tris.as_ptr(), tris.len(), 3)
    };
    assert!(!robust.is_null());
    assert_eq!(unsafe { manifold_rs_status(robust) }, 0, "robust: NoError");
    unsafe { manifold_rs_destroy(robust) };
}

#[test]
fn robust_import_rejects_open_mesh_with_not_closed() {
    let (verts, mut tris) = cube_soup([0.0; 3], 2.0);
    tris.truncate(tris.len() - 3); // drop one triangle → open
    let handle = unsafe {
        manifold_rs_from_mesh_robust(verts.as_ptr(), verts.len(), tris.as_ptr(), tris.len(), 3)
    };
    assert!(!handle.is_null());
    assert_eq!(unsafe { manifold_rs_status(handle) }, 15, "NotClosed");
    unsafe { manifold_rs_destroy(handle) };
}

// The engine selector is process-global state, so everything touching it
// lives in ONE test — the parallel test runner would otherwise race the
// setter against tests reading the default.
#[test]
fn engine_selector_and_auto_boolean_on_soup() {
    assert_eq!(manifold_rs_get_boolean_engine(), 0, "default is Exact");
    assert_eq!(manifold_rs_set_boolean_engine(2), 0);
    assert_eq!(manifold_rs_get_boolean_engine(), 2);
    assert_eq!(manifold_rs_set_boolean_engine(7), -1, "unknown engine");
    assert_eq!(manifold_rs_get_boolean_engine(), 2, "unchanged after error");
    assert_eq!(manifold_rs_set_boolean_engine(0), 0);
    assert_eq!(manifold_rs_get_boolean_engine(), 0);

    auto_engine_boolean_on_soup_through_batch();
}

fn auto_engine_boolean_on_soup_through_batch() {
    // Edge-kissing cubes as one soup + a manifold cutter; under the Auto
    // engine the difference must succeed with the right volume-bearing mesh
    // (checked via triangle count > 0 and status 0 — volume isn't exposed
    // through this ABI).
    let (mut verts, _) = cube_soup([0.0; 3], 2.0);
    let (verts2, _) = cube_soup([2.0, 2.0, 0.0], 2.0);
    verts.extend(verts2);
    let tris: Vec<u32> = (0..verts.len() as u32 / 3).collect();
    let soup = unsafe {
        manifold_rs_from_mesh_robust(verts.as_ptr(), verts.len(), tris.as_ptr(), tris.len(), 3)
    };
    assert_eq!(unsafe { manifold_rs_status(soup) }, 0);

    let cutter = super::tests::ffi_cube([0.5, 0.5, 0.5], 1.0);

    assert_eq!(manifold_rs_set_boolean_engine(2), 0); // Auto
    let inputs = [soup as *const ManifoldRs, cutter as *const ManifoldRs];
    let diff = unsafe { manifold_rs_batch_boolean(inputs.as_ptr(), 2, 1) };
    assert_eq!(manifold_rs_set_boolean_engine(0), 0); // restore before asserts

    assert!(!diff.is_null());
    assert_eq!(unsafe { manifold_rs_status(diff) }, 0);
    let mesh = unsafe { manifold_rs_get_meshgl(diff) };
    assert!(!mesh.is_null());
    let mut len = 0usize;
    let tri_ptr = unsafe { manifold_rs_meshgl_tri_verts(mesh, &mut len) };
    assert!(!tri_ptr.is_null() && len > 0, "difference produced no triangles");
    unsafe {
        manifold_rs_meshgl_destroy(mesh);
        manifold_rs_destroy(diff);
        manifold_rs_destroy(cutter);
        manifold_rs_destroy(soup);
    }
}
