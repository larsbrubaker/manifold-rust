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

// Unit tests for the double-precision mesh path (src/meshgl64.rs), called the
// way a C caller does: raw pointers, manual destroy.

use crate::meshgl64::*;
use crate::tests::{cube_mesh, export, ffi_cube, read_array};
use crate::*;

/// The f32 cube fixture widened to the f64/u64 element types the 64-bit entry
/// point takes. Keeping one source of geometry makes the f32-vs-f64 comparison
/// tests meaningful.
fn cube_mesh64(origin: [f32; 3], size: f32) -> (Vec<f64>, Vec<u64>) {
    let (verts, tris) = cube_mesh(origin, size);
    (
        verts.into_iter().map(f64::from).collect(),
        tris.into_iter().map(u64::from).collect(),
    )
}

fn ffi_cube64(origin: [f32; 3], size: f32) -> *mut ManifoldRs {
    let (verts, tris) = cube_mesh64(origin, size);
    let handle = unsafe {
        manifold_rs_from_mesh64(verts.as_ptr(), verts.len(), tris.as_ptr(), tris.len(), 3)
    };
    assert!(!handle.is_null(), "cube64 import returned NULL");
    handle
}

struct ExportedMesh64 {
    num_prop: u32,
    vert_properties: Vec<f64>,
    tri_verts: Vec<u64>,
    run_index: Vec<u64>,
    run_original_id: Vec<u32>,
    face_id: Vec<u64>,
}

fn export64(m: *const ManifoldRs) -> ExportedMesh64 {
    let g = unsafe { manifold_rs_get_meshgl64(m) };
    assert!(!g.is_null(), "get_meshgl64 returned NULL");
    let exported = unsafe {
        ExportedMesh64 {
            num_prop: manifold_rs_meshgl64_num_prop(g),
            vert_properties: read_array(manifold_rs_meshgl64_vert_properties, g),
            tri_verts: read_array(manifold_rs_meshgl64_tri_verts, g),
            run_index: read_array(manifold_rs_meshgl64_run_index, g),
            run_original_id: read_array(manifold_rs_meshgl64_run_original_id, g),
            face_id: read_array(manifold_rs_meshgl64_face_id, g),
        }
    };
    unsafe { manifold_rs_meshgl64_destroy(g) };
    exported
}

#[test]
fn cube_round_trips_through_the_f64_path() {
    let cube = unsafe { manifold_rs_as_original(ffi_cube64([0.0, 0.0, 0.0], 1.0)) };
    assert_eq!(unsafe { manifold_rs_status(cube) }, 0);

    let mesh = export64(cube);
    assert_eq!(mesh.num_prop, 3);
    assert_eq!(mesh.tri_verts.len(), 36, "12 tris x 3 indices");
    assert_eq!(mesh.vert_properties.len() % 3, 0);
    assert_eq!(mesh.face_id.len(), 12);
    assert_eq!(mesh.run_original_id.len(), 1);
    assert_eq!(mesh.run_index.len(), 2, "one run plus the end sentinel");

    unsafe { manifold_rs_destroy(cube) };
}

#[test]
fn f64_and_f32_paths_describe_the_same_solid() {
    // The cube's coordinates are exactly representable in f32, so even though
    // the f64 path is lossless (no narrowing anywhere), both paths must
    // produce bit-identical geometry — not merely to some tolerance.
    let a64 = ffi_cube64([0.0, 0.0, 0.0], 1.0);
    let a32 = ffi_cube([0.0, 0.0, 0.0], 1.0);

    let m64 = export64(a64);
    let m32 = export(a32);

    assert_eq!(m64.num_prop, m32.num_prop);
    assert_eq!(
        m64.tri_verts,
        m32.tri_verts
            .iter()
            .map(|&v| u64::from(v))
            .collect::<Vec<u64>>()
    );
    assert_eq!(
        m64.face_id,
        m32.face_id
            .iter()
            .map(|&v| u64::from(v))
            .collect::<Vec<u64>>()
    );
    assert_eq!(
        m64.vert_properties,
        m32.vert_properties
            .iter()
            .map(|&v| f64::from(v))
            .collect::<Vec<f64>>()
    );

    unsafe {
        manifold_rs_destroy(a64);
        manifold_rs_destroy(a32);
    }
}

#[test]
fn f64_coordinates_survive_the_ffi_boundary_losslessly() {
    // A coordinate with more significand bits than f32 holds must come back
    // bit-identical through from_mesh64 → export. This is the FFI-level
    // tripwire for any regression back to the old narrow-through-f32 import.
    let third = 1.0_f64 / 3.0;
    assert_ne!(
        third as f32 as f64, third,
        "premise: 1/3 is not f32-representable"
    );

    let (mut verts, tris) = cube_mesh64([0.0, 0.0, 0.0], 1.0);
    for v in verts.iter_mut() {
        *v *= third; // 0 and 1 scale exactly; the far corner becomes 1/3
    }
    let cube = unsafe {
        manifold_rs_from_mesh64(verts.as_ptr(), verts.len(), tris.as_ptr(), tris.len(), 3)
    };
    assert!(!cube.is_null());
    assert_eq!(unsafe { manifold_rs_status(cube) }, 0);

    let mesh = export64(cube);
    assert!(
        mesh.vert_properties.iter().any(|&v| v == third),
        "1/3 did not survive the f64 round trip bit-exactly"
    );

    unsafe { manifold_rs_destroy(cube) };
}

#[test]
fn f64_import_feeds_the_boolean_pipeline() {
    let a = unsafe { manifold_rs_as_original(ffi_cube64([0.0, 0.0, 0.0], 1.0)) };
    let b = unsafe { manifold_rs_as_original(ffi_cube64([0.5, 0.0, 0.0], 1.0)) };
    let inputs = [a as *const ManifoldRs, b as *const ManifoldRs];

    let result = unsafe { manifold_rs_batch_boolean(inputs.as_ptr(), inputs.len(), 0) };
    assert!(!result.is_null());
    assert_eq!(unsafe { manifold_rs_status(result) }, 0);

    let mesh = export64(result);
    assert!(!mesh.tri_verts.is_empty(), "union should not be empty");
    assert_eq!(mesh.run_original_id.len(), 2, "one run per source mesh");

    unsafe {
        manifold_rs_destroy(result);
        manifold_rs_destroy(a);
        manifold_rs_destroy(b);
    }
}

#[test]
fn error_status_manifold_exports_as_an_empty_f64_mesh() {
    // Same wound-wrong tetrahedron as the f32 suite: clears the size checks,
    // fails the topology check.
    let verts: Vec<f64> = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    let tris: Vec<u64> = vec![0, 1, 2, 0, 1, 3, 0, 3, 2, 1, 2, 3];
    let bad = unsafe {
        manifold_rs_from_mesh64(verts.as_ptr(), verts.len(), tris.as_ptr(), tris.len(), 3)
    };
    assert!(!bad.is_null(), "a validation failure still yields a handle");
    assert_eq!(
        unsafe { manifold_rs_status(bad) },
        2,
        "expected NotManifold"
    );

    let mesh = export64(bad);
    assert!(mesh.tri_verts.is_empty());
    assert!(mesh.vert_properties.is_empty());

    unsafe { manifold_rs_destroy(bad) };
}

#[test]
fn invalid_from_mesh64_arguments_return_null() {
    let (verts, tris) = cube_mesh64([0.0, 0.0, 0.0], 1.0);
    unsafe {
        // num_prop below the three position slots.
        assert!(
            manifold_rs_from_mesh64(verts.as_ptr(), verts.len(), tris.as_ptr(), tris.len(), 2)
                .is_null()
        );
        // vert_properties not a whole number of vertices.
        assert!(manifold_rs_from_mesh64(
            verts.as_ptr(),
            verts.len() - 1,
            tris.as_ptr(),
            tris.len(),
            3
        )
        .is_null());
        // tri_verts not a whole number of triangles.
        assert!(manifold_rs_from_mesh64(
            verts.as_ptr(),
            verts.len(),
            tris.as_ptr(),
            tris.len() - 1,
            3
        )
        .is_null());
        // Null array with a non-zero length.
        assert!(manifold_rs_from_mesh64(
            std::ptr::null(),
            verts.len(),
            tris.as_ptr(),
            tris.len(),
            3
        )
        .is_null());
        // A length no allocation could ever have — a negative int widened to
        // size_t. usize::MAX divides by both 3 and num_prop, so it reaches the
        // isize::MAX guard rather than the earlier checks.
        assert!(
            manifold_rs_from_mesh64(verts.as_ptr(), usize::MAX, tris.as_ptr(), tris.len(), 3)
                .is_null()
        );
        assert!(
            manifold_rs_from_mesh64(verts.as_ptr(), verts.len(), tris.as_ptr(), usize::MAX, 3)
                .is_null()
        );
    }
}

#[test]
fn tri_index_beyond_u32_is_rejected_rather_than_wrapped() {
    // The kernel indexes with 32 bits and the core's from_mesh_gl64 truncates
    // indices to u32, which wraps. 1 << 32 would become 0 — an in-range index
    // pointing at the wrong vertex, producing wrong geometry that reports
    // status 0. That has to be a NULL plus a message, not a silent success.
    let (verts, mut tris) = cube_mesh64([0.0, 0.0, 0.0], 1.0);
    tris[0] = 1u64 << 32;
    let handle = unsafe {
        manifold_rs_from_mesh64(verts.as_ptr(), verts.len(), tris.as_ptr(), tris.len(), 3)
    };
    assert!(
        handle.is_null(),
        "an index that wraps to 0 must be rejected"
    );

    let mut buf = [0u8; 256];
    let len = unsafe { manifold_rs_last_error(buf.as_mut_ptr(), buf.len()) };
    let message = std::str::from_utf8(&buf[..len]).unwrap();
    assert!(
        message.contains("exceeds u32::MAX"),
        "unexpected message: {message}"
    );

    // u32::MAX itself is representable, so it is not a NULL case — it is an
    // ordinary out-of-bounds vertex index, which the kernel reports as a
    // non-zero status on a real handle.
    tris[0] = u64::from(u32::MAX);
    let handle = unsafe {
        manifold_rs_from_mesh64(verts.as_ptr(), verts.len(), tris.as_ptr(), tris.len(), 3)
    };
    assert!(
        !handle.is_null(),
        "u32::MAX is representable, so not a NULL case"
    );
    assert_ne!(
        unsafe { manifold_rs_status(handle) },
        0,
        "an out-of-bounds vertex index should surface as a non-zero status"
    );
    unsafe { manifold_rs_destroy(handle) };
}

#[test]
fn null_handles_return_sentinels_on_the_f64_path() {
    unsafe {
        assert!(manifold_rs_get_meshgl64(std::ptr::null()).is_null());
        assert_eq!(manifold_rs_meshgl64_num_prop(std::ptr::null()), 0);
        manifold_rs_meshgl64_destroy(std::ptr::null_mut());

        let mut len = usize::MAX;
        assert!(manifold_rs_meshgl64_vert_properties(std::ptr::null(), &mut len).is_null());
        assert_eq!(len, 0);
        len = usize::MAX;
        assert!(manifold_rs_meshgl64_tri_verts(std::ptr::null(), &mut len).is_null());
        assert_eq!(len, 0);
        len = usize::MAX;
        assert!(manifold_rs_meshgl64_run_index(std::ptr::null(), &mut len).is_null());
        assert_eq!(len, 0);
        len = usize::MAX;
        assert!(manifold_rs_meshgl64_run_original_id(std::ptr::null(), &mut len).is_null());
        assert_eq!(len, 0);
        len = usize::MAX;
        assert!(manifold_rs_meshgl64_face_id(std::ptr::null(), &mut len).is_null());
        assert_eq!(len, 0);

        // out_len itself may be NULL, for a null and for a live handle.
        assert!(manifold_rs_meshgl64_tri_verts(std::ptr::null(), std::ptr::null_mut()).is_null());
        let cube = ffi_cube64([0.0, 0.0, 0.0], 1.0);
        let g = manifold_rs_get_meshgl64(cube);
        assert!(!manifold_rs_meshgl64_tri_verts(g, std::ptr::null_mut()).is_null());
        manifold_rs_meshgl64_destroy(g);
        manifold_rs_destroy(cube);
    }
}
