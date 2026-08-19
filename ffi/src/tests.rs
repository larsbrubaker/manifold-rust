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

// Unit tests for the C ABI. They call the exported functions the way a C
// caller does — raw pointers, manual destroy — so the tests exercise the real
// entry points rather than the Rust API underneath them.

use super::*;

/// Axis-aligned cube as interleaved positions plus CCW-outward triangles.
/// Written out by hand so the tests do not depend on the core crate's own
/// cube constructor for their input.
pub(crate) fn cube_mesh(origin: [f32; 3], size: f32) -> (Vec<f32>, Vec<u32>) {
    let [x, y, z] = origin;
    let s = size;
    let verts = vec![
        x,
        y,
        z,
        x + s,
        y,
        z,
        x + s,
        y + s,
        z,
        x,
        y + s,
        z,
        x,
        y,
        z + s,
        x + s,
        y,
        z + s,
        x + s,
        y + s,
        z + s,
        x,
        y + s,
        z + s,
    ];
    let tris = vec![
        0, 2, 1, 0, 3, 2, // -Z
        4, 5, 6, 4, 6, 7, // +Z
        0, 1, 5, 0, 5, 4, // -Y
        1, 2, 6, 1, 6, 5, // +X
        2, 3, 7, 2, 7, 6, // +Y
        3, 0, 4, 3, 4, 7, // -X
    ];
    (verts, tris)
}

/// Build a cube through the FFI entry point. Panics on failure so a broken
/// import shows up as a test failure at the point of construction.
pub(crate) fn ffi_cube(origin: [f32; 3], size: f32) -> *mut ManifoldRs {
    let (verts, tris) = cube_mesh(origin, size);
    let handle =
        unsafe { manifold_rs_from_mesh(verts.as_ptr(), verts.len(), tris.as_ptr(), tris.len(), 3) };
    assert!(!handle.is_null(), "cube import returned NULL");
    handle
}

/// Copy an exported mesh's arrays out of the FFI so assertions can compare
/// them without juggling raw pointers.
pub(crate) struct ExportedMesh {
    pub num_prop: u32,
    pub vert_properties: Vec<f32>,
    pub tri_verts: Vec<u32>,
    pub run_index: Vec<u32>,
    pub run_original_id: Vec<u32>,
    pub face_id: Vec<u32>,
}

/// Copy one array accessor's data out. Generic over the handle type so the
/// f32 and f64 mesh test suites share it.
///
/// # Safety
/// `g` must be a live handle of the type `accessor` expects.
pub(crate) unsafe fn read_array<H, T: Copy>(
    accessor: unsafe extern "C" fn(*const H, *mut usize) -> *const T,
    g: *const H,
) -> Vec<T> {
    let mut len = usize::MAX;
    // SAFETY: caller guarantees g is live, and len is writable.
    let ptr = unsafe { accessor(g, &mut len) };
    assert!(!ptr.is_null(), "accessor returned NULL for a live handle");
    assert_ne!(len, usize::MAX, "accessor did not write out_len");
    if len == 0 {
        Vec::new()
    } else {
        // SAFETY: the accessor just reported len elements at ptr, and the data
        // outlives the copy because g is still alive.
        unsafe { std::slice::from_raw_parts(ptr, len) }.to_vec()
    }
}

pub(crate) fn export(m: *const ManifoldRs) -> ExportedMesh {
    let g = unsafe { manifold_rs_get_meshgl(m) };
    assert!(!g.is_null(), "get_meshgl returned NULL");

    let exported = unsafe {
        ExportedMesh {
            num_prop: manifold_rs_meshgl_num_prop(g),
            vert_properties: read_array(manifold_rs_meshgl_vert_properties, g),
            tri_verts: read_array(manifold_rs_meshgl_tri_verts, g),
            run_index: read_array(manifold_rs_meshgl_run_index, g),
            run_original_id: read_array(manifold_rs_meshgl_run_original_id, g),
            face_id: read_array(manifold_rs_meshgl_face_id, g),
        }
    };
    unsafe { manifold_rs_meshgl_destroy(g) };
    exported
}

#[test]
fn version_string_is_nul_terminated() {
    let ptr = manifold_rs_version();
    assert!(!ptr.is_null());
    let text = unsafe { std::ffi::CStr::from_ptr(ptr) }.to_str().unwrap();
    assert!(
        text.starts_with("manifold-ffi "),
        "unexpected version: {text}"
    );
    assert!(
        text.contains("manifold-rust "),
        "unexpected version: {text}"
    );
}

#[test]
fn cube_round_trips_through_the_ffi() {
    let cube = ffi_cube([0.0, 0.0, 0.0], 1.0);
    let original = unsafe { manifold_rs_as_original(cube) };
    assert!(!original.is_null());

    assert_eq!(unsafe { manifold_rs_status(original) }, 0);
    assert!(
        unsafe { manifold_rs_original_id(original) } >= 0,
        "as_original should assign a mesh ID"
    );

    let mesh = export(original);
    assert_eq!(mesh.num_prop, 3);
    assert_eq!(mesh.tri_verts.len(), 36, "12 tris x 3 indices");
    assert!(!mesh.vert_properties.is_empty());
    assert_eq!(mesh.vert_properties.len() % 3, 0);
    assert_eq!(mesh.face_id.len(), 12);

    unsafe {
        manifold_rs_destroy(original);
        manifold_rs_destroy(cube);
    }
}

#[test]
fn union_of_two_cubes_keeps_both_original_ids() {
    let a = unsafe { manifold_rs_as_original(ffi_cube([0.0, 0.0, 0.0], 1.0)) };
    let b = unsafe { manifold_rs_as_original(ffi_cube([0.5, 0.0, 0.0], 1.0)) };
    let id_a = unsafe { manifold_rs_original_id(a) } as u32;
    let id_b = unsafe { manifold_rs_original_id(b) } as u32;

    let inputs = [a as *const ManifoldRs, b as *const ManifoldRs];
    let result = unsafe { manifold_rs_batch_boolean(inputs.as_ptr(), inputs.len(), 0) };
    assert!(!result.is_null());
    assert_eq!(unsafe { manifold_rs_status(result) }, 0);

    let mesh = export(result);
    assert!(!mesh.tri_verts.is_empty(), "union should not be empty");
    assert_eq!(mesh.run_original_id.len(), 2, "one run per source mesh");
    assert_eq!(mesh.run_index.len(), 3, "run starts plus the end sentinel");

    let mut ids = mesh.run_original_id.clone();
    ids.sort_unstable();
    let mut expected = vec![id_a, id_b];
    expected.sort_unstable();
    assert_eq!(ids, expected);

    unsafe {
        manifold_rs_destroy(result);
        manifold_rs_destroy(a);
        manifold_rs_destroy(b);
    }
}

#[test]
fn three_operand_subtract_removes_both_cutters() {
    // A unit cube with two small cubes carved out of opposite corners. The
    // cutters are disjoint from each other, so (a - b) - c and a - (b u c)
    // describe the same solid — the point of the test is that the FFI applies
    // both cutters to the *first* operand rather than folding them together
    // in some other order.
    let a = unsafe { manifold_rs_as_original(ffi_cube([0.0, 0.0, 0.0], 1.0)) };
    let b = unsafe { manifold_rs_as_original(ffi_cube([-0.1, -0.1, -0.1], 0.3)) };
    let c = unsafe { manifold_rs_as_original(ffi_cube([0.8, 0.8, 0.8], 0.3)) };

    let inputs = [
        a as *const ManifoldRs,
        b as *const ManifoldRs,
        c as *const ManifoldRs,
    ];
    let result = unsafe { manifold_rs_batch_boolean(inputs.as_ptr(), inputs.len(), 1) };
    assert!(!result.is_null());
    assert_eq!(unsafe { manifold_rs_status(result) }, 0);
    let ffi_mesh = export(result);

    // Same operation via the core crate's CSG tree: bit-exact equality pins
    // that the FFI hands the operands to the tree in the documented order.
    let leaf = |m: *const ManifoldRs| CsgNode::leaf(unsafe { &*m }.inner.as_impl().clone());
    let expected_impl = CsgNode::op_n(OpType::Subtract, vec![leaf(a), leaf(b), leaf(c)]).evaluate();
    let expected = Manifold::from_impl(expected_impl);
    let expected_mesh = expected.get_mesh_gl(-1);
    assert_eq!(ffi_mesh.tri_verts, expected_mesh.tri_verts);
    assert_eq!(ffi_mesh.vert_properties, expected_mesh.vert_properties);

    // And the solid itself matches the sequential (a - b) - c, which is the
    // semantics C++ BatchBoolean promises for difference.
    let sequential = unsafe { &*a }
        .inner
        .difference(&unsafe { &*b }.inner)
        .difference(&unsafe { &*c }.inner);
    assert!(
        (expected.volume() - sequential.volume()).abs() < 1e-9,
        "batch subtract volume {} != sequential {}",
        expected.volume(),
        sequential.volume()
    );
    assert!(
        expected.volume() < 1.0 && expected.volume() > 0.9,
        "two small corners should be missing, got volume {}",
        expected.volume()
    );

    unsafe {
        manifold_rs_destroy(result);
        manifold_rs_destroy(a);
        manifold_rs_destroy(b);
        manifold_rs_destroy(c);
    }
}

#[test]
fn intersect_of_overlapping_cubes_is_non_empty() {
    let a = ffi_cube([0.0, 0.0, 0.0], 1.0);
    let b = ffi_cube([0.5, 0.0, 0.0], 1.0);
    let inputs = [a as *const ManifoldRs, b as *const ManifoldRs];

    let result = unsafe { manifold_rs_batch_boolean(inputs.as_ptr(), inputs.len(), 2) };
    assert!(!result.is_null());
    assert_eq!(unsafe { manifold_rs_status(result) }, 0);

    let mesh = export(result);
    assert!(
        !mesh.tri_verts.is_empty(),
        "intersection should not be empty"
    );

    unsafe {
        manifold_rs_destroy(result);
        manifold_rs_destroy(a);
        manifold_rs_destroy(b);
    }
}

#[test]
fn single_operand_batch_boolean_clones() {
    let a = unsafe { manifold_rs_as_original(ffi_cube([0.0, 0.0, 0.0], 1.0)) };
    let inputs = [a as *const ManifoldRs];
    let result = unsafe { manifold_rs_batch_boolean(inputs.as_ptr(), 1, 0) };
    assert!(!result.is_null());
    assert_eq!(
        unsafe { manifold_rs_original_id(result) },
        unsafe { manifold_rs_original_id(a) },
        "a one-operand boolean must not re-tag the mesh"
    );
    assert_eq!(export(result).tri_verts, export(a).tri_verts);

    unsafe {
        manifold_rs_destroy(result);
        manifold_rs_destroy(a);
    }
}

#[test]
fn non_manifold_input_reports_status_without_panicking() {
    // A tetrahedron with one face wound the wrong way: four verts and four
    // tris, so it clears the size checks, but two half-edges then run in the
    // same direction and the topology check rejects it.
    let verts: Vec<f32> = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    let tris: Vec<u32> = vec![0, 1, 2, 0, 1, 3, 0, 3, 2, 1, 2, 3];

    let bad =
        unsafe { manifold_rs_from_mesh(verts.as_ptr(), verts.len(), tris.as_ptr(), tris.len(), 3) };
    assert!(!bad.is_null(), "a validation failure still yields a handle");
    assert_eq!(
        unsafe { manifold_rs_status(bad) },
        2,
        "expected NotManifold (code 2)"
    );

    // The C++ binding crashes when exporting an error-status manifold; here it
    // must simply produce an empty mesh.
    let mesh = export(bad);
    assert!(mesh.tri_verts.is_empty());
    assert!(mesh.vert_properties.is_empty());

    unsafe { manifold_rs_destroy(bad) };
}

#[test]
fn invalid_from_mesh_arguments_return_null() {
    let (verts, tris) = cube_mesh([0.0, 0.0, 0.0], 1.0);
    unsafe {
        // num_prop below the three position slots.
        assert!(
            manifold_rs_from_mesh(verts.as_ptr(), verts.len(), tris.as_ptr(), tris.len(), 2)
                .is_null()
        );
        // vert_properties not a whole number of vertices.
        assert!(manifold_rs_from_mesh(
            verts.as_ptr(),
            verts.len() - 1,
            tris.as_ptr(),
            tris.len(),
            3
        )
        .is_null());
        // tri_verts not a whole number of triangles.
        assert!(manifold_rs_from_mesh(
            verts.as_ptr(),
            verts.len(),
            tris.as_ptr(),
            tris.len() - 1,
            3
        )
        .is_null());
        // Null array with a non-zero length.
        assert!(
            manifold_rs_from_mesh(std::ptr::null(), verts.len(), tris.as_ptr(), tris.len(), 3)
                .is_null()
        );
        // A length no allocation could ever have — what a negative int widened
        // to size_t looks like. usize::MAX happens to divide by both 3 and
        // num_prop, so it reaches the size check rather than the earlier ones.
        assert!(
            manifold_rs_from_mesh(verts.as_ptr(), usize::MAX, tris.as_ptr(), tris.len(), 3)
                .is_null()
        );
        assert!(
            manifold_rs_from_mesh(verts.as_ptr(), verts.len(), tris.as_ptr(), usize::MAX, 3)
                .is_null()
        );

        // Same guard on the operand count.
        let a = ffi_cube([0.0, 0.0, 0.0], 1.0);
        let inputs = [a as *const ManifoldRs];
        assert!(manifold_rs_batch_boolean(inputs.as_ptr(), usize::MAX, 0).is_null());
        manifold_rs_destroy(a);
    }
}

#[test]
fn null_handles_return_sentinels() {
    unsafe {
        assert!(manifold_rs_as_original(std::ptr::null()).is_null());
        assert_eq!(manifold_rs_original_id(std::ptr::null()), -1);
        assert_eq!(manifold_rs_status(std::ptr::null()), -1);
        assert!(manifold_rs_get_meshgl(std::ptr::null()).is_null());
        assert_eq!(manifold_rs_meshgl_num_prop(std::ptr::null()), 0);

        // Destroys are null-safe.
        manifold_rs_destroy(std::ptr::null_mut());
        manifold_rs_meshgl_destroy(std::ptr::null_mut());

        // A null manifold list, and a null entry inside a valid list.
        assert!(manifold_rs_batch_boolean(std::ptr::null(), 2, 0).is_null());
        let a = ffi_cube([0.0, 0.0, 0.0], 1.0);
        let with_null = [a as *const ManifoldRs, std::ptr::null()];
        assert!(manifold_rs_batch_boolean(with_null.as_ptr(), 2, 0).is_null());
        let valid = [a as *const ManifoldRs];
        assert!(manifold_rs_batch_boolean(valid.as_ptr(), 0, 0).is_null());
        manifold_rs_destroy(a);

        // Array accessors: NULL result, and out_len zeroed. out_len itself may
        // be null.
        let mut len = usize::MAX;
        assert!(manifold_rs_meshgl_vert_properties(std::ptr::null(), &mut len).is_null());
        assert_eq!(len, 0);
        len = usize::MAX;
        assert!(manifold_rs_meshgl_tri_verts(std::ptr::null(), &mut len).is_null());
        assert_eq!(len, 0);
        len = usize::MAX;
        assert!(manifold_rs_meshgl_run_index(std::ptr::null(), &mut len).is_null());
        assert_eq!(len, 0);
        len = usize::MAX;
        assert!(manifold_rs_meshgl_run_original_id(std::ptr::null(), &mut len).is_null());
        assert_eq!(len, 0);
        len = usize::MAX;
        assert!(manifold_rs_meshgl_face_id(std::ptr::null(), &mut len).is_null());
        assert_eq!(len, 0);
        assert!(
            manifold_rs_meshgl_vert_properties(std::ptr::null(), std::ptr::null_mut()).is_null()
        );
    }
}

#[test]
fn out_len_may_be_null_for_live_handles() {
    let cube = ffi_cube([0.0, 0.0, 0.0], 1.0);
    let g = unsafe { manifold_rs_get_meshgl(cube) };
    unsafe {
        assert!(!manifold_rs_meshgl_tri_verts(g, std::ptr::null_mut()).is_null());
        manifold_rs_meshgl_destroy(g);
        manifold_rs_destroy(cube);
    }
}

#[test]
fn a_panic_becomes_the_fallback_value_and_a_message() {
    // The geometry kernel can assert on degenerate input, and unwinding out of
    // an extern "C" function is undefined behaviour, so exercise the guard the
    // exported functions all wrap themselves in. (The panic message printed by
    // the default hook during this test is expected.)
    let value = crate::error::guard(-1i32, || panic!("synthetic kernel failure"));
    assert_eq!(value, -1);

    let mut buf = [0u8; 256];
    let len = unsafe { manifold_rs_last_error(buf.as_mut_ptr(), buf.len()) };
    let message = std::str::from_utf8(&buf[..len]).unwrap();
    assert!(
        message.contains("synthetic kernel failure"),
        "unexpected message: {message}"
    );
}

#[test]
fn last_error_is_empty_on_a_thread_that_has_not_failed() {
    // Checked on a freshly spawned thread rather than the test's own: the slot
    // is thread-local, and with --test-threads=1 every test shares one thread,
    // so this thread may already carry a message from an earlier test.
    std::thread::spawn(|| {
        let mut buf = [0u8; 64];
        assert_eq!(
            unsafe { manifold_rs_last_error(buf.as_mut_ptr(), buf.len()) },
            0
        );
        assert_eq!(
            unsafe { manifold_rs_last_error(std::ptr::null_mut(), 0) },
            0
        );
    })
    .join()
    .unwrap();
}

#[test]
fn last_error_reports_the_most_recent_failure() {
    let mut buf = [0u8; 256];
    let a = ffi_cube([0.0, 0.0, 0.0], 1.0);
    let inputs = [a as *const ManifoldRs];

    // An unknown op both fails and records why.
    assert!(unsafe { manifold_rs_batch_boolean(inputs.as_ptr(), 1, 99) }.is_null());

    let len = unsafe { manifold_rs_last_error(buf.as_mut_ptr(), buf.len()) };
    assert!(len > 0, "unknown op should record a message");
    let message = std::str::from_utf8(&buf[..len]).unwrap().to_string();
    assert!(
        message.contains("unknown op 99"),
        "unexpected message: {message}"
    );

    // Passing NULL asks only for the length, and never writes.
    assert_eq!(
        unsafe { manifold_rs_last_error(std::ptr::null_mut(), 0) },
        len
    );

    // A short buffer truncates but still reports the full length.
    let mut small = [0u8; 4];
    assert_eq!(
        unsafe { manifold_rs_last_error(small.as_mut_ptr(), small.len()) },
        len
    );
    assert_eq!(&small, &message.as_bytes()[..4]);

    // The next failure replaces the message rather than appending to it, so a
    // caller always reads the failure it just saw.
    assert!(
        unsafe { manifold_rs_from_mesh(std::ptr::null(), 0, std::ptr::null(), 0, 1) }.is_null()
    );
    let len = unsafe { manifold_rs_last_error(buf.as_mut_ptr(), buf.len()) };
    let second = std::str::from_utf8(&buf[..len]).unwrap();
    assert!(
        second.contains("num_prop 1 < 3"),
        "unexpected message: {second}"
    );
    assert!(
        !second.contains("unknown op 99"),
        "message was appended, not replaced"
    );

    unsafe { manifold_rs_destroy(a) };
}
