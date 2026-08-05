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

// Unit tests for the cancellation C ABI (src/cancel.rs) and the cancellable
// batch boolean, called the way a C caller does.

use std::time::{Duration, Instant};

use crate::cancel::*;
use crate::tests::{export, ffi_cube};
use crate::*;

/// Status code for `Error::Cancelled`; the value the header documents.
const CANCELLED: i32 = 14;

/// Raw handles are not `Send`, but the C contract explicitly allows a manifold
/// to be used on one thread while its token is cancelled from another. This
/// wrapper carries a handle across the boundary for exactly that scenario.
struct SendPtr<T>(*const T);
// SAFETY: the pointed-to handle is only read (never freed) on the worker
// thread, and the main thread does not touch it until after `join()`.
unsafe impl<T> Send for SendPtr<T> {}

/// A pair of heavily overlapping high-resolution spheres, imported through the
/// FFI. Big enough that a union takes long enough to be cancelled mid-flight.
fn ffi_sphere_pair() -> (*mut ManifoldRs, *mut ManifoldRs) {
    use manifold_rust::linalg::Vec3;

    let build = |offset: f64| {
        let mesh = Manifold::sphere(1.0, 256)
            .translate(Vec3::new(offset, 0.0, 0.0))
            .get_mesh_gl(-1);
        let handle = unsafe {
            manifold_rs_from_mesh(
                mesh.vert_properties.as_ptr(),
                mesh.vert_properties.len(),
                mesh.tri_verts.as_ptr(),
                mesh.tri_verts.len(),
                mesh.num_prop,
            )
        };
        assert!(!handle.is_null(), "sphere import returned NULL");
        assert_eq!(unsafe { manifold_rs_status(handle) }, 0);
        handle
    };
    (build(0.0), build(0.5))
}

#[test]
fn token_lifecycle_and_null_safety() {
    let t = manifold_rs_cancel_token_new();
    assert!(!t.is_null(), "token allocation should not fail");

    unsafe {
        assert_eq!(manifold_rs_cancel_token_is_cancelled(t), 0);
        manifold_rs_cancel_token_cancel(t);
        assert_eq!(manifold_rs_cancel_token_is_cancelled(t), 1);
        // Sticky: cancelling twice is a no-op, not an error.
        manifold_rs_cancel_token_cancel(t);
        assert_eq!(manifold_rs_cancel_token_is_cancelled(t), 1);
        manifold_rs_cancel_token_destroy(t);

        // Every entry point tolerates NULL.
        manifold_rs_cancel_token_cancel(std::ptr::null());
        assert_eq!(manifold_rs_cancel_token_is_cancelled(std::ptr::null()), 0);
        manifold_rs_cancel_token_destroy(std::ptr::null_mut());
    }
}

#[test]
fn pre_cancelled_token_yields_a_handle_with_the_cancelled_status() {
    let a = ffi_cube([0.0, 0.0, 0.0], 1.0);
    let b = ffi_cube([0.5, 0.0, 0.0], 1.0);
    let inputs = [a as *const ManifoldRs, b as *const ManifoldRs];

    let t = manifold_rs_cancel_token_new();
    unsafe { manifold_rs_cancel_token_cancel(t) };

    let result = unsafe { manifold_rs_batch_boolean_ct(inputs.as_ptr(), inputs.len(), 0, t) };
    // NOT NULL: a caller must be able to tell cancellation from a panic or a
    // bad argument, both of which are NULL.
    assert!(!result.is_null(), "cancellation must not look like a failure");
    assert_eq!(unsafe { manifold_rs_status(result) }, CANCELLED);
    assert!(export(result).tri_verts.is_empty());

    unsafe {
        manifold_rs_destroy(result);
        manifold_rs_destroy(a);
        manifold_rs_destroy(b);
        manifold_rs_cancel_token_destroy(t);
    }
}

#[test]
fn a_pre_cancelled_token_beats_the_single_operand_shortcut() {
    let a = ffi_cube([0.0, 0.0, 0.0], 1.0);
    let inputs = [a as *const ManifoldRs];
    let t = manifold_rs_cancel_token_new();
    unsafe { manifold_rs_cancel_token_cancel(t) };

    let result = unsafe { manifold_rs_batch_boolean_ct(inputs.as_ptr(), 1, 0, t) };
    assert!(!result.is_null());
    assert_eq!(unsafe { manifold_rs_status(result) }, CANCELLED);

    unsafe {
        manifold_rs_destroy(result);
        manifold_rs_destroy(a);
        manifold_rs_cancel_token_destroy(t);
    }
}

#[test]
fn null_token_is_identical_to_the_uncancellable_entry_point() {
    let a = unsafe { manifold_rs_as_original(ffi_cube([0.0, 0.0, 0.0], 1.0)) };
    let b = unsafe { manifold_rs_as_original(ffi_cube([0.5, 0.0, 0.0], 1.0)) };
    let inputs = [a as *const ManifoldRs, b as *const ManifoldRs];

    let old = unsafe { manifold_rs_batch_boolean(inputs.as_ptr(), inputs.len(), 0) };
    let new =
        unsafe { manifold_rs_batch_boolean_ct(inputs.as_ptr(), inputs.len(), 0, std::ptr::null()) };
    assert!(!old.is_null() && !new.is_null());
    assert_eq!(unsafe { manifold_rs_status(old) }, 0);
    assert_eq!(unsafe { manifold_rs_status(new) }, 0);

    let (old_mesh, new_mesh) = (export(old), export(new));
    assert_eq!(old_mesh.tri_verts, new_mesh.tri_verts);
    assert_eq!(old_mesh.vert_properties, new_mesh.vert_properties);
    assert_eq!(old_mesh.run_index, new_mesh.run_index);

    // And the argument-validation sentinels are shared, since one delegates to
    // the other.
    unsafe {
        assert!(manifold_rs_batch_boolean_ct(inputs.as_ptr(), 0, 0, std::ptr::null()).is_null());
        assert!(manifold_rs_batch_boolean_ct(std::ptr::null(), 2, 0, std::ptr::null()).is_null());
        assert!(manifold_rs_batch_boolean_ct(inputs.as_ptr(), 2, 99, std::ptr::null()).is_null());
        let with_null = [a as *const ManifoldRs, std::ptr::null()];
        assert!(manifold_rs_batch_boolean_ct(with_null.as_ptr(), 2, 0, std::ptr::null()).is_null());
    }

    unsafe {
        manifold_rs_destroy(old);
        manifold_rs_destroy(new);
        manifold_rs_destroy(a);
        manifold_rs_destroy(b);
    }
}

#[test]
fn cross_thread_cancel_interrupts_a_slow_boolean() {
    let (a, b) = ffi_sphere_pair();
    let inputs = [a as *const ManifoldRs, b as *const ManifoldRs];

    // Baseline: the same op, uncancelled, timed in-test so the assertion is a
    // ratio rather than an absolute millisecond count.
    let start = Instant::now();
    let baseline = unsafe { manifold_rs_batch_boolean(inputs.as_ptr(), inputs.len(), 0) };
    let uncancelled = start.elapsed();
    assert!(!baseline.is_null());
    assert_eq!(unsafe { manifold_rs_status(baseline) }, 0);
    assert!(
        uncancelled > Duration::from_millis(20),
        "test input is too fast ({uncancelled:?}) to be a meaningful cancel target"
    );
    unsafe { manifold_rs_destroy(baseline) };

    let t = manifold_rs_cancel_token_new();
    let worker_inputs = SendPtr(inputs.as_ptr());
    let worker_token = SendPtr(t as *const CancelTokenRs);
    // The handshake pins the start: the main thread does not begin its delay
    // until the worker is at the call, so the cancel lands on work in flight
    // rather than racing the thread spawn.
    let (started_tx, started_rx) = std::sync::mpsc::channel();
    let worker = std::thread::spawn(move || {
        let inputs = worker_inputs;
        let token = worker_token;
        started_tx.send(()).expect("main thread went away");
        let start = Instant::now();
        // SAFETY: the handles and the token outlive the thread — the main
        // thread joins before destroying any of them.
        let result = unsafe { manifold_rs_batch_boolean_ct(inputs.0, 2, 0, token.0) };
        let elapsed = start.elapsed();
        assert!(!result.is_null(), "cancellation must not return NULL");
        let status = unsafe { manifold_rs_status(result) };
        unsafe { manifold_rs_destroy(result) };
        (status, elapsed)
    });
    started_rx.recv().expect("worker never started");

    // Small fraction of the runtime, scaled off the baseline so the ratio
    // stays stable on a loaded machine where both numbers inflate.
    std::thread::sleep((uncancelled / 16).max(Duration::from_millis(1)));
    // The point of the whole feature: cancel from a different thread while the
    // kernel is running.
    unsafe { manifold_rs_cancel_token_cancel(t) };
    let (status, cancelled_elapsed) = worker.join().expect("worker panicked");

    assert_eq!(status, CANCELLED);
    assert!(
        cancelled_elapsed * 2 < uncancelled,
        "cancelled boolean took {cancelled_elapsed:?}, which is not well under \
         the uncancelled {uncancelled:?}"
    );

    // Destroying the token only after the call using it has returned — the
    // ordering the header makes the caller's responsibility.
    unsafe {
        manifold_rs_cancel_token_destroy(t);
        manifold_rs_destroy(a);
        manifold_rs_destroy(b);
    }
}
