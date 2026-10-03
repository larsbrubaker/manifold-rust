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

use std::os::raw::c_void;

use crate::cancel::*;
use crate::progress::manifold_rs_boolean_progress;
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

/// A pair of heavily overlapping spheres, imported through the FFI. Their
/// union runs the whole robust pipeline, which is all the cross-thread test
/// needs: it counts phases, so the input does not have to be slow.
fn ffi_sphere_pair() -> (*mut ManifoldRs, *mut ManifoldRs) {
    use manifold_rust::linalg::Vec3;

    let build = |offset: f64| {
        let mesh = Manifold::sphere(1.0, 32)
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
    assert!(
        !result.is_null(),
        "cancellation must not look like a failure"
    );
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

/// What the progress callback's `user` word points at: the phases entered,
/// consecutive repeats folded, plus an optional pair of rendezvous that park
/// the kernel in its first report.
struct PhaseLog {
    phases: std::sync::Mutex<Vec<u32>>,
    park: Option<(std::sync::Barrier, std::sync::Barrier)>,
}

extern "C" fn log_phase(phase_id: u32, _fraction: f64, user: *mut c_void) {
    // SAFETY: `user` is the &PhaseLog the test passed in, which outlives the
    // boolean call that drives this callback.
    let log = unsafe { &*(user as *const PhaseLog) };
    let first = {
        let mut seen = log.phases.lock().expect("phase list poisoned");
        let first = seen.is_empty();
        if first || seen.last() != Some(&phase_id) {
            seen.push(phase_id);
        }
        first
    };
    if let (true, Some((inside, sent))) = (first, &log.park) {
        inside.wait();
        sent.wait();
    }
}

/// A cancel sent from another thread through the C ABI stops the work before
/// it finishes. Measured in work, not time, like the crate's
/// `cancel_from_another_thread_interrupts_a_boolean_in_flight` (see its doc):
/// the earlier wall-time ratio compared two tens-of-milliseconds runs and
/// flaked, its `uncancelled > 20ms` precondition failing outright on a fast
/// machine. The robust engine's progress phases are driven by work, so the
/// cancelled run, parked in its first report until the token is cancelled
/// from this thread, must enter fewer phases than the full run.
#[test]
fn cross_thread_cancel_interrupts_a_slow_boolean() {
    let (a, b) = ffi_sphere_pair();
    // Captures nothing, so the worker thread can share it; the handles are
    // passed in, wrapped in `SendPtr` on the way across.
    let run = |a: *const ManifoldRs,
               b: *const ManifoldRs,
               token: *const CancelTokenRs,
               log: &PhaseLog|
     -> i32 {
        // SAFETY: live handles; `log` outlives the call. 0 = union, 1 = Robust.
        let result = unsafe {
            manifold_rs_boolean_progress(
                a,
                b,
                0,
                1,
                token,
                Some(log_phase),
                log as *const PhaseLog as *mut c_void,
            )
        };
        assert!(!result.is_null(), "cancellation must not return NULL");
        let status = unsafe { manifold_rs_status(result) };
        unsafe { manifold_rs_destroy(result) };
        status
    };

    let full = PhaseLog {
        phases: Default::default(),
        park: None,
    };
    assert_eq!(run(a, b, std::ptr::null(), &full), 0);
    let full_run = full.phases.into_inner().expect("phase list poisoned");
    assert!(full_run.len() > 1, "no later work to skip: {full_run:?}");

    let t = manifold_rs_cancel_token_new();
    let worker_args = (
        SendPtr(a as *const ManifoldRs),
        SendPtr(b as *const ManifoldRs),
        SendPtr(t as *const CancelTokenRs),
    );
    let parked = PhaseLog {
        phases: Default::default(),
        park: Some((std::sync::Barrier::new(2), std::sync::Barrier::new(2))),
    };
    let status = std::thread::scope(|s| {
        let parked = &parked;
        let worker = s.spawn(move || {
            let (a, b, token) = worker_args;
            run(a.0, b.0, token.0, parked)
        });
        let (inside, sent) = parked.park.as_ref().expect("parked log");
        // Not a delay: the boolean has reported its first phase and is
        // parked there.
        inside.wait();
        // The point of the whole feature: cancel from a different thread
        // while the kernel is running.
        unsafe { manifold_rs_cancel_token_cancel(t) };
        sent.wait();
        worker.join().expect("worker panicked")
    });

    assert_eq!(status, CANCELLED);
    let cancelled_run = parked.phases.into_inner().expect("phase list poisoned");
    assert!(
        cancelled_run.len() < full_run.len(),
        "the cancelled boolean went on to report {cancelled_run:?}, every phase of the \
         full run {full_run:?} - the cancel is being ignored until the work finishes"
    );

    // Destroying the token only after the call using it has returned — the
    // ordering the header makes the caller's responsibility.
    unsafe {
        manifold_rs_cancel_token_destroy(t);
        manifold_rs_destroy(a);
        manifold_rs_destroy(b);
    }
}
