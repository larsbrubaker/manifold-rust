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

// Unit tests for the progress C ABI (src/progress.rs), called the way a C
// caller does: a plain `extern "C"` function pointer plus an opaque user word.

use std::ffi::CStr;
use std::os::raw::c_void;
use std::ptr;
use std::sync::Mutex;

use crate::cancel::{
    manifold_rs_cancel_token_cancel, manifold_rs_cancel_token_destroy, manifold_rs_cancel_token_new,
};
use crate::progress::*;
use crate::tests::{export, ffi_cube};
use crate::*;

/// What a C caller would put behind its `user` pointer. A mutex because the
/// header promises only that the callback is never re-entered *concurrently*,
/// not that it runs on the calling thread.
#[derive(Default)]
struct Collected {
    events: Mutex<Vec<(u32, f64)>>,
}

extern "C" fn collect(phase_id: u32, fraction: f64, user: *mut c_void) {
    assert!(!user.is_null(), "user word must round-trip unchanged");
    // SAFETY: `user` is the &Collected the test passed in, which outlives the
    // boolean call that drives this callback.
    let sink = unsafe { &*(user as *const Collected) };
    sink.events.lock().unwrap().push((phase_id, fraction));
}

fn run(engine: i32, progress: ManifoldRsProgressFn, user: *mut c_void) -> *mut ManifoldRs {
    let a = ffi_cube([0.0, 0.0, 0.0], 1.0);
    let b = ffi_cube([0.5, 0.25, 0.25], 1.0);
    let out = unsafe {
        manifold_rs_boolean_progress(a, b, 0, engine, ptr::null(), progress, user)
    };
    unsafe {
        manifold_rs_destroy(a);
        manifold_rs_destroy(b);
    }
    out
}

#[test]
fn every_phase_id_has_a_name_and_out_of_range_is_null() {
    let count = manifold_rs_progress_phase_count();
    assert!(count >= 8, "expected the documented phase table, got {count}");
    for id in 0..count {
        let ptr = manifold_rs_progress_phase_name(id);
        assert!(!ptr.is_null(), "phase {id} has no name");
        let name = unsafe { CStr::from_ptr(ptr) }.to_str().unwrap();
        assert!(!name.is_empty());
    }
    assert!(manifold_rs_progress_phase_name(count).is_null());
    assert!(manifold_rs_progress_phase_name(u32::MAX).is_null());
    // Pointers are static: two calls agree, so a caller may cache them.
    assert_eq!(
        manifold_rs_progress_phase_name(0),
        manifold_rs_progress_phase_name(0)
    );
}

#[test]
fn a_robust_boolean_reports_named_phases_with_valid_fractions() {
    let sink = Collected::default();
    let result = run(1, Some(collect), &sink as *const _ as *mut c_void);
    assert!(!result.is_null());
    assert_eq!(unsafe { manifold_rs_status(result) }, 0);

    let events = sink.events.lock().unwrap().clone();
    assert!(!events.is_empty(), "robust boolean reported nothing");
    let mut last = 0;
    for &(id, fraction) in &events {
        assert!(
            !manifold_rs_progress_phase_name(id).is_null(),
            "reported unnamed phase {id}"
        );
        assert!(id >= last, "phase {id} went backwards from {last}");
        last = id;
        assert!(
            fraction < 0.0 || (0.0..=1.0).contains(&fraction),
            "fraction {fraction} is neither indeterminate nor in [0,1]"
        );
    }
    unsafe { manifold_rs_destroy(result) };
}

#[test]
fn a_null_callback_matches_the_uninstrumented_result() {
    let plain = run(1, None, ptr::null_mut());
    let sink = Collected::default();
    let watched = run(1, Some(collect), &sink as *const _ as *mut c_void);
    assert!(!plain.is_null() && !watched.is_null());

    let (p, w) = (export(plain), export(watched));
    assert_eq!(p.tri_verts, w.tri_verts);
    assert_eq!(p.vert_properties, w.vert_properties);
    assert_eq!(p.run_index, w.run_index);
    assert!(!sink.events.lock().unwrap().is_empty());

    unsafe {
        manifold_rs_destroy(plain);
        manifold_rs_destroy(watched);
    }
}

#[test]
fn the_exact_engine_reports_a_single_indeterminate_phase() {
    let sink = Collected::default();
    let result = run(0, Some(collect), &sink as *const _ as *mut c_void);
    assert!(!result.is_null());
    let events = sink.events.lock().unwrap().clone();
    assert_eq!(events.len(), 1, "exact engine should report exactly one phase");
    assert!(events[0].1 < 0.0, "an indeterminate phase reports a negative fraction");
    let name = unsafe { CStr::from_ptr(manifold_rs_progress_phase_name(events[0].0)) };
    assert_eq!(name.to_str().unwrap(), "exact boolean");
    unsafe { manifold_rs_destroy(result) };
}

#[test]
fn argument_errors_return_null() {
    let a = ffi_cube([0.0, 0.0, 0.0], 1.0);
    unsafe {
        // Unknown op / engine.
        assert!(manifold_rs_boolean_progress(a, a, 99, 0, ptr::null(), None, ptr::null_mut())
            .is_null());
        assert!(manifold_rs_boolean_progress(a, a, 0, 99, ptr::null(), None, ptr::null_mut())
            .is_null());
        // NULL operands.
        assert!(manifold_rs_boolean_progress(
            ptr::null(),
            a,
            0,
            0,
            ptr::null(),
            None,
            ptr::null_mut()
        )
        .is_null());
        assert!(manifold_rs_boolean_progress(
            a,
            ptr::null(),
            0,
            0,
            ptr::null(),
            None,
            ptr::null_mut()
        )
        .is_null());
        manifold_rs_destroy(a);
    }
}

#[test]
fn a_pre_cancelled_token_still_wins_with_a_reporter_attached() {
    let a = ffi_cube([0.0, 0.0, 0.0], 1.0);
    let b = ffi_cube([0.5, 0.0, 0.0], 1.0);
    let t = manifold_rs_cancel_token_new();
    unsafe { manifold_rs_cancel_token_cancel(t) };

    let sink = Collected::default();
    let result = unsafe {
        manifold_rs_boolean_progress(
            a,
            b,
            0,
            1,
            t,
            Some(collect),
            &sink as *const _ as *mut c_void,
        )
    };
    assert!(!result.is_null(), "cancellation must not look like a failure");
    // 14 = Error::Cancelled, the code the header documents.
    assert_eq!(unsafe { manifold_rs_status(result) }, 14);

    unsafe {
        manifold_rs_destroy(result);
        manifold_rs_destroy(a);
        manifold_rs_destroy(b);
        manifold_rs_cancel_token_destroy(t);
    }
}
