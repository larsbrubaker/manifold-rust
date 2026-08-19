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

//! C ABI for the core crate's progress reporting
//! ([`manifold_rust::progress`]).
//!
//! Shaped like `cancel.rs`: an optional extra argument on a boolean entry
//! point, with a NULL callback meaning "no reporting" and costing nothing.
//! The callback is passed as a plain function pointer plus an opaque `user`
//! word, so a C or C# consumer can marshal it however it likes.
//!
//! Two contracts the header repeats for callers:
//! - the callback is **never re-entered concurrently** (the core reporter
//!   serializes it), but it *may* run on a worker thread rather than the one
//!   that made the call, because the boolean pipeline is internally parallel;
//! - a negative `fraction` means "indeterminate" — the phase is running but has
//!   no meaningful completion ratio. Otherwise it is in `[0, 1]`.

use std::ffi::CString;
use std::os::raw::{c_char, c_void};
use std::ptr;
use std::sync::OnceLock;

use manifold_rust::cancel::CancelToken;
use manifold_rust::manifold::Manifold;
use manifold_rust::progress::{Phase, ProgressReporter};
use manifold_rust::types::{BooleanEngine, OpType, WindingRule};

use crate::cancel::CancelTokenRs;
use crate::error::{guard, set_last_error};
use crate::{into_handle, ManifoldRs};

/// Progress callback: `phase_id` indexes the phase table (see
/// [`manifold_rs_progress_phase_name`]), `fraction` is in `[0, 1]` or negative
/// for indeterminate, and `user` is the pointer handed to the boolean call.
pub type ManifoldRsProgressFn =
    Option<extern "C" fn(phase_id: u32, fraction: f64, user: *mut c_void)>;

/// Sentinel `fraction` for a phase with no completion ratio.
const INDETERMINATE: f64 = -1.0;

/// The C callback plus its user word, packaged so it can be moved into the
/// core crate's `Send + Sync` reporter closure.
struct CCallback {
    f: extern "C" fn(u32, f64, *mut c_void),
    user: *mut c_void,
}

// SAFETY: the raw `user` word is only ever handed straight back to the
// caller's own function, never dereferenced here. The core reporter holds a
// mutex across every invocation, so the callback is never re-entered
// concurrently; what crossing threads costs the caller is that the callback
// may run on a pipeline worker thread, which `manifold_rs.h` makes part of the
// contract. A caller whose `user` cannot survive that must marshal inside its
// own callback.
unsafe impl Send for CCallback {}
unsafe impl Sync for CCallback {}

/// NUL-terminated copies of every [`Phase`] name, built once so the exported
/// accessor can hand out stable pointers without the header duplicating the
/// name table (and drifting from it).
fn phase_names() -> &'static [CString] {
    static NAMES: OnceLock<Vec<CString>> = OnceLock::new();
    NAMES.get_or_init(|| {
        Phase::ALL
            .iter()
            .map(|p| CString::new(p.name()).expect("phase names contain no NUL"))
            .collect()
    })
}

/// Static NUL-terminated display name for a phase id, or NULL for an id this
/// build does not know (ids are append-only, so a newer library can report a
/// phase an older caller has no name for — display the id then).
///
/// The returned pointer is valid for the process lifetime and must not be
/// freed.
#[no_mangle]
pub extern "C" fn manifold_rs_progress_phase_name(phase_id: u32) -> *const c_char {
    guard(ptr::null(), || match phase_names().get(phase_id as usize) {
        Some(name) => name.as_ptr(),
        None => ptr::null(),
    })
}

/// Number of phase ids this build defines; ids `0..count` all have names.
#[no_mangle]
pub extern "C" fn manifold_rs_progress_phase_count() -> u32 {
    Phase::ALL.len() as u32
}

/// Binary boolean with optional cancellation *and* optional progress
/// reporting.
///
/// `op` is 0 = union, 1 = difference, 2 = intersection (the
/// `manifold_rs_batch_boolean` ordering). `engine` is 0 = Exact, 1 = Robust,
/// 2 = Auto, matching `manifold_rs_set_boolean_engine`; unlike the batch entry
/// points this one takes the engine explicitly rather than reading the
/// process-global default, since a caller that wants progress detail usually
/// also wants to pin the engine.
///
/// A NULL `token` is uncancellable and a NULL `progress` reports nothing —
/// both are then exactly the un-instrumented pipeline. Returns NULL only for
/// an argument error or a caught panic; a cancelled operation returns a valid
/// handle whose status is the cancelled code (14).
///
/// # Safety
/// `a` and `b` must be live handles. `token` must be NULL or a live handle at
/// the moment of the call. `progress` must be NULL or a function pointer that
/// stays valid — and whose `user` word stays usable, from any thread — until
/// this call returns.
#[no_mangle]
pub unsafe extern "C" fn manifold_rs_boolean_progress(
    a: *const ManifoldRs,
    b: *const ManifoldRs,
    op: i32,
    engine: i32,
    token: *const CancelTokenRs,
    progress: ManifoldRsProgressFn,
    user: *mut c_void,
) -> *mut ManifoldRs {
    // SAFETY: identical contract, forwarded unchanged.
    unsafe { manifold_rs_boolean_progress_rule(a, b, op, engine, 0, token, progress, user) }
}

/// [`manifold_rs_boolean_progress`] with an explicit winding rule:
/// 0 = positive (`w >= 1`, the default and the historical behavior),
/// 1 = nonzero (`w != 0`, inside-out geometry stays solid).
///
/// The rule is a robust-engine semantic — the exact engine ignores it, and
/// `engine = MANIFOLD_RS_ENGINE_AUTO` resolves to the robust engine whenever
/// the rule is nonzero. Returns NULL for an unknown rule value.
///
/// # Safety
/// Same contract as [`manifold_rs_boolean_progress`].
#[no_mangle]
pub unsafe extern "C" fn manifold_rs_boolean_progress_rule(
    a: *const ManifoldRs,
    b: *const ManifoldRs,
    op: i32,
    engine: i32,
    winding_rule: i32,
    token: *const CancelTokenRs,
    progress: ManifoldRsProgressFn,
    user: *mut c_void,
) -> *mut ManifoldRs {
    guard(ptr::null_mut(), || {
        let op = match op {
            0 => OpType::Add,
            1 => OpType::Subtract,
            2 => OpType::Intersect,
            other => {
                set_last_error(format!("manifold_rs_boolean_progress: unknown op {other}"));
                return ptr::null_mut();
            }
        };
        let engine = match engine {
            0 => BooleanEngine::Exact,
            1 => BooleanEngine::Robust,
            2 => BooleanEngine::Auto,
            other => {
                set_last_error(format!(
                    "manifold_rs_boolean_progress: unknown engine {other}"
                ));
                return ptr::null_mut();
            }
        };
        let rule = match winding_rule {
            0 => WindingRule::Positive,
            1 => WindingRule::Nonzero,
            other => {
                set_last_error(format!(
                    "manifold_rs_boolean_progress: unknown winding rule {other}"
                ));
                return ptr::null_mut();
            }
        };
        // SAFETY: caller contract; as_ref() handles the NULL case.
        let (Some(a), Some(b)) = (unsafe { a.as_ref() }, unsafe { b.as_ref() }) else {
            set_last_error("manifold_rs_boolean_progress: null manifold");
            return ptr::null_mut();
        };

        // Same ownership reasoning as manifold_rs_batch_boolean_ct: clone the
        // shared flag so the operation no longer depends on the token *handle*
        // outliving the call.
        // SAFETY: caller contract; as_ref() handles the NULL case.
        let token: Option<CancelToken> = unsafe { token.as_ref() }.map(|t| t.inner.clone());
        let token: Option<&CancelToken> = token.as_ref();

        let reporter = progress.map(|f| {
            let cb = CCallback { f, user };
            ProgressReporter::new(move |phase: Phase, fraction| {
                // Bind the whole struct first: edition-2021 disjoint capture
                // would otherwise capture the bare `*mut c_void` field, losing
                // the `Send + Sync` the wrapper exists to provide.
                let cb = &cb;
                (cb.f)(phase.id(), fraction.unwrap_or(INDETERMINATE), cb.user);
            })
        });

        let result: Manifold = a.inner.boolean_with_engine_rule_and_progress(
            &b.inner,
            op,
            engine,
            rule,
            token,
            reporter.as_ref(),
        );
        into_handle(result)
    })
}
