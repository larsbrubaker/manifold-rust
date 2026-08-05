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

//! C ABI for the core crate's cooperative cancellation
//! ([`manifold_rust::cancel`]).
//!
//! The token is the one handle in this library that is deliberately meant to be
//! touched from two threads at once: one thread sits inside
//! `manifold_rs_batch_boolean_ct`, another calls
//! `manifold_rs_cancel_token_cancel`. That is the whole point of the type, and
//! it is sound because the underlying flag is an atomic.
//!
//! Destroy is memory-safe even against a running operation: the operation
//! clones the shared flag on entry (see `manifold_rs_batch_boolean_ct`), so the
//! flag outlives the handle. Destroying early is still a caller bug — after it
//! there is nothing left to cancel with, so the operation runs to completion —
//! but it is a lost cancellation, not a use-after-free. See `manifold_rs.h`.

use std::ptr;

use manifold_rust::cancel::CancelToken;

use crate::error::{guard, set_last_error};

/// Opaque handle wrapping a [`CancelToken`]. Created by
/// [`manifold_rs_cancel_token_new`], released by
/// [`manifold_rs_cancel_token_destroy`].
pub struct CancelTokenRs {
    pub(crate) inner: CancelToken,
}

/// A fresh, uncancelled token. NULL only if the allocation panics.
#[no_mangle]
pub extern "C" fn manifold_rs_cancel_token_new() -> *mut CancelTokenRs {
    guard(ptr::null_mut(), || {
        Box::into_raw(Box::new(CancelTokenRs {
            inner: CancelToken::new(),
        }))
    })
}

/// Request cancellation. Safe to call from any thread, including while another
/// thread is inside an operation that was handed this token — that is the
/// intended use. NULL is ignored (with a last-error message).
///
/// Cancellation is sticky: there is no "uncancel". Create a new token for work
/// that should be allowed to run.
///
/// # Safety
/// `t` must be NULL or a live handle from this library.
#[no_mangle]
pub unsafe extern "C" fn manifold_rs_cancel_token_cancel(t: *const CancelTokenRs) {
    guard((), || {
        // SAFETY: caller contract; as_ref() handles the NULL case. A shared
        // reference is enough — the flag is an atomic, so concurrent readers
        // inside a running operation are fine.
        match unsafe { t.as_ref() } {
            Some(handle) => handle.inner.cancel(),
            None => set_last_error("manifold_rs_cancel_token_cancel: null token"),
        }
    })
}

/// Whether cancellation has been requested. 0 for a NULL handle, which is also
/// the "not cancelled" answer — a NULL token is uncancellable by definition.
///
/// # Safety
/// `t` must be NULL or a live handle from this library.
#[no_mangle]
pub unsafe extern "C" fn manifold_rs_cancel_token_is_cancelled(t: *const CancelTokenRs) -> i32 {
    guard(0, || {
        // SAFETY: caller contract; as_ref() handles the NULL case.
        match unsafe { t.as_ref() } {
            Some(handle) => i32::from(handle.inner.is_cancelled()),
            None => {
                set_last_error("manifold_rs_cancel_token_is_cancelled: null token");
                0
            }
        }
    })
}

/// Release a token handle. NULL is ignored.
///
/// Destroying a token that a running operation was given is safe — the
/// operation holds its own share of the shared flag — but pointless: nothing
/// can cancel that operation afterwards. Keep the handle until the call
/// returns if you want the cancel to be able to land.
///
/// # Safety
/// `t` must be NULL or a handle from this library that has not already been
/// destroyed, and no other thread may be calling
/// [`manifold_rs_cancel_token_cancel`] on it concurrently (that call reads the
/// handle itself, unlike an in-flight operation, which does not).
#[no_mangle]
pub unsafe extern "C" fn manifold_rs_cancel_token_destroy(t: *mut CancelTokenRs) {
    guard((), || {
        if !t.is_null() {
            // SAFETY: caller guarantees single ownership of a live handle that
            // is not concurrently passed to another token entry point.
            drop(unsafe { Box::from_raw(t) });
        }
    })
}
