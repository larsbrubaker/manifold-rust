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

// Panic containment and the thread-local "last error" slot used by lib.rs.
//
// Unwinding across an `extern "C"` boundary is undefined behaviour, and the
// geometry kernel can panic on degenerate input (e.g. the assert in
// boolean_result.rs). Every exported function therefore runs its body through
// `guard`, which converts a panic into the function's failure sentinel and
// records the payload message for `manifold_rs_last_error`.

use std::any::Any;
use std::cell::RefCell;
use std::panic::{catch_unwind, AssertUnwindSafe};

thread_local! {
    /// Overwritten on every failure, never cleared on success: a caller that
    /// sees a failure sentinel can still read the message afterwards.
    static LAST_ERROR: RefCell<String> = const { RefCell::new(String::new()) };
}

/// Record a message describing why an exported function is about to fail.
///
/// `try_with`, not `with`: once a thread is being torn down the thread-local is
/// destroyed and `with` panics. That panic would fire in `guard`'s error arm —
/// outside `catch_unwind` — and abort the process, contradicting the header's
/// promise that no function aborts. Losing a message on a dying thread is the
/// better trade.
pub(crate) fn set_last_error(message: impl Into<String>) {
    let message = message.into();
    let _ = LAST_ERROR.try_with(|slot| *slot.borrow_mut() = message);
}

/// Copy of the current thread's last error message. Empty when there was no
/// failure, or when the thread-local is already gone (see `set_last_error`).
pub(crate) fn last_error_message() -> String {
    LAST_ERROR
        .try_with(|slot| slot.borrow().clone())
        .unwrap_or_default()
}

/// Run `body`, returning `fallback` if it panics.
///
/// `AssertUnwindSafe` is sound here because nothing observable survives a
/// panic: on the unwind path we drop every value the body created and return
/// the sentinel, so a caller can never see a half-built handle.
pub(crate) fn guard<T>(fallback: T, body: impl FnOnce() -> T) -> T {
    match catch_unwind(AssertUnwindSafe(body)) {
        Ok(value) => value,
        Err(payload) => {
            set_last_error(panic_message(payload.as_ref()));
            fallback
        }
    }
}

/// Extract a readable message from a panic payload. `panic!` with a literal
/// yields `&str`, formatted panics yield `String`; anything else is opaque.
fn panic_message(payload: &(dyn Any + Send)) -> String {
    if let Some(s) = payload.downcast_ref::<&str>() {
        format!("panic: {s}")
    } else if let Some(s) = payload.downcast_ref::<String>() {
        format!("panic: {s}")
    } else {
        "panic: unknown panic".to_string()
    }
}
