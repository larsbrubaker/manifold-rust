// Progress bridge between the core crate's ProgressReporter and a JS callback.
//
// The core reporter's callback must be `Send + Sync` (rayon workers drive it
// in native parallel builds), but a `js_sys::Function` is neither — and this
// crate is compiled for the host as well as for wasm32, so a blanket
// `unsafe impl Send` would be a lie on the host.
//
// The way out needs no unsafe at all: park the JS function in a thread-local
// for the duration of the call and hand the reporter a closure that captures
// *nothing* (so it is trivially `Send + Sync`) and looks the function up on
// whatever thread it runs on. On wasm32 there is exactly one thread, so the
// lookup always finds it. On any other build a foreign worker thread simply
// finds no function and drops the update — the only correct thing to do with
// a JS value it could not touch anyway.
//
// The JS side (demo/src/boolean-worker.ts) forwards these to the main thread
// with postMessage, which is exactly why this must work from *inside* a
// blocking wasm call: the message queues while the worker is still busy.

use std::cell::RefCell;

use manifold_rust::progress::{Phase, ProgressReporter};
use wasm_bindgen::prelude::*;

thread_local! {
    /// The JS callback for the boolean currently running on this thread.
    static SINK: RefCell<Option<js_sys::Function>> = const { RefCell::new(None) };
}

/// Installs `f` as this thread's progress sink and removes it on drop, so an
/// early return or a panic cannot leave a stale function behind.
struct SinkGuard;

impl SinkGuard {
    fn install(f: js_sys::Function) -> Self {
        SINK.with(|s| *s.borrow_mut() = Some(f));
        SinkGuard
    }
}

impl Drop for SinkGuard {
    fn drop(&mut self) {
        SINK.with(|s| *s.borrow_mut() = None);
    }
}

/// Call `body` with a reporter that forwards to `callback`, if one was given.
///
/// `callback` is invoked as `callback(phaseName: string, fraction: number |
/// null)` — the same shape `BusyIndicator.setPhase` takes. A callback that
/// throws is ignored: a broken progress sink must not fail the geometry
/// operation.
pub(crate) fn with_reporter<T>(
    callback: Option<js_sys::Function>,
    body: impl FnOnce(Option<&ProgressReporter>) -> T,
) -> T {
    let Some(callback) = callback else {
        return body(None);
    };
    let _guard = SinkGuard::install(callback);
    // Captures nothing, so it satisfies the reporter's `Send + Sync` bound
    // without any unsafe assertion about JsValue.
    let reporter = ProgressReporter::new(|phase: Phase, fraction| {
        // Clone the function out before calling it: holding the RefCell borrow
        // across a call into JS would turn any re-entrancy into a panic.
        let f = SINK.with(|s| s.borrow().clone());
        if let Some(f) = f {
            let frac = match fraction {
                Some(v) => JsValue::from_f64(v),
                None => JsValue::NULL,
            };
            let _ = f.call2(&JsValue::NULL, &JsValue::from_str(phase.name()), &frac);
        }
    });
    body(Some(&reporter))
}
