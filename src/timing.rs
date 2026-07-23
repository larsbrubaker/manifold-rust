// Lightweight, env-gated phase timing for the boolean pipeline, mirroring the
// C++ MANIFOLD_TIMING instrumentation (Timer::Print in boolean3.cpp /
// boolean_result.cpp). Enabled by setting the MANIFOLD_TIMING environment
// variable to any non-empty value; otherwise every call is a no-op so release
// performance is unaffected. Used to compare per-stage wall-clock against the
// C++ reference when hunting performance gaps (see CLAUDE.md "Instrumentation
// Strategy").

use std::sync::OnceLock;
use std::time::Instant;

fn enabled() -> bool {
    static ENABLED: OnceLock<bool> = OnceLock::new();
    *ENABLED.get_or_init(|| {
        std::env::var("MANIFOLD_TIMING").map_or(false, |v| !v.is_empty())
    })
}

/// Start a stage timer. Returns None (and times nothing) unless the
/// MANIFOLD_TIMING environment variable is set.
pub(crate) fn start() -> Option<Instant> {
    if enabled() {
        Some(Instant::now())
    } else {
        None
    }
}

/// Print the elapsed time for a stage started with `start`, matching the C++
/// Timer::Print format ("label: N sec") on stderr.
pub(crate) fn print(label: &str, t0: Option<Instant>) {
    if let Some(t0) = t0 {
        eprintln!("{}: {} sec", label, t0.elapsed().as_secs_f64());
    }
}
