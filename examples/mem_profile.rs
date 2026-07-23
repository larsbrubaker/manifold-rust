// Heap profiling driver for the boolean pipeline: runs one sphere-minus-sphere
// boolean at the perf_test tessellation for round `i` (default 6 -> 2M input
// tris) under a counting global allocator. With MANIFOLD_TIMING set, each
// pipeline stage line also reports current heap and that stage's peak, via the
// timing::set_mem_hook bridge. Used to chase the peak-memory gap vs the C++
// reference (see PORTING_PLAN.md "Performance vs C++").
//
// Run with: MANIFOLD_TIMING=1 cargo run --release --example mem_profile [round]

use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::atomic::{AtomicUsize, Ordering};

use manifold_rust::linalg::Vec3;
use manifold_rust::manifold::Manifold;

static CURRENT: AtomicUsize = AtomicUsize::new(0);
static PEAK: AtomicUsize = AtomicUsize::new(0);

struct CountingAlloc;

// SAFETY: delegates all allocation to System; only maintains atomic counters.
unsafe impl GlobalAlloc for CountingAlloc {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        let p = System.alloc(layout);
        if !p.is_null() {
            let c = CURRENT.fetch_add(layout.size(), Ordering::Relaxed) + layout.size();
            PEAK.fetch_max(c, Ordering::Relaxed);
        }
        p
    }
    unsafe fn dealloc(&self, p: *mut u8, layout: Layout) {
        CURRENT.fetch_sub(layout.size(), Ordering::Relaxed);
        System.dealloc(p, layout);
    }
    unsafe fn realloc(&self, p: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        let q = System.realloc(p, layout, new_size);
        if !q.is_null() {
            if new_size >= layout.size() {
                let c = CURRENT.fetch_add(new_size - layout.size(), Ordering::Relaxed)
                    + (new_size - layout.size());
                PEAK.fetch_max(c, Ordering::Relaxed);
            } else {
                CURRENT.fetch_sub(layout.size() - new_size, Ordering::Relaxed);
            }
        }
        q
    }
}

#[global_allocator]
static ALLOC: CountingAlloc = CountingAlloc;

/// timing::MemHook — returns (current, peak since last call), resetting the
/// watermark so each stage reports its own peak.
fn mem_stats() -> (usize, usize) {
    let current = CURRENT.load(Ordering::Relaxed);
    let peak = PEAK.swap(current, Ordering::Relaxed);
    (current, peak)
}

fn mb(bytes: usize) -> f64 {
    bytes as f64 / 1048576.0
}

fn main() {
    let round: usize = std::env::args()
        .nth(1)
        .and_then(|a| a.parse().ok())
        .unwrap_or(6);

    manifold_rust::timing::set_mem_hook(mem_stats);

    let sphere = Manifold::sphere(1.0, ((8 << round) * 4) as i32);
    let sphere2 = sphere.translate(Vec3::splat(0.5));
    let (current, _) = mem_stats();
    eprintln!(
        "inputs ready: nTri = {} each, current = {:.1} MB",
        sphere.num_tri(),
        mb(current)
    );

    let diff = sphere.difference(&sphere2);

    let (current, peak) = mem_stats();
    eprintln!(
        "done: nTri(out) = {}, current = {:.1} MB, tail peak = {:.1} MB",
        diff.num_tri(),
        mb(current),
        mb(peak)
    );
}
