// Rust port of cpp-reference/manifold/extras/perf_test.cpp.
//
// Benchmarks a single sphere-minus-sphere boolean at doubling tessellation
// levels, printing the same "nTri = N, time = S sec" lines as the C++ driver
// so the two outputs can be diffed side by side. An optional CLI argument
// caps the number of doubling rounds (C++ hard-codes 8).
//
// Run with: cargo run --release --example perf_test [rounds]

use std::time::Instant;

use manifold_rust::linalg::Vec3;
use manifold_rust::manifold::Manifold;

fn main() {
    let rounds: usize = std::env::args()
        .nth(1)
        .and_then(|a| a.parse().ok())
        .unwrap_or(8);

    for i in 0..rounds {
        let sphere = Manifold::sphere(1.0, (8 << i) * 4);
        let sphere2 = sphere.translate(Vec3::splat(0.5));
        let start = Instant::now();
        let diff = sphere.difference(&sphere2);
        let n_tri_result = diff.num_tri();
        let elapsed = start.elapsed().as_secs_f64();
        // Keep result live so the boolean cannot be optimized away.
        std::hint::black_box(n_tri_result);
        println!("nTri = {}, time = {} sec", sphere.num_tri(), elapsed);
    }
}
