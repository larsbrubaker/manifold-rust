// Timing driver for large booleans. `menger N`: the C++ sample
// MengerSponge(n) (samples/src/menger_sponge.cpp). `drill S`: a sphere of S
// segments minus the union of 49 tilted cylinders. Prints the union and
// difference times, the triangle count and an FNV-1a hash of the MeshGL64
// (every field but `run_original_id`, as in src/par_tests.rs), to compare
// builds and thread counts.
//
// Run with: cargo run --release [--features parallel] --example boolean_perf [menger N | drill S] [repeats]

use std::time::Instant;

use manifold_rust::linalg::{Vec2, Vec3};
use manifold_rust::manifold::Manifold;
use manifold_rust::types::OpType;

fn fractal(holes: &mut Vec<Manifold>, hole: &Manifold, w: f64, p: Vec2, depth: i32, max: i32) {
    let w = w / 3.0;
    holes.push(
        hole.scale(Vec3::new(w, w, 1.0))
            .translate(Vec3::new(p.x, p.y, 0.0)),
    );
    if depth == max {
        return;
    }
    let offsets = [
        (-w, -w),
        (-w, 0.0),
        (-w, w),
        (0.0, w),
        (w, w),
        (w, 0.0),
        (w, -w),
        (0.0, -w),
    ];
    for (x, y) in offsets {
        fractal(holes, hole, w, Vec2::new(p.x + x, p.y + y), depth + 1, max);
    }
}

/// MengerSponge(n), with the times of the hole union and the differences.
fn menger(n: i32) -> (Manifold, f64, f64) {
    let cube = Manifold::cube(Vec3::splat(1.0), true);
    let mut holes = Vec::new();
    fractal(&mut holes, &cube, 1.0, Vec2::new(0.0, 0.0), 1, n);
    let start = Instant::now();
    let hole = Manifold::batch_boolean(&holes, OpType::Add);
    let union = start.elapsed().as_secs_f64();
    let start = Instant::now();
    let r = cube.difference(&hole);
    let hole = hole.rotate(90.0, 0.0, 0.0);
    let r = r.difference(&hole);
    let hole = hole.rotate(0.0, 0.0, 90.0);
    let r = r.difference(&hole);
    (r, union, start.elapsed().as_secs_f64())
}

/// A sphere of `segments` segments minus 49 tilted cylinders.
fn drill(segments: i32) -> (Manifold, f64, f64) {
    let mut cutters = Vec::new();
    for i in -3..=3 {
        for j in -3..=3 {
            cutters.push(
                Manifold::cylinder_centered(50.0, 2.0, 2.0, 64, true)
                    .rotate(f64::from(i) * 9.0, f64::from(j) * 7.0, 0.0)
                    .translate(Vec3::new(f64::from(i) * 5.5, f64::from(j) * 5.5, 0.0)),
            );
        }
    }
    let sphere = Manifold::sphere(20.0, segments);
    let start = Instant::now();
    let cutters = Manifold::batch_boolean(&cutters, OpType::Add);
    let union = start.elapsed().as_secs_f64();
    let start = Instant::now();
    let r = sphere.difference(&cutters);
    (r, union, start.elapsed().as_secs_f64())
}

/// FNV-1a (64-bit) over every field of the output `MeshGL64` but
/// `run_original_id`, which comes from the process-wide ID counter. In struct
/// order: `num_prop`; then `vert_properties`, `tri_verts`, `merge_from_vert`,
/// `merge_to_vert`, `run_index`, `run_transform`, `face_id`,
/// `halfedge_tangent` and `run_flags`, each as its length then its elements;
/// then `tolerance`. Lengths and integers are u64 little-endian, floats their
/// IEEE bits as u64 little-endian, `run_flags` one byte each.
/// `src/par_tests.rs` hashes the same way.
fn fingerprint(m: &Manifold) -> u64 {
    let gl = m.get_mesh_gl64(-1);
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    let mut eat = |bytes: &[u8]| {
        for &b in bytes {
            h ^= u64::from(b);
            h = h.wrapping_mul(0x0100_0000_01b3);
        }
    };
    let ints = |eat: &mut dyn FnMut(&[u8]), list: &[u64]| {
        eat(&(list.len() as u64).to_le_bytes());
        for x in list {
            eat(&x.to_le_bytes());
        }
    };
    let floats = |eat: &mut dyn FnMut(&[u8]), list: &[f64]| {
        eat(&(list.len() as u64).to_le_bytes());
        for x in list {
            eat(&x.to_bits().to_le_bytes());
        }
    };
    eat(&gl.num_prop.to_le_bytes());
    floats(&mut eat, &gl.vert_properties);
    ints(&mut eat, &gl.tri_verts);
    ints(&mut eat, &gl.merge_from_vert);
    ints(&mut eat, &gl.merge_to_vert);
    ints(&mut eat, &gl.run_index);
    floats(&mut eat, &gl.run_transform);
    ints(&mut eat, &gl.face_id);
    floats(&mut eat, &gl.halfedge_tangent);
    eat(&(gl.run_flags.len() as u64).to_le_bytes());
    eat(&gl.run_flags);
    eat(&gl.tolerance.to_bits().to_le_bytes());
    h
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let arg = |i: usize, default: usize| -> usize {
        args.get(i).and_then(|a| a.parse().ok()).unwrap_or(default)
    };
    let model = args.first().map(String::as_str).unwrap_or("menger");
    let (size, repeats) = match model {
        "drill" => (arg(1, 224), arg(2, 3).max(1)),
        // MengerSponge(0) would never reach its recursion's base case.
        _ => (arg(1, 4).max(1), arg(2, 3).max(1)),
    };
    let (mut best_union, mut best_diff) = (f64::INFINITY, f64::INFINITY);
    let mut summary = String::new();
    for _ in 0..repeats {
        // Drop each result before the next run, so peak memory is one run's.
        let (result, union, diff) = match model {
            "drill" => drill(size as i32),
            _ => menger(size as i32),
        };
        best_union = best_union.min(union);
        best_diff = best_diff.min(diff);
        summary = format!(
            "{} tris, hash {:#018x}",
            result.num_tri(),
            fingerprint(&result)
        );
    }
    println!(
        "{model}({size}): {summary}, best of {repeats}: union {best_union:.3} s, differences {best_diff:.3} s"
    );
}
