// Timing driver for keyholing (CutKeyhole / FindCloserBridge) in the ear
// clipper. `holes`: one square with an N x N grid of 32-gon holes. `glyphs`:
// R x C square outer rings, each with one or two octagonal holes. Prints the
// triangle count and an FNV-1a hash of the triangle list next to the time.
//
// Run with: cargo run --release --example keyhole_perf [holes N | glyphs R C] [repeats]

use std::time::Instant;

use manifold_rust::linalg::{IVec3, Vec2};
use manifold_rust::polygon::triangulate;
use manifold_rust::types::Polygons;

fn fnv(tris: &[IVec3]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for t in tris {
        for c in [t.x, t.y, t.z] {
            for b in c.to_le_bytes() {
                h ^= u64::from(b);
                h = h.wrapping_mul(0x0100_0000_01b3);
            }
        }
    }
    h
}

fn holes(n: usize) -> Polygons {
    let size = 3.0 * n as f64 + 1.0;
    let mut polys = vec![vec![
        Vec2::new(0.0, 0.0),
        Vec2::new(size, 0.0),
        Vec2::new(size, size),
        Vec2::new(0.0, size),
    ]];
    for i in 0..n {
        for j in 0..n {
            let (cx, cy) = (3.0 * i as f64 + 2.0, 3.0 * j as f64 + 2.0);
            // Clockwise, so each ring is a hole.
            polys.push(
                (0..32)
                    .map(|k| {
                        let a = -(k as f64) * std::f64::consts::TAU / 32.0;
                        Vec2::new(cx + a.cos(), cy + a.sin())
                    })
                    .collect(),
            );
        }
    }
    polys
}

fn glyphs(rows: usize, cols: usize) -> Polygons {
    let octagon = |cx: f64, cy: f64, r: f64| -> Vec<Vec2> {
        [
            (1.0, -0.5),
            (0.5, -1.0),
            (-0.5, -1.0),
            (-1.0, -0.5),
            (-1.0, 0.5),
            (-0.5, 1.0),
            (0.5, 1.0),
            (1.0, 0.5),
        ]
        .iter()
        .map(|&(x, y)| Vec2::new(cx + r * x, cy + r * y))
        .collect()
    };
    let mut polys = Vec::new();
    for row in 0..rows {
        for col in 0..cols {
            let x0 = 5.0 * col as f64;
            let y0 = 7.0 * row as f64;
            polys.push(vec![
                Vec2::new(x0, y0),
                Vec2::new(x0 + 4.0, y0),
                Vec2::new(x0 + 4.0, y0 + 6.0),
                Vec2::new(x0, y0 + 6.0),
            ]);
            if (row + col) % 2 == 0 {
                polys.push(octagon(x0 + 2.0, y0 + 3.0, 1.0));
            } else {
                polys.push(octagon(x0 + 2.0, y0 + 1.5, 0.75));
                polys.push(octagon(x0 + 2.0, y0 + 4.5, 0.75));
            }
        }
    }
    polys
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let arg = |i: usize, default: usize| -> usize {
        args.get(i).and_then(|a| a.parse().ok()).unwrap_or(default)
    };
    let (name, polys, repeats) = match args.first().map(String::as_str) {
        Some("glyphs") => ("glyphs", glyphs(arg(1, 100), arg(2, 100)), arg(3, 3)),
        _ => ("holes", holes(arg(1, 71)), arg(2, 3)),
    };
    let verts: usize = polys.iter().map(Vec::len).sum();
    let mut best = f64::INFINITY;
    let mut tris = Vec::new();
    for _ in 0..repeats.max(1) {
        let start = Instant::now();
        tris = triangulate(&polys, -1.0, true);
        best = best.min(start.elapsed().as_secs_f64());
    }
    println!(
        "{name}: {} rings, {verts} verts -> {} tris, hash {:#018x}, best of {} {:.3} s",
        polys.len(),
        tris.len(),
        fnv(&tris),
        repeats.max(1),
        best
    );
}
