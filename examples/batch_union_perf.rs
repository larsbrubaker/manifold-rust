// Timing driver for batch unions through the CSG tree (`CsgNode::op_n`): N
// overlapping spheres of S segments in a 3D grid. Prints the time, triangle
// count and an FNV-1a hash of the MeshGL64 (less `run_original_id`, which
// comes from the process-wide ID counter), to compare builds and threads.
//
// Run with: cargo run --release [--features parallel] --example batch_union_perf [N S] [repeats]

use std::time::Instant;

use manifold_rust::csg_tree::CsgNode;
use manifold_rust::linalg::Vec3;
use manifold_rust::manifold::Manifold;
use manifold_rust::types::OpType;

fn fingerprint(m: &Manifold) -> u64 {
    let gl = m.get_mesh_gl64(-1);
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    let mut eat = |bytes: &[u8]| {
        for &b in bytes {
            h ^= u64::from(b);
            h = h.wrapping_mul(0x0100_0000_01b3);
        }
    };
    for x in &gl.vert_properties {
        eat(&x.to_bits().to_le_bytes());
    }
    for x in gl.tri_verts.iter().chain(&gl.run_index) {
        eat(&x.to_le_bytes());
    }
    h
}

fn main() {
    let arg = |i: usize, default: usize| -> usize {
        std::env::args()
            .nth(i)
            .and_then(|a| a.parse().ok())
            .unwrap_or(default)
    };
    let (n, segments, repeats) = (arg(1, 125), arg(2, 32) as i32, arg(3, 3).max(1));
    let sphere = Manifold::sphere(1.0, segments);
    let side = (n as f64).cbrt().ceil() as usize;
    let mut best = f64::INFINITY;
    let mut summary = String::new();
    for _ in 0..repeats {
        let leaves: Vec<CsgNode> = (0..n)
            .map(|i| {
                let (x, y, z) = (i % side, (i / side) % side, i / (side * side));
                let at = Vec3::new(x as f64, y as f64, z as f64) * 1.5;
                CsgNode::leaf(sphere.translate(at).as_impl().clone())
            })
            .collect();
        let start = Instant::now();
        let result = Manifold::from_impl(CsgNode::op_n(OpType::Add, leaves).evaluate());
        best = best.min(start.elapsed().as_secs_f64());
        summary = format!(
            "{} tris, hash {:#018x}",
            result.num_tri(),
            fingerprint(&result)
        );
    }
    println!("{n} spheres of {segments} segments: {summary}, best of {repeats} {best:.3} s");
}
