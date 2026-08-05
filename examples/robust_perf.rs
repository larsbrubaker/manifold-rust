// Stage-by-stage timing driver for the robust boolean engine.
//
// Runs the demo's spiky-dodecahedron-vs-spiky-dodecahedron boolean (the
// Boolean Gallery default) through both engines and prints per-stage wall
// times for the robust pipeline, so optimization work has a baseline and a
// regression check. Usage:
//
//   cargo run --release --example robust_perf [-- <rot_x> <rot_y> <rot_z>]
//
// The rotation (degrees, default 30/45/60) matches an animated gallery frame.

use std::time::Instant;

use manifold_rust::linalg::Vec3;
use manifold_rust::manifold::Manifold;
use manifold_rust::robust::{classify, intersection_graph, propagate, soup};
use manifold_rust::types::{BooleanEngine, MeshGL};

fn make_spiky_dodecahedron(spike_height: f64) -> Manifold {
    let phi: f64 = (1.0 + 5.0_f64.sqrt()) / 2.0;
    let inv_phi = 1.0 / phi;
    let scale = 0.5;
    let raw_verts: [(f64, f64, f64); 20] = [
        (1.0, 1.0, 1.0), (1.0, 1.0, -1.0), (1.0, -1.0, 1.0), (1.0, -1.0, -1.0),
        (-1.0, 1.0, 1.0), (-1.0, 1.0, -1.0), (-1.0, -1.0, 1.0), (-1.0, -1.0, -1.0),
        (0.0, inv_phi, phi), (0.0, inv_phi, -phi), (0.0, -inv_phi, phi), (0.0, -inv_phi, -phi),
        (inv_phi, phi, 0.0), (-inv_phi, phi, 0.0), (inv_phi, -phi, 0.0), (-inv_phi, -phi, 0.0),
        (phi, 0.0, inv_phi), (phi, 0.0, -inv_phi), (-phi, 0.0, inv_phi), (-phi, 0.0, -inv_phi),
    ];
    let faces: [[usize; 5]; 12] = [
        [0, 8, 10, 2, 16], [0, 16, 17, 1, 12], [0, 12, 13, 4, 8],
        [1, 17, 3, 11, 9], [1, 9, 5, 13, 12], [2, 10, 6, 15, 14],
        [2, 14, 3, 17, 16], [4, 13, 5, 19, 18], [4, 18, 6, 10, 8],
        [5, 9, 11, 7, 19], [6, 18, 19, 7, 15], [3, 14, 15, 7, 11],
    ];
    let verts: Vec<(f64, f64, f64)> =
        raw_verts.iter().map(|&(x, y, z)| (x * scale, y * scale, z * scale)).collect();
    let mut positions: Vec<f32> = Vec::with_capacity(32 * 3);
    let mut tri_verts: Vec<u32> = Vec::with_capacity(60 * 3);
    for &(x, y, z) in &verts {
        positions.extend([x as f32, y as f32, z as f32]);
    }
    for face in &faces {
        let cx: f64 = face.iter().map(|&i| verts[i].0).sum::<f64>() / 5.0;
        let cy: f64 = face.iter().map(|&i| verts[i].1).sum::<f64>() / 5.0;
        let cz: f64 = face.iter().map(|&i| verts[i].2).sum::<f64>() / 5.0;
        let len = (cx * cx + cy * cy + cz * cz).sqrt();
        let (nx, ny, nz) = (cx / len, cy / len, cz / len);
        let spike_idx = (positions.len() / 3) as u32;
        positions.extend([
            (cx + nx * spike_height) as f32,
            (cy + ny * spike_height) as f32,
            (cz + nz * spike_height) as f32,
        ]);
        for j in 0..5 {
            tri_verts.extend([spike_idx, face[j] as u32, face[(j + 1) % 5] as u32]);
        }
    }
    let mut mesh = MeshGL::default();
    mesh.num_prop = 3;
    mesh.vert_properties = positions;
    mesh.tri_verts = tri_verts;
    Manifold::from_mesh_gl(&mesh)
}

fn main() {
    let args: Vec<f64> = std::env::args()
        .skip(1)
        .filter_map(|a| a.parse().ok())
        .collect();
    let (rx, ry, rz) = if args.len() == 3 {
        (args[0], args[1], args[2])
    } else {
        (30.0, 45.0, 60.0)
    };

    let a = make_spiky_dodecahedron(0.4);
    let b = make_spiky_dodecahedron(0.4)
        .rotate(rx, ry, rz)
        .translate(Vec3::new(0.3, 0.0, 0.0));
    println!(
        "spiky vs spiky: {} + {} tris, rot=({rx},{ry},{rz})",
        a.num_tri(),
        b.num_tri()
    );

    // Whole-op comparison, best of 3.
    for (name, engine) in [("exact ", BooleanEngine::Exact), ("robust", BooleanEngine::Robust)] {
        let mut best = f64::INFINITY;
        let mut out_tris = 0;
        for _ in 0..3 {
            let t0 = Instant::now();
            let out = a.union_with_engine(&b, engine);
            let dt = t0.elapsed().as_secs_f64();
            best = best.min(dt);
            out_tris = out.num_tri();
        }
        println!("{name} union: {:8.2} ms   ({out_tris} tris out)", best * 1e3);
    }

    // Dense-mesh case (sparse intersections relative to size): sphere pair.
    // This is where the BVH broad phase matters — brute-force box pruning is
    // O(|P|·|Q|) (e.g. 8k×8k tris = 68M box tests) while the BVH visits only
    // the overlapping region.
    let s1 = Manifold::sphere(1.0, 64);
    let s2 = Manifold::sphere(1.0, 64).translate(Vec3::new(1.7, 0.0, 0.0));
    for (name, engine) in [("exact ", BooleanEngine::Exact), ("robust", BooleanEngine::Robust)] {
        let t0 = Instant::now();
        let out = s1.union_with_engine(&s2, engine);
        println!(
            "{name} sphere64 union: {:8.2} ms   ({} tris in, {} out)",
            t0.elapsed().as_secs_f64() * 1e3,
            s1.num_tri() + s2.num_tri(),
            out.num_tri()
        );
    }

    // Robust pipeline stage timings (single run; stages mirror robust::boolean).
    let p_tris = soup::impl_to_tris(a.as_impl());
    let q_tris = soup::impl_to_tris(b.as_impl());

    let t0 = Instant::now();
    let graph = intersection_graph::build_graph(&p_tris, &q_tris);
    let t_graph = t0.elapsed().as_secs_f64();

    let t0 = Instant::now();
    let cls = classify::classify_rings(&graph);
    let t_classify = t0.elapsed().as_secs_f64();

    let t0 = Instant::now();
    let prop = propagate::propagate(&graph, &cls);
    let t_prop = t0.elapsed().as_secs_f64();

    let t0 = Instant::now();
    let out = manifold_rust::robust::assemble::assemble(
        &graph.pieces,
        |pi| !cls.discarded[pi] && prop.tags[pi] == Some(classify::Tag::Union),
        None,
    );
    let t_assemble = t0.elapsed().as_secs_f64();
    let _ = out;

    println!("robust stages:");
    println!("  build_graph    {:8.2} ms   ({} pieces)", t_graph * 1e3, graph.pieces.len());
    println!("  classify_rings {:8.2} ms", t_classify * 1e3);
    println!("  propagate      {:8.2} ms   ({} untagged roots)", t_prop * 1e3, prop.untagged.len());
    println!("  assemble       {:8.2} ms", t_assemble * 1e3);
}
