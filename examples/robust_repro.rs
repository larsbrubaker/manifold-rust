// One-off reproduction driver for Boolean Gallery debug records: replays a
// logged frame (two STL files, op, engine, rotation, offset) through the
// same import pipeline as the demo. Usage:
//
//   cargo run --release --example robust_repro -- a.stl b.stl \
//       [op] [engine] [rx ry rz] [ox oy oz]
//
// op: union|intersection|difference (default union)
// engine: exact|robust|auto (default robust)

use manifold_rust::linalg::Vec3;
use manifold_rust::manifold::Manifold;
use manifold_rust::types::{BooleanEngine, MeshGL, OpType};

fn import_stl(path: &str) -> Manifold {
    let data = std::fs::read(path).expect("STL file");
    let head = String::from_utf8_lossy(&data[..data.len().min(512)]).to_string();
    let mut positions: Vec<f32> = Vec::new();
    if head.trim_start().starts_with("solid") && head.contains("facet") {
        for line in String::from_utf8_lossy(&data).lines() {
            if let Some(rest) = line.trim_start().strip_prefix("vertex") {
                for tok in rest.split_whitespace().take(3) {
                    positions.push(tok.parse::<f64>().unwrap() as f32);
                }
            }
        }
    } else {
        let n = u32::from_le_bytes(data[80..84].try_into().unwrap()) as usize;
        let n = n.min((data.len() - 84) / 50);
        for f in 0..n {
            let base = 84 + f * 50 + 12;
            for v in 0..9 {
                let o = base + v * 4;
                positions.push(f32::from_le_bytes(data[o..o + 4].try_into().unwrap()));
            }
        }
    }
    let nv = positions.len() / 3;
    let (mut min, mut max) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for i in 0..nv {
        for k in 0..3 {
            let v = positions[i * 3 + k] as f64;
            min[k] = min[k].min(v);
            max[k] = max[k].max(v);
        }
    }
    let c = [(min[0] + max[0]) / 2.0, (min[1] + max[1]) / 2.0, (min[2] + max[2]) / 2.0];
    let side = (max[0] - min[0]).max(max[1] - min[1]).max(max[2] - min[2]);
    let s = if side > 0.0 { 2.0 / side } else { 1.0 };
    for i in 0..nv {
        for k in 0..3 {
            positions[i * 3 + k] = ((positions[i * 3 + k] as f64 - c[k]) * s) as f32;
        }
    }
    let mut mesh = MeshGL::default();
    mesh.num_prop = 3;
    mesh.tri_verts = (0..nv as u32).collect();
    mesh.vert_properties = positions;
    mesh.merge();
    Manifold::from_mesh_gl_robust(&mesh)
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.len() < 2 {
        eprintln!("usage: robust_repro a.stl b.stl [op] [engine] [rx ry rz] [ox oy oz]");
        std::process::exit(2);
    }
    let op = match args.get(2).map(String::as_str).unwrap_or("union") {
        "intersection" => OpType::Intersect,
        "difference" => OpType::Subtract,
        _ => OpType::Add,
    };
    let engine = match args.get(3).map(String::as_str).unwrap_or("robust") {
        "exact" => BooleanEngine::Exact,
        "auto" => BooleanEngine::Auto,
        _ => BooleanEngine::Robust,
    };
    let num = |i: usize| args.get(i).and_then(|a| a.parse::<f64>().ok()).unwrap_or(0.0);
    let (rx, ry, rz) = (num(4), num(5), num(6));
    let (ox, oy, oz) = (num(7), num(8), num(9));

    // The demo colorizes non-soup operands (4 RGBA property channels) before
    // the boolean, so property interpolation is part of the reproduced path.
    let colorize = |m: Manifold, r: f64, g: f64, b: f64, al: f64| -> Manifold {
        if m.as_impl().is_soup {
            return m;
        }
        m.set_properties(4, move |new_prop: &mut [f64], _pos, _old: &[f64]| {
            new_prop[0] = r;
            new_prop[1] = g;
            new_prop[2] = b;
            new_prop[3] = al;
        })
    };
    let a0 = import_stl(&args[0]);
    let b0 = import_stl(&args[1]);
    let both_solid = !a0.as_impl().is_soup && !b0.as_impl().is_soup;
    let (a, b) = if both_solid {
        (colorize(a0, 0.27, 0.53, 0.80, 1.0), colorize(b0, 0.85, 0.25, 0.25, 0.6))
    } else {
        (a0, b0)
    };
    let b = b.rotate(rx, ry, rz).translate(Vec3::new(ox, oy, oz));
    println!(
        "A: {} tris (status {:?}), B: {} tris (status {:?}), op {:?}, engine {:?}, rot ({rx},{ry},{rz}), off ({ox},{oy},{oz})",
        a.num_tri(), a.status(), b.num_tri(), b.status(), op, engine,
    );
    let t0 = std::time::Instant::now();
    let out = a.boolean_with_engine(&b, op, engine);
    println!(
        "result: {} tris / {} verts, volume {:.7}, area {:.7}, status {:?}, {:.0} ms",
        out.num_tri(), out.num_vert(), out.volume(), out.surface_area(), out.status(),
        t0.elapsed().as_secs_f64() * 1e3,
    );
}
