// One-off reproduction driver for Boolean Gallery debug records: replays a
// logged frame (two STL files, op, engine, rotation, offset) through the
// same import pipeline as the demo. Usage:
//
//   cargo run --release --example robust_repro -- a.stl b.stl \
//       [op] [engine] [rx ry rz] [ox oy oz]
//
// op: union|intersection|difference (default union)
// engine: exact|robust|auto (default robust)
//
// REPAIR=1 in the environment applies `Manifold::repair_orientation` to both
// operands after import (rewinding inside-out shells), so the effect of the
// repair on a logged frame can be compared against the same frame without it.
//
// WINDING=nonzero runs the boolean under the {w != 0} solid rule (inside-out
// geometry stays material) instead of the default {w >= 1}; REFEREE_RULE
// takes the same values and selects which rule the Monte-Carlo referee below
// measures the ground truth under.

use manifold_rust::linalg::Vec3;
use manifold_rust::manifold::Manifold;
use manifold_rust::types::{BooleanEngine, MeshGL, OpType, WindingRule};

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
    // WINDING=nonzero selects the {w != 0} solid rule (inside-out geometry
    // counts as material); anything else keeps the default {w >= 1}.
    let rule = if std::env::var("WINDING").is_ok_and(|v| v == "nonzero") {
        WindingRule::Nonzero
    } else {
        WindingRule::Positive
    };
    let repair = std::env::var("REPAIR").is_ok_and(|v| v == "1");
    let mut a0 = import_stl(&args[0]);
    let mut b0 = import_stl(&args[1]);
    if repair {
        a0 = a0.repair_orientation();
        b0 = b0.repair_orientation();
        println!("repair_orientation applied to both operands");
    }
    let both_solid = !a0.as_impl().is_soup && !b0.as_impl().is_soup;
    let (a, b) = if both_solid {
        (colorize(a0, 0.27, 0.53, 0.80, 1.0), colorize(b0, 0.85, 0.25, 0.25, 0.6))
    } else {
        (a0, b0)
    };
    let b = b.rotate(rx, ry, rz).translate(Vec3::new(ox, oy, oz));
    println!(
        "A: {} tris (status {:?}), B: {} tris (status {:?}), op {:?}, engine {:?}, winding {:?}, rot ({rx},{ry},{rz}), off ({ox},{oy},{oz})",
        a.num_tri(), a.status(), b.num_tri(), b.status(), op, engine, rule,
    );
    let t0 = std::time::Instant::now();
    let out = a.boolean_with_engine_and_rule(&b, op, engine, rule);
    println!(
        "result: {} tris / {} verts, volume {:.7}, area {:.7}, status {:?}, {:.0} ms",
        out.num_tri(), out.num_vert(), out.volume(), out.surface_area(), out.status(),
        t0.elapsed().as_secs_f64() * 1e3,
    );

    // REFEREE_SAMPLES=N arbitrates the union volume independently of both
    // engines: stratified Monte Carlo over the joint bounding box with the
    // exact winding query, plus the arrangement's inconsistent-wall count
    // (same approach as examples/volume_referee.rs, generalized to a
    // two-mesh frame). Union only.
    let samples: usize = std::env::var("REFEREE_SAMPLES")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(0);
    if samples > 0 && op == OpType::Add {
        use manifold_rust::robust::exact::rational::R3;
        use manifold_rust::robust::ray_shoot::{winding_number_indexed, WindingIndex};
        use manifold_rust::robust::{cells, intersection_graph, soup};
        let p = soup::impl_to_tris(a.as_impl());
        let q = soup::impl_to_tris(b.as_impl());
        let (mut min, mut max) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
        for t in p.iter().chain(q.iter()) {
            for v in t {
                let c = [v.x, v.y, v.z];
                for k in 0..3 {
                    min[k] = min[k].min(c[k]);
                    max[k] = max[k].max(c[k]);
                }
            }
        }
        let ext = [max[0] - min[0], max[1] - min[1], max[2] - min[2]];
        let idx_p = WindingIndex::new(&p);
        let idx_q = WindingIndex::new(&q);
        // REFEREE_RULE=nonzero measures the ground truth under {w != 0} so the
        // referee can arbitrate either winding rule; default is {w >= 1}.
        let ref_rule = if std::env::var("REFEREE_RULE").is_ok_and(|v| v == "nonzero") {
            WindingRule::Nonzero
        } else {
            WindingRule::Positive
        };
        let solid = move |w: i32| match ref_rule {
            WindingRule::Positive => w >= 1,
            WindingRule::Nonzero => w != 0,
        };
        // SplitMix64, fixed seed.
        let mut state = 0x5EED_CAFE_F00D_D1CEu64;
        let mut rand = move || {
            state = state.wrapping_add(0x9E3779B97F4A7C15);
            let mut z = state;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
            ((z ^ (z >> 31)) >> 11) as f64 / (1u64 << 53) as f64
        };
        let mut hits = 0usize;
        for _ in 0..samples {
            let pt = R3::from_vec3(Vec3::new(
                min[0] + ext[0] * rand(),
                min[1] + ext[1] * rand(),
                min[2] + ext[2] * rand(),
            ));
            if solid(winding_number_indexed(&pt, &p, &idx_p))
                || solid(winding_number_indexed(&pt, &q, &idx_q))
            {
                hits += 1;
            }
        }
        let vol_box = ext[0] * ext[1] * ext[2];
        let frac = hits as f64 / samples as f64;
        let (est, sigma) = (
            vol_box * frac,
            vol_box * (frac * (1.0 - frac) / samples as f64).sqrt(),
        );
        let bad = intersection_graph::build_graph_with_token(&p, &q, None)
            .map(|graph| {
                let complex = cells::build_cells(&graph);
                let wind = cells::windings(&graph, &complex, [&p, &q]);
                cells::inconsistent_walls(&complex, &wind).len()
            })
            .unwrap_or(usize::MAX);
        println!("referee ({ref_rule:?}): volume {est:.7} ± {sigma:.7} ({samples} samples), inconsistent walls {bad}");
        // Per-operand self-overlap: divergence-theorem volume (counts a
        // doubly-wound region twice) over sampled once-counted material.
        // A ratio away from 1.0 identifies the operand whose geometry
        // self-overlaps despite clean manifold topology.
        for (name, mani, tris, idx) in [("A", &a, &p, &idx_p), ("B", &b, &q, &idx_q)] {
            let (mut omin, mut omax) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
            for t in tris.iter() {
                for v in t {
                    let c = [v.x, v.y, v.z];
                    for k in 0..3 {
                        omin[k] = omin[k].min(c[k]);
                        omax[k] = omax[k].max(c[k]);
                    }
                }
            }
            let oext = [omax[0] - omin[0], omax[1] - omin[1], omax[2] - omin[2]];
            let n = samples / 2;
            let mut ohits = 0usize;
            for _ in 0..n {
                let pt = R3::from_vec3(Vec3::new(
                    omin[0] + oext[0] * rand(),
                    omin[1] + oext[1] * rand(),
                    omin[2] + oext[2] * rand(),
                ));
                if solid(winding_number_indexed(&pt, tris, idx)) {
                    ohits += 1;
                }
            }
            let sampled = oext[0] * oext[1] * oext[2] * ohits as f64 / n as f64;
            println!(
                "operand {name} ({ref_rule:?}): divergence volume {:.7}, sampled solid {:.7}, overlap ratio {:.3}",
                mani.volume(),
                sampled,
                mani.volume() / sampled
            );
        }
    }
}
