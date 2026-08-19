// Third-opinion arbiter for exact-vs-robust volume disagreements found by
// examples/thingi_sweep.rs.
//
// Both engines claim to compute the same solid: for A ∪ B, the region
// {w_A ≥ 1} ∪ {w_B ≥ 1} of the *input* soups. This referee measures that
// region directly — stratified Monte Carlo point sampling over the joint
// bounding box, with each point classified by the exact rational winding
// query in robust::ray_shoot — touching neither engine's boolean pipeline.
// The estimate carries a standard error, so a disagreement of tens of
// percent is arbitrated decisively with modest sample counts.
//
// Per mesh it also reports the robust arrangement's self-consistency
// (cells::inconsistent_walls) so "robust is wrong" can be split into
// "arrangement is broken" vs "arrangement fine, classification wrong".
//
// Usage:
//   cargo run --release --example volume_referee -- --ids A,B,C [flags]
//   cargo run --release --example volume_referee -- --from-run N [flags]
//     --root DIR      mesh repo root (default C:\Development\rust-apps\Thingi10K\meshes)
//     --db FILE       sqlite db for --from-run (default thingi_sweep.db)
//     --from-run N    referee every match_ok=0 mesh of sweep run N
//     --limit K       cap the --from-run list
//     --samples N     Monte Carlo points (default 100000)
//
// The pass replicates thingi_sweep exactly: import (normalize to side 2,
// weld, robust import), B = A rotated (30,45,60)° and offset (0.3,0,0),
// operation A ∪ B.

use std::io::Read;
use std::panic::{catch_unwind, AssertUnwindSafe};

use manifold_rust::linalg::Vec3;
use manifold_rust::manifold::Manifold;
use manifold_rust::robust::exact::rational::R3;
use manifold_rust::robust::ray_shoot::{winding_number_indexed, WindingIndex};
use manifold_rust::robust::{cells, intersection_graph, soup};
use manifold_rust::types::{BooleanEngine, Error, MeshGL, OpType};

struct Args {
    root: String,
    db: String,
    ids: Vec<u64>,
    from_run: Option<i64>,
    limit: usize,
    samples: usize,
}

fn parse_args() -> Args {
    let mut args = Args {
        root: r"C:\Development\rust-apps\Thingi10K\meshes".to_string(),
        db: "thingi_sweep.db".to_string(),
        ids: Vec::new(),
        from_run: None,
        limit: usize::MAX,
        samples: 100_000,
    };
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    while i < argv.len() {
        let val = |i: &mut usize| -> String {
            *i += 1;
            argv.get(*i).cloned().unwrap_or_default()
        };
        match argv[i].as_str() {
            "--root" => args.root = val(&mut i),
            "--db" => args.db = val(&mut i),
            "--from-run" => args.from_run = val(&mut i).parse().ok(),
            "--limit" => args.limit = val(&mut i).parse().unwrap_or(usize::MAX),
            "--samples" => args.samples = val(&mut i).parse().unwrap_or(100_000),
            "--ids" => {
                args.ids = val(&mut i)
                    .split(',')
                    .filter_map(|s| s.trim().parse().ok())
                    .collect()
            }
            other => {
                eprintln!("unknown flag {other}");
                std::process::exit(2);
            }
        }
        i += 1;
    }
    args
}

/// Locate mesh files by id under the Thingi10K repo layout (same walk as
/// thingi_sweep's discover, filtered to the requested ids).
fn locate(root: &str, ids: &[u64]) -> Vec<(u64, std::path::PathBuf)> {
    let mut out = Vec::new();
    let Ok(repos) = std::fs::read_dir(root) else {
        eprintln!("cannot read mesh root {root}");
        std::process::exit(2);
    };
    for repo in repos.flatten() {
        let meshes = repo.path().join("meshes");
        let Ok(entries) = std::fs::read_dir(&meshes) else {
            continue;
        };
        for f in entries.flatten() {
            let name = f.file_name().to_string_lossy().to_string();
            if let Some(stem) = name.strip_suffix(".stl.zip") {
                if let Ok(id) = stem.parse::<u64>() {
                    if ids.contains(&id) {
                        out.push((id, f.path()));
                    }
                }
            }
        }
    }
    out.sort();
    out.dedup_by_key(|(id, _)| *id);
    out
}

fn read_zipped_stl(path: &std::path::Path) -> Result<Vec<u8>, String> {
    let file = std::fs::File::open(path).map_err(|e| e.to_string())?;
    let mut zip = zip::ZipArchive::new(file).map_err(|e| e.to_string())?;
    for i in 0..zip.len() {
        let mut entry = zip.by_index(i).map_err(|e| e.to_string())?;
        if entry.name().to_ascii_lowercase().ends_with(".stl") {
            let mut buf = Vec::with_capacity(entry.size() as usize);
            entry.read_to_end(&mut buf).map_err(|e| e.to_string())?;
            return Ok(buf);
        }
    }
    Err("no .stl entry in zip".to_string())
}

/// STL bytes → robust-imported Manifold. Byte-for-byte the thingi_sweep
/// pipeline — the referee must judge the same inputs the sweep judged.
fn import_stl_bytes(data: &[u8]) -> Manifold {
    let head = String::from_utf8_lossy(&data[..data.len().min(512)]).to_string();
    let mut positions: Vec<f32> = Vec::new();
    if head.trim_start().starts_with("solid") && head.contains("facet") {
        for line in String::from_utf8_lossy(data).lines() {
            if let Some(rest) = line.trim_start().strip_prefix("vertex") {
                for tok in rest.split_whitespace().take(3) {
                    positions.push(tok.parse::<f64>().unwrap_or(0.0) as f32);
                }
            }
        }
    } else if data.len() >= 84 {
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
    let c = [
        (min[0] + max[0]) / 2.0,
        (min[1] + max[1]) / 2.0,
        (min[2] + max[2]) / 2.0,
    ];
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

/// SplitMix64: deterministic, seedable, and good enough for stratification
/// jitter. Avoids pulling a rand dependency into the examples.
struct Rng(u64);
impl Rng {
    fn next_f64(&mut self) -> f64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z = z ^ (z >> 31);
        (z >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// Monte Carlo volume of {w ≥ 1} of a single soup — the once-counted
/// material of one operand. Compared against the divergence-theorem volume
/// (which integrates winding, counting a doubly-wound region twice), the gap
/// measures how self-overlapping the input is: the root cause when the exact
/// engine reports more volume than physically exists.
fn operand_volume(tris: &[[Vec3; 3]], samples: usize) -> f64 {
    let (mut min, mut max) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for t in tris {
        for v in t {
            let c = [v.x, v.y, v.z];
            for k in 0..3 {
                min[k] = min[k].min(c[k]);
                max[k] = max[k].max(c[k]);
            }
        }
    }
    let ext = [max[0] - min[0], max[1] - min[1], max[2] - min[2]];
    let idx = WindingIndex::new(tris);
    let mut rng = Rng(0xA11_0F_5A);
    let mut hits = 0usize;
    for _ in 0..samples {
        let pt = R3::from_vec3(Vec3::new(
            min[0] + ext[0] * rng.next_f64(),
            min[1] + ext[1] * rng.next_f64(),
            min[2] + ext[2] * rng.next_f64(),
        ));
        if winding_number_indexed(&pt, tris, &idx) >= 1 {
            hits += 1;
        }
    }
    ext[0] * ext[1] * ext[2] * hits as f64 / samples as f64
}

/// Monte Carlo volume of {w_P ≥ 1} ∪ {w_Q ≥ 1} over the joint bounding box.
/// Returns (estimate, standard error). Jittered-grid stratification cuts the
/// variance well below the iid binomial bound the error term reports, so the
/// reported σ is conservative.
fn referee_volume(p: &[[Vec3; 3]], q: &[[Vec3; 3]], samples: usize) -> (f64, f64) {
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
    let vol_box = ext[0] * ext[1] * ext[2];
    let idx_p = WindingIndex::new(p);
    let idx_q = WindingIndex::new(q);

    // Jittered grid: n^3 cells (n = cube root of samples), one point each,
    // remainder drawn uniformly.
    let n = (samples as f64).cbrt().floor() as usize;
    let mut rng = Rng(0x5EED_CAFE_F00D_D1CE);
    let mut hits = 0usize;
    let mut total = 0usize;
    let classify = |x: f64, y: f64, z: f64| -> bool {
        let pt = R3::from_vec3(Vec3::new(x, y, z));
        winding_number_indexed(&pt, p, &idx_p) >= 1 || winding_number_indexed(&pt, q, &idx_q) >= 1
    };
    for ix in 0..n {
        for iy in 0..n {
            for iz in 0..n {
                let x = min[0] + ext[0] * ((ix as f64 + rng.next_f64()) / n as f64);
                let y = min[1] + ext[1] * ((iy as f64 + rng.next_f64()) / n as f64);
                let z = min[2] + ext[2] * ((iz as f64 + rng.next_f64()) / n as f64);
                total += 1;
                if classify(x, y, z) {
                    hits += 1;
                }
            }
        }
    }
    for _ in total..samples {
        let x = min[0] + ext[0] * rng.next_f64();
        let y = min[1] + ext[1] * rng.next_f64();
        let z = min[2] + ext[2] * rng.next_f64();
        total += 1;
        if classify(x, y, z) {
            hits += 1;
        }
    }
    let frac = hits as f64 / total as f64;
    (
        vol_box * frac,
        vol_box * (frac * (1.0 - frac) / total as f64).sqrt(),
    )
}

/// One engine pass with panic capture; (status, volume, tris).
fn run_pass(a: &Manifold, b: &Manifold, engine: BooleanEngine) -> (String, f64, usize) {
    match catch_unwind(AssertUnwindSafe(|| {
        a.boolean_with_engine_and_token(b, OpType::Add, engine, None)
    })) {
        Ok(out) => (format!("{:?}", out.status()), out.volume(), out.num_tri()),
        Err(_) => ("PANIC".to_string(), f64::NAN, 0),
    }
}

fn main() {
    let args = parse_args();
    let mut ids = args.ids.clone();
    if let Some(run) = args.from_run {
        let conn = rusqlite::Connection::open(&args.db).expect("open sqlite db");
        let mut stmt = conn
            .prepare(
                "SELECT mesh_id FROM results
                 WHERE run_id = ?1 AND match_ok = 0
                   AND exact_status = 'NoError' AND robust_status = 'NoError'
                 ORDER BY mesh_id",
            )
            .expect("query mismatches");
        let got: Vec<u64> = stmt
            .query_map([run], |r| r.get::<_, i64>(0))
            .expect("query")
            .filter_map(|r| r.ok())
            .map(|v| v as u64)
            .take(args.limit)
            .collect();
        ids.extend(got);
    }
    if ids.is_empty() {
        eprintln!("no mesh ids: pass --ids or --from-run");
        std::process::exit(2);
    }
    let files = locate(&args.root, &ids);
    println!(
        "{} of {} requested meshes found; {} samples per referee",
        files.len(),
        ids.len(),
        args.samples
    );
    println!(
        "{:>9} {:>13} {:>10} {:>13} {:>13} {:>7} {:>8} {:>8}  verdict",
        "mesh", "referee", "±σ", "exact", "robust", "badwall", "soup", "overlapA"
    );

    for (mesh_id, path) in &files {
        let stl = match read_zipped_stl(path) {
            Ok(b) => b,
            Err(e) => {
                println!("{mesh_id:>9} ZIP error: {e}");
                continue;
            }
        };
        let Ok(a) = catch_unwind(AssertUnwindSafe(|| import_stl_bytes(&stl))) else {
            println!("{mesh_id:>9} PANIC on import");
            continue;
        };
        if a.status() != Error::NoError || a.num_tri() == 0 {
            println!("{mesh_id:>9} import status {:?}", a.status());
            continue;
        }
        let is_soup = a.as_impl().is_soup;
        let b = a
            .rotate(30.0, 45.0, 60.0)
            .translate(Vec3::new(0.3, 0.0, 0.0));
        let p_tris = soup::impl_to_tris(a.as_impl());
        let q_tris = soup::impl_to_tris(b.as_impl());

        let (ref_vol, sigma) = referee_volume(&p_tris, &q_tris, args.samples);
        let (e_status, e_vol, _) = if is_soup {
            ("n/a".to_string(), f64::NAN, 0)
        } else {
            run_pass(&a, &b, BooleanEngine::Exact)
        };
        let (r_status, r_vol, _) = run_pass(&a, &b, BooleanEngine::Robust);

        // Arrangement self-consistency: walls whose winding step disagrees
        // with the resolved cell windings. Non-zero means the cell complex
        // itself is broken, not just the final classification.
        let bad_walls = catch_unwind(AssertUnwindSafe(|| {
            let graph = intersection_graph::build_graph(&p_tris, &q_tris);
            let complex = cells::build_cells(&graph);
            let wind = cells::windings(&graph, &complex, [&p_tris, &q_tris]);
            cells::inconsistent_walls(&complex, &wind).len()
        }))
        .unwrap_or(usize::MAX);

        // Self-overlap of operand A: divergence-theorem volume vs sampled
        // once-counted material. div/sampled > 1 means the input winds some
        // region more than once — exactly what makes the exact engine's
        // volume exceed physical material.
        let a_div = a.volume();
        let a_sampled = operand_volume(&p_tris, args.samples / 2);
        let overlap = a_div / a_sampled;

        // A zero-hit sample makes σ collapse to 0; fall back to an absolute
        // band so "referee saw no material at all" can still judge an empty
        // result as correct.
        let z = |v: f64| {
            if sigma > 0.0 {
                (v - ref_vol).abs() / sigma
            } else if (v - ref_vol).abs() < 1e-6 {
                0.0
            } else {
                f64::INFINITY
            }
        };
        let verdict = match (e_vol.is_finite(), r_vol.is_finite()) {
            (true, true) => match (z(e_vol) <= 4.0, z(r_vol) <= 4.0) {
                (true, true) => "both agree with referee",
                (true, false) => "EXACT right, robust wrong",
                (false, true) => "ROBUST right, exact wrong",
                (false, false) => "neither matches referee",
            },
            (false, true) => {
                if z(r_vol) <= 4.0 {
                    "robust agrees with referee"
                } else {
                    "robust disagrees with referee"
                }
            }
            _ => "no engine volume",
        };
        println!(
            "{mesh_id:>9} {ref_vol:>13.6} {sigma:>10.6} {e_vol:>13.6} {r_vol:>13.6} {bad_walls:>7} {is_soup:>8} {overlap:>8.3}  {verdict} (exact {e_status}, robust {r_status})"
        );
    }
}
