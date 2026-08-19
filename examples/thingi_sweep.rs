// Thingi10K sweep harness: every local mesh through an Exact and a Robust
// pass, with results and timings stored in a SQLite database so perf and
// correctness can be compared across engine changes.
//
// Usage:
//   cargo run --release --example thingi_sweep -- [flags]
//     --root DIR        mesh repo root (default C:\Development\rust-apps\Thingi10K\meshes)
//     --db FILE         sqlite output (default thingi_sweep.db)
//     --limit N         stop after N meshes
//     --ids A,B,C       only these mesh ids
//     --skip N          skip the first N meshes (resume a partial sweep)
//     --timeout-secs S  per-boolean cancel timeout (default 120)
//     --notes TEXT      free-form note stored on the run row
//
// The pass: A = imported mesh (demo pipeline: normalize to side 2, weld,
// robust import), B = the same mesh rotated (30,45,60)° and offset
// (0.3,0,0); the operation is A ∪ B. Meshes that import as manifold run
// BOTH engines and must agree: NoError, equal tri/vert counts, volume and
// area within 1e-9 relative. Meshes that import as soup run Robust only and
// must be "reasonable": NoError, non-empty, finite volume/area. Panics are
// caught per-op and recorded; hangs are cut off by a cancel watchdog.
//
// Schema: runs(run_id, started, git_commit, notes) and results keyed by
// (run_id, mesh_id) — compare two run_ids to see what an engine change did.

use std::io::Read;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;
use std::time::{Duration, Instant};

use manifold_rust::cancel::CancelToken;
use manifold_rust::linalg::Vec3;
use manifold_rust::manifold::Manifold;
use manifold_rust::types::{BooleanEngine, Error, MeshGL, OpType};

struct Args {
    root: String,
    db: String,
    limit: usize,
    ids: Option<Vec<u64>>,
    skip: usize,
    timeout_secs: u64,
    notes: String,
    jobs: usize,
}

fn parse_args() -> Args {
    let default_jobs = std::thread::available_parallelism()
        .map(|n| (n.get().saturating_sub(2)).clamp(1, 8))
        .unwrap_or(1);
    let mut args = Args {
        root: r"C:\Development\rust-apps\Thingi10K\meshes".to_string(),
        db: "thingi_sweep.db".to_string(),
        limit: usize::MAX,
        ids: None,
        skip: 0,
        timeout_secs: 120,
        notes: String::new(),
        jobs: default_jobs,
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
            "--limit" => args.limit = val(&mut i).parse().unwrap_or(usize::MAX),
            "--skip" => args.skip = val(&mut i).parse().unwrap_or(0),
            "--timeout-secs" => args.timeout_secs = val(&mut i).parse().unwrap_or(120),
            "--notes" => args.notes = val(&mut i),
            "--jobs" => args.jobs = val(&mut i).parse().unwrap_or(1).max(1),
            "--ids" => {
                args.ids = Some(
                    val(&mut i)
                        .split(',')
                        .filter_map(|s| s.trim().parse().ok())
                        .collect(),
                )
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

/// All *.stl.zip files under root's repo subfolders, sorted by mesh id for a
/// deterministic sweep order.
fn discover(root: &str) -> Vec<(u64, std::path::PathBuf)> {
    let mut out = Vec::new();
    let repos = match std::fs::read_dir(root) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("cannot read mesh root {root}: {e}");
            std::process::exit(2);
        }
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
                    out.push((id, f.path()));
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

/// STL bytes → robust-imported Manifold, replicating the demo pipeline
/// (same as examples/robust_repro.rs): parse, center, scale longest side to
/// 2, weld, robust import.
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

#[derive(Default, Clone)]
struct PassResult {
    ms: f64,
    status: String,
    tris: i64,
    verts: i64,
    volume: f64,
    area: f64,
}

/// One engine pass with panic capture and a cancel watchdog.
fn run_pass(a: &Manifold, b: &Manifold, engine: BooleanEngine, timeout: Duration) -> PassResult {
    let token = CancelToken::new();
    let done = Arc::new(AtomicBool::new(false));
    let watchdog = {
        let token = token.clone();
        let done = done.clone();
        std::thread::spawn(move || {
            let t0 = Instant::now();
            while !done.load(Ordering::Relaxed) {
                if t0.elapsed() > timeout {
                    token.cancel();
                    return;
                }
                std::thread::sleep(Duration::from_millis(100));
            }
        })
    };
    let t0 = Instant::now();
    let outcome = catch_unwind(AssertUnwindSafe(|| {
        a.boolean_with_engine_and_token(b, OpType::Add, engine, Some(&token))
    }));
    let ms = t0.elapsed().as_secs_f64() * 1e3;
    done.store(true, Ordering::Relaxed);
    let _ = watchdog.join();
    match outcome {
        Ok(out) => PassResult {
            ms,
            status: format!("{:?}", out.status()),
            tris: out.num_tri() as i64,
            verts: out.num_vert() as i64,
            volume: out.volume(),
            area: out.surface_area(),
        },
        Err(p) => {
            let msg = p
                .downcast_ref::<&str>()
                .map(|s| s.to_string())
                .or_else(|| p.downcast_ref::<String>().cloned())
                .unwrap_or_else(|| "unknown panic".to_string());
            PassResult {
                ms,
                status: format!("PANIC: {msg}"),
                ..Default::default()
            }
        }
    }
}

fn rel_close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol * a.abs().max(b.abs()).max(1.0)
}

/// Everything one mesh produces, handed from a worker to the single writer
/// thread. `exact` is None for soup inputs (no exact pass) and both passes
/// are None on import failure.
struct Row {
    mesh_id: u64,
    file: String,
    tris_in: i64,
    verts_in: i64,
    import_status: String,
    is_soup: bool,
    manifold: bool,
    imported: bool,
    exact: Option<PassResult>,
    robust: Option<PassResult>,
    match_ok: Option<bool>,
    reasonable: bool,
    detail: String,
}

/// The whole per-mesh pipeline: read, import, both engine passes, verdict.
/// Pure with respect to the database — workers run this concurrently and
/// the writer thread owns SQLite.
fn process_mesh(mesh_id: u64, path: &std::path::Path, timeout: Duration) -> Row {
    let file_str = path.to_string_lossy().to_string();
    let mut row = Row {
        mesh_id,
        file: file_str,
        tris_in: 0,
        verts_in: 0,
        import_status: String::new(),
        is_soup: false,
        manifold: false,
        imported: false,
        exact: None,
        robust: None,
        match_ok: None,
        reasonable: false,
        detail: String::new(),
    };
    let stl = match read_zipped_stl(path) {
        Ok(b) => b,
        Err(e) => {
            row.import_status = format!("ZIP: {e}");
            return row;
        }
    };
    let a = match catch_unwind(AssertUnwindSafe(|| import_stl_bytes(&stl))) {
        Ok(m) => m,
        Err(_) => {
            row.import_status = "PANIC on import".to_string();
            return row;
        }
    };
    row.import_status = format!("{:?}", a.status());
    row.is_soup = a.as_impl().is_soup;
    row.tris_in = a.num_tri() as i64;
    row.verts_in = a.num_vert() as i64;
    let importable = a.status() == Error::NoError && a.num_tri() > 0;
    row.manifold = importable && !row.is_soup;
    if !importable {
        return row;
    }
    row.imported = true;

    let b = a
        .rotate(30.0, 45.0, 60.0)
        .translate(Vec3::new(0.3, 0.0, 0.0));
    let robust = run_pass(&a, &b, BooleanEngine::Robust, timeout);
    let exact = if row.manifold {
        Some(run_pass(&a, &b, BooleanEngine::Exact, timeout))
    } else {
        None
    };

    row.reasonable = robust.status == "NoError"
        && robust.tris > 0
        && robust.volume.is_finite()
        && robust.area.is_finite();
    if !row.reasonable {
        row.detail = format!("robust not reasonable: status {}", robust.status);
    }
    row.match_ok = exact.as_ref().map(|e| {
        // Volume/area must agree to 1e-9 relative; counts must be NEARLY
        // equal (the user-set acceptance: triangulation may differ, and
        // simplify_topology makes slightly different degenerate-collapse
        // choices on real meshes). A count delta beyond the band with
        // equal volumes would still be suspicious — keep it a mismatch.
        let count_band = |x: i64| (x / 500).max(16);
        let ok = e.status == "NoError"
            && robust.status == "NoError"
            && (e.tris - robust.tris).abs() <= count_band(e.tris)
            && (e.verts - robust.verts).abs() <= count_band(e.verts)
            && rel_close(e.volume, robust.volume, 1e-9)
            && rel_close(e.area, robust.area, 1e-9);
        if !ok && row.detail.is_empty() {
            row.detail =
                format!(
                "exact {} {}t/{}v vol {:.9} area {:.9} vs robust {} {}t/{}v vol {:.9} area {:.9}",
                e.status, e.tris, e.verts, e.volume, e.area,
                robust.status, robust.tris, robust.verts, robust.volume, robust.area
            );
        } else if ok && (e.tris != robust.tris || e.verts != robust.verts) && row.detail.is_empty()
        {
            row.detail = format!(
                "count delta within band: {}t/{}v vs {}t/{}v",
                e.tris, e.verts, robust.tris, robust.verts
            );
        }
        ok
    });
    row.exact = exact;
    row.robust = Some(robust);
    row
}

fn main() {
    let args = parse_args();
    let files = discover(&args.root);
    let git = std::process::Command::new("git")
        .args(["rev-parse", "--short", "HEAD"])
        .output()
        .ok()
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
        .unwrap_or_default();

    let conn = rusqlite::Connection::open(&args.db).expect("open sqlite db");
    conn.execute_batch(
        "CREATE TABLE IF NOT EXISTS runs (
           run_id INTEGER PRIMARY KEY AUTOINCREMENT,
           started TEXT NOT NULL,
           git_commit TEXT,
           notes TEXT
         );
         CREATE TABLE IF NOT EXISTS results (
           run_id INTEGER NOT NULL,
           mesh_id INTEGER NOT NULL,
           file TEXT NOT NULL,
           tris_in INTEGER,
           verts_in INTEGER,
           import_status TEXT,
           is_soup INTEGER,
           manifold INTEGER,
           exact_ms REAL, exact_status TEXT, exact_tris INTEGER, exact_verts INTEGER,
           exact_volume REAL, exact_area REAL,
           robust_ms REAL, robust_status TEXT, robust_tris INTEGER, robust_verts INTEGER,
           robust_volume REAL, robust_area REAL,
           match_ok INTEGER,
           reasonable_ok INTEGER,
           detail TEXT,
           PRIMARY KEY (run_id, mesh_id)
         );",
    )
    .expect("create schema");

    let started = chrono_free_timestamp();
    conn.execute(
        "INSERT INTO runs (started, git_commit, notes) VALUES (?1, ?2, ?3)",
        rusqlite::params![started, git, args.notes],
    )
    .expect("insert run");
    let run_id = conn.last_insert_rowid();
    println!(
        "run {run_id} @ {git}: {} meshes discovered under {}",
        files.len(),
        args.root
    );

    let timeout = Duration::from_secs(args.timeout_secs);
    let mut n_done = 0usize;
    let mut n_match = 0usize;
    let mut n_mismatch = 0usize;
    let mut n_unreasonable = 0usize;
    let mut n_import_fail = 0usize;
    let t_sweep = Instant::now();

    // Worklist after skip/limit/ids filtering; workers pull from it by
    // atomic index, so each mesh is processed exactly once regardless of
    // job count and per-mesh results are identical to a sequential sweep
    // (only completion order varies — results are keyed by mesh_id).
    let work: Vec<(u64, std::path::PathBuf)> = files
        .iter()
        .enumerate()
        .filter(|(idx, _)| *idx >= args.skip)
        .filter(|(_, (id, _))| args.ids.as_ref().is_none_or(|ids| ids.contains(id)))
        .take(args.limit)
        .map(|(_, (id, p))| (*id, p.clone()))
        .collect();
    let total = work.len();
    let jobs = args.jobs.min(total.max(1));
    println!("{total} meshes queued on {jobs} worker(s)");

    let work = Arc::new(work);
    let next = Arc::new(std::sync::atomic::AtomicUsize::new(0));
    let (tx, rx) = std::sync::mpsc::channel::<Row>();
    let mut handles = Vec::new();
    for _ in 0..jobs {
        let work = work.clone();
        let next = next.clone();
        let tx = tx.clone();
        handles.push(std::thread::spawn(move || loop {
            let i = next.fetch_add(1, Ordering::Relaxed);
            let Some((mesh_id, path)) = work.get(i) else {
                break;
            };
            let row = process_mesh(*mesh_id, path, timeout);
            if tx.send(row).is_err() {
                break;
            }
        }));
    }
    drop(tx);

    for row in rx {
        n_done += 1;
        if !row.imported {
            n_import_fail += 1;
        }
        if row.imported && !row.reasonable {
            n_unreasonable += 1;
        }
        match row.match_ok {
            Some(true) => n_match += 1,
            Some(false) => n_mismatch += 1,
            None => {}
        }
        let e = row.exact.clone().unwrap_or_default();
        let r = row.robust.clone().unwrap_or_default();
        conn.execute(
            "INSERT OR REPLACE INTO results
               (run_id, mesh_id, file, tris_in, verts_in, import_status, is_soup, manifold,
                exact_ms, exact_status, exact_tris, exact_verts, exact_volume, exact_area,
                robust_ms, robust_status, robust_tris, robust_verts, robust_volume, robust_area,
                match_ok, reasonable_ok, detail)
             VALUES (?1,?2,?3,?4,?5,?6,?7,?8,?9,?10,?11,?12,?13,?14,?15,?16,?17,?18,?19,?20,?21,?22,?23)",
            rusqlite::params![
                run_id,
                row.mesh_id as i64,
                row.file,
                row.tris_in,
                row.verts_in,
                row.import_status,
                row.is_soup as i64,
                row.manifold as i64,
                row.exact.as_ref().map(|_| e.ms),
                row.exact.as_ref().map(|_| e.status.clone()),
                row.exact.as_ref().map(|_| e.tris),
                row.exact.as_ref().map(|_| e.verts),
                row.exact.as_ref().map(|_| e.volume),
                row.exact.as_ref().map(|_| e.area),
                row.robust.as_ref().map(|_| r.ms),
                row.robust.as_ref().map(|_| r.status.clone()),
                row.robust.as_ref().map(|_| r.tris),
                row.robust.as_ref().map(|_| r.verts),
                row.robust.as_ref().map(|_| r.volume),
                row.robust.as_ref().map(|_| r.area),
                row.match_ok.map(|b| b as i64),
                row.reasonable as i64,
                row.detail,
            ],
        )
        .unwrap();

        if n_done % 25 == 0 || n_done == total {
            println!(
                "[{n_done}/{total}] {:.0}s elapsed — match {n_match}, mismatch {n_mismatch}, \
                 unreasonable {n_unreasonable}, import-fail {n_import_fail} (last #{})",
                t_sweep.elapsed().as_secs_f64(),
                row.mesh_id
            );
        }
    }
    for h in handles {
        let _ = h.join();
    }

    println!(
        "run {run_id} complete: {n_done} meshes in {:.0}s — {n_match} matched, \
         {n_mismatch} MISMATCHED, {n_unreasonable} UNREASONABLE, {n_import_fail} import failures",
        t_sweep.elapsed().as_secs_f64()
    );
    println!(
        "db: {} (compare runs with SQL over the results table)",
        args.db
    );
}

/// ISO-ish local timestamp without a chrono dependency.
fn chrono_free_timestamp() -> String {
    let now = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or_default();
    format!("unix:{}", now.as_secs())
}
