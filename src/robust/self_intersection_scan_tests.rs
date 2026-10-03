// Copyright 2026 Lars Brubaker
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// Tests for the Auto engine's self-intersection pre-check
// (robust/soup.rs `has_self_intersections` and the operand order in
// boolean3.rs `boolean_dispatch_full`). They pin the verdicts the scan
// gives, so the speedups it carries (smaller operand first, the parallel
// per-triangle loop, the symmetric vertex-neighbour shortcut in
// `genuine_contact`) are held to "same answer". Shared 1:1 with
// manifold-sharp's SelfIntersectionScanTests (its commit 1d9162c), fixtures
// included; where sharp flips its runtime parallelism switch, these assert
// whichever build they run in (run the suite with and without
// `--features parallel`).
//
// The selfisect-*.txt fixtures are MatterCAD scene parts dumped from their
// .3mf through f32, the way agg-sharp's Mesh stores them, and imported the
// way its ManifoldKernel imports them: from_mesh_gl64_robust,
// repair_orientation, as_original. The two "fold" parts carry ulp-scale
// folds, some between triangles sharing a single vertex.

use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

use super::super::intersection_graph::{real_self_contact, SelfCutStats};
use super::SELF_INTERSECT_PAR_THRESHOLD;
use super::{genuine_contact, has_self_intersections, has_self_intersections_with_token};
use crate::cancel::CancelToken;
use crate::linalg::Vec3;
use crate::manifold::Manifold;
use crate::progress::{Phase, ProgressReporter};
use crate::types::{BooleanEngine, Error, MeshGL64, OpType};

const FOLD_F594: &str = include_str!("testdata/selfisect-fold-F594AC8BB1542BE9.txt");
const FOLD_BC6E: &str = include_str!("testdata/selfisect-fold-BC6E6DD496D38966.txt");
const CLEAN_63F7: &str = include_str!("testdata/selfisect-clean-63F7A7FFD6A5F7C6.txt");

fn v(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z)
}

/// Reads a "numVert numTri" fixture (one "x y z" line per vertex, one
/// "v0 v1 v2" line per triangle) and imports it the way agg-sharp's
/// ManifoldKernel does.
fn load_fixture(text: &str) -> Manifold {
    let mut lines = text.lines();
    let header: Vec<usize> = lines
        .next()
        .expect("fixture header")
        .split_whitespace()
        .map(|s| s.parse().expect("header count"))
        .collect();
    let mut mesh = MeshGL64 {
        num_prop: 3,
        ..Default::default()
    };
    for line in lines.by_ref().take(header[0]) {
        mesh.vert_properties.extend(
            line.split_whitespace()
                .map(|s| s.parse::<f64>().expect("coordinate")),
        );
    }
    for line in lines.take(header[1]) {
        mesh.tri_verts.extend(
            line.split_whitespace()
                .map(|s| s.parse::<u64>().expect("index")),
        );
    }
    let imported = Manifold::from_mesh_gl64_robust(&mesh);
    let repaired = imported.repair_orientation();
    let original = repaired.as_original();
    if original.status() == Error::NoError {
        original
    } else {
        repaired
    }
}

/// Whether [`crate::par::maybe_par_any_ct`] ran its predicate on a rayon
/// worker at this size — the parallel loop rather than the plain one.
fn any_ran_in_a_parallel_loop(n: usize, threshold: usize) -> bool {
    let inside = AtomicBool::new(false);
    crate::par::maybe_par_any_ct(
        n,
        threshold,
        None,
        || (),
        |_, _| {
            #[cfg(feature = "parallel")]
            if rayon::current_thread_index().is_some() {
                inside.store(true, Ordering::Relaxed);
            }
            false
        },
    );
    inside.load(Ordering::Relaxed)
}

/// The fixtures keep their verdicts. Each run imports afresh, because the
/// verdict is cached per impl.
#[test]
fn imported_parts_keep_their_verdicts() {
    let parallel = cfg!(feature = "parallel");
    let cases = [
        ("selfisect-fold-F594AC8BB1542BE9", FOLD_F594, true),
        ("selfisect-fold-BC6E6DD496D38966", FOLD_BC6E, true),
        ("selfisect-clean-63F7A7FFD6A5F7C6", CLEAN_63F7, false),
    ];
    for (name, text, expected) in cases {
        let m = load_fixture(text);
        assert_eq!(m.status(), Error::NoError, "{name}");
        assert_eq!(
            m.has_self_intersections(),
            expected,
            "{name}, parallel={parallel}"
        );
    }
    // A clean part big enough that the parallel build really runs in parallel.
    let sphere = Manifold::sphere(1.0, 64);
    assert!(
        !sphere.has_self_intersections(),
        "sphere, parallel={parallel}"
    );

    // Anti-vacuity: the scan's own helper, at each part's triangle count and
    // the scan's threshold, really runs on rayon workers in the parallel
    // build and the plain loop otherwise, so the parallel verdicts above came
    // from the parallel loop. The parts that go parallel are the two folds and
    // the sphere.
    for n in [
        load_fixture(FOLD_F594).num_tri(),
        load_fixture(FOLD_BC6E).num_tri(),
        sphere.num_tri(),
    ] {
        assert_eq!(
            any_ran_in_a_parallel_loop(n, SELF_INTERSECT_PAR_THRESHOLD),
            parallel,
            "n={n}"
        );
    }
}

/// A cancelled scan answers true (route to robust) and caches nothing, so the
/// next uncancelled call computes the real verdict.
#[test]
fn cancelled_scan_answers_true_and_caches_nothing() {
    let sphere = Manifold::sphere(1.0, 64);
    let imp = sphere.as_impl();
    let token = CancelToken::new();
    token.cancel();
    assert!(has_self_intersections_with_token(imp, Some(&token)));
    assert_eq!(imp.self_intersects.get(), None);
    assert!(!has_self_intersections(imp));
    assert_eq!(imp.self_intersects.get(), Some(false));
}

/// Two triangles sharing one vertex, where t2's corners straddle t1's plane
/// but t1's two other corners are strictly above t2's plane: they meet only
/// at the shared vertex. In this orientation real_self_contact's own shortcut
/// (t2's corners against t1's plane) cannot decide and pays for the full
/// tri-tri test; the mirror shortcut in genuine_contact decides it before
/// real_self_contact is reached, with the same verdict.
#[test]
fn vertex_neighbour_with_only_its_own_corners_one_sided_is_benign() {
    let t1 = [v(0.0, 0.0, 0.0), v(1.0, 0.0, 1.0), v(0.0, -1.0, 1.0)];
    let t2 = [v(0.0, 0.0, 0.0), v(2.0, 0.0, 0.0), v(0.0, 2.0, 0.0)];

    // The old path alone: its shortcut misses, tri-tri finds the point contact.
    let mut old = SelfCutStats::default();
    assert!(real_self_contact(t1, t2, &mut old).is_none());
    assert_eq!(
        old.vert_benign, 0,
        "real_self_contact's shortcut cannot see this orientation"
    );
    assert_eq!(old.full, 1);
    assert_eq!(old.full_point, 1);

    // genuine_contact: the same verdict, from the new shortcut, without tri-tri.
    let mut stats = SelfCutStats::default();
    assert!(!genuine_contact(t1, t2, &mut stats));
    assert_eq!(stats.vert_benign, 1, "the new shortcut counted it");
    assert_eq!(stats.full, 0, "the new shortcut, not tri-tri, decided");
}

/// A vertex neighbour with one of t1's other corners exactly on t2's plane
/// (sign Zero) is not "strictly on one side", so the new shortcut must fall
/// through to real_self_contact and give its verdict, by the same path.
#[test]
fn vertex_neighbour_with_a_corner_on_the_plane_falls_through() {
    // (-1, -1, 0) lies on t2's plane z = 0, outside t2; (0, -1, 1) is above it.
    let t1 = [v(0.0, 0.0, 0.0), v(-1.0, -1.0, 0.0), v(0.0, -1.0, 1.0)];
    let t2 = [v(0.0, 0.0, 0.0), v(2.0, 0.0, 0.0), v(0.0, 2.0, 0.0)];

    let mut old = SelfCutStats::default();
    let old_verdict = real_self_contact(t1, t2, &mut old).is_some();

    let mut stats = SelfCutStats::default();
    assert_eq!(genuine_contact(t1, t2, &mut stats), old_verdict);
    assert!(
        !old_verdict,
        "t1 meets z = 0 outside t2 except at the shared vertex"
    );
    assert_eq!(
        stats.vert_benign, old.vert_benign,
        "the new shortcut did not count it"
    );
    assert_eq!(
        stats.full, 1,
        "it fell through to real_self_contact's full test"
    );
    assert_eq!(stats.full, old.full);
}

/// A hit and a cancel together answer true: the hit is a genuine witness, so
/// the true stands whether or not other workers then saw the cancel. The first
/// predicate call cancels the token and reports a hit, so the cancel is always
/// in flight with the hit, on both loops.
#[test]
fn a_hit_racing_a_cancel_answers_true() {
    let token = CancelToken::new();
    let claimed = AtomicUsize::new(0);
    let answer = crate::par::maybe_par_any_ct(
        100_000,
        SELF_INTERSECT_PAR_THRESHOLD,
        Some(&token),
        || (),
        |_, _| {
            if claimed.swap(1, Ordering::SeqCst) != 0 {
                return false;
            }
            token.cancel();
            true
        },
    );
    assert!(token.is_cancelled());
    assert_eq!(answer, Some(true));
    // And in the parallel build the call above really ran the parallel loop.
    assert_eq!(
        any_ran_in_a_parallel_loop(100_000, SELF_INTERSECT_PAR_THRESHOLD),
        cfg!(feature = "parallel")
    );
}

/// Two triangles sharing one vertex whose other corners straddle each other's
/// planes and cross through each other's interior: still a genuine contact.
#[test]
fn vertex_neighbour_that_crosses_is_a_contact() {
    let t1 = [v(0.0, 0.0, 0.0), v(1.0, 0.5, 1.0), v(0.5, 1.0, -1.0)];
    let t2 = [v(0.0, 0.0, 0.0), v(2.0, 0.0, 0.0), v(0.0, 2.0, 0.0)];
    assert!(genuine_contact(t1, t2, &mut SelfCutStats::default()));
    assert!(genuine_contact(t2, t1, &mut SelfCutStats::default()));
}

/// Did an Auto union run the exact engine (its one "exact boolean" report)?
fn ran_exact_engine(a: &Manifold, b: &Manifold) -> bool {
    let phases = Arc::new(Mutex::new(Vec::<&'static str>::new()));
    let sink = Arc::clone(&phases);
    let reporter = ProgressReporter::new(move |phase: Phase, _| {
        sink.lock().expect("sink poisoned").push(phase.name());
    });
    let r = a.boolean_with_engine_and_progress(
        b,
        OpType::Add,
        BooleanEngine::Auto,
        None,
        Some(&reporter),
    );
    assert_eq!(r.status(), Error::NoError, "Auto union failed");
    let seen = phases.lock().expect("sink poisoned");
    seen.contains(&Phase::ExactBoolean.name())
}

/// Auto picks the same engine whichever operand comes first, and scans the
/// smaller operand first: when that one self-intersects, the larger is never
/// scanned.
#[test]
fn auto_picks_the_same_engine_whichever_operand_comes_first() {
    let fold = load_fixture(FOLD_F594);
    let clean = Manifold::sphere(1.0, 64);
    assert!(clean.num_tri() > fold.num_tri());

    assert!(!ran_exact_engine(&clean, &fold), "clean, fold");
    assert_eq!(
        clean.as_impl().self_intersects.get(),
        None,
        "the smaller, self-intersecting operand decided alone"
    );
    assert!(!ran_exact_engine(&fold, &clean), "fold, clean");

    let small = Manifold::cube(v(1.0, 1.0, 1.0), false).translate(v(0.5, 0.0, 0.0));
    assert!(ran_exact_engine(&clean, &small), "clean, small");
    assert!(ran_exact_engine(&small, &clean), "small, clean");
}
