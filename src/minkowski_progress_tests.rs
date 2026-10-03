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

// Tests for minkowski.rs's cancellation and progress parameters
// (`minkowski_with_progress`), shared 1:1 with the last six tests of
// manifold-sharp's MinkowskiTests (its commit cf2c170).
//
// ORDERING ASSERTIONS ONLY HOLD SEQUENTIALLY. Either the fixture keeps the
// per-face hull map below its parallel threshold of 8 (a tetrahedron's four
// faces), or the assertion is skipped under `--features parallel`: two
// workers can cross the throttle together and both report, so a report count
// over a parallel map would assert something the reporter does not promise.
// The closing 1.0 is asserted in both builds: it is emitted on the calling
// thread after the workers have joined. (manifold-sharp switches parallelism
// at runtime and runs both modes in one test; here the build decides.)

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

use super::*;
use crate::linalg::{mat4_to_mat3x4, scaling_matrix, translation_matrix, Mat3x4};
use crate::manifold::Manifold;
use crate::types::Error;

type Event = (&'static str, Option<f64>);

/// Collects every callback a run emits (the reporter may call back from any
/// thread).
#[derive(Default)]
struct Sink {
    events: Arc<Mutex<Vec<Event>>>,
}

impl Sink {
    fn reporter(&self) -> ProgressReporter {
        let events = Arc::clone(&self.events);
        ProgressReporter::new(move |phase: Phase, fraction| {
            events
                .lock()
                .expect("sink poisoned")
                .push((phase.name(), fraction));
        })
    }

    fn events(&self) -> Vec<Event> {
        self.events.lock().expect("sink poisoned").clone()
    }
}

/// Every reported fraction is in range and no smaller than the one before
/// it, the first is 0.0 and the last exactly 1.0.
fn assert_rises_to_one(sink: &Sink) {
    let events = sink.events();
    assert!(events.len() > 1, "begin_phase alone is not progress");
    let mut previous = -1.0;
    for (name, fraction) in &events {
        assert_eq!(*name, "minkowski");
        let f = fraction.expect("Minkowski knows its work total, so no report is indeterminate");
        assert!(
            f >= previous,
            "progress went backwards, {previous} then {f}"
        );
        assert!(f <= 1.0);
        previous = f;
    }
    assert_eq!(events[0].1, Some(0.0));
    assert_eq!(
        events[events.len() - 1].1,
        Some(1.0),
        "the closing merge is reported by complete_phase, so a finished run lands on 1.0"
    );
}

/// Two unit cubes unioned into an L — the cheapest non-convex solid the crate
/// can build out of its own primitives.
fn non_convex_pair() -> ManifoldImpl {
    let cube = Manifold::cube(Vec3::splat(1.0), true);
    cube.union(&cube.translate(Vec3::new(0.6, 0.6, 0.0)))
        .as_impl()
        .clone()
}

/// FNV-1a over the raw bit patterns of every vertex coordinate and halfedge
/// index — a bit-exact fingerprint. Deliberately does NOT read mesh IDs,
/// which move between two sequential runs because every hull mints fresh
/// ones from a process-global counter.
fn geometry_hash(mesh: &ManifoldImpl) -> u64 {
    let mut hash: u64 = 14695981039346656037;
    let mut mix = |value: u64| {
        for shift in (0..64).step_by(8) {
            hash ^= (value >> shift) & 0xFF;
            hash = hash.wrapping_mul(1099511628211);
        }
    };
    for v in &mesh.vert_pos {
        mix(v.x.to_bits());
        mix(v.y.to_bits());
        mix(v.z.to_bits());
    }
    for e in &mesh.halfedge {
        mix(e.start_vert as u32 as u64);
        mix(e.end_vert as u32 as u64);
        mix(e.paired_halfedge as u32 as u64);
    }
    hash
}

/// The convex×convex fast path reports the whole bar: it does one hull and
/// one merge, so it has only 0 and 1 to say, and it has to say both.
#[test]
fn sum_reports_progress_that_only_rises_and_ends_at_one() {
    let a = ManifoldImpl::cube(&Mat3x4::identity());
    let b = ManifoldImpl::cube(&Mat3x4::identity());
    let sink = Sink::default();
    let sum = minkowski_sum_with_progress(&a, &b, None, Some(&sink.reporter()));
    assert!(sum.num_tri() > 0);
    assert_rises_to_one(&sink);
}

/// The per-face path — every erosion takes it, convex operands or not —
/// reports one unit per hull, one per batch reduction and one for the
/// closing merge. A tetrahedron because its four faces keep the hull map
/// under the parallel threshold of 8.
#[test]
fn difference_reports_progress_that_only_rises_and_ends_at_one() {
    let a = ManifoldImpl::tetrahedron(&mat4_to_mat3x4(scaling_matrix(Vec3::splat(4.0))));
    let b = ManifoldImpl::cube(&mat4_to_mat3x4(
        translation_matrix(Vec3::splat(-0.25)) * scaling_matrix(Vec3::splat(0.5)),
    ));
    let sink = Sink::default();
    let diff = minkowski_difference_with_progress(&a, &b, None, Some(&sink.reporter()));
    assert!(diff.num_tri() > 0);
    assert_rises_to_one(&sink);
    // Four hulls, one batch reduction, one closing merge, plus begin_phase's
    // own report: the unit model is asserted rather than only its endpoints,
    // because it is the thing a consumer's bar is scaled against.
    assert_eq!(sink.events().len(), 7);
}

/// A cancel that lands mid-computation stops the work almost at once,
/// answers with the empty `Cancelled` mesh the boolean pipeline answers with,
/// and leaves both operands exactly as it found them. The token is tripped
/// from the progress callback — the earliest in-kernel moment a test can
/// reach — so this measures cancellation latency in work units, not time.
#[test]
fn cancel_mid_computation_aborts_promptly_and_leaves_the_inputs_untouched() {
    // 128 triangles, so an uncancelled erosion reports 130 units; the
    // assertion below is that a cancel after the first hull costs a small
    // fraction of that.
    let a = Manifold::sphere(1.0, 16).as_impl().clone();
    let b = ManifoldImpl::cube(&mat4_to_mat3x4(scaling_matrix(Vec3::splat(0.1))));
    assert!(
        a.num_tri() > 100,
        "the fixture has to be big enough that stopping early is visible"
    );
    let (a_before, b_before) = (geometry_hash(&a), geometry_hash(&b));

    let token = CancelToken::new();
    let reports = Arc::new(AtomicUsize::new(0));
    // Cancel on the second callback: the first is begin_phase's
    // unconditional 0.0, which arrives before any hull runs, so tripping on
    // that would only re-test the entry gate.
    let (trip, count) = (token.clone(), Arc::clone(&reports));
    let reporter = ProgressReporter::new(move |_, _| {
        if count.fetch_add(1, Ordering::SeqCst) + 1 >= 2 {
            trip.cancel();
        }
    });

    let result = minkowski_difference_with_progress(&a, &b, Some(&token), Some(&reporter));

    assert_eq!(result.status, Error::Cancelled);
    assert_eq!(
        result.num_tri(),
        0,
        "a cancelled result must be empty, as it is for a boolean"
    );
    // Well under the 130 an uncancelled run reports. Loose because the
    // parallel map cannot recall iterations already in flight.
    let n = reports.load(Ordering::SeqCst);
    assert!(n < 40, "cancel was ignored for {n} of 130 work units");
    assert_eq!(
        geometry_hash(&a),
        a_before,
        "Minkowski must not mutate the solid it was handed"
    );
    assert_eq!(
        geometry_hash(&b),
        b_before,
        "Minkowski must not mutate the structuring element it was handed"
    );
}

/// The whole point of the additive design: the plain entry points produce
/// the bit-identical mesh an instrumented, tokened run produces. The
/// live-token half is the load-bearing one: a token routes the hull maps
/// through a structurally different collect inside `maybe_par_map_ct`, and a
/// reporter wraps the map's closure. Both branches are covered — the convex
/// sum and the per-face erosion.
#[test]
fn default_parameters_are_bit_identical_to_an_instrumented_run() {
    let cube = ManifoldImpl::cube(&Mat3x4::identity());
    let tetra = ManifoldImpl::tetrahedron(&mat4_to_mat3x4(scaling_matrix(Vec3::splat(4.0))));
    let tool = ManifoldImpl::cube(&mat4_to_mat3x4(scaling_matrix(Vec3::splat(0.5))));

    let live = CancelToken::new();
    let sink = Sink::default();
    assert_eq!(
        geometry_hash(&minkowski_sum_with_progress(
            &cube,
            &cube,
            Some(&live),
            Some(&sink.reporter())
        )),
        geometry_hash(&minkowski_sum(&cube, &cube)),
        "instrumenting the convex sum moved a bit"
    );
    assert_eq!(
        geometry_hash(&minkowski_difference_with_progress(
            &tetra,
            &tool,
            Some(&live),
            Some(&sink.reporter())
        )),
        geometry_hash(&minkowski_difference(&tetra, &tool)),
        "instrumenting the per-face erosion moved a bit"
    );
    assert!(!live.is_cancelled());
    assert!(
        !sink.events().is_empty(),
        "a run that reported nothing would prove nothing"
    );
}

/// The bar has to land on 1.0 even when the work total is past the
/// throttle's 100 reports per phase — the regime where `step` exceeds 1 and
/// every unit after the last step boundary is swallowed. Without
/// `complete_phase` a 290-unit erosion stopped reporting at 288/290.
#[test]
fn a_large_run_ends_at_exactly_one_in_both_parallel_modes() {
    // 288 triangles: 288 hulls + one batch reduction + the closing merge is
    // 290 units, so the throttle's step is 2 and roughly half the units never
    // report.
    let solid = Manifold::sphere(1.0, 24).as_impl().clone();
    let tool = ManifoldImpl::cube(&mat4_to_mat3x4(scaling_matrix(Vec3::splat(0.1))));
    assert_eq!(
        solid.num_tri(),
        288,
        "the unit total below is computed from this count"
    );
    let parallel = cfg!(feature = "parallel");

    let sink = Sink::default();
    minkowski_difference_with_progress(&solid, &tool, None, Some(&sink.reporter()));
    let events = sink.events();
    assert!(
        events.len() < 290,
        "parallel={parallel}: fewer reports than units is what makes this the throttled regime"
    );
    assert_eq!(
        events[events.len() - 1].1,
        Some(1.0),
        "parallel={parallel}: a finished run must leave the bar full, not at 288/290"
    );
    for (name, fraction) in &events {
        assert_eq!(*name, "minkowski");
        assert!(fraction.expect("determinate") <= 1.0);
    }
}

/// The non-convex × non-convex branch — a hull per face *pair* — reports the
/// total its work model predicts, and stops promptly when cancelled. Two
/// overlapping cubes unioned into an L, 28 triangles each: 28 × 28 hulls, 28
/// per-face reductions and the closing merge is 813 units, so `step` is 8 and
/// the throttle emits 101 times (at 8, 16 … 808) plus the opening and closing
/// reports. That exact count is asserted in the sequential build only: the
/// per-B-face map is 28 wide and goes parallel under `--features parallel`.
#[test]
fn non_convex_pair_reports_its_exact_total_and_cancels_promptly() {
    let l_shape = non_convex_pair();
    assert!(
        !l_shape.is_convex(),
        "both operands have to be non-convex to reach the per-face-pair branch"
    );
    assert_eq!(
        l_shape.num_tri(),
        28,
        "the 813-unit total below is computed from this count"
    );

    let sink = Sink::default();
    let sum = minkowski_sum_with_progress(&l_shape, &l_shape, None, Some(&sink.reporter()));
    assert!(sum.num_tri() > 0);
    assert_rises_to_one(&sink);
    if !cfg!(feature = "parallel") {
        assert_eq!(
            sink.events().len(),
            103,
            "one opening report, 101 throttled ones and the closing 1.0"
        );
    }

    // And a cancel lands within a face pair or two rather than at the end.
    let token = CancelToken::new();
    let reports = Arc::new(AtomicUsize::new(0));
    let (trip, count) = (token.clone(), Arc::clone(&reports));
    let reporter = ProgressReporter::new(move |_, _| {
        if count.fetch_add(1, Ordering::SeqCst) + 1 >= 2 {
            trip.cancel();
        }
    });
    let cancelled = minkowski_sum_with_progress(&l_shape, &l_shape, Some(&token), Some(&reporter));
    assert_eq!(cancelled.status, Error::Cancelled);
    assert_eq!(cancelled.num_tri(), 0);
    let n = reports.load(Ordering::SeqCst);
    assert!(
        n < 10,
        "cancel was ignored for {n} of the run's 103 reports"
    );
}
