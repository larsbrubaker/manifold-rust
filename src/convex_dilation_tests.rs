// Tests for convex_dilation.rs, ported 1:1 from manifold-sharp's
// ConvexDilationTests.cs (same names in snake_case, values and tolerances;
// `[Arguments]` rows become one test per row). The oracle is the ported
// `minkowski_sum`; volumes agree at 1e-9 relative, not bit for bit, because
// the same hulls unioned in a different order round differently.
//
// manifold-sharp toggles `ManifoldParallel.Enabled` at runtime inside a test;
// here the `parallel` feature decides per build, so the `(bool parallel)` rows
// collapse into one test that runs in whichever mode the build has.
// The fixtures (`l_shape`, `drilled_part`) are shared with the sibling test
// files of this module.

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

use super::*;
use crate::linalg::{length, Vec3};
use crate::manifold::Manifold;
use crate::types::Error;

const VOLUME_TOLERANCE: f64 = 1e-9;

pub(crate) fn l_shape() -> Manifold {
    Manifold::cube(Vec3::new(4.0, 1.0, 1.0), false)
        .union(&Manifold::cube(Vec3::new(1.0, 3.0, 1.0), false))
}

pub(crate) fn drilled_part(segments: i32) -> Manifold {
    Manifold::cube(Vec3::splat(4.0), true)
        .difference(&Manifold::sphere(1.5, segments).translate(Vec3::new(2.0, 2.0, 2.0)))
        .difference(&Manifold::cylinder_centered(6.0, 0.8, -1.0, segments, true))
}

fn shape_named(shape: &str) -> Manifold {
    if shape == "L-shape" {
        l_shape()
    } else {
        drilled_part(16)
    }
}

fn matches_the_minkowski_sum_on_a_non_convex_solid(shape: &str) {
    let solid = shape_named(shape);
    let ball = Manifold::sphere(0.3, 8);

    let reference = solid.minkowski_sum(&ball);
    let tree = solid
        .try_dilate_by_convex(&ball, None, None)
        .expect("non-convex ⊕ convex applies");

    assert_eq!(tree.status(), Error::NoError);
    assert_eq!(tree.genus(), reference.genus());
    let relative = (tree.volume() - reference.volume()).abs() / reference.volume();
    assert!(
        relative <= VOLUME_TOLERANCE,
        "tree {} against the ported sum's {}",
        tree.volume(),
        reference.volume()
    );
}

#[test]
fn matches_the_minkowski_sum_on_a_non_convex_solid_l_shape() {
    matches_the_minkowski_sum_on_a_non_convex_solid("L-shape");
}

#[test]
fn matches_the_minkowski_sum_on_a_non_convex_solid_drilled() {
    matches_the_minkowski_sum_on_a_non_convex_solid("drilled");
}

fn an_asymmetric_off_centre_tool_matches(tool: &str, shape: &str) {
    let solid = shape_named(shape);
    let b = if tool == "box" {
        Manifold::cube(Vec3::new(0.3, 0.2, 0.1), false).translate(Vec3::new(-0.05, -0.05, -0.03))
    } else {
        Manifold::tetrahedron()
            .scale(Vec3::new(0.3, 0.2, 0.1))
            .translate(Vec3::new(0.05, 0.0, 0.0))
    };

    let reference = solid.minkowski_sum(&b);
    let tree = solid.try_dilate_by_convex(&b, None, None).expect("applies");

    assert_eq!(tree.status(), Error::NoError);
    assert_eq!(tree.genus(), reference.genus());
    let relative = (tree.volume() - reference.volume()).abs() / reference.volume();
    assert!(
        relative <= VOLUME_TOLERANCE,
        "tree {} against the ported path's {}",
        tree.volume(),
        reference.volume()
    );
    let tree_box = tree.bounding_box();
    let reference_box = reference.bounding_box();
    assert!(length(tree_box.min - reference_box.min) <= 1e-9);
    assert!(length(tree_box.max - reference_box.max) <= 1e-9);
}

#[test]
fn an_asymmetric_off_centre_tool_matches_box_l_shape() {
    an_asymmetric_off_centre_tool_matches("box", "L-shape");
}

#[test]
fn an_asymmetric_off_centre_tool_matches_box_drilled() {
    an_asymmetric_off_centre_tool_matches("box", "drilled");
}

#[test]
fn an_asymmetric_off_centre_tool_matches_tetrahedron_l_shape() {
    an_asymmetric_off_centre_tool_matches("tetrahedron", "L-shape");
}

#[test]
fn an_asymmetric_off_centre_tool_matches_tetrahedron_drilled() {
    an_asymmetric_off_centre_tool_matches("tetrahedron", "drilled");
}

fn a_convex_solid_is_declined(shape: &str) {
    let solid = if shape == "cube" {
        Manifold::cube(Vec3::splat(2.0), true)
    } else {
        Manifold::sphere(1.0, 16)
    };

    assert!(solid
        .try_dilate_by_convex(&Manifold::sphere(0.3, 8), None, None)
        .is_none());
}

#[test]
fn a_convex_solid_is_declined_cube() {
    a_convex_solid_is_declined("cube");
}

#[test]
fn a_convex_solid_is_declined_sphere() {
    a_convex_solid_is_declined("sphere");
}

#[test]
fn a_non_convex_tool_is_declined() {
    assert!(drilled_part(8)
        .try_dilate_by_convex(&l_shape(), None, None)
        .is_none());
}

#[test]
fn a_cancel_mid_run_returns_cancelled() {
    let token = Arc::new(crate::cancel::CancelToken::new());
    let reports = Arc::new(AtomicUsize::new(0));
    let cb_token = Arc::clone(&token);
    let cb_reports = Arc::clone(&reports);

    // Cancel on the second report, after begin_phase's opening one.
    let reporter = ProgressReporter::new(move |_: Phase, _| {
        if cb_reports.fetch_add(1, Ordering::SeqCst) + 1 == 2 {
            cb_token.cancel();
        }
    });

    let result = drilled_part(16)
        .try_dilate_by_convex(&Manifold::sphere(0.3, 8), Some(&token), Some(&reporter))
        .expect("a cancelled run is an applied run");

    assert_eq!(result.status(), Error::Cancelled);
    assert!(result.is_empty());
}

#[test]
fn a_pre_cancelled_token_returns_cancelled() {
    let token = crate::cancel::CancelToken::new();
    token.cancel();

    let result = l_shape()
        .try_dilate_by_convex(&Manifold::sphere(0.3, 8), Some(&token), None)
        .expect("cancelled is applied");
    assert_eq!(result.status(), Error::Cancelled);
}

type Event = (&'static str, Option<f64>);

#[test]
fn progress_is_monotonic_and_ends_at_one() {
    let events: Arc<Mutex<Vec<Event>>> = Arc::new(Mutex::new(Vec::new()));
    let sink = Arc::clone(&events);
    let reporter = ProgressReporter::new(move |phase: Phase, fraction| {
        sink.lock().expect("sink").push((phase.name(), fraction));
    });

    assert!(drilled_part(16)
        .try_dilate_by_convex(&Manifold::sphere(0.3, 8), None, Some(&reporter))
        .is_some());

    let seen = events.lock().expect("sink").clone();
    assert!(seen.len() > 2);
    let mut previous = -1.0;
    for &(name, fraction) in seen.iter() {
        assert_eq!(name, "minkowski");
        let f = fraction.expect("determinate");
        assert!(f >= previous, "the bar went backwards, {previous} then {f}");
        previous = f;
    }

    assert_eq!(seen[seen.len() - 1].1, Some(1.0));

    events.lock().expect("sink").clear();
    assert!(Manifold::cube(Vec3::splat(2.0), true)
        .try_dilate_by_convex(&Manifold::sphere(0.3, 8), None, Some(&reporter))
        .is_none());
    assert_eq!(
        events.lock().expect("sink").len(),
        0,
        "a decline happens before the phase opens"
    );
}

#[test]
fn top_unions_report_fractional_progress_from_inside() {
    let solid = drilled_part(16);

    let fractions: Arc<Mutex<Vec<f64>>> = Arc::new(Mutex::new(Vec::new()));
    let sink = Arc::clone(&fractions);
    let reporter = ProgressReporter::new(move |_: Phase, fraction: Option<f64>| {
        sink.lock()
            .expect("sink")
            .push(fraction.expect("determinate"));
    });

    assert!(solid
        .try_dilate_by_convex(&Manifold::sphere(0.3, 8), None, Some(&reporter))
        .is_some());
    let fractions = fractions.lock().expect("sink").clone();
    assert_eq!(fractions[fractions.len() - 1], 1.0);

    // One unit per hull, then per leaf, per union and the closing pass.
    let num_hulls = LAST_HULL_COUNT.with(|c| c.get()) as f64;
    let num_leaves = ((num_hulls as usize + 15) / 16 + 1) as f64;
    let total = num_hulls + num_leaves + (num_leaves - 1.0) + 1.0;

    let is_fractional = |fraction: f64| {
        let units = fraction * total;
        (units - units.round()).abs() > 1e-6
    };

    assert!(
        fractions.iter().any(|&f| is_fractional(f)),
        "the top unions must report between their unit boundaries"
    );

    // The top union (the last node before the closing unit) must be heard inside.
    let top_unit_start = (total - 2.0) / total;
    let top_unit_end = (total - 1.0) / total;
    let inside_top: Vec<f64> = fractions
        .iter()
        .copied()
        .filter(|&f| f > top_unit_start && f < top_unit_end)
        .collect();
    assert!(inside_top.len() >= 2);
    for i in 1..inside_top.len() {
        assert!(inside_top[i] > inside_top[i - 1]);
    }

    for i in 1..fractions.len() {
        assert!(
            fractions[i] >= fractions[i - 1],
            "the bar went backwards at report {i}, {} then {}",
            fractions[i - 1],
            fractions[i]
        );
    }
}
