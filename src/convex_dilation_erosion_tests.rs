// Tests for convex_dilation.rs's erosion path (`try_erode_by_convex`), ported
// 1:1 from manifold-sharp's ConvexDilationTests.Erosion.cs
// (ConvexDilationErosionTests). The oracle is the ported sweep,
// `minkowski_difference`; volumes agree at 1e-9 relative. The C# pins the
// sequential mode for its progress-count tests; here those assertions hold in
// either build (sub-unit reports are monotone and the closing reports run on
// the calling thread after the maps join).

use std::sync::{Arc, Mutex};

use super::tests::{drilled_part, l_shape};
use super::*;
use crate::linalg::{length, Vec3};
use crate::manifold::Manifold;
use crate::types::Error;

const VOLUME_TOLERANCE: f64 = 1e-9;

pub(crate) fn frame() -> Manifold {
    Manifold::cube(Vec3::new(4.0, 4.0, 1.0), true)
        .difference(&Manifold::cube(Vec3::new(2.0, 2.0, 2.0), true))
}

pub(crate) fn thin_wall_dumbbell() -> Manifold {
    Manifold::cube(Vec3::new(2.0, 2.0, 2.0), false)
        .union(&Manifold::cube(Vec3::new(2.0, 2.0, 2.0), false).translate(Vec3::new(4.0, 0.0, 0.0)))
        .union(&Manifold::cube(Vec3::new(3.0, 0.2, 2.0), false).translate(Vec3::new(1.5, 0.9, 0.0)))
}

pub(crate) fn thin_l() -> Manifold {
    Manifold::cube(Vec3::new(3.0, 0.2, 2.0), false)
        .union(&Manifold::cube(Vec3::new(0.2, 3.0, 2.0), false))
}

fn shape(name: &str) -> Manifold {
    match name {
        "L-shape" => l_shape().scale(Vec3::splat(2.0)),
        "drilled" => drilled_part(16),
        "frame" => frame(),
        "thin-wall" => thin_wall_dumbbell(),
        "thin-L" => thin_l(),
        "cube" => Manifold::cube(Vec3::splat(2.0), true),
        _ => panic!("unknown fixture {name}"),
    }
}

fn matches_the_minkowski_difference(name: &str) {
    let solid = shape(name);
    let ball = Manifold::sphere(0.3, 8);

    let reference = solid.minkowski_difference(&ball);
    let tree = solid
        .try_erode_by_convex(&ball, None, None)
        .expect("applies");

    assert_eq!(tree.status(), Error::NoError);
    assert_eq!(tree.genus(), reference.genus());
    assert_eq!(tree.decompose().len(), reference.decompose().len());
    if name == "thin-L" {
        // Anti-vacuity for the vanishing case: the sweep really leaves nothing.
        assert!(reference.is_empty());
        assert!(tree.is_empty());
        return;
    }

    if name == "thin-wall" {
        // Anti-vacuity for the split case: the wall really is gone.
        assert_eq!(reference.decompose().len(), 2);
    }

    assert!(reference.volume() > 0.0);
    let relative = (tree.volume() - reference.volume()).abs() / reference.volume();
    assert!(
        relative <= VOLUME_TOLERANCE,
        "tree {} against the ported sweep's {}",
        tree.volume(),
        reference.volume()
    );
}

#[test]
fn matches_the_minkowski_difference_l_shape() {
    matches_the_minkowski_difference("L-shape");
}

#[test]
fn matches_the_minkowski_difference_drilled() {
    matches_the_minkowski_difference("drilled");
}

#[test]
fn matches_the_minkowski_difference_frame() {
    matches_the_minkowski_difference("frame");
}

#[test]
fn matches_the_minkowski_difference_thin_wall() {
    matches_the_minkowski_difference("thin-wall");
}

#[test]
fn matches_the_minkowski_difference_thin_l() {
    matches_the_minkowski_difference("thin-L");
}

#[test]
fn matches_the_minkowski_difference_cube() {
    matches_the_minkowski_difference("cube");
}

fn an_asymmetric_off_centre_tool_matches(tool: &str, name: &str) {
    let solid = shape(name);
    let b = if tool == "box" {
        Manifold::cube(Vec3::new(0.3, 0.2, 0.1), false).translate(Vec3::new(-0.05, -0.05, -0.03))
    } else {
        Manifold::tetrahedron()
            .scale(Vec3::new(0.3, 0.2, 0.1))
            .translate(Vec3::new(0.05, 0.0, 0.0))
    };

    let reference = solid.minkowski_difference(&b);
    let tree = solid.try_erode_by_convex(&b, None, None).expect("applies");

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

#[test]
fn non_convex_tools_and_empty_operands_are_declined() {
    let ball = Manifold::sphere(0.3, 8);
    assert!(drilled_part(8)
        .try_erode_by_convex(&l_shape(), None, None)
        .is_none());
    assert!(Manifold::empty()
        .try_erode_by_convex(&ball, None, None)
        .is_none());
    assert!(l_shape()
        .try_erode_by_convex(&Manifold::empty(), None, None)
        .is_none());
}

#[test]
fn a_pre_cancelled_token_returns_cancelled() {
    let token = crate::cancel::CancelToken::new();
    token.cancel();

    let result = l_shape()
        .try_erode_by_convex(&Manifold::sphere(0.3, 8), Some(&token), None)
        .expect("cancelled is applied");
    assert_eq!(result.status(), Error::Cancelled);
}

fn collecting_reporter() -> (ProgressReporter, Arc<Mutex<Vec<f64>>>) {
    let fractions: Arc<Mutex<Vec<f64>>> = Arc::new(Mutex::new(Vec::new()));
    let sink = Arc::clone(&fractions);
    let reporter = ProgressReporter::new(move |_: Phase, fraction: Option<f64>| {
        sink.lock()
            .expect("sink")
            .push(fraction.expect("determinate"));
    });
    (reporter, fractions)
}

#[test]
fn progress_is_monotonic_and_ends_at_one() {
    let (reporter, fractions) = collecting_reporter();
    let solid = l_shape();

    // One unit per hull, per leaf (no solid leaf when eroding), per tree
    // union, the subtraction, and the unit complete_phase spends.
    let num_tri = solid.num_tri();
    let num_leaves = num_tri.div_ceil(16);
    let total = (num_tri + num_leaves + (num_leaves - 1) + 1 + 1) as f64;
    assert!(
        total <= 100.0,
        "above 100 the throttle skips units and the penultimate report proves nothing"
    );

    assert!(solid
        .try_erode_by_convex(&Manifold::sphere(0.3, 8), None, Some(&reporter))
        .is_some());

    let fractions = fractions.lock().expect("sink").clone();
    assert!(fractions.len() > 2);
    for i in 1..fractions.len() {
        assert!(fractions[i] >= fractions[i - 1]);
    }

    assert_eq!(fractions[fractions.len() - 1], 1.0);
    assert_eq!(
        fractions[fractions.len() - 2],
        (total - 1.0) / total,
        "{num_tri} triangles in {num_leaves} leaves should cost {total} units"
    );
}

#[test]
fn the_closing_subtraction_reports_fractional_progress() {
    let solid = drilled_part(16);
    let num_tri = solid.num_tri();
    let num_leaves = num_tri.div_ceil(16);

    // Hulls, leaves, tree nodes, the subtraction and the normals pass.
    let total = (num_tri + num_leaves + (num_leaves - 1) + 1 + 1) as f64;

    let (reporter, fractions) = collecting_reporter();
    assert!(solid
        .try_erode_by_convex(&Manifold::sphere(0.3, 8), None, Some(&reporter))
        .is_some());

    let fractions = fractions.lock().expect("sink").clone();
    assert_eq!(fractions[fractions.len() - 1], 1.0);

    let subtraction_start = (total - 2.0) / total;
    let subtraction_end = (total - 1.0) / total;
    let inside: Vec<f64> = fractions
        .iter()
        .copied()
        .filter(|&f| f > subtraction_start && f < subtraction_end)
        .collect();
    assert!(
        inside.len() >= 2,
        "the closing subtraction must report between its unit boundaries"
    );
    for i in 1..inside.len() {
        assert!(inside[i] > inside[i - 1]);
    }
}
