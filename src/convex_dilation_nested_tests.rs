// Tests for convex_dilation.rs's nested/crossing-shell handling, ported 1:1
// from manifold-sharp's ConvexDilationTests.Nested.cs
// (ConvexDilationNestedTests). The oracle is the sweep of the rebuilt union
// (`rebuild_solid(Positive)`), since the raw erosion sweep carves an inner
// shell's boundary out of kept material. `REBUILDS_RUN` is the thread-local
// mirror of the C# [ThreadStatic] counter; the rebuild runs on the calling
// thread before any parallel map.

use std::sync::{Arc, Mutex};

use super::tests::l_shape;
use super::*;
use crate::linalg::Vec3;
use crate::manifold::Manifold;
use crate::types::{Error, WindingRule};

const VOLUME_TOLERANCE: f64 = 1e-9;

fn outer() -> Manifold {
    l_shape().scale(Vec3::new(2.0, 2.0, 2.0))
}

fn inner() -> Manifold {
    Manifold::cube(Vec3::splat(0.5), false).translate(Vec3::new(0.75, 0.75, 0.75))
}

fn fixture(shape: &str) -> Manifold {
    match shape {
        "hollow" => outer().difference(&inner()),
        "nested" => Manifold::compose(&[outer(), inner()]),
        "interlocked" => Manifold::compose(&[
            outer(),
            Manifold::cube(Vec3::splat(1.0), false).translate(Vec3::new(4.0, 3.5, 0.5)),
        ]),
        "crossing" => Manifold::compose(&[
            outer(),
            Manifold::cube(Vec3::splat(1.0), false).translate(Vec3::new(1.5, 3.0, 0.5)),
        ]),
        "side-by-side" => {
            Manifold::compose(&[outer(), outer().translate(Vec3::new(1.0, 1.0, 0.5))])
        }
        _ => Manifold::compose(&[
            Manifold::sphere(3.0, 16).difference(&Manifold::sphere(2.0, 16)),
            Manifold::sphere(1.0, 16),
        ]),
    }
}

fn sweep_operand(shape: &str) -> Manifold {
    fixture(shape).rebuild_solid(WindingRule::Positive)
}

fn rebuilds() -> usize {
    REBUILDS_RUN.with(|c| c.get())
}

fn a_nesting_solid_dilates_through_the_tree(shape: &str) {
    let ball = Manifold::sphere(0.3, 8);
    let reference = sweep_operand(shape).minkowski_sum(&ball);

    let before = rebuilds();
    let tree = fixture(shape).try_dilate_by_convex(&ball, None, None);
    let ran = rebuilds() - before;
    let tree = tree.expect("overlapping boxes no longer decline by themselves");
    assert_eq!(ran, 1, "every fixture has overlapping component boxes");
    assert_eq!(tree.status(), Error::NoError);
    assert_eq!(tree.genus(), reference.genus());
    let relative = (tree.volume() - reference.volume()).abs() / reference.volume();
    assert!(
        relative <= VOLUME_TOLERANCE,
        "tree {} against the sweep's {}",
        tree.volume(),
        reference.volume()
    );
}

fn a_nesting_solid_erodes_through_the_tree(shape: &str) {
    let ball = Manifold::sphere(0.3, 8);
    let reference = sweep_operand(shape).minkowski_difference(&ball);

    let before = rebuilds();
    let tree = fixture(shape).try_erode_by_convex(&ball, None, None);
    let ran = rebuilds() - before;
    let tree = tree.expect("overlapping boxes no longer decline by themselves");
    assert_eq!(ran, 1, "every fixture has overlapping component boxes");
    assert_eq!(tree.status(), Error::NoError);
    assert_eq!(tree.genus(), reference.genus());
    let relative = (tree.volume() - reference.volume()).abs() / reference.volume();
    assert!(
        relative <= VOLUME_TOLERANCE,
        "tree {} against the sweep's {}",
        tree.volume(),
        reference.volume()
    );
}

macro_rules! nesting_rows {
    ($($name:ident => $shape:expr),* $(,)?) => {
        mod dilates {
            $(#[test] fn $name() { super::a_nesting_solid_dilates_through_the_tree($shape); })*
        }
        mod erodes {
            $(#[test] fn $name() { super::a_nesting_solid_erodes_through_the_tree($shape); })*
        }
    };
}

nesting_rows! {
    hollow => "hollow",
    nested => "nested",
    interlocked => "interlocked",
    crossing => "crossing",
    side_by_side => "side-by-side",
    ball_cavity_ball => "ball-cavity-ball",
}

#[test]
fn the_raw_sweep_of_same_orientation_nesting_carves_the_inner_shell() {
    let ball = Manifold::sphere(0.3, 8);
    let raw = fixture("nested").minkowski_difference(&ball).volume();
    assert!(raw > 0.0);
    let union = sweep_operand("nested").minkowski_difference(&ball).volume();

    assert!(raw < union, "raw sweep {raw} against the union's {union}");
}

#[test]
fn the_rebuild_runs_only_when_component_boxes_overlap() {
    let ball = Manifold::sphere(0.3, 8);
    let apart = Manifold::compose(&[outer(), inner().translate(Vec3::new(20.0, 0.0, 0.0))]);

    let before = rebuilds();
    let apart_applied = apart.try_dilate_by_convex(&ball, None, None).is_some();
    let single_applied = outer().try_erode_by_convex(&ball, None, None).is_some();
    let skipped = rebuilds() - before;

    let before = rebuilds();
    let overlap_applied = fixture("side-by-side")
        .try_erode_by_convex(&ball, None, None)
        .is_some();
    let ran = rebuilds() - before;

    assert!(apart_applied && single_applied && overlap_applied);
    assert_eq!(skipped, 0, "disjoint boxes cannot overlap in winding");
    assert_eq!(ran, 1, "overlapping boxes rebuild the union once");
}

#[test]
fn the_rebuild_does_not_report_a_finished_bar() {
    let fractions: Arc<Mutex<Vec<f64>>> = Arc::new(Mutex::new(Vec::new()));
    let sink = Arc::clone(&fractions);
    let reporter = ProgressReporter::new(move |_: Phase, fraction: Option<f64>| {
        sink.lock().expect("sink").push(fraction.unwrap_or(-1.0));
    });

    let before = rebuilds();
    let applied = fixture("side-by-side")
        .try_dilate_by_convex(&Manifold::sphere(0.3, 8), None, Some(&reporter))
        .is_some();
    let ran = rebuilds() - before;

    assert!(applied);
    assert_eq!(ran, 1);
    let fractions = fractions.lock().expect("sink").clone();
    assert!(fractions.len() > 2);
    for i in 0..fractions.len() - 1 {
        assert!(
            fractions[i] < 1.0,
            "report {i} of {} already read finished",
            fractions.len()
        );
        assert!(
            fractions[i + 1] >= fractions[i],
            "the bar went backwards, {} then {}",
            fractions[i],
            fractions[i + 1]
        );
    }

    assert_eq!(fractions[fractions.len() - 1], 1.0);
}
