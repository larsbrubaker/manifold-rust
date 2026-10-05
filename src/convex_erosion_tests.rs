// Tests for convex_erosion.rs, ported 1:1 from manifold-sharp's
// ConvexErosionTests.cs and ConvexErosionTests.Contract.cs (same test names in
// snake_case, same values, same tolerances; `[Arguments]` rows become one test
// per row). The subject has no C++ counterpart (docs/CPP_DIVERGENCES.md entry
// 13), so no C++ expected value appears here.
//
// TWO ORACLES, because the closed form is a second answer to a question the
// library already answers:
//  1. An INDEPENDENT one, `halfspace_intersection_by_triples`: every triple of
//     offset face planes, filtered by the DEFINITION of erosion (x - b inside
//     the solid for every tool vertex b), then hulled. O(faces^3), test only.
//  2. The GENERAL SWEEP, `Manifold::minkowski_difference`, which the fast path
//     promises to agree with.
// Volumes are compared relatively at 1e-6: the sweep's vertices carry its
// union's rounding, the closed form solves each vertex from three planes.
//
// The contract half (routing pin, every decline, progress and cancellation)
// follows the geometry half below.

// `!(x > 0.0)` is deliberate throughout: a NaN length or slack must decline,
// exactly as the C# `!(x > 0.0)` it mirrors does.
#![allow(clippy::neg_cmp_op_on_partial_ord)]

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

use crate::cancel::CancelToken;
use crate::linalg::{cross, dot, length, Vec3};
use crate::manifold::Manifold;
use crate::progress::{Phase, ProgressReporter};
use crate::quickhull::convex_hull;
use crate::types::Error;

/// Relative slop between the closed form and either oracle.
const VOLUME_TOLERANCE: f64 = 1e-6;

fn try_erode(solid: &Manifold, tool: &Manifold) -> Option<Manifold> {
    solid.try_convex_erosion(tool, None, None)
}

/// A cube eroded by the kernel's ball is a smaller cube, exactly: the ball's
/// support along an axis is exactly the radius, so every plane offset and
/// corner solve is whole-number arithmetic.
#[test]
fn a_cube_erodes_to_an_exactly_smaller_cube() {
    let cube = Manifold::cube(Vec3::splat(20.0), true);
    let ball = Manifold::sphere(1.0, 12);

    let eroded = try_erode(&cube, &ball)
        .expect("a cube and a ball are both convex, which is the whole gate");

    assert_eq!(eroded.volume(), 5832.0, "18^3 exactly");
    assert_eq!(eroded.num_tri(), 12);

    let bounds = eroded.bounding_box();
    assert_eq!(bounds.min.x, -9.0);
    assert_eq!(bounds.min.y, -9.0);
    assert_eq!(bounds.min.z, -9.0);
    assert_eq!(bounds.max.x, 9.0);
    assert_eq!(bounds.max.y, 9.0);
    assert_eq!(bounds.max.z, 9.0);
}

fn the_closed_form_is_the_halfspace_intersection(label: &str, shape: i32) {
    let solid = polyhedron(shape);
    let ball = Manifold::sphere(1.0, 12);

    let eroded = try_erode(&solid, &ball).unwrap_or_else(|| panic!("{label} is convex"));
    let expected = halfspace_intersection_by_triples(&solid, &ball);

    assert!(
        relative(eroded.volume(), expected.volume()) < VOLUME_TOLERANCE,
        "{label}: closed form {} vs brute-force halfspace intersection {}",
        eroded.volume(),
        expected.volume()
    );
    assert!(
        relative(eroded.surface_area(), expected.surface_area()) < VOLUME_TOLERANCE,
        "{label}: the surfaces have to agree too, or the volumes agreeing was luck"
    );
}

#[test]
fn the_closed_form_is_the_halfspace_intersection_cube() {
    the_closed_form_is_the_halfspace_intersection("cube", 0);
}

#[test]
fn the_closed_form_is_the_halfspace_intersection_tetrahedron() {
    the_closed_form_is_the_halfspace_intersection("tetrahedron", 1);
}

#[test]
fn the_closed_form_is_the_halfspace_intersection_icosahedron() {
    the_closed_form_is_the_halfspace_intersection("icosahedron", 2);
}

fn every_eroded_vertex_is_somewhere_the_tool_still_fits(label: &str, shape: i32) {
    let solid = polyhedron(shape);
    let ball = Manifold::sphere(1.0, 12);
    let eroded = try_erode(&solid, &ball).expect("convex");

    let (normals, offsets) = face_planes(&solid);
    let slop = length(solid.bounding_box().size()) * 1e-9;

    let verts = &eroded.as_impl().vert_pos;
    assert!(
        verts.len() > 3,
        "{label}: a solid needs four vertices before this proves anything"
    );

    for &vert in verts.iter() {
        assert!(
            tool_fits_at(vert, &normals, &offsets, &ball, slop),
            "{label}: the ball does not fit inside the solid when centred on {},{},{}",
            vert.x,
            vert.y,
            vert.z
        );
    }
}

#[test]
fn every_eroded_vertex_is_somewhere_the_tool_still_fits_cube() {
    every_eroded_vertex_is_somewhere_the_tool_still_fits("cube", 0);
}

#[test]
fn every_eroded_vertex_is_somewhere_the_tool_still_fits_tetrahedron() {
    every_eroded_vertex_is_somewhere_the_tool_still_fits("tetrahedron", 1);
}

#[test]
fn every_eroded_vertex_is_somewhere_the_tool_still_fits_icosahedron() {
    every_eroded_vertex_is_somewhere_the_tool_still_fits("icosahedron", 2);
}

/// Erode by r then dilate by r is a uniform fillet; an opening built on the
/// closed form must be the same rounded solid as one built on the sweep.
fn an_opening_on_the_closed_form_matches_an_opening_on_the_sweep(label: &str, shape: i32) {
    let solid = polyhedron(shape);
    let ball = Manifold::sphere(1.0, 12);

    let fast_eroded = try_erode(&solid, &ball).expect("convex");
    let swept_eroded = solid.minkowski_difference(&ball);

    assert!(
        relative(fast_eroded.volume(), swept_eroded.volume()) < VOLUME_TOLERANCE,
        "{label}: erosions differ — closed form {}, sweep {}",
        fast_eroded.volume(),
        swept_eroded.volume()
    );

    let fast_opening = fast_eroded.minkowski_sum(&ball);
    let swept_opening = swept_eroded.minkowski_sum(&ball);

    assert!(
        relative(fast_opening.volume(), swept_opening.volume()) < VOLUME_TOLERANCE,
        "{label}: openings differ — closed form {}, sweep {}",
        fast_opening.volume(),
        swept_opening.volume()
    );
    assert!(
        relative(fast_opening.surface_area(), swept_opening.surface_area()) < VOLUME_TOLERANCE,
        "{label}: the rounded surfaces have to agree, not just the volume they enclose"
    );
}

#[test]
fn an_opening_on_the_closed_form_matches_an_opening_on_the_sweep_cube() {
    an_opening_on_the_closed_form_matches_an_opening_on_the_sweep("cube", 0);
}

#[test]
fn an_opening_on_the_closed_form_matches_an_opening_on_the_sweep_tetrahedron() {
    an_opening_on_the_closed_form_matches_an_opening_on_the_sweep("tetrahedron", 1);
}

#[test]
fn an_opening_on_the_closed_form_matches_an_opening_on_the_sweep_icosahedron() {
    an_opening_on_the_closed_form_matches_an_opening_on_the_sweep("icosahedron", 2);
}

/// The support is taken against -n, not +n: this tool reaches 1.0 in +x and
/// only 0.3 in -x (0.8 / 0.4 in z), so a flipped sign gives the mirrored solid.
/// Asserted bit-exactly and against the sweep.
#[test]
fn an_asymmetric_tool_erodes_the_way_the_sweep_does() {
    let tool = asymmetric_tool();
    let tool_bounds = tool.bounding_box();
    assert_eq!(tool_bounds.min.x, -0.3);
    assert_eq!(tool_bounds.max.x, 1.0);
    assert_eq!(tool_bounds.min.z, -0.4);
    assert_eq!(tool_bounds.max.z, 0.8);

    let bx = Manifold::cube(Vec3::splat(40.0), true);
    let eroded = try_erode(&bx, &tool).expect("convex box, convex tool");

    let bounds = eroded.bounding_box();
    assert_eq!(
        bounds.max.x, 19.7,
        "+x face moves in by the tool's -x reach, 0.3"
    );
    assert_eq!(
        bounds.min.x, -19.0,
        "-x face moves in by the tool's +x reach, 1.0"
    );
    assert_eq!(
        bounds.max.z, 19.6,
        "+z face moves in by the tool's -z reach, 0.4"
    );
    assert_eq!(
        bounds.min.z, -19.2,
        "-z face moves in by the tool's +z reach, 0.8"
    );

    let swept = bx.minkowski_difference(&tool).bounding_box();
    assert_eq!(bounds.min.x, swept.min.x);
    assert_eq!(bounds.max.x, swept.max.x);
    assert_eq!(bounds.min.y, swept.min.y);
    assert_eq!(bounds.max.y, swept.max.y);
    assert_eq!(bounds.min.z, swept.min.z);
    assert_eq!(bounds.max.z, swept.max.z);
}

/// On a densely tessellated solid the two paths enclose the same volume
/// (measured at 8e-15) though their triangulations part company.
#[test]
fn a_dense_solid_agrees_with_the_sweep_in_volume() {
    let sphere = Manifold::sphere(10.0, 64);
    let ball = Manifold::sphere(1.0, 12);

    assert_eq!(
        sphere.num_tri(),
        2048,
        "the measured figures are for this tessellation"
    );

    let eroded = try_erode(&sphere, &ball).expect("convex");
    let swept = sphere.minkowski_difference(&ball);

    let rel = relative(eroded.volume(), swept.volume());
    assert!(
        rel < 1e-12,
        "measured at 8e-15 ({} against {})",
        eroded.volume(),
        swept.volume()
    );
}

// ── Contract half (ConvexErosionTests.Contract.cs) ──────────────────────────

/// `minkowski_difference` is not rerouted: on a 40x20x10 box the sweep leaves
/// 36 triangles, the closed form 12, with equal volumes.
#[test]
fn the_general_difference_still_runs_the_sweep() {
    let bx = Manifold::cube(Vec3::new(40.0, 20.0, 10.0), true);
    let ball = Manifold::sphere(1.0, 12);

    let swept = bx.minkowski_difference(&ball);
    let closed = try_erode(&bx, &ball).expect("convex");

    assert_eq!(swept.volume(), closed.volume(), "38 x 18 x 8 either way");
    assert_eq!(swept.num_tri(), 36, "the sweep's triangulation, unchanged");
    assert_eq!(
        closed.num_tri(),
        12,
        "the closed form hulls the eight corners"
    );
}

fn l_shape() -> Manifold {
    let cube = Manifold::cube(Vec3::splat(20.0), true);
    let notch = Manifold::cube(Vec3::splat(10.0), true).translate(Vec3::new(10.0, 10.0, 10.0));
    cube.difference(&notch)
}

#[test]
fn a_non_convex_solid_is_declined() {
    let l = l_shape();
    assert!(
        !l.as_impl().is_convex(),
        "the fixture has to actually be non-convex for the decline to mean anything"
    );

    assert!(
        try_erode(&l, &Manifold::sphere(1.0, 12)).is_none(),
        "no closed form for a solid that is not an intersection of its face halfspaces"
    );
}

/// The mathematics would allow it, but the sweep swaps its operands for a
/// convex solid and a non-convex tool, so the promise to agree does not.
#[test]
fn a_non_convex_tool_is_declined() {
    let cube = Manifold::cube(Vec3::splat(1.0), true);
    let l_shaped_tool = cube.union(&cube.translate(Vec3::new(0.6, 0.6, 0.0)));
    assert!(!l_shaped_tool.as_impl().is_convex());

    assert!(try_erode(&Manifold::cube(Vec3::splat(20.0), true), &l_shaped_tool).is_none());
}

/// A tool that does not contain the origin is where the closed form and the
/// sweep stop being the same function; (0.8, 0.8) is the regression that a
/// gate on the solid's own normals let through.
fn a_tool_that_misses_the_origin_is_declined(x: f64, y: f64, swept_volume: f64) {
    let cube = Manifold::cube(Vec3::splat(20.0), true);
    let offset_ball = Manifold::sphere(1.0, 12).translate(Vec3::new(x, y, 0.0));

    assert!(
        offset_ball.as_impl().is_convex(),
        "convexity must not be what rejects this - the origin is the point"
    );

    assert!(
        try_erode(&cube, &offset_ball).is_none(),
        "tool at ({x},{y},0): the sweep answers {swept_volume} here"
    );

    let swept = cube.minkowski_difference(&offset_ball);
    assert!(
        relative(swept.volume(), swept_volume) < 1e-9,
        "the sweep's own answer moved: {}",
        swept.volume()
    );
    assert_ne!(
        swept.volume(),
        5832.0,
        "if the sweep agreed with the erosion here there would be nothing to gate"
    );
}

#[test]
fn a_tool_that_misses_the_origin_is_declined_2_0() {
    a_tool_that_misses_the_origin_is_declined(2.0, 0.0, 5908.0);
}

#[test]
fn a_tool_that_misses_the_origin_is_declined_0_8_0_8() {
    a_tool_that_misses_the_origin_is_declined(0.8, 0.8, 5832.547441);
}

#[test]
fn a_solid_too_thin_to_hold_the_tool_is_declined() {
    let thin = Manifold::cube(Vec3::splat(1.5), true);
    let ball = Manifold::sphere(1.0, 12);

    assert!(try_erode(&thin, &ball).is_none());
    assert!(
        thin.minkowski_difference(&ball).is_empty(),
        "the sweep is what the decline hands the question to, and it answers empty"
    );
}

#[test]
fn a_pre_cancelled_token_answers_cancelled_rather_than_declining() {
    let token = CancelToken::new();
    token.cancel();

    let result = Manifold::cube(Vec3::splat(20.0), true)
        .try_convex_erosion(&Manifold::sphere(1.0, 12), Some(&token), None)
        .expect("None would send a cancelled caller off to run the minutes-long path");

    assert_eq!(result.status(), Error::Cancelled);
    assert!(result.is_empty());
}

type Event = (&'static str, Option<f64>);

#[test]
fn progress_is_reported_on_success_and_not_on_a_decline_before_any_work() {
    let events: Arc<Mutex<Vec<Event>>> = Arc::new(Mutex::new(Vec::new()));
    let sink = Arc::clone(&events);
    let reporter = ProgressReporter::new(move |phase: Phase, fraction| {
        sink.lock().expect("sink").push((phase.name(), fraction));
    });

    assert!(Manifold::cube(Vec3::splat(20.0), true)
        .try_convex_erosion(&Manifold::sphere(1.0, 12), None, Some(&reporter))
        .is_some());

    let seen = events.lock().expect("sink").clone();
    assert!(!seen.is_empty());
    let mut previous = -1.0;
    for &(name, fraction) in seen.iter() {
        assert_eq!(name, "minkowski");
        let f = fraction.expect("determinate phase");
        assert!(f >= previous, "the bar went backwards, {previous} then {f}");
        previous = f;
    }

    assert_eq!(
        seen[seen.len() - 1].1,
        Some(1.0),
        "a finished operation leaves a full bar, the rule complete_phase exists for"
    );

    events.lock().expect("sink").clear();
    assert!(l_shape()
        .try_convex_erosion(&Manifold::sphere(1.0, 12), None, Some(&reporter))
        .is_none());
    assert_eq!(
        events.lock().expect("sink").len(),
        0,
        "the sweep the caller is about to run opens this phase itself"
    );
}

/// A cancel raised from the second progress report (around face 10 of 512)
/// is observed by the support pass's next 64-face poll.
#[test]
fn a_cancel_raised_during_the_support_pass_is_observed_there() {
    let solid = Manifold::sphere(10.0, 32);
    let ball = Manifold::sphere(1.0, 12);
    assert_eq!(
        solid.num_tri(),
        512,
        "the report arithmetic is computed from this count"
    );

    let token = Arc::new(CancelToken::new());
    let reports = Arc::new(AtomicUsize::new(0));
    let cb_token = Arc::clone(&token);
    let cb_reports = Arc::clone(&reports);
    let reporter = ProgressReporter::new(move |_: Phase, _| {
        // Not the first: that is begin_phase, before any face was measured.
        if cb_reports.fetch_add(1, Ordering::SeqCst) + 1 >= 2 {
            cb_token.cancel();
        }
    });

    let result = solid
        .try_convex_erosion(&ball, Some(&token), Some(&reporter))
        .expect("a cancelled run is an applied run, or the caller falls into the sweep");

    assert_eq!(result.status(), Error::Cancelled);
    assert!(result.is_empty());
    let n = reports.load(Ordering::SeqCst);
    assert!(
        n < 30,
        "cancel was ignored for {n} of the ~103 reports a completed pass emits"
    );
}

#[test]
fn neither_operand_is_mutated() {
    let solid = polyhedron(2);
    let ball = Manifold::sphere(1.0, 12);
    let solid_before = geometry_hash(&solid);
    let ball_before = geometry_hash(&ball);

    assert!(try_erode(&solid, &ball).is_some());

    assert_eq!(geometry_hash(&solid), solid_before);
    assert_eq!(geometry_hash(&ball), ball_before);
}

// ── Helpers ─────────────────────────────────────────────────────────────────

/// 0 cube, 1 tetrahedron, 2 icosahedron — sized so a unit ball fits well inside.
fn polyhedron(shape: i32) -> Manifold {
    if shape == 0 {
        return Manifold::cube(Vec3::splat(20.0), true);
    }

    if shape == 1 {
        return Manifold::tetrahedron().scale(Vec3::splat(8.0));
    }

    // The regular icosahedron: cyclic permutations of (0, +-1, +-phi).
    let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
    let mut points: Vec<Vec3> = Vec::with_capacity(12);
    for a in [-1.0, 1.0] {
        for b in [-phi, phi] {
            points.push(Vec3::new(0.0, a, b) * 6.0);
            points.push(Vec3::new(a, b, 0.0) * 6.0);
            points.push(Vec3::new(b, 0.0, a) * 6.0);
        }
    }

    Manifold::from_impl(convex_hull(&points))
}

/// A skewed octahedron: x in [-0.3, 1], y in [-0.5, 0.6], z in [-0.4, 0.8],
/// containing the origin strictly.
fn asymmetric_tool() -> Manifold {
    let points = vec![
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(-0.3, 0.0, 0.0),
        Vec3::new(0.0, 0.6, 0.0),
        Vec3::new(0.0, -0.5, 0.0),
        Vec3::new(0.0, 0.0, 0.8),
        Vec3::new(0.0, 0.0, -0.4),
    ];
    Manifold::from_impl(convex_hull(&points))
}

/// The solid's outward unit face normals and their plane offsets.
fn face_planes(solid: &Manifold) -> (Vec<Vec3>, Vec<f64>) {
    let imp = solid.as_impl();
    let mut normals = Vec::new();
    let mut offsets = Vec::new();

    for tri in 0..imp.num_tri() {
        let mut normal = imp.face_normal[tri];
        let len = length(normal);
        if !(len > 0.0) {
            continue;
        }

        normal /= len;
        normals.push(normal);
        offsets.push(dot(
            normal,
            imp.vert_pos[imp.halfedge[tri * 3].start_vert as usize],
        ));
    }

    (normals, offsets)
}

/// The erosion by brute force: every triple of offset planes is a candidate
/// vertex, kept when the tool centred there is still inside the solid.
fn halfspace_intersection_by_triples(solid: &Manifold, tool: &Manifold) -> Manifold {
    let (normals, offsets) = face_planes(solid);
    let tool_verts = &tool.as_impl().vert_pos;
    let slop = length(solid.bounding_box().size()) * 1e-9;

    let mut pushed: Vec<f64> = Vec::with_capacity(normals.len());
    for i in 0..normals.len() {
        let mut support = f64::NEG_INFINITY;
        for &b in tool_verts.iter() {
            support = support.max(-dot(normals[i], b));
        }
        pushed.push(offsets[i] - support);
    }

    let mut survivors: Vec<Vec3> = Vec::new();
    for a in 0..normals.len() {
        for b in (a + 1)..normals.len() {
            for c in (b + 1)..normals.len() {
                let cr = cross(normals[b], normals[c]);
                let det = dot(normals[a], cr);
                if det.abs() < 1e-9 {
                    continue;
                }

                let point = ((cr * pushed[a])
                    + (cross(normals[c], normals[a]) * pushed[b])
                    + (cross(normals[a], normals[b]) * pushed[c]))
                    / det;

                if tool_fits_at(point, &normals, &offsets, tool, slop) {
                    survivors.push(point);
                }
            }
        }
    }

    Manifold::from_impl(convex_hull(&survivors))
}

/// Whether the tool placed as the sweep places it (`centre - b`) lies inside
/// the solid the planes describe.
fn tool_fits_at(
    centre: Vec3,
    normals: &[Vec3],
    offsets: &[f64],
    tool: &Manifold,
    slop: f64,
) -> bool {
    for &b in tool.as_impl().vert_pos.iter() {
        let placed = centre - b;
        for i in 0..normals.len() {
            if dot(normals[i], placed) > offsets[i] + slop {
                return false;
            }
        }
    }
    true
}

fn relative(actual: f64, expected: f64) -> f64 {
    (actual - expected).abs() / expected.abs().max(1e-30)
}

/// FNV-1a over the raw bits of every vertex coordinate and halfedge index.
fn geometry_hash(manifold: &Manifold) -> u64 {
    let mesh = manifold.as_impl();
    let mut hash: u64 = 14695981039346656037;
    let mut mix = |value: u64| {
        for shift in (0..64).step_by(8) {
            hash ^= (value >> shift) & 0xFF;
            hash = hash.wrapping_mul(1099511628211);
        }
    };

    for v in mesh.vert_pos.iter() {
        mix(v.x.to_bits());
        mix(v.y.to_bits());
        mix(v.z.to_bits());
    }

    for e in mesh.halfedge.iter() {
        mix(e.start_vert as u32 as u64);
        mix(e.end_vert as u32 as u64);
        mix(e.paired_halfedge as u32 as u64);
    }

    hash
}
