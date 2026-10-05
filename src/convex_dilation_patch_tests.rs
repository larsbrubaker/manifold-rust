// Tests for convex_patches.rs through the dilation tree, ported 1:1 from
// manifold-sharp's ConvexDilationTests.Patches.cs (ConvexDilationPatchTests).
// The thread-local hooks (`PATCH_SIZE_OVERRIDE`, `SKIP_GUARD_FOR_TESTS`,
// `LAST_HULL_COUNT`) mirror the C# [ThreadStatic] ones and are read on the
// calling thread before any parallel map starts.

use super::erosion_tests::{frame, thin_l, thin_wall_dumbbell};
use super::tests::{drilled_part, l_shape};
use super::*;
use crate::convex_patches::{faces_support, SKIP_GUARD_FOR_TESTS};
use crate::linalg::Vec3;
use crate::manifold::Manifold;
use crate::types::Error;

const VOLUME_TOLERANCE: f64 = 1e-9;

fn drilled_plate() -> Manifold {
    Manifold::cube(Vec3::new(4.0, 4.0, 1.0), true)
        .difference(&Manifold::cylinder_centered(3.0, 0.8, -1.0, 16, true))
}

fn shape(name: &str) -> Manifold {
    match name {
        "L-shape" => l_shape(),
        "drilled" => drilled_part(16),
        "drilled-plate" => drilled_plate(),
        "frame" => frame(),
        "thin-wall" => thin_wall_dumbbell(),
        "thin-L" => thin_l(),
        "hollow" => l_shape().scale(Vec3::splat(2.0)).difference(
            &Manifold::cube(Vec3::splat(0.5), false).translate(Vec3::new(0.75, 0.75, 0.75)),
        ),
        _ => panic!("unknown fixture {name}"),
    }
}

/// Restores the thread-local hooks even when the dilation panics.
struct HookGuard;

impl Drop for HookGuard {
    fn drop(&mut self) {
        PATCH_SIZE_OVERRIDE.with(|c| c.set(None));
        SKIP_GUARD_FOR_TESTS.with(|c| c.set(false));
    }
}

fn dilate(
    solid: &Manifold,
    tool: &Manifold,
    patch_size: Option<usize>,
    skip_guard: bool,
) -> (Manifold, usize) {
    PATCH_SIZE_OVERRIDE.with(|c| c.set(patch_size));
    SKIP_GUARD_FOR_TESTS.with(|c| c.set(skip_guard));
    let _guard = HookGuard;
    let result = solid
        .try_dilate_by_convex(tool, None, None)
        .expect("declined");
    (result, LAST_HULL_COUNT.with(|c| c.get()))
}

fn patches_match_the_per_triangle_tree(name: &str) {
    let solid = shape(name);
    let ball = Manifold::sphere(0.3, 8);
    let (reference, triangles) = dilate(&solid, &ball, Some(1), false);
    let (patched, hulls) = dilate(&solid, &ball, None, false);

    assert!(
        hulls < triangles,
        "no patch formed, so the comparison proves nothing"
    );
    assert_eq!(patched.status(), Error::NoError);
    assert_eq!(patched.genus(), reference.genus());
    let relative = (patched.volume() - reference.volume()).abs() / reference.volume();
    assert!(
        relative <= VOLUME_TOLERANCE,
        "patched {} against per-triangle {} ({hulls} hulls for {triangles} triangles)",
        patched.volume(),
        reference.volume()
    );
}

#[test]
fn patches_match_the_per_triangle_tree_l_shape() {
    patches_match_the_per_triangle_tree("L-shape");
}

#[test]
fn patches_match_the_per_triangle_tree_drilled() {
    patches_match_the_per_triangle_tree("drilled");
}

#[test]
fn patches_match_the_per_triangle_tree_drilled_plate() {
    patches_match_the_per_triangle_tree("drilled-plate");
}

#[test]
fn patches_match_the_per_triangle_tree_frame() {
    patches_match_the_per_triangle_tree("frame");
}

#[test]
fn patches_match_the_per_triangle_tree_thin_wall() {
    patches_match_the_per_triangle_tree("thin-wall");
}

#[test]
fn patches_match_the_per_triangle_tree_thin_l() {
    patches_match_the_per_triangle_tree("thin-L");
}

#[test]
fn patches_match_the_per_triangle_tree_hollow() {
    patches_match_the_per_triangle_tree("hollow");
}

#[test]
fn a_zero_area_hull_face_is_refused() {
    let o = Vec3::new(0.0, 0.0, 0.0);
    let x = Vec3::new(1.0, 0.0, 0.0);
    let y = Vec3::new(0.0, 1.0, 0.0);
    let z = Vec3::new(0.0, 0.0, 1.0);
    let points = vec![o, x, y, z];

    // The unit tetrahedron, outward counterclockwise.
    let tetra = [o, y, x, o, x, z, o, z, y, x, y, z];
    assert!(faces_support(&tetra, 4, &points));

    // The same faces plus one collapsed onto a segment.
    let with_sliver = [o, y, x, o, x, z, o, z, y, x, y, z, o, x, x];
    assert!(!faces_support(&with_sliver, 5, &points));
}

#[test]
fn the_guard_keeps_a_drilled_hole_open() {
    let solid = drilled_plate();
    let ball = Manifold::sphere(0.1, 8);

    let (guarded, _) = dilate(&solid, &ball, None, false);
    assert_eq!(guarded.genus(), 1);

    let (unguarded, _) = dilate(&solid, &ball, None, true);
    assert_eq!(
        unguarded.genus(),
        0,
        "without the guard the fixture must close the bore, or it does not test the guard"
    );
}
