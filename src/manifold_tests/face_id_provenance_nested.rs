// Face-ID provenance through nested booleans, modelled on a design that lost
// its colours: a thin-walled open box (a shell cut from a block) with eight
// cylinders standing on its inner floor, re-imported with the colours the
// boolean recovered, drilled with holes and given small blocks, then
// re-imported and cut once more; and one flat face carrying two colours, cut
// across the colour seam. Every level checks the invariant from
// `face_id_provenance.rs`, the inner floor must keep the shell's colour, and
// no triangle may straddle a colour seam. Both cases take the boolean engine;
// `face_id_provenance_robust.rs` runs them on the robust one.

use super::face_id_provenance::*;
use super::*;
use crate::linalg::Vec2;
use crate::types::BooleanEngine;
use crate::types::{MeshIndex, MeshPrecision};

const FLOOR: f64 = 1.0;
/// The inner floor's footprint, inside the walls.
const INNER: (Vec2, Vec2) = (Vec2::new(2.0, 2.0), Vec2::new(38.0, 28.0));

pub(super) fn shell_with_posts<K: GlKind>(tol: f64, engine: BooleanEngine) {
    let mut reg = Registry::default();

    // Level 1: shell (face 0) unioned with posts (face 1) on the inner floor.
    let block = Manifold::cube(Vec3::new(40.0, 30.0, 15.0), false);
    let cavity =
        Manifold::cube(Vec3::new(36.0, 26.0, 15.0), false).translate(Vec3::new(2.0, 2.0, FLOOR));
    let shell = block.difference_with_engine(&cavity, engine);
    let (shell, shell_id) = reg.import::<K>(&shell, solid(0));
    let mut parts = vec![shell];
    for i in 0..4 {
        for j in 0..2 {
            let post = Manifold::cylinder(8.0, 2.5, 2.5, 32).translate(Vec3::new(
                7.0 + 8.5 * i as f64,
                9.0 + 12.0 * j as f64,
                FLOOR,
            ));
            parts.push(reg.import::<K>(&post, solid(1)).0);
        }
    }
    let level1 = Manifold::batch_boolean_with_engine(&parts, OpType::Add, engine);
    let mesh1 = reg.check::<K>(&level1, tol, "level 1 shell + posts");
    assert_floor_face(&mesh1, FLOOR, INNER, tol, shell_id, 0);

    // Level 2: the result, re-imported with its recovered colours, drilled
    // through every post and the floor (face 2), plus blocks standing on the
    // floor in the corners (face 3).
    let (level1, l1_id) = reg.import_mesh::<K>(mesh1);
    let mut cut = level1;
    for i in 0..4 {
        for j in 0..2 {
            let hole = Manifold::cylinder(20.0, 1.0, 1.0, 16).translate(Vec3::new(
                7.0 + 8.5 * i as f64,
                9.0 + 12.0 * j as f64,
                -2.0,
            ));
            let (hole, _) = reg.import::<K>(&hole, solid(2));
            cut = cut.difference_with_engine(&hole, engine);
        }
    }
    for (x, y) in [(2.0, 2.0), (35.0, 2.0), (2.0, 25.0), (35.0, 25.0)] {
        let pad = Manifold::cube(Vec3::splat(3.0), false).translate(Vec3::new(x, y, FLOOR));
        let (pad, _) = reg.import::<K>(&pad, solid(3));
        cut = cut.union_with_engine(&pad, engine);
    }
    let mesh2 = reg.check::<K>(&cut, tol, "level 2 holes + pads");
    assert_floor_face(&mesh2, FLOOR, INNER, tol, l1_id, 0);

    // Level 3: re-imported again, a slot cut through one wall flush with the
    // floor and a rotated rib unioned across the cavity.
    let (level2, l2_id) = reg.import_mesh::<K>(mesh2);
    let slot =
        Manifold::cube(Vec3::new(6.0, 4.0, 5.0), false).translate(Vec3::new(17.0, -1.0, FLOOR));
    let (slot, _) = reg.import::<K>(&slot, box_palette);
    let rib = Manifold::cube(Vec3::new(30.0, 1.0, 4.0), true)
        .rotate(0.0, 0.0, 20.0)
        .translate(Vec3::new(20.0, 15.0, FLOOR + 2.0));
    let (rib, _) = reg.import::<K>(&rib, solid(4));
    let level3 = level2
        .difference_with_engine(&slot, engine)
        .union_with_engine(&rib, engine);
    let mesh3 = reg.check::<K>(&level3, tol, "level 3 slot + rib");
    assert_floor_face(&mesh3, FLOOR, INNER, tol, l2_id, 0);
}

#[test]
fn face_id_provenance_nested_shell_64() {
    shell_with_posts::<Gl64>(1e-9, EXACT);
}

#[test]
fn face_id_provenance_nested_shell_32() {
    shell_with_posts::<Gl32>(1e-4, EXACT);
}

/// One flat face, two colours: two blocks with flush tops are unioned, the
/// result re-imported with its recovered colours (so one coplanar top spans
/// two face IDs), then cut across the seam. Simplification must not merge
/// triangles across the colour seam even though they are coplanar.
pub(super) fn two_colour_plane<K: GlKind>(tol: f64, engine: BooleanEngine) {
    let mut reg = Registry::default();
    let (a, _) = reg.import::<K>(&Manifold::cube(Vec3::new(10.0, 10.0, 4.0), false), solid(0));
    let b = Manifold::cube(Vec3::new(10.0, 10.0, 4.0), false).translate(Vec3::new(10.0, 0.0, 0.0));
    let (b, _) = reg.import::<K>(&b, solid(1));
    let joined = a.union_with_engine(&b, engine);
    let mesh = reg.check::<K>(&joined, tol, "flush blocks");
    let (joined, _) = reg.import_mesh::<K>(mesh);
    let mut cut = joined;
    for (x, y) in [(10.0, 8.0), (5.0, 5.0)] {
        let hole = Manifold::cylinder(10.0, 1.5, 1.5, 24).translate(Vec3::new(x, y, -2.0));
        let (hole, _) = reg.import::<K>(&hole, solid(2));
        cut = cut.difference_with_engine(&hole, engine);
    }
    let notch =
        Manifold::cube(Vec3::new(8.0, 3.0, 2.0), false).translate(Vec3::new(6.0, -1.0, 2.0));
    let (notch, _) = reg.import::<K>(&notch, box_palette);
    cut = cut.difference_with_engine(&notch, engine);
    // A groove across the seam: its walls meet the top at new vertices on the
    // seam that lie on straight lines, which simplification would remove if
    // it ignored the colours.
    let groove =
        Manifold::cube(Vec3::new(10.0, 1.5, 3.0), false).translate(Vec3::new(5.0, 4.0, 3.0));
    let (groove, _) = reg.import::<K>(&groove, box_palette);
    cut = cut.difference_with_engine(&groove, engine);
    let mesh = reg.check::<K>(&cut, tol, "cut across the colour seam");
    // Every top triangle left of the seam is colour 0, right of it colour 1.
    let pos: Vec<Vec3> = (0..mesh.vert_properties.len() / mesh.num_prop.to_u64() as usize)
        .map(|i| {
            let np = mesh.num_prop.to_u64() as usize;
            let p = |k: usize| mesh.vert_properties[np * i + k].to_f64();
            Vec3::new(p(0), p(1), p(2))
        })
        .collect();
    for tri in 0..mesh.tri_verts.len() / 3 {
        let v: Vec<Vec3> = (0..3)
            .map(|j| pos[mesh.tri_verts[3 * tri + j].to_u64() as usize])
            .collect();
        if !v.iter().all(|p| (p.z - 4.0).abs() <= tol) {
            continue;
        }
        let want = if v.iter().all(|p| p.x <= 10.0 + tol) {
            0
        } else {
            1
        };
        assert!(
            v.iter().all(|p| p.x <= 10.0 + tol) || v.iter().all(|p| p.x >= 10.0 - tol),
            "top triangle {tri} {v:?} straddles the colour seam"
        );
        assert_eq!(mesh.face_id[tri].to_u64(), want, "top triangle {tri} {v:?}");
    }
}

#[test]
fn face_id_provenance_two_colour_plane_64() {
    two_colour_plane::<Gl64>(1e-9, EXACT);
}

#[test]
fn face_id_provenance_two_colour_plane_32() {
    two_colour_plane::<Gl32>(1e-4, EXACT);
}
