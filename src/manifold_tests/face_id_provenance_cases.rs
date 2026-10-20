// Face-ID provenance on single booleans: box and cylinder operands with
// palette face IDs under union, difference and intersection, coplanar
// touching cases (shared planes, flush cuts, a cylinder standing on a floor),
// and batch booleans of translated and rotated operands, each on both the
// f64 MeshGL64 and the f32 MeshGL path. Each case takes the boolean engine;
// the tests here run the exact engine and `face_id_provenance_robust.rs` runs
// the same cases on the robust one. The harness and the invariant live in
// `face_id_provenance.rs`.

use super::face_id_provenance::*;
use super::*;
use crate::linalg::Vec2;
use crate::types::BooleanEngine;

const TOL64: f64 = 1e-9;
const TOL32: f64 = 1e-4;

fn box_and_cylinder<K: GlKind>(reg: &mut Registry) -> (Manifold, Manifold) {
    let (b, _) = reg.import::<K>(&Manifold::cube(Vec3::splat(10.0), false), box_palette);
    let cyl = Manifold::cylinder(14.0, 3.0, 3.0, 32).translate(Vec3::new(5.0, 5.0, -2.0));
    let (c, _) = reg.import::<K>(&cyl, cylinder_palette);
    (b, c)
}

pub(super) fn box_cylinder_ops<K: GlKind>(tol: f64, engine: BooleanEngine) {
    let mut reg = Registry::default();
    let (b, c) = box_and_cylinder::<K>(&mut reg);
    for op in [OpType::Add, OpType::Subtract, OpType::Intersect] {
        reg.check::<K>(
            &b.boolean_with_engine(&c, op, engine),
            tol,
            &format!("box {op:?} cyl"),
        );
        reg.check::<K>(
            &c.boolean_with_engine(&b, op, engine),
            tol,
            &format!("cyl {op:?} box"),
        );
    }
    // A cylinder partly outside the box, so curved side meets flat faces.
    let side = c.translate(Vec3::new(4.0, 0.0, 0.0)).rotate(0.0, 0.0, 15.0);
    for op in [OpType::Add, OpType::Subtract, OpType::Intersect] {
        reg.check::<K>(
            &b.boolean_with_engine(&side, op, engine),
            tol,
            &format!("box {op:?} offset cyl"),
        );
    }
}

#[test]
fn face_id_provenance_box_cylinder_ops_64() {
    box_cylinder_ops::<Gl64>(TOL64, EXACT);
}

#[test]
fn face_id_provenance_box_cylinder_ops_32() {
    box_cylinder_ops::<Gl32>(TOL32, EXACT);
}

pub(super) fn coplanar_touching<K: GlKind>(tol: f64, engine: BooleanEngine) {
    let mut reg = Registry::default();
    // A box resting on a larger box: the top of one is the bottom of the other.
    let (base, _) = reg.import::<K>(
        &Manifold::cube(Vec3::new(10.0, 10.0, 5.0), false),
        box_palette,
    );
    let top = Manifold::cube(Vec3::splat(4.0), false).translate(Vec3::new(3.0, 3.0, 5.0));
    let (top, _) = reg.import::<K>(&top, solid(7));
    reg.check::<K>(
        &base.boolean_with_engine(&top, OpType::Add, engine),
        tol,
        "box on box",
    );
    // The same box hanging over an edge: shared plane and shared side planes.
    let over = top.translate(Vec3::new(5.0, -3.0, 0.0));
    reg.check::<K>(
        &base.boolean_with_engine(&over, OpType::Add, engine),
        tol,
        "box over edge",
    );

    // A notch whose cutter top is flush with the base's top and whose ends
    // are flush with the base's sides.
    let cutter =
        Manifold::cube(Vec3::new(4.0, 10.0, 3.0), false).translate(Vec3::new(3.0, 0.0, 2.0));
    let (cutter, _) = reg.import::<K>(&cutter, box_palette);
    reg.check::<K>(
        &base.boolean_with_engine(&cutter, OpType::Subtract, engine),
        tol,
        "flush notch",
    );
    reg.check::<K>(
        &base.boolean_with_engine(&cutter, OpType::Intersect, engine),
        tol,
        "flush overlap",
    );

    // A cylinder standing exactly on a box floor.
    let (floor, floor_id) =
        reg.import::<K>(&Manifold::cube(Vec3::new(20.0, 20.0, 2.0), false), solid(0));
    let cyl = Manifold::cylinder(5.0, 3.0, 3.0, 32).translate(Vec3::new(10.0, 10.0, 2.0));
    let (cyl, _) = reg.import::<K>(&cyl, solid(1));
    let out = floor.boolean_with_engine(&cyl, OpType::Add, engine);
    let mesh = reg.check::<K>(&out, tol, "cylinder on floor");
    assert_floor_face(
        &mesh,
        2.0,
        (Vec2::new(0.0, 0.0), Vec2::new(20.0, 20.0)),
        tol,
        floor_id,
        0,
    );
}

#[test]
fn face_id_provenance_coplanar_touching_64() {
    coplanar_touching::<Gl64>(TOL64, EXACT);
}

#[test]
fn face_id_provenance_coplanar_touching_32() {
    coplanar_touching::<Gl32>(TOL32, EXACT);
}

pub(super) fn batch<K: GlKind>(tol: f64, engine: BooleanEngine) {
    let mut reg = Registry::default();
    let mut parts = Vec::new();
    let (b, _) = reg.import::<K>(
        &Manifold::cube(Vec3::new(30.0, 30.0, 4.0), false),
        box_palette,
    );
    parts.push(b);
    for i in 0..6 {
        let x = 4.0 + 4.5 * i as f64;
        let cyl = Manifold::cylinder(6.0, 1.5, 1.5, 24);
        let (cyl, _) = reg.import::<K>(&cyl, cylinder_palette);
        // Moved and turned after import, so the run transforms are not identity.
        parts.push(
            cyl.rotate(0.0, 0.0, 10.0 * i as f64)
                .translate(Vec3::new(x, 8.0, 4.0)),
        );
        let cube = Manifold::cube(Vec3::splat(3.0), true);
        let (cube, _) = reg.import::<K>(&cube, box_palette);
        parts.push(
            cube.rotate(15.0 * i as f64, 30.0, 0.0)
                .translate(Vec3::new(x, 20.0, 4.0)),
        );
    }
    let union = Manifold::batch_boolean_with_engine(&parts, OpType::Add, engine);
    reg.check::<K>(&union, tol, "batch union");
    let cut = Manifold::batch_boolean_with_engine(&parts, OpType::Subtract, engine);
    reg.check::<K>(&cut, tol, "batch difference");
}

#[test]
fn face_id_provenance_batch_64() {
    batch::<Gl64>(TOL64, EXACT);
}

#[test]
fn face_id_provenance_batch_32() {
    batch::<Gl32>(TOL32, EXACT);
}

/// The import route keeps the input face IDs and reports the reserved
/// original ID per run; `original_id()` stays -1 because an imported mesh is
/// not itself an original, as in the reference.
#[test]
fn face_id_import_keeps_face_ids_and_run_original_id() {
    let mut reg = Registry::default();
    let (b, id) = reg.import::<Gl64>(&Manifold::cube(Vec3::splat(1.0), false), box_palette);
    assert_eq!(b.original_id(), -1);
    let mesh = b.get_mesh_gl64(-1);
    assert_eq!(mesh.run_original_id, vec![id]);
    let mut faces: Vec<u64> = mesh.face_id.clone();
    faces.sort_unstable();
    faces.dedup();
    assert_eq!(faces, vec![0, 1, 2]);
    reg.check::<Gl64>(&b, TOL64, "plain import");
}
