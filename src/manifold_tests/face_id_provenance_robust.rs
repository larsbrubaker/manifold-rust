// Face-ID provenance under the robust (rational-arithmetic) boolean engine:
// the single-boolean and nested cases of `face_id_provenance_cases.rs` and
// `face_id_provenance_nested.rs` forced onto `BooleanEngine::Robust`; `Auto`
// with operands that send it to the robust engine on its own (two boxes
// sharing an edge, imported as a soup, and two overlapping boxes in one
// self-intersecting mesh, both painted with face IDs that span many
// non-coplanar triangles); the subtracted operand coming back as backside
// runs; and a cross-check that for manifold operands the fast and robust
// engines paint every point of the result with the same colour. Each runs on
// both the f64 MeshGL64 and the f32 MeshGL path, through the harness in
// `face_id_provenance.rs`.

use super::face_id_provenance::*;
use super::face_id_provenance_cases::{batch, box_cylinder_ops, coplanar_touching};
use super::face_id_provenance_nested::{shell_with_posts, two_colour_plane};
use super::*;
use crate::linalg::{cross, dot, length, normalize, Vec2};
use crate::types::{BooleanEngine, MeshGLP, MeshIndex, MeshPrecision};

const TOL64: f64 = 1e-9;
const TOL32: f64 = 1e-4;
const AUTO: BooleanEngine = BooleanEngine::Auto;

#[test]
fn face_id_provenance_robust_box_cylinder_ops_64() {
    box_cylinder_ops::<Gl64>(TOL64, ROBUST);
}

#[test]
fn face_id_provenance_robust_box_cylinder_ops_32() {
    box_cylinder_ops::<Gl32>(TOL32, ROBUST);
}

#[test]
fn face_id_provenance_robust_coplanar_touching_64() {
    coplanar_touching::<Gl64>(TOL64, ROBUST);
}

#[test]
fn face_id_provenance_robust_coplanar_touching_32() {
    coplanar_touching::<Gl32>(TOL32, ROBUST);
}

#[test]
fn face_id_provenance_robust_batch_64() {
    batch::<Gl64>(TOL64, ROBUST);
}

#[test]
fn face_id_provenance_robust_batch_32() {
    batch::<Gl32>(TOL32, ROBUST);
}

#[test]
fn face_id_provenance_robust_nested_shell_64() {
    shell_with_posts::<Gl64>(TOL64, ROBUST);
}

#[test]
fn face_id_provenance_robust_nested_shell_32() {
    shell_with_posts::<Gl32>(TOL32, ROBUST);
}

#[test]
fn face_id_provenance_robust_two_colour_plane_64() {
    two_colour_plane::<Gl64>(TOL64, ROBUST);
}

#[test]
fn face_id_provenance_robust_two_colour_plane_32() {
    two_colour_plane::<Gl32>(TOL32, ROBUST);
}

/// Concatenate the exported meshes of `parts` into one positions-only mesh.
/// With `unindexed`, every triangle corner is its own vertex — a soup no
/// edge of which pairs by index, the shape an STL import hands over — else
/// every part keeps its own shared vertices (closed shells that may cross).
fn concat<K: GlKind>(parts: &[Manifold], unindexed: bool) -> MeshGLP<K::P, K::I> {
    let mut out = MeshGLP::<K::P, K::I> {
        num_prop: K::I::from_usize(3),
        ..Default::default()
    };
    let push = |out: &mut MeshGLP<K::P, K::I>, p: Vec3| -> usize {
        out.vert_properties
            .extend([p.x, p.y, p.z].map(<K::P as MeshPrecision>::from_f64));
        out.vert_properties.len() / 3 - 1
    };
    for part in parts {
        let mesh = K::export(part);
        let pos = positions(&mesh);
        if unindexed {
            for &v in &mesh.tri_verts {
                let id = push(&mut out, pos[v.to_u64() as usize]);
                out.tri_verts.push(K::I::from_usize(id));
            }
        } else {
            let ids: Vec<usize> = pos.iter().map(|&p| push(&mut out, p)).collect();
            out.tri_verts.extend(
                mesh.tri_verts
                    .iter()
                    .map(|&v| K::I::from_usize(ids[v.to_u64() as usize])),
            );
        }
    }
    out
}

/// Paint every triangle of `mesh` with `face_of(unit normal, centroid)`.
fn paint<K: GlKind>(
    mut mesh: MeshGLP<K::P, K::I>,
    face_of: impl Fn(Vec3, Vec3) -> u64,
) -> MeshGLP<K::P, K::I> {
    let pos = positions(&mesh);
    mesh.face_id = (0..mesh.tri_verts.len() / 3)
        .map(|t| {
            let v = tri_verts(&mesh, &pos, t);
            let c = (v[0] + v[1] + v[2]) / 3.0;
            K::I::from_usize(face_of(normalize(tri_normal(&v)), c) as usize)
        })
        .collect();
    mesh
}

/// A palette whose IDs span many non-coplanar triangles, the way a
/// MatterCAD part's colour indexes cover whole regions: left of `x = split`
/// is colour 5, right of it colour 6, except downward faces (colour 7).
fn region_palette(split: f64) -> impl Fn(Vec3, Vec3) -> u64 {
    move |n, c| {
        if n.z < -0.5 {
            7
        } else if c.x < split {
            5
        } else {
            6
        }
    }
}

/// Assert `mesh` has a run of original `orig` and every such run is a
/// backside run (run_flags bit 0): the subtracted operand's faces are flipped.
fn assert_backside<P: MeshPrecision, I: MeshIndex>(mesh: &MeshGLP<P, I>, orig: u32, what: &str) {
    let runs: Vec<usize> = (0..mesh.run_original_id.len())
        .filter(|&r| {
            mesh.run_original_id[r] == orig
                && mesh.run_index[r].to_u64() != mesh.run_index[r + 1].to_u64()
        })
        .collect();
    assert!(!runs.is_empty(), "{what}: no run of original {orig}");
    for r in runs {
        assert!(
            mesh.run_flags.get(r).is_some_and(|f| f & 1 != 0),
            "{what}: run {r} of the subtracted original {orig} is not a backside run"
        );
    }
}

fn soup_operand<K: GlKind>(tol: f64) {
    let mut reg = Registry::default();
    // Two boxes touching along one edge, as an unindexed soup: four faces
    // meet on that edge and no edge pairs by index, so only the soup-keeping
    // import accepts the mesh.
    let a = Manifold::cube(Vec3::splat(5.0), false);
    let b = Manifold::cube(Vec3::splat(5.0), false).translate(Vec3::new(5.0, 5.0, 0.0));
    let mesh = paint::<K>(concat::<K>(&[a, b], true), region_palette(5.0));
    assert!(
        K::import(&mesh).status() != Error::NoError,
        "the operand should be non-manifold"
    );
    let (soup, _) = reg.import_mesh_with::<K>(mesh, K::import_robust);
    assert!(
        soup.as_impl().is_soup,
        "the operand should import as a soup"
    );
    let cyl = Manifold::cylinder(9.0, 2.0, 2.0, 24).translate(Vec3::new(5.0, 5.0, -2.0));
    let (cyl, cyl_id) = reg.import::<K>(&cyl, cylinder_palette);
    for op in [OpType::Add, OpType::Subtract, OpType::Intersect] {
        let what = format!("soup {op:?} cyl");
        let mesh = reg.check::<K>(&soup.boolean_with_engine(&cyl, op, AUTO), tol, &what);
        if op == OpType::Subtract {
            assert_backside(&mesh, cyl_id, &what);
        }
        reg.check::<K>(
            &cyl.boolean_with_engine(&soup, op, AUTO),
            tol,
            &format!("cyl {op:?} soup"),
        );
    }
}

#[test]
fn face_id_provenance_auto_soup_operand_64() {
    soup_operand::<Gl64>(TOL64);
}

#[test]
fn face_id_provenance_auto_soup_operand_32() {
    soup_operand::<Gl32>(TOL32);
}

fn self_intersecting_operand<K: GlKind>(tol: f64) {
    let mut reg = Registry::default();
    // Two overlapping slabs in one mesh: closed shells that cross, which
    // `Auto` detects and hands to the robust engine.
    let a = Manifold::cube(Vec3::new(12.0, 12.0, 2.0), false);
    let b = Manifold::cube(Vec3::new(14.0, 8.0, 3.0), false)
        .rotate(0.0, 0.0, 10.0)
        .translate(Vec3::new(7.0, 2.0, -1.0));
    let mesh = paint::<K>(concat::<K>(&[a, b], false), region_palette(9.0));
    let (body, body_id) = reg.import_mesh::<K>(mesh);
    let cyl = Manifold::cylinder(9.0, 2.0, 2.0, 24).translate(Vec3::new(9.0, 6.0, -3.0));
    let (cyl, cyl_id) = reg.import::<K>(&cyl, cylinder_palette);
    for op in [OpType::Add, OpType::Subtract, OpType::Intersect] {
        let what = format!("self-intersecting {op:?} cyl");
        let mesh = reg.check::<K>(&body.boolean_with_engine(&cyl, op, AUTO), tol, &what);
        if op == OpType::Subtract {
            assert_backside(&mesh, cyl_id, &what);
        }
    }
    // A post standing on the body's top where the slabs do not overlap: the
    // floor around it keeps the body's colour.
    let post = Manifold::cylinder(4.0, 1.5, 1.5, 24).translate(Vec3::new(3.0, 3.0, 2.0));
    let (post, _) = reg.import::<K>(&post, solid(1));
    let out = body.boolean_with_engine(&post, OpType::Add, AUTO);
    let mesh = reg.check::<K>(&out, tol, "post on self-intersecting floor");
    assert_floor_face(
        &mesh,
        2.0,
        (Vec2::new(0.0, 0.0), Vec2::new(6.0, 6.0)),
        tol,
        body_id,
        5,
    );
}

#[test]
fn face_id_provenance_auto_self_intersecting_operand_64() {
    self_intersecting_operand::<Gl64>(TOL64);
}

#[test]
fn face_id_provenance_auto_self_intersecting_operand_32() {
    self_intersecting_operand::<Gl32>(TOL32);
}

/// The robust engine's subtracted operand comes back as backside runs, as
/// the fast engine's does.
#[test]
fn face_id_provenance_robust_difference_is_backside() {
    let mut reg = Registry::default();
    let (b, _) = reg.import::<Gl64>(&Manifold::cube(Vec3::splat(10.0), false), box_palette);
    let cyl = Manifold::cylinder(14.0, 3.0, 3.0, 32).translate(Vec3::new(5.0, 5.0, -2.0));
    let (c, cyl_id) = reg.import::<Gl64>(&cyl, cylinder_palette);
    for engine in [EXACT, ROBUST] {
        let mesh = b
            .boolean_with_engine(&c, OpType::Subtract, engine)
            .get_mesh_gl64(-1);
        assert_backside(&mesh, cyl_id, &format!("{engine:?} box - cyl"));
    }
}

/// One triangle of an exported mesh with its colour: (run original ID, face
/// ID, backside).
struct Painted {
    v: [Vec3; 3],
    n: Vec3,
    colour: (u32, u64, bool),
    height: f64,
}

fn painted<P: MeshPrecision, I: MeshIndex>(mesh: &MeshGLP<P, I>) -> Vec<Painted> {
    let pos = positions(mesh);
    let mut out = Vec::new();
    for run in 0..mesh.run_original_id.len() {
        let start = mesh.run_index[run].to_u64() as usize / 3;
        let end = mesh.run_index[run + 1].to_u64() as usize / 3;
        let back = mesh.run_flags.get(run).is_some_and(|f| f & 1 != 0);
        for tri in start..end {
            let v = tri_verts(mesh, &pos, tri);
            let n = tri_normal(&v);
            let longest = length(v[1] - v[0])
                .max(length(v[2] - v[1]))
                .max(length(v[0] - v[2]));
            let height = if longest > 0.0 {
                length(n) / longest
            } else {
                0.0
            };
            let colour = (mesh.run_original_id[run], mesh.face_id[tri].to_u64(), back);
            out.push(Painted {
                v,
                n,
                colour,
                height,
            });
        }
    }
    out
}

/// Does triangle `t` (same facing as `n`) cover point `p`, within `tol`?
fn covers(t: &Painted, p: Vec3, n: Vec3, tol: f64) -> bool {
    if length(t.n) == 0.0 || dot(normalize(t.n), n) < 0.999 {
        return false;
    }
    let tn = normalize(t.n);
    if dot(tn, p - t.v[0]).abs() > tol {
        return false;
    }
    (0..3).all(|i| {
        let (a, b) = (t.v[i], t.v[(i + 1) % 3]);
        let e = b - a;
        dot(cross(e, p - a), tn) >= -tol * length(e).max(1.0)
    })
}

/// The colours (run original ID, face ID) of `mesh`'s triangles covering `p`
/// with facing `n`.
fn colours_at(mesh: &[Painted], p: Vec3, n: Vec3, tol: f64) -> Vec<(u32, u64, bool)> {
    mesh.iter()
        .filter(|u| covers(u, p, n, tol))
        .map(|u| u.colour)
        .collect()
}

/// Every interior sample point of every well-shaped triangle of `a` must be
/// covered in `b` by a triangle of the same colour — except, for a union,
/// where both operands' faces cover the point with the same facing: there
/// the fast engine's choice of surviving face follows its symbolic
/// perturbation, which differs with the face's direction and position (even
/// within one face), while the robust engine always keeps the first
/// operand's. Either operand's colour is then correct; `ins` are the two
/// operands' own triangles, to recognise such points.
fn same_colours_one_way(
    a: &[Painted],
    b: &[Painted],
    ins: Option<[&[Painted]; 2]>,
    tol: f64,
    what: &str,
) {
    let mut failures = Vec::new();
    for (i, t) in a.iter().enumerate() {
        if t.height <= 1000.0 * tol {
            continue; // a sliver's facing is noise
        }
        let n = normalize(t.n);
        let c = (t.v[0] + t.v[1] + t.v[2]) / 3.0;
        let samples = [
            c,
            (t.v[0] + c * 2.0) / 3.0,
            (t.v[1] + c * 2.0) / 3.0,
            (t.v[2] + c * 2.0) / 3.0,
        ];
        for p in samples {
            let found = colours_at(b, p, n, tol);
            if found.contains(&t.colour) {
                continue;
            }
            if let Some([in_a, in_b]) = ins {
                let (ca, cb) = (colours_at(in_a, p, n, tol), colours_at(in_b, p, n, tol));
                let either =
                    |c: &(u32, u64, bool)| ca.iter().chain(&cb).any(|x| (x.0, x.1) == (c.0, c.1));
                if !ca.is_empty() && !cb.is_empty() && either(&t.colour) && found.iter().any(either)
                {
                    continue;
                }
            }
            failures.push(format!(
                "tri {i} {:?} at {p:?}: colour {:?}, other has {found:?}",
                t.v, t.colour
            ));
        }
    }
    assert!(
        failures.is_empty(),
        "{what}: {} samples disagree:\n{}",
        failures.len(),
        failures
            .iter()
            .take(20)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
}

/// The fast and robust results of one boolean paint the same surface with
/// the same colours, up to triangulation (and, for a union, up to which
/// operand's face survives where both coincide — see
/// [`same_colours_one_way`]).
fn engines_agree<K: GlKind>(a: &Manifold, b: &Manifold, op: OpType, tol: f64, what: &str) {
    let fast = painted(&K::export(&a.boolean_with_engine(b, op, EXACT)));
    let robust = painted(&K::export(&a.boolean_with_engine(b, op, ROBUST)));
    assert!(
        !fast.is_empty() && !robust.is_empty(),
        "{what}: empty result"
    );
    let (in_a, in_b) = (painted(&K::export(a)), painted(&K::export(b)));
    let ins = (op == OpType::Add).then_some([&in_a[..], &in_b[..]]);
    same_colours_one_way(&robust, &fast, ins, tol, &format!("{what}: robust vs fast"));
    same_colours_one_way(&fast, &robust, ins, tol, &format!("{what}: fast vs robust"));
}

fn cross_engine<K: GlKind>(tol: f64) {
    let mut reg = Registry::default();
    let (b, _) = reg.import::<K>(&Manifold::cube(Vec3::splat(10.0), false), box_palette);
    let cyl = Manifold::cylinder(14.0, 3.0, 3.0, 32).translate(Vec3::new(5.0, 5.0, -2.0));
    let (c, _) = reg.import::<K>(&cyl, cylinder_palette);
    let side = c.translate(Vec3::new(4.0, 0.0, 0.0)).rotate(0.0, 0.0, 15.0);
    for op in [OpType::Add, OpType::Subtract, OpType::Intersect] {
        engines_agree::<K>(&b, &c, op, tol, &format!("box {op:?} cyl"));
        engines_agree::<K>(&c, &b, op, tol, &format!("cyl {op:?} box"));
        engines_agree::<K>(&b, &side, op, tol, &format!("box {op:?} offset cyl"));
    }
    // Coplanar contacts: a box on a box, a flush notch, a post on a floor.
    let (base, _) = reg.import::<K>(
        &Manifold::cube(Vec3::new(10.0, 10.0, 5.0), false),
        box_palette,
    );
    let top = Manifold::cube(Vec3::splat(4.0), false).translate(Vec3::new(3.0, 3.0, 5.0));
    let (top, _) = reg.import::<K>(&top, solid(7));
    engines_agree::<K>(&base, &top, OpType::Add, tol, "box on box");
    let over = top.translate(Vec3::new(5.0, -3.0, 0.0));
    engines_agree::<K>(&base, &over, OpType::Add, tol, "box over edge");
    let cutter =
        Manifold::cube(Vec3::new(4.0, 10.0, 3.0), false).translate(Vec3::new(3.0, 0.0, 2.0));
    let (cutter, _) = reg.import::<K>(&cutter, box_palette);
    for op in [OpType::Subtract, OpType::Intersect, OpType::Add] {
        engines_agree::<K>(&base, &cutter, op, tol, &format!("flush notch {op:?}"));
    }
    // Same-facing coplanar overlap: both tops at z = 5.
    let (twin, _) = reg.import::<K>(
        &Manifold::cube(Vec3::new(6.0, 6.0, 5.0), false).translate(Vec3::new(5.0, 5.0, 0.0)),
        solid(9),
    );
    for op in [OpType::Add, OpType::Subtract, OpType::Intersect] {
        engines_agree::<K>(&base, &twin, op, tol, &format!("shared top {op:?}"));
        engines_agree::<K>(
            &twin,
            &base,
            op,
            tol,
            &format!("shared top reversed {op:?}"),
        );
    }
}

#[test]
fn face_id_provenance_fast_and_robust_agree_64() {
    cross_engine::<Gl64>(TOL64);
}

#[test]
fn face_id_provenance_fast_and_robust_agree_32() {
    cross_engine::<Gl32>(TOL32);
}

/// Rebuilding a soup keeps each triangle's source face: a self-intersecting
/// painted body rebuilt under the positive rule, and a body with one
/// inside-out shell rebuilt under the nonzero rule, whose rewound triangles
/// come back as a backside run of the same original.
fn rebuild_keeps_faces<K: GlKind>(tol: f64) {
    use crate::types::WindingRule;
    let mut reg = Registry::default();
    let a = Manifold::cube(Vec3::new(12.0, 12.0, 2.0), false);
    let b = Manifold::cube(Vec3::new(14.0, 8.0, 3.0), false)
        .rotate(0.0, 0.0, 10.0)
        .translate(Vec3::new(7.0, 2.0, -1.0));
    let mesh = paint::<K>(concat::<K>(&[a, b], false), region_palette(9.0));
    let (body, _) = reg.import_mesh::<K>(mesh);
    reg.check::<K>(
        &body.rebuild_solid(WindingRule::Positive),
        tol,
        "rebuilt crossing slabs",
    );

    // A box beside an inside-out box (its triangles wound inward).
    let c = Manifold::cube(Vec3::splat(4.0), false);
    let d = Manifold::cube(Vec3::splat(4.0), false).translate(Vec3::new(8.0, 0.0, 0.0));
    let mut mesh = concat::<K>(&[c, d], false);
    let half = mesh.tri_verts.len() / 2;
    for t in mesh.tri_verts[half..].chunks_mut(3) {
        t.swap(1, 2);
    }
    let mesh = paint::<K>(mesh, region_palette(6.0));
    let (inverted, id) = reg.import_mesh_with::<K>(mesh, K::import_robust);
    let out = reg.check::<K>(
        &inverted.rebuild_solid(WindingRule::Nonzero),
        tol,
        "rebuilt inside-out shell",
    );
    let flags: Vec<bool> = (0..out.run_original_id.len())
        .filter(|&r| out.run_original_id[r] == id)
        .map(|r| out.run_flags.get(r).is_some_and(|f| f & 1 != 0))
        .collect();
    assert!(
        flags.contains(&true) && flags.contains(&false),
        "the rewound shell should be a backside run beside a front run: {flags:?}"
    );
}

#[test]
fn face_id_provenance_robust_rebuild_keeps_faces_64() {
    rebuild_keeps_faces::<Gl64>(TOL64);
}

#[test]
fn face_id_provenance_robust_rebuild_keeps_faces_32() {
    rebuild_keeps_faces::<Gl32>(TOL32);
}

/// Where same-facing coplanar faces of both operands coincide, the robust
/// engine keeps the first operand's face, for every operation that keeps one.
#[test]
fn face_id_provenance_robust_coplanar_overlap_keeps_first_operand() {
    let mut reg = Registry::default();
    let (base, base_id) =
        reg.import::<Gl64>(&Manifold::cube(Vec3::new(10.0, 10.0, 5.0), false), solid(0));
    let twin = Manifold::cube(Vec3::splat(5.0), false).translate(Vec3::new(5.0, 5.0, 0.0));
    let (twin, twin_id) = reg.import::<Gl64>(&twin, solid(9));
    let probes = [
        (Vec3::new(7.5, 7.5, 5.0), Vec3::new(0.0, 0.0, 1.0)),
        (Vec3::new(7.5, 7.5, 0.0), Vec3::new(0.0, 0.0, -1.0)),
        (Vec3::new(10.0, 7.5, 2.5), Vec3::new(1.0, 0.0, 0.0)),
        (Vec3::new(7.5, 10.0, 2.5), Vec3::new(0.0, 1.0, 0.0)),
    ];
    for op in [OpType::Add, OpType::Intersect] {
        for (a, b, want) in [(&base, &twin, (base_id, 0)), (&twin, &base, (twin_id, 9))] {
            let out = painted(&a.boolean_with_engine(b, op, ROBUST).get_mesh_gl64(-1));
            for (p, n) in probes {
                let got: Vec<(u32, u64)> = colours_at(&out, p, n, TOL64)
                    .iter()
                    .map(|c| (c.0, c.1))
                    .collect();
                assert!(
                    !got.is_empty() && got.iter().all(|&c| c == want),
                    "{op:?} at {p:?}: {got:?}, expected {want:?}"
                );
            }
        }
    }
}
