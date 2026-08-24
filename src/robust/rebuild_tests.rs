// robust/rebuild_tests.rs — Unit tests for the single-operand "rebuild solid"
// entry (`robust::rebuild_with_rule` / `Manifold::rebuild_solid`).
//
// Where robust/repair_tests.rs covers the cheap winding-only repair, these
// cover the full arrangement pipeline run on one mesh: the fixtures are
// deliberately *not* fixable by rewinding shells (doubled faces, mutually
// overlapping bodies, redundant interior shells) and can only come out
// 2-manifold if the mesh is genuinely re-derived from the winding numbers.

use crate::linalg::Vec3;
use crate::manifold::Manifold;
use crate::types::{Error, MeshGL, WindingRule};

fn v(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z)
}

/// Axis-aligned cube [lo,hi]³ as 12 outward-wound triangles.
fn cube_tris(lo: f64, hi: f64) -> Vec<[Vec3; 3]> {
    let quads: [([f64; 3], [f64; 3], [f64; 3], [f64; 3]); 6] = [
        ([0., 0., 0.], [0., 1., 0.], [1., 1., 0.], [1., 0., 0.]), // -z
        ([0., 0., 1.], [1., 0., 1.], [1., 1., 1.], [0., 1., 1.]), // +z
        ([0., 0., 0.], [1., 0., 0.], [1., 0., 1.], [0., 0., 1.]), // -y
        ([0., 1., 0.], [0., 1., 1.], [1., 1., 1.], [1., 1., 0.]), // +y
        ([0., 0., 0.], [0., 0., 1.], [0., 1., 1.], [0., 1., 0.]), // -x
        ([1., 0., 0.], [1., 1., 0.], [1., 1., 1.], [1., 0., 1.]), // +x
    ];
    let s = hi - lo;
    let m = |q: [f64; 3]| v(lo + q[0] * s, lo + q[1] * s, lo + q[2] * s);
    let mut out = Vec::new();
    for (a, b, c, d) in quads {
        out.push([m(a), m(b), m(c)]);
        out.push([m(a), m(c), m(d)]);
    }
    out
}

/// Box [lo,hi] with independent per-axis extents, outward-wound.
fn box_tris(lo: Vec3, hi: Vec3) -> Vec<[Vec3; 3]> {
    let unit = cube_tris(0.0, 1.0);
    let m = |p: Vec3| {
        v(
            lo.x + p.x * (hi.x - lo.x),
            lo.y + p.y * (hi.y - lo.y),
            lo.z + p.z * (hi.z - lo.z),
        )
    };
    unit.iter().map(|t| [m(t[0]), m(t[1]), m(t[2])]).collect()
}

fn flipped(tris: &[[Vec3; 3]]) -> Vec<[Vec3; 3]> {
    tris.iter().map(|t| [t[0], t[2], t[1]]).collect()
}

/// Import a raw soup exactly as given — no welding away of duplicate faces,
/// no orientation repair — so the rebuild is what has to clean it up.
fn mesh_from_tris(tris: &[[Vec3; 3]]) -> Manifold {
    let mut mesh = MeshGL::default();
    mesh.num_prop = 3;
    for t in tris {
        for p in t {
            mesh.vert_properties
                .extend([p.x as f32, p.y as f32, p.z as f32]);
        }
    }
    mesh.tri_verts = (0..(tris.len() * 3) as u32).collect();
    mesh.merge();
    Manifold::from_mesh_gl_robust(&mesh)
}

/// Signed (divergence-theorem) volume — `volume()` takes the absolute value,
/// which would hide an inside-out result.
fn signed_volume(m: &Manifold) -> f64 {
    let tris = crate::robust::soup::impl_to_tris(m.as_impl());
    tris.iter()
        .map(|t| crate::linalg::dot(t[0], crate::linalg::cross(t[1], t[2])) / 6.0)
        .sum()
}

/// Every assertion a "properly paired 2-manifold" has to satisfy: real
/// halfedge pairing (not a soup fallback), a valid import status, and no
/// duplicated face.
fn assert_clean_manifold(m: &Manifold) {
    assert_eq!(m.status(), Error::NoError, "rebuilt mesh must import clean");
    assert!(
        !m.as_impl().is_soup,
        "rebuilt mesh must not fall back to soup"
    );
    assert!(m.as_impl().is_manifold(), "rebuilt mesh must be 2-manifold");
    let tris = crate::robust::soup::impl_to_tris(m.as_impl());
    let mut keys: Vec<[[u64; 3]; 3]> = tris
        .iter()
        .map(|t| {
            let mut k = [
                [t[0].x.to_bits(), t[0].y.to_bits(), t[0].z.to_bits()],
                [t[1].x.to_bits(), t[1].y.to_bits(), t[1].z.to_bits()],
                [t[2].x.to_bits(), t[2].y.to_bits(), t[2].z.to_bits()],
            ];
            k.sort();
            k
        })
        .collect();
    let before = keys.len();
    keys.sort();
    keys.dedup();
    assert_eq!(
        before,
        keys.len(),
        "rebuilt mesh must have no doubled faces"
    );
}

#[test]
fn duplicated_faces_collapse_to_one_shell() {
    // Every triangle present twice: >2 faces per edge everywhere, so nothing
    // pairs. The winding jumps 0 → 2 across the surface, which {w >= 1} still
    // reads as one solid boundary.
    let mut tris = cube_tris(0.0, 2.0);
    tris.extend(cube_tris(0.0, 2.0));
    let m = mesh_from_tris(&tris);
    // The import keeps both copies (it never dedups faces), so the fixture
    // really does carry the doubled sheet the rebuild has to collapse — its
    // divergence-theorem volume double-counts at 16.
    assert_eq!(m.num_tri(), 24, "fixture must keep both copies");
    assert!((signed_volume(&m) - 16.0).abs() < 1e-9);

    let out = m.rebuild_solid(WindingRule::Positive);
    assert_clean_manifold(&out);
    assert!(
        (signed_volume(&out) - 8.0).abs() < 1e-9,
        "expected 8, got {}",
        signed_volume(&out)
    );
}

#[test]
fn overlapping_bodies_in_one_soup_become_their_union() {
    // Two mutually penetrating cubes concatenated into a single mesh: the
    // interior walls have material on both sides and must dissolve.
    let mut tris = cube_tris(0.0, 2.0);
    tris.extend(box_tris(v(1.0, 1.0, 1.0), v(3.0, 3.0, 3.0)));
    let soup = mesh_from_tris(&tris);
    assert!(
        soup.has_self_intersections(),
        "fixture must genuinely self-intersect"
    );
    let out = soup.rebuild_solid(WindingRule::Positive);
    assert_clean_manifold(&out);

    let a = Manifold::cube(v(2.0, 2.0, 2.0), false);
    let b = Manifold::cube(v(2.0, 2.0, 2.0), false).translate(v(1.0, 1.0, 1.0));
    let reference = a.union(&b);
    // Signed, not `volume()`: an inside-out union would have the same absolute
    // volume and would sail through.
    let expected = signed_volume(&reference);
    assert!(expected > 0.0, "reference union must be outward-wound");
    assert!(
        (signed_volume(&out) - expected).abs() < 1e-9,
        "expected {expected}, got {}",
        signed_volume(&out)
    );
}

#[test]
fn genuine_cavity_survives_rebuild() {
    // Outer [0,6] outward + inner [2,4] inward: a solid with a void. The
    // inner shell is a real boundary (material outside it, none inside), so
    // the rebuild must keep it.
    let mut tris = cube_tris(0.0, 6.0);
    tris.extend(flipped(&cube_tris(2.0, 4.0)));
    let out = mesh_from_tris(&tris).rebuild_solid(WindingRule::Positive);
    assert_clean_manifold(&out);
    assert!(
        (signed_volume(&out) - (216.0 - 8.0)).abs() < 1e-9,
        "expected 208, got {}",
        signed_volume(&out)
    );
}

#[test]
fn redundant_interior_shell_is_removed() {
    // Outer [0,6] and inner [2,4] both wound outward: the inner shell is
    // material on both sides (winding 2 inside, 1 outside), so it is not a
    // boundary at all and must vanish.
    let mut tris = cube_tris(0.0, 6.0);
    tris.extend(cube_tris(2.0, 4.0));
    let m = mesh_from_tris(&tris);
    assert_eq!(m.num_tri(), 24, "fixture must carry both shells");
    let out = m.rebuild_solid(WindingRule::Positive);
    assert_clean_manifold(&out);
    assert!(
        (signed_volume(&out) - 216.0).abs() < 1e-9,
        "expected 216, got {}",
        signed_volume(&out)
    );
    assert_eq!(
        out.num_tri(),
        12,
        "only the outer cube's faces should remain"
    );
}

#[test]
fn inside_out_cube_depends_on_the_winding_rule() {
    // Winding -1 inside. {w != 0} calls that solid and rewinds it outward;
    // {w >= 1} calls it nothing at all.
    let m = mesh_from_tris(&flipped(&cube_tris(0.0, 2.0)));
    assert!(signed_volume(&m) < 0.0, "fixture must import inverted");

    let nonzero = m.rebuild_solid(WindingRule::Nonzero);
    assert_clean_manifold(&nonzero);
    assert!(
        (signed_volume(&nonzero) - 8.0).abs() < 1e-9,
        "expected +8, got {}",
        signed_volume(&nonzero)
    );

    let positive = m.rebuild_solid(WindingRule::Positive);
    assert!(
        positive.is_empty(),
        "no positive material, got {} tris",
        positive.num_tri()
    );
    assert_eq!(positive.status(), Error::NoError);
}

/// The operand occupies mesh slot 0, so `assemble::PropCtx` interpolates its
/// corner properties onto the rebuilt triangles exactly as a two-operand
/// boolean would. An affine-in-position property must therefore reproduce
/// itself at every output vertex, including the ones the arrangement invents.
#[test]
fn corner_properties_survive_the_rebuild() {
    let prop = |pos: Vec3| 0.25 + pos.x * 0.5 + pos.y * 0.125 - pos.z * 0.0625;
    // Two overlapping cubes in one body, so the rebuild really does
    // re-triangulate and must interpolate rather than copy.
    let a = Manifold::cube(v(2.0, 2.0, 2.0), false);
    let b = Manifold::cube(v(2.0, 2.0, 2.0), false).translate(v(1.0, 1.0, 1.0));
    let joined = a
        .union(&b)
        .set_properties(1, |new_prop, pos, _old| new_prop[0] = prop(pos));
    let out = joined.rebuild_solid(WindingRule::Positive);
    assert_eq!(out.status(), Error::NoError);
    assert_eq!(out.as_impl().num_prop, 1, "the property must survive");

    let np = 4usize; // xyz + 1
    let gl = out.get_mesh_gl(-1);
    assert_eq!(gl.num_prop as usize, np);
    for i in 0..gl.vert_properties.len() / np {
        let pos = v(
            gl.vert_properties[i * np] as f64,
            gl.vert_properties[i * np + 1] as f64,
            gl.vert_properties[i * np + 2] as f64,
        );
        let got = gl.vert_properties[i * np + 3] as f64;
        assert!(
            (got - prop(pos)).abs() < 1e-6,
            "property at {pos:?}: expected {}, got {got}",
            prop(pos)
        );
    }
}

/// The no-op case has to actually be a no-op. Nothing in a strictly manifold,
/// self-intersection-free mesh needs cutting, so every original triangle must
/// survive whole: same count, and a volume that matches bit-for-bit rather
/// than merely to a tolerance. A rebuild that quietly perturbed clean input
/// would be invisible to every other test here, which all work in tolerances.
#[test]
fn already_manifold_input_round_trips_unchanged() {
    for input in [
        Manifold::cube(v(2.0, 3.0, 5.0), false),
        Manifold::sphere(1.5, 32),
    ] {
        let out = input.rebuild_solid(WindingRule::Positive);
        assert_clean_manifold(&out);
        assert_eq!(
            out.num_tri(),
            input.num_tri(),
            "clean input must not be re-triangulated"
        );
        assert_eq!(
            out.volume().to_bits(),
            input.volume().to_bits(),
            "expected {}, got {}",
            input.volume(),
            out.volume()
        );
    }
}

#[test]
fn empty_input_rebuilds_to_empty() {
    let empty = Manifold::new();
    let out = empty.rebuild_solid(WindingRule::Positive);
    assert!(out.is_empty());
    assert_eq!(out.status(), Error::NoError);
}

#[test]
fn rebuild_respects_a_cancel_token() {
    let token = crate::cancel::CancelToken::new();
    token.cancel();
    let mut tris = cube_tris(0.0, 2.0);
    tris.extend(cube_tris(0.0, 2.0));
    let m = mesh_from_tris(&tris);
    let out = super::rebuild_with_rule(m.as_impl(), WindingRule::Positive, Some(&token), None);
    assert!(out.is_empty());
    assert_eq!(out.status, Error::Cancelled);
}
