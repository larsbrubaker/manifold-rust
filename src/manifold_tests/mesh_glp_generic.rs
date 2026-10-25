// Tests for the generic MeshGLP entry points (`Manifold::from_mesh_glp`,
// `from_mesh_glp_robust`, `get_mesh_glp`) — GitHub issue #11: an f64-vertex /
// u32-index `MeshGLP<f64, u32>` must import and export without the caller
// converting buffers to MeshGL or MeshGL64.
//
// The C++ templates these bodies as `Impl(const MeshGLP<Precision, I>&)` /
// `GetMeshGLImpl<Precision, I>` and only names the MeshGL and MeshGL64
// instantiations. Here, precision behaviour must follow `P` alone and the
// index width must not change any value, so every mixed instantiation is
// checked bit-for-bit against the named one sharing its precision:
// f64/u32 against `from_mesh_gl64` / `get_mesh_gl64`, f32/u32 against
// `from_mesh_gl` / `get_mesh_gl`.
//
// The samples include meshes with extra properties (normals, so property
// dedup and the merge vectors are populated), halfedge tangents, and legacy
// per-mesh-frame normals on transformed runs, and every export is checked
// with both the auto (-1) and an explicit (0) normal index;
// `samples_populate_every_export_field` pins that those fields really are
// non-empty so the comparisons cannot silently compare empty vectors.

use super::*;
use crate::types::{Error, MeshGL64, MeshGLP, MeshIndex, MeshPrecision};

/// Normal indices every export comparison runs with: -1 (auto-substitutes
/// slot 0 when every run has world-frame normals) and an explicit slot 0
/// (normalizes the slot, transforming legacy per-mesh-frame runs).
const NORMAL_IDXS: [i32; 2] = [-1, 0];

/// Shapes that exercise curved faces, many runs, transforms, a boolean,
/// extra properties with normals, tangents, and legacy normals.
fn sample_manifolds() -> Vec<(&'static str, Manifold)> {
    let sphere = Manifold::sphere(1.0, 48);
    let cylinder = Manifold::cylinder(2.0, 0.7, 0.4, 37).rotate(10.0, 20.0, 30.0);
    let cube = Manifold::cube(Vec3::new(1.5, 1.5, 1.5), true);
    let boolean = cube
        .difference(&sphere.translate(Vec3::new(0.3, -0.2, 0.1)))
        .union(&cylinder.translate(Vec3::new(0.5, 0.5, 0.5)));
    assert_eq!(boolean.status(), Error::NoError);
    // Sharp edges split the normals, so vertices carry several property
    // vertices and the export fills merge_from_vert / merge_to_vert.
    let normals = boolean.calculate_normals(0, 40.0);
    assert_eq!(normals.status(), Error::NoError);
    let tangents = Manifold::cube(Vec3::new(1.0, 1.0, 1.0), true)
        .difference(&Manifold::sphere(0.65, 24))
        .smooth_out(60.0, 0.0);
    assert_eq!(tangents.status(), Error::NoError);
    vec![
        ("sphere", sphere),
        ("cylinder", cylinder),
        ("boolean", boolean),
        ("normals", normals),
        ("tangents", tangents),
        ("legacy_normals", legacy_normals()),
    ]
}

/// Normals in property slot 0..2 on runs *without* the hasNormals bit, then
/// rotated and combined so the result is non-original: exporting with
/// normal index 0 re-normalizes them through each run's transform
/// (`normal_transform_for_run`).
fn legacy_normals() -> Manifold {
    let mut mesh = Manifold::sphere(1.0, 32)
        .calculate_normals(0, 60.0)
        .get_mesh_gl64(-1);
    mesh.run_flags.clear();
    let imported = Manifold::from_mesh_gl64(&mesh);
    assert_eq!(imported.status(), Error::NoError);
    let m = imported.rotate(15.0, 25.0, 35.0).union(
        &Manifold::cube(Vec3::new(1.0, 1.0, 1.0), false).translate(Vec3::new(0.5, 0.0, 0.0)),
    );
    assert_eq!(m.status(), Error::NoError);
    m
}

/// Distinct vertex positions as sorted bit triples, independent of vertex
/// order and of how many property vertices share a position.
fn sorted_positions(props: &[f64], num_prop: usize) -> Vec<[u64; 3]> {
    let mut out: Vec<[u64; 3]> = props
        .chunks(num_prop)
        .map(|p| [p[0].to_bits(), p[1].to_bits(), p[2].to_bits()])
        .collect();
    out.sort();
    out.dedup();
    out
}

fn widen(v: &[u32]) -> Vec<u64> {
    v.iter().map(|&x| x as u64).collect()
}

fn prec_bits<P: MeshPrecision>(v: &[P]) -> Vec<u64> {
    // f32 -> f64 widening is exact, so this compares f32 data bit for bit too.
    v.iter().map(|x| x.to_f64().to_bits()).collect()
}

fn idx<I: MeshIndex>(v: &[I]) -> Vec<u64> {
    v.iter().map(|x| x.to_u64()).collect()
}

/// Every field of two exports with the same precision must match: float
/// buffers bit for bit, index buffers after widening to u64.
fn assert_mesh_eq<P: MeshPrecision, IA: MeshIndex, IB: MeshIndex>(
    a: &MeshGLP<P, IA>,
    b: &MeshGLP<P, IB>,
    what: &str,
) {
    assert_eq!(a.num_prop.to_u64(), b.num_prop.to_u64(), "{what}: num_prop");
    assert_eq!(
        prec_bits(&a.vert_properties),
        prec_bits(&b.vert_properties),
        "{what}: vert_properties"
    );
    assert_eq!(idx(&a.tri_verts), idx(&b.tri_verts), "{what}: tri_verts");
    assert_eq!(
        idx(&a.merge_from_vert),
        idx(&b.merge_from_vert),
        "{what}: merge_from_vert"
    );
    assert_eq!(
        idx(&a.merge_to_vert),
        idx(&b.merge_to_vert),
        "{what}: merge_to_vert"
    );
    assert_eq!(idx(&a.run_index), idx(&b.run_index), "{what}: run_index");
    assert_eq!(
        a.run_original_id, b.run_original_id,
        "{what}: run_original_id"
    );
    assert_eq!(
        prec_bits(&a.run_transform),
        prec_bits(&b.run_transform),
        "{what}: run_transform"
    );
    assert_eq!(idx(&a.face_id), idx(&b.face_id), "{what}: face_id");
    assert_eq!(
        prec_bits(&a.halfedge_tangent),
        prec_bits(&b.halfedge_tangent),
        "{what}: halfedge_tangent"
    );
    assert_eq!(a.run_flags, b.run_flags, "{what}: run_flags");
    assert_eq!(
        a.tolerance.to_f64().to_bits(),
        b.tolerance.to_f64().to_bits(),
        "{what}: tolerance"
    );
}

/// Two manifolds are identical when their status, counts, f64 volume and
/// surface area, and full f64 export agree bit for bit.
fn assert_same_manifold(a: &Manifold, b: &Manifold, what: &str) {
    assert_eq!(a.status(), b.status(), "{what}: status");
    assert_eq!(a.num_tri(), b.num_tri(), "{what}: num_tri");
    assert_eq!(a.num_vert(), b.num_vert(), "{what}: num_vert");
    assert_eq!(a.volume().to_bits(), b.volume().to_bits(), "{what}: volume");
    assert_eq!(
        a.surface_area().to_bits(),
        b.surface_area().to_bits(),
        "{what}: surface_area"
    );
    assert_eq!(a.imp.is_soup, b.imp.is_soup, "{what}: is_soup");
    for n in NORMAL_IDXS {
        let what = format!("{what} (exported, normal_idx {n})");
        assert_mesh_eq(&a.get_mesh_gl64(n), &b.get_mesh_gl64(n), &what);
    }
}

#[test]
fn get_mesh_glp_f64_u32_matches_get_mesh_gl64() {
    for (name, m) in sample_manifolds() {
        for n in NORMAL_IDXS {
            let what = format!("{name} (normal_idx {n})");
            let glp = m.get_mesh_glp::<f64, u32>(n);
            let gl64 = m.get_mesh_gl64(n);
            assert!(glp.num_tri() > 0, "{what}: empty export");
            assert_mesh_eq(&glp, &gl64, &what);
            // The named instantiation is the generic one with P = f64, I = u64.
            assert_mesh_eq(&m.get_mesh_glp::<f64, u64>(n), &gl64, &what);
        }
    }
}

#[test]
fn get_mesh_glp_f32_u32_matches_get_mesh_gl() {
    for (name, m) in sample_manifolds() {
        for n in NORMAL_IDXS {
            let what = format!("{name} (normal_idx {n})");
            // Includes the f32 tolerance floor, exactly like get_mesh_gl.
            assert_mesh_eq(&m.get_mesh_glp::<f32, u32>(n), &m.get_mesh_gl(n), &what);
        }
    }
}

/// The samples must actually populate the fields the comparisons above
/// cover; otherwise those comparisons would pass on empty vectors.
#[test]
fn samples_populate_every_export_field() {
    let samples: std::collections::BTreeMap<_, _> = sample_manifolds().into_iter().collect();
    for n in NORMAL_IDXS {
        let normals = samples["normals"].get_mesh_glp::<f64, u32>(n);
        assert_eq!(
            normals.num_prop, 6,
            "normals: extra properties (normal_idx {n})"
        );
        assert!(
            !normals.merge_from_vert.is_empty(),
            "normals: merge_from_vert (normal_idx {n})"
        );
        assert!(
            !normals.merge_to_vert.is_empty(),
            "normals: merge_to_vert (normal_idx {n})"
        );
        assert!(
            normals.run_flags.iter().all(|f| f & 2 != 0),
            "normals: hasNormals runs"
        );
        let normals32 = samples["normals"].get_mesh_glp::<f32, u32>(n);
        assert!(
            !normals32.merge_from_vert.is_empty(),
            "normals f32: merge_from_vert"
        );

        let tangents = samples["tangents"].get_mesh_glp::<f64, u32>(n);
        assert!(
            !tangents.halfedge_tangent.is_empty(),
            "tangents: halfedge_tangent (normal_idx {n})"
        );
        let tangents32 = samples["tangents"].get_mesh_glp::<f32, u32>(n);
        assert!(
            !tangents32.halfedge_tangent.is_empty(),
            "tangents f32: halfedge_tangent"
        );

        // Legacy normals: transformed runs (run_transform present) that lack
        // the hasNormals bit, so slot 0 goes through normal_transform_for_run.
        let legacy = samples["legacy_normals"].get_mesh_glp::<f64, u32>(n);
        assert_eq!(legacy.num_prop, 6, "legacy_normals: extra properties");
        assert!(
            !legacy.run_transform.is_empty(),
            "legacy_normals: run_transform"
        );
        assert!(
            legacy.run_flags.iter().any(|f| f & 2 == 0),
            "legacy_normals: run without hasNormals"
        );
    }
    // An explicit normal index re-normalizes the legacy normals, so the
    // exported slot differs from the raw one (the -1 export leaves it as is).
    let legacy = &samples["legacy_normals"];
    assert_ne!(
        prec_bits(&legacy.get_mesh_glp::<f64, u32>(0).vert_properties),
        prec_bits(&legacy.get_mesh_glp::<f64, u32>(-1).vert_properties),
        "legacy_normals: normal_idx 0 must transform the normals"
    );
}

#[test]
fn from_mesh_glp_f64_u32_matches_from_mesh_gl64() {
    for (name, m) in sample_manifolds() {
        let glp = m.get_mesh_glp::<f64, u32>(-1);
        let gl64 = m.get_mesh_gl64(-1);
        let a = Manifold::from_mesh_glp(&glp);
        let b = Manifold::from_mesh_gl64(&gl64);
        assert_eq!(a.status(), Error::NoError, "{name}");
        assert_same_manifold(&a, &b, name);
        // Lossless: every source vertex position comes back bit for bit
        // (the import re-sorts vertices and may regroup the property vertices
        // sharing a position, so compare distinct positions as sorted sets).
        let back = a.get_mesh_glp::<f64, u32>(-1);
        assert_eq!(
            sorted_positions(&back.vert_properties, back.num_prop as usize),
            sorted_positions(&gl64.vert_properties, gl64.num_prop as usize),
            "{name}: round-trip positions"
        );
        assert_mesh_eq(&back, &a.get_mesh_gl64(-1), name);
    }
}

#[test]
fn from_mesh_glp_f32_u32_matches_from_mesh_gl() {
    for (name, m) in sample_manifolds() {
        let gl = m.get_mesh_gl(-1);
        let a = Manifold::from_mesh_glp::<f32, u32>(&gl);
        let b = Manifold::from_mesh_gl(&gl);
        assert_eq!(a.status(), Error::NoError, "{name}");
        assert_same_manifold(&a, &b, name);
    }
}

/// Triangle soup of a [0,2]³ cube: three unshared vertices per triangle, so
/// strict halfedge pairing fails and only the robust import keeps it.
fn soup_cube_f64_u32() -> MeshGLP<f64, u32> {
    let cube = Manifold::cube(Vec3::new(2.0, 2.0, 2.0), false).get_mesh_gl64(-1);
    let np = cube.num_prop as usize;
    let mut mesh = MeshGLP::<f64, u32> {
        num_prop: 3,
        run_original_id: vec![Manifold::reserve_ids(1)],
        ..Default::default()
    };
    for &v in &cube.tri_verts {
        let o = v as usize * np;
        mesh.vert_properties
            .extend_from_slice(&cube.vert_properties[o..o + 3]);
    }
    mesh.tri_verts = (0..cube.tri_verts.len() as u32).collect();
    mesh
}

fn to_gl64(mesh: &MeshGLP<f64, u32>) -> MeshGL64 {
    MeshGL64 {
        num_prop: mesh.num_prop as u64,
        vert_properties: mesh.vert_properties.clone(),
        tri_verts: widen(&mesh.tri_verts),
        merge_from_vert: widen(&mesh.merge_from_vert),
        merge_to_vert: widen(&mesh.merge_to_vert),
        run_index: widen(&mesh.run_index),
        run_original_id: mesh.run_original_id.clone(),
        run_transform: mesh.run_transform.clone(),
        face_id: widen(&mesh.face_id),
        halfedge_tangent: mesh.halfedge_tangent.clone(),
        run_flags: mesh.run_flags.clone(),
        tolerance: mesh.tolerance,
    }
}

#[test]
fn from_mesh_glp_robust_f64_u32_matches_from_mesh_gl64_robust() {
    // Manifold input: the robust import behaves like the strict one.
    for (name, m) in sample_manifolds() {
        let glp = m.get_mesh_glp::<f64, u32>(-1);
        let a = Manifold::from_mesh_glp_robust(&glp);
        let b = Manifold::from_mesh_gl64_robust(&m.get_mesh_gl64(-1));
        assert_eq!(a.status(), Error::NoError, "{name}");
        assert_same_manifold(&a, &b, name);
    }
    // Non-manifold input: kept as a soup on both paths, identically.
    let soup = soup_cube_f64_u32();
    let a = Manifold::from_mesh_glp_robust(&soup);
    let b = Manifold::from_mesh_gl64_robust(&to_gl64(&soup));
    assert_eq!(a.status(), Error::NoError);
    assert!(a.imp.is_soup, "disconnected cube must import as a soup");
    assert_eq!(a.volume(), 8.0);
    assert_same_manifold(&a, &b, "soup cube");
    // The strict generic import rejects it, as from_mesh_gl64 does.
    let strict = Manifold::from_mesh_glp(&soup);
    assert_eq!(
        strict.status(),
        Manifold::from_mesh_gl64(&to_gl64(&soup)).status()
    );
    assert_ne!(strict.status(), Error::NoError);
}

/// Issue #11's use case: buffers the caller built as f64 vertices with u32
/// indices import directly, with the coordinates kept at full f64 precision.
#[test]
fn user_built_f64_u32_mesh_imports_directly() {
    // Unit cube scaled per axis by values f32 cannot represent.
    let s = [1.0_f64 / 3.0, 2.0_f64 / 3.0, 5.0_f64 / 7.0];
    let mut cube = MeshGLP::<f64, u32> {
        num_prop: 3,
        ..Default::default()
    };
    for i in 0..8u32 {
        let corner = [(i & 1) as f64, ((i >> 1) & 1) as f64, ((i >> 2) & 1) as f64];
        cube.vert_properties
            .extend((0..3).map(|j| corner[j] * s[j]));
    }
    // Outward-facing (CCW from outside) triangles; vertex i = x + 2y + 4z.
    cube.tri_verts = vec![
        0, 2, 1, 1, 2, 3, // z = 0
        4, 5, 6, 5, 7, 6, // z = 1
        0, 1, 4, 1, 5, 4, // y = 0
        2, 6, 3, 3, 6, 7, // y = 1
        0, 4, 2, 2, 4, 6, // x = 0
        1, 3, 5, 3, 7, 5, // x = 1
    ];
    let m = Manifold::from_mesh_glp(&cube);
    assert_eq!(m.status(), Error::NoError);
    assert_eq!(m.num_tri(), 12);
    assert_eq!(m.num_vert(), 8);
    let expected = s[0] * s[1] * s[2];
    assert!(
        (m.volume() - expected).abs() <= 1e-15,
        "volume {} != {expected}",
        m.volume()
    );
    // Round trip: the far corner survives export at full f64 precision.
    let out = m.get_mesh_glp::<f64, u32>(-1);
    let np = out.num_prop as usize;
    assert!(
        (0..out.num_vert()).any(|i| (0..3).all(|j| out.vert_properties[i * np + j] == s[j])),
        "far corner lost precision on the f64/u32 round trip"
    );

    // A hand-written tetrahedron, again straight from f64/u32 buffers.
    let tet = MeshGLP::<f64, u32> {
        num_prop: 3,
        vert_properties: vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
        tri_verts: vec![0, 2, 1, 0, 1, 3, 0, 3, 2, 1, 2, 3],
        ..Default::default()
    };
    let t = Manifold::from_mesh_glp(&tet);
    assert_eq!(t.status(), Error::NoError);
    assert_eq!(t.num_tri(), 4);
    assert!(
        (t.volume() - 1.0 / 6.0).abs() <= 1e-15,
        "tet volume {}",
        t.volume()
    );
}
