// Face-ID provenance harness: imports operands as originals whose triangles
// carry palette-style face IDs (many triangles per ID, IDs spanning
// non-coplanar triangles), records every input triangle per original ID, and
// checks the provenance invariant on boolean results: each output triangle's
// (run original ID, face ID) names an input triangle whose plane contains the
// output triangle with the same orientation. The cases live in
// `face_id_provenance_cases.rs` and `face_id_provenance_nested.rs`; this file
// is shared by both and works for both the f64 `MeshGL64` and f32 `MeshGL`
// import/export paths.

use super::*;
use crate::linalg::{cross, dot, length, normalize};
use crate::types::{BooleanEngine, MeshGL64, MeshGLP, MeshIndex, MeshPrecision};
use std::collections::BTreeMap;

/// One MeshGL instantiation: how to export a Manifold to it and import it back.
pub(super) trait GlKind {
    type P: MeshPrecision;
    type I: MeshIndex;
    fn export(m: &Manifold) -> MeshGLP<Self::P, Self::I>;
    fn import(mesh: &MeshGLP<Self::P, Self::I>) -> Manifold;
}

pub(super) struct Gl64;
pub(super) struct Gl32;

impl GlKind for Gl64 {
    type P = f64;
    type I = u64;
    fn export(m: &Manifold) -> MeshGL64 {
        m.get_mesh_gl64(-1)
    }
    fn import(mesh: &MeshGL64) -> Manifold {
        Manifold::from_mesh_gl64(mesh)
    }
}

impl GlKind for Gl32 {
    type P = f32;
    type I = u32;
    fn export(m: &Manifold) -> MeshGL {
        m.get_mesh_gl(-1)
    }
    fn import(mesh: &MeshGL) -> Manifold {
        Manifold::from_mesh_gl(mesh)
    }
}

pub(super) const EXACT: BooleanEngine = BooleanEngine::Exact;

/// A distinct plane of one face ID of one input mesh: three points of a
/// representative triangle, in the input mesh's own coordinates.
struct SourcePlane {
    face: u64,
    v: [Vec3; 3],
}

/// Every registered input mesh, keyed by the original ID it was imported with.
#[derive(Default)]
pub(super) struct Registry {
    sources: BTreeMap<u32, Vec<SourcePlane>>,
}

fn tri_normal(v: &[Vec3; 3]) -> Vec3 {
    cross(v[1] - v[0], v[2] - v[0])
}

fn positions<P: MeshPrecision, I: MeshIndex>(mesh: &MeshGLP<P, I>) -> Vec<Vec3> {
    let np = mesh.num_prop.to_u64() as usize;
    (0..mesh.vert_properties.len() / np)
        .map(|i| {
            Vec3::new(
                mesh.vert_properties[np * i].to_f64(),
                mesh.vert_properties[np * i + 1].to_f64(),
                mesh.vert_properties[np * i + 2].to_f64(),
            )
        })
        .collect()
}

fn tri_verts<P: MeshPrecision, I: MeshIndex>(
    mesh: &MeshGLP<P, I>,
    pos: &[Vec3],
    tri: usize,
) -> [Vec3; 3] {
    let at = |j: usize| pos[mesh.tri_verts[3 * tri + j].to_u64() as usize];
    [at(0), at(1), at(2)]
}

impl Registry {
    /// Export `m`, colour each triangle with `face_of(unit normal, centroid)`,
    /// and import it back as a fresh original with those face IDs: the import
    /// path a consumer uses (MeshGL with `face_id` set and one reserved
    /// `run_original_id`).
    pub fn import<K: GlKind>(
        &mut self,
        m: &Manifold,
        face_of: impl Fn(Vec3, Vec3) -> u64,
    ) -> (Manifold, u32) {
        let mut mesh = K::export(m);
        let pos = positions(&mesh);
        let n = mesh.tri_verts.len() / 3;
        mesh.face_id = (0..n)
            .map(|t| {
                let v = tri_verts(&mesh, &pos, t);
                let c = (v[0] + v[1] + v[2]) / 3.0;
                K::I::from_usize(face_of(normalize(tri_normal(&v)), c) as usize)
            })
            .collect();
        self.import_mesh::<K>(mesh)
    }

    /// Import `mesh` (whose `face_id` is already set, e.g. the colours a
    /// previous boolean recovered) as a fresh original, dropping its runs.
    pub fn import_mesh<K: GlKind>(&mut self, mut mesh: MeshGLP<K::P, K::I>) -> (Manifold, u32) {
        let id = Manifold::reserve_ids(1);
        mesh.run_index.clear();
        mesh.run_transform.clear();
        mesh.run_flags.clear();
        mesh.halfedge_tangent.clear();
        mesh.run_original_id = vec![id];
        let pos = positions(&mesh);
        let mut planes: Vec<SourcePlane> = Vec::new();
        for t in 0..mesh.tri_verts.len() / 3 {
            let v = tri_verts(&mesh, &pos, t);
            let face = mesh.face_id[t].to_u64();
            let n = tri_normal(&v);
            if length(n) == 0.0 {
                continue;
            }
            let nu = normalize(n);
            let dup = planes.iter().any(|p| {
                let pn = normalize(tri_normal(&p.v));
                p.face == face && dot(pn, nu) > 1.0 - 1e-12 && dot(pn, v[0] - p.v[0]).abs() < 1e-9
            });
            if !dup {
                planes.push(SourcePlane { face, v });
            }
        }
        self.sources.insert(id, planes);
        let out = K::import(&mesh);
        assert_eq!(out.status(), Error::NoError, "import failed");
        (out, id)
    }

    /// Check the provenance invariant on `m`'s exported mesh; returns it so
    /// a caller can re-import it.
    pub fn check<K: GlKind>(&self, m: &Manifold, tol: f64, what: &str) -> MeshGLP<K::P, K::I> {
        assert_eq!(m.status(), Error::NoError, "{what}: boolean failed");
        assert!(m.num_tri() > 0, "{what}: empty result");
        let mesh = K::export(m);
        let pos = positions(&mesh);
        let num_run = mesh.run_original_id.len();
        let mut failures = Vec::new();
        for run in 0..num_run {
            let orig = mesh.run_original_id[run];
            let start = mesh.run_index[run].to_u64() as usize / 3;
            let end = mesh.run_index[run + 1].to_u64() as usize / 3;
            if start == end {
                continue;
            }
            let Some(planes) = self.sources.get(&orig) else {
                failures.push(format!(
                    "run {run}: original id {orig} was never registered"
                ));
                continue;
            };
            let xf = |p: Vec3| -> Vec3 {
                if mesh.run_transform.len() < 12 * (run + 1) {
                    return p;
                }
                let m = |k: usize| mesh.run_transform[12 * run + k].to_f64();
                Vec3::new(
                    m(0) * p.x + m(3) * p.y + m(6) * p.z + m(9),
                    m(1) * p.x + m(4) * p.y + m(7) * p.z + m(10),
                    m(2) * p.x + m(5) * p.y + m(8) * p.z + m(11),
                )
            };
            let backside = mesh.run_flags.get(run).is_some_and(|f| f & 1 != 0);
            // Source planes in world space, oriented as the run's triangles must be.
            let world: Vec<(u64, Vec3, Vec3)> = planes
                .iter()
                .map(|p| {
                    let w = [xf(p.v[0]), xf(p.v[1]), xf(p.v[2])];
                    let n = normalize(tri_normal(&w));
                    (p.face, if backside { -n } else { n }, w[0])
                })
                .collect();
            for tri in start..end {
                let face = mesh.face_id[tri].to_u64();
                let v = tri_verts(&mesh, &pos, tri);
                let n = tri_normal(&v);
                let longest = length(v[1] - v[0])
                    .max(length(v[2] - v[1]))
                    .max(length(v[0] - v[2]));
                // Height over the longest edge: a sliver's normal is noise.
                let height = if longest > 0.0 {
                    length(n) / longest
                } else {
                    0.0
                };
                let ok = world.iter().any(|&(f, wn, w0)| {
                    if f != face || !v.iter().all(|&p| dot(wn, p - w0).abs() <= tol) {
                        return false;
                    }
                    if length(n) == 0.0 {
                        return true;
                    }
                    let d = dot(normalize(n), wn);
                    if height > 100.0 * tol {
                        d > 0.999
                    } else {
                        d > 0.0
                    }
                });
                if !ok {
                    failures.push(format!(
                        "tri {tri} (run {run}, original {orig}, face {face}): {:?} {:?} {:?} lies on no plane of that face",
                        v[0], v[1], v[2]
                    ));
                }
            }
        }
        assert!(
            failures.is_empty(),
            "{what}: {} of {} triangles break provenance:\n{}",
            failures.len(),
            mesh.tri_verts.len() / 3,
            failures
                .iter()
                .take(20)
                .cloned()
                .collect::<Vec<_>>()
                .join("\n")
        );
        mesh
    }
}

/// Assert every exported triangle lying on the plane `z = height` facing +Z
/// with its centroid inside the XY rectangle `region` carries `face` from
/// original `orig`, and that at least one does.
pub(super) fn assert_floor_face<P: MeshPrecision, I: MeshIndex>(
    mesh: &MeshGLP<P, I>,
    height: f64,
    region: (Vec2, Vec2),
    tol: f64,
    orig: u32,
    face: u64,
) {
    let pos = positions(mesh);
    let mut found = 0;
    for run in 0..mesh.run_original_id.len() {
        let start = mesh.run_index[run].to_u64() as usize / 3;
        let end = mesh.run_index[run + 1].to_u64() as usize / 3;
        for tri in start..end {
            let v = tri_verts(mesh, &pos, tri);
            let n = tri_normal(&v);
            if length(n) == 0.0 || normalize(n).z < 0.999 {
                continue;
            }
            if !v.iter().all(|p| (p.z - height).abs() <= tol) {
                continue;
            }
            let c = (v[0] + v[1] + v[2]) / 3.0;
            if c.x < region.0.x || c.y < region.0.y || c.x > region.1.x || c.y > region.1.y {
                continue;
            }
            found += 1;
            assert_eq!(
                (mesh.run_original_id[run], mesh.face_id[tri].to_u64()),
                (orig, face),
                "floor triangle {tri} {:?} {:?} {:?} carries the wrong colour",
                v[0],
                v[1],
                v[2]
            );
        }
    }
    assert!(found > 0, "no triangle lies on the floor z = {height}");
}

/// Box palette: top and bottom share 0, the X sides 1, the Y sides 2.
pub(super) fn box_palette(n: Vec3, _c: Vec3) -> u64 {
    if n.z.abs() > 0.5 {
        0
    } else if n.x.abs() > 0.5 {
        1
    } else {
        2
    }
}

/// Cylinder palette: the curved side and the top share 0, the bottom is 1.
pub(super) fn cylinder_palette(n: Vec3, _c: Vec3) -> u64 {
    if n.z < -0.5 {
        1
    } else {
        0
    }
}

/// One colour for the whole mesh.
pub(super) fn solid(face: u64) -> impl Fn(Vec3, Vec3) -> u64 {
    move |_, _| face
}
