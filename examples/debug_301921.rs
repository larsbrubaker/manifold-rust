// Scratch diagnostic for Thingi10K #301921 ∪ rotated self: robust loses ~7%
// volume with a consistent arrangement. Locates points inside the exact
// result but outside the robust result, then reports the input windings at
// those points and the resolved cell windings of the nearest pieces.

use manifold_rust::linalg::Vec3;
use manifold_rust::manifold::Manifold;
use manifold_rust::robust::exact::rational::R3;
use manifold_rust::robust::ray_shoot::{winding_number_indexed, WindingIndex};
use manifold_rust::robust::{cells, intersection_graph, soup};
use manifold_rust::types::{BooleanEngine, MeshGL, OpType, WindingRule};

fn import(path: &str) -> Manifold {
    let data = std::fs::read(path).expect("read stl");
    let n = u32::from_le_bytes(data[80..84].try_into().unwrap()) as usize;
    let mut positions: Vec<f32> = Vec::new();
    for f in 0..n.min((data.len() - 84) / 50) {
        let base = 84 + f * 50 + 12;
        for v in 0..9 {
            let o = base + v * 4;
            positions.push(f32::from_le_bytes(data[o..o + 4].try_into().unwrap()));
        }
    }
    let nv = positions.len() / 3;
    let (mut min, mut max) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for i in 0..nv {
        for k in 0..3 {
            let v = positions[i * 3 + k] as f64;
            min[k] = min[k].min(v);
            max[k] = max[k].max(v);
        }
    }
    let c = [(min[0] + max[0]) / 2.0, (min[1] + max[1]) / 2.0, (min[2] + max[2]) / 2.0];
    let side = (max[0] - min[0]).max(max[1] - min[1]).max(max[2] - min[2]);
    let s = if side > 0.0 { 2.0 / side } else { 1.0 };
    for i in 0..nv {
        for k in 0..3 {
            positions[i * 3 + k] = ((positions[i * 3 + k] as f64 - c[k]) * s) as f32;
        }
    }
    let mut mesh = MeshGL::default();
    mesh.num_prop = 3;
    mesh.tri_verts = (0..nv as u32).collect();
    mesh.vert_properties = positions;
    mesh.merge();
    Manifold::from_mesh_gl_robust(&mesh)
}

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> f64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        ((z ^ (z >> 31)) >> 11) as f64 / (1u64 << 53) as f64
    }
}

fn main() {
    let a = import("src/robust/testdata/301921.stl");
    let b = a.rotate(30.0, 45.0, 60.0).translate(Vec3::new(0.3, 0.0, 0.0));
    let robust = a.union_with_engine(&b, BooleanEngine::Robust);
    let exact = a.union_with_engine(&b, BooleanEngine::Exact);
    println!(
        "robust vol {:.6} ({} tris), exact vol {:.6} ({} tris)",
        robust.volume(),
        robust.num_tri(),
        exact.volume(),
        exact.num_tri()
    );

    let p = soup::impl_to_tris(a.as_impl());
    let q = soup::impl_to_tris(b.as_impl());
    let r_tris = soup::impl_to_tris(robust.as_impl());
    let e_tris = soup::impl_to_tris(exact.as_impl());
    let idx_p = WindingIndex::new(&p);
    let idx_q = WindingIndex::new(&q);
    let idx_r = WindingIndex::new(&r_tris);
    let idx_e = WindingIndex::new(&e_tris);

    let (mut min, mut max) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for t in p.iter().chain(q.iter()) {
        for v in t {
            let c = [v.x, v.y, v.z];
            for k in 0..3 {
                min[k] = min[k].min(c[k]);
                max[k] = max[k].max(c[k]);
            }
        }
    }
    let ext = [max[0] - min[0], max[1] - min[1], max[2] - min[2]];

    let graph = intersection_graph::build_graph(&p, &q);
    let complex = cells::build_cells(&graph);
    let wind = cells::windings(&graph, &complex, [&p, &q]);
    println!(
        "{} pieces, {} cells, complete {}",
        graph.pieces.len(),
        complex.num_cells,
        wind.complete()
    );

    // Volume of the raw extraction, before assembly and simplification:
    // isolates whether the material is lost at classification or later.
    let pieces = cells::extract(&graph, &complex, &wind, OpType::Add, WindingRule::Positive);
    let mut vol6 = 0.0f64;
    for pc in &pieces {
        let v = [
            graph.verts_f64[pc.vi[0] as usize],
            graph.verts_f64[pc.vi[1] as usize],
            graph.verts_f64[pc.vi[2] as usize],
        ];
        vol6 += v[0].x * (v[1].y * v[2].z - v[2].y * v[1].z)
            - v[1].x * (v[0].y * v[2].z - v[2].y * v[0].z)
            + v[2].x * (v[0].y * v[1].z - v[1].y * v[0].z);
    }
    println!(
        "extraction: {} pieces, signed volume {:.6}",
        pieces.len(),
        vol6 / 6.0
    );

    // Replicate assemble's stages with a volume probe after each: weld to
    // MeshGL64, robust import, then simplify_topology.
    {
        use manifold_rust::types::MeshGL64;
        let mut vert_index: std::collections::HashMap<u32, u64> = std::collections::HashMap::new();
        let mut vert_order: Vec<u32> = Vec::new();
        let mut tri_verts: Vec<u64> = Vec::new();
        for pc in &pieces {
            for &vid in &pc.vi {
                let next = vert_order.len() as u64;
                let id = *vert_index.entry(vid).or_insert_with(|| {
                    vert_order.push(vid);
                    next
                });
                tri_verts.push(id);
            }
        }
        let mut mesh = MeshGL64::default();
        mesh.num_prop = 3;
        mesh.vert_properties = Vec::with_capacity(3 * vert_order.len());
        for vid in &vert_order {
            let p = graph.verts_f64[*vid as usize];
            mesh.vert_properties.extend([p.x, p.y, p.z]);
        }
        mesh.tri_verts = tri_verts;
        let imported = Manifold::from_mesh_gl64_robust(&mesh);
        println!(
            "after import: status {:?}, soup {}, {} tris, volume {:.6}",
            imported.status(),
            imported.as_impl().is_soup,
            imported.num_tri(),
            imported.volume()
        );
        let mut imp = imported.into_impl();
        manifold_rust::edge_op::simplify_topology(&mut imp, 0);
        imp.remove_unreferenced_verts();
        imp.calculate_bbox();
        imp.sort_geometry();
        let simplified = Manifold::from_impl(imp);
        println!(
            "after simplify: {} tris, volume {:.6}",
            simplified.num_tri(),
            simplified.volume()
        );
    }

    // Points inside exact result but outside robust result.
    let mut rng = Rng(42);
    let mut missing: Vec<Vec3> = Vec::new();
    let mut n_probe = 0;
    while missing.len() < 12 && n_probe < 400_000 {
        n_probe += 1;
        let v = Vec3::new(
            min[0] + ext[0] * rng.next(),
            min[1] + ext[1] * rng.next(),
            min[2] + ext[2] * rng.next(),
        );
        let pt = R3::from_vec3(v);
        let in_e = winding_number_indexed(&pt, &e_tris, &idx_e) >= 1;
        let in_r = winding_number_indexed(&pt, &r_tris, &idx_r) >= 1;
        if in_e && !in_r {
            let wa = winding_number_indexed(&pt, &p, &idx_p);
            let wb = winding_number_indexed(&pt, &q, &idx_q);
            let wr = winding_number_indexed(&pt, &r_tris, &idx_r);
            println!(
                "missing point ({:.5},{:.5},{:.5}): input windings A={wa} B={wb}, robust output winding {wr}",
                v.x, v.y, v.z,
            );
            missing.push(v);
        }
    }
    println!("({} probes for {} missing points)", n_probe, missing.len());

    // For each missing point, the nearest piece centroid and its cells.
    for v in &missing {
        let mut best = (f64::INFINITY, usize::MAX);
        for pi in 0..graph.pieces.len() {
            let pv = graph.pieces[pi].vi;
            let c = [
                graph.verts_f64[pv[0] as usize],
                graph.verts_f64[pv[1] as usize],
                graph.verts_f64[pv[2] as usize],
            ];
            let cx = Vec3::new(
                (c[0].x + c[1].x + c[2].x) / 3.0,
                (c[0].y + c[1].y + c[2].y) / 3.0,
                (c[0].z + c[1].z + c[2].z) / 3.0,
            );
            let dv = cx - *v;
            let d = (dv.x * dv.x + dv.y * dv.y + dv.z * dv.z).sqrt();
            if d < best.0 {
                best = (d, pi);
            }
        }
        let pi = best.1;
        let (cn, ca) = (
            complex.cell(pi, cells::NORMAL),
            complex.cell(pi, cells::ANTI),
        );
        println!(
            "point ({:.4},{:.4},{:.4}): nearest piece {} (mesh {} tri {}, dist {:.5}) cells n{}={:?} a{}={:?} in_result n={} a={}",
            v.x, v.y, v.z,
            pi,
            graph.pieces[pi].mesh,
            graph.pieces[pi].tri,
            best.0,
            cn,
            wind.w[cn],
            ca,
            wind.w[ca],
            cells::in_result(OpType::Add, WindingRule::Positive, wind.w[cn]),
            cells::in_result(OpType::Add, WindingRule::Positive, wind.w[ca]),
        );
    }
}
