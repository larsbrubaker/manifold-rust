// robust/assemble.rs — From tagged pieces to a Manifold (paper §7.5, output
// browse).
//
// The selected pieces are welded on their exact rational coordinates first
// (so identical points are identical regardless of construction path), then
// each unique vertex rounds once to the nearest f64
// (robust/exact/rational.rs) and the result re-enters the library through
// the robust MeshGL64 import: manifold results get the full strict pipeline
// (normals, degenerate removal, sorting — the same post-processing the
// exact engine's outputs receive), while legitimately non-manifold results
// (booleans of non-manifold inputs) are retained as soup impls, ready for
// chained robust operations.

use std::collections::BTreeMap;

use crate::manifold::Manifold;
use crate::types::MeshGL64;

use super::exact::rational::R3;
use super::intersection_graph::Piece;

/// Build the output manifold from every piece whose index passes `select`.
pub fn assemble<F: Fn(usize) -> bool>(pieces: &[Piece], select: F) -> Manifold {
    let mut vert_index: BTreeMap<R3, u64> = BTreeMap::new();
    let mut vert_order: Vec<R3> = Vec::new();
    let mut tri_verts: Vec<u64> = Vec::new();

    for (pi, piece) in pieces.iter().enumerate() {
        if !select(pi) {
            continue;
        }
        for v in &piece.v {
            let next = vert_order.len() as u64;
            let id = *vert_index.entry(v.clone()).or_insert_with(|| {
                vert_order.push(v.clone());
                next
            });
            tri_verts.push(id);
        }
    }
    if tri_verts.is_empty() {
        return Manifold::empty();
    }

    let mut mesh = MeshGL64::default();
    mesh.num_prop = 3;
    mesh.vert_properties = Vec::with_capacity(3 * vert_order.len());
    for v in &vert_order {
        let p = v.to_vec3_rounded();
        mesh.vert_properties.extend([p.x, p.y, p.z]);
    }
    mesh.tri_verts = tri_verts;

    // The robust import handles everything rounding can produce: verts that
    // collapsed to identical f64 positions, exactly-degenerate triangles,
    // and non-manifold connectivity (kept as a soup impl).
    let out = Manifold::from_mesh_gl64_robust(&mesh);

    // Manifold results get the same topology simplification the exact
    // engine's boolean_result applies: without it the CDT's coplanar
    // subdivision vertices survive and the output carries more (redundant)
    // vertices than the exact engine produces for the same inputs.
    if out.status() == crate::types::Error::NoError && !out.as_impl().is_soup && !out.is_empty() {
        let mut imp = out.into_impl();
        crate::edge_op::simplify_topology(&mut imp, 0);
        imp.remove_unreferenced_verts();
        imp.calculate_bbox();
        imp.sort_geometry();
        return Manifold::from_impl(imp);
    }
    out
}
