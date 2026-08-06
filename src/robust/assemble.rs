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
//
// When either operand carries vertex properties (colors, UVs, …), each
// output vertex's properties are barycentrically interpolated from its
// originating input triangle — exact rational barycentrics, one f64
// rounding — so constant per-operand properties survive exactly and
// interpolated ones agree with the exact engine to double precision.
// Coincident vertices with different properties stay separate property
// vertices linked by merge vectors, mirroring the exact engine's MeshGL
// output shape.

use num_rational::BigRational;
use num_traits::One;

use crate::linalg::Vec3;
use crate::manifold::Manifold;
use crate::types::MeshGL64;

use super::exact::rational::{rat_to_f64, R3};
use super::intersection_graph::Piece;
use super::tri_tri::dominant_axis;

/// Per-operand property data for interpolation. `props[m]` is flattened as
/// `props[m][(3*tri + corner) * num_prop[m] + channel]`, aligned with the
/// operand's soup triangle order (the `Piece::tri` indexing).
pub struct PropCtx<'a> {
    pub num_prop: [usize; 2],
    pub tris: [&'a [[Vec3; 3]]; 2],
    pub props: [&'a [f64]; 2],
}

impl<'a> PropCtx<'a> {
    pub fn out_num_prop(&self) -> usize {
        self.num_prop[0].max(self.num_prop[1])
    }
}

/// Exact barycentric coordinates of `p` on triangle `tri` (p must lie on the
/// triangle's plane), computed in the dominant-axis projection. The three
/// weights sum to exactly 1.
fn barycentric_r(p: &R3, tri: &[R3; 3]) -> [BigRational; 3] {
    use super::exact::predicates::tri_normal_r;
    let n = tri_normal_r(&tri[0], &tri[1], &tri[2]);
    let axis = dominant_axis(&n);
    let p2 = p.project_drop(axis);
    let a = tri[0].project_drop(axis);
    let b = tri[1].project_drop(axis);
    let c = tri[2].project_drop(axis);
    let total = b.sub(&a).cross(&c.sub(&a));
    let w0 = b.sub(&p2).cross(&c.sub(&p2)) / &total;
    let w1 = c.sub(&p2).cross(&a.sub(&p2)) / &total;
    let w2 = BigRational::one() - &w0 - &w1;
    [w0, w1, w2]
}

/// Interpolated properties (padded to `out` channels) for piece vertex `v`.
fn interpolate_props(ctx: &PropCtx, piece: &Piece, v: &R3, out: usize) -> Vec<f64> {
    let m = piece.mesh as usize;
    let np = ctx.num_prop[m];
    let mut result = vec![0.0f64; out];
    if np == 0 {
        return result;
    }
    let base = 3 * piece.tri * np;
    let corner = |i: usize| &ctx.props[m][base + i * np..base + (i + 1) * np];
    let (c0, c1, c2) = (corner(0), corner(1), corner(2));

    // Constant-per-face channels pass through exactly, no arithmetic.
    let all_const = (0..np).all(|k| c0[k] == c1[k] && c0[k] == c2[k]);
    if all_const {
        result[..np].copy_from_slice(c0);
        return result;
    }

    let t = ctx.tris[m][piece.tri];
    let corners = [
        R3::from_vec3(t[0]),
        R3::from_vec3(t[1]),
        R3::from_vec3(t[2]),
    ];
    let w = barycentric_r(v, &corners);
    let wf = [rat_to_f64(&w[0]), rat_to_f64(&w[1]), rat_to_f64(&w[2])];
    for k in 0..np {
        result[k] = if c0[k] == c1[k] && c0[k] == c2[k] {
            c0[k]
        } else {
            wf[0] * c0[k] + wf[1] * c1[k] + wf[2] * c2[k]
        };
    }
    result
}

/// Build the output manifold from every piece whose index passes `select`.
/// `verts` / `verts_f64` are the graph's interned tables: exact coordinates
/// for property interpolation, cached correctly rounded positions for the
/// output — no per-vertex rational rounding here.
/// With a `PropCtx` whose operands carry properties, output vertices get
/// interpolated properties; otherwise the output is positions-only and
/// byte-identical to the pre-property behavior.
pub fn assemble<F: Fn(usize) -> bool>(
    pieces: &[Piece],
    verts: &[R3],
    verts_f64: &[Vec3],
    select: F,
    props: Option<&PropCtx>,
) -> Manifold {
    let out_prop = props.map_or(0, |p| p.out_num_prop());

    // Property-vertex identity: interned position id + property bit pattern
    // (id equality is exact geometric identity — see VertInterner).
    type Key = (u32, Vec<u64>);
    let mut vert_index: std::collections::HashMap<Key, u64> = std::collections::HashMap::new();
    let mut vert_order: Vec<(u32, Vec<f64>)> = Vec::new();
    let mut tri_verts: Vec<u64> = Vec::new();

    for (pi, piece) in pieces.iter().enumerate() {
        if !select(pi) {
            continue;
        }
        for &vid in &piece.vi {
            let pvals = match props {
                Some(ctx) if out_prop > 0 => {
                    interpolate_props(ctx, piece, &verts[vid as usize], out_prop)
                }
                _ => Vec::new(),
            };
            let key = (vid, pvals.iter().map(|x| x.to_bits()).collect());
            let next = vert_order.len() as u64;
            let id = *vert_index.entry(key).or_insert_with(|| {
                vert_order.push((vid, pvals));
                next
            });
            tri_verts.push(id);
        }
    }
    if tri_verts.is_empty() {
        return Manifold::empty();
    }

    let stride = 3 + out_prop;
    let mut mesh = MeshGL64::default();
    mesh.num_prop = stride as u64;
    mesh.vert_properties = Vec::with_capacity(stride * vert_order.len());
    for (vid, pvals) in &vert_order {
        let p = verts_f64[*vid as usize];
        mesh.vert_properties.extend([p.x, p.y, p.z]);
        mesh.vert_properties.extend(pvals.iter());
    }
    mesh.tri_verts = tri_verts;

    // Coincident positions with different properties are distinct property
    // vertices; merge vectors tell the import they are topologically one.
    if out_prop > 0 {
        let mut by_pos: std::collections::HashMap<u32, u64> = std::collections::HashMap::new();
        for (i, (vid, _)) in vert_order.iter().enumerate() {
            match by_pos.get(vid) {
                Some(&first) => {
                    mesh.merge_from_vert.push(i as u64);
                    mesh.merge_to_vert.push(first);
                }
                None => {
                    by_pos.insert(*vid, i as u64);
                }
            }
        }
    }

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
