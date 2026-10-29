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
//
// With a provenance context (robust/provenance.rs), each output triangle
// also carries the `TriRef` of the input triangle its piece lies on, so the
// simplification below merges only triangles of one source face and the
// export reports run original IDs and face IDs as for an exact-engine result.
//
// The assembled MeshGL64 carries no tolerance, so the import alone would
// reset the result to the automatic `kPrecision * bbox.scale()`. Callers
// pass the operands' [`CarriedTolerance`] instead, and the result keeps
// `max` of their epsilon and tolerance — what the exact engine's
// boolean_result (C++ boolean_result.cpp:832) and C++ `CsgLeafNode::Compose`
// give theirs. They are applied after the simplification below, which keeps
// running at the automatic values: the operands' tolerance is carried as
// metadata for what comes next (`simplify`, later booleans, export), not
// used to coarsen the robust result itself.

use super::exact::backend::{rat_one, Rational};

use crate::impl_mesh::ManifoldImpl;
use crate::linalg::Vec3;
use crate::manifold::Manifold;
use crate::types::MeshGL64;

use super::cells::VertTables;
use super::exact::rational::{rat_to_f64, R3};
use super::intersection_graph::Piece;
use super::provenance::RefCtx;
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

/// The epsilon and tolerance a robust result inherits from its operands:
/// the largest of each, as the exact engine and C++ `Compose` take them.
/// The result's own automatic floor (`kPrecision * bbox.scale()`, set by the
/// import) still applies, so a carried value only ever raises it.
///
/// They are carried as metadata: the assembly's own simplification runs at
/// the automatic values, so the robust engine's geometry does not depend on
/// them (see the note in [`assemble_with_tolerance`]). What they govern is
/// everything downstream — `simplify`, later booleans, the exported
/// `MeshGL.tolerance`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarriedTolerance {
    pub epsilon: f64,
    pub tolerance: f64,
}

impl CarriedTolerance {
    /// `max` over `operands`' epsilon and tolerance. Unset values (the -1 of
    /// a fresh impl) never win over a real one.
    pub fn of<'a, I: IntoIterator<Item = &'a ManifoldImpl>>(operands: I) -> Self {
        operands.into_iter().fold(
            Self {
                epsilon: -1.0,
                tolerance: -1.0,
            },
            |acc, m| Self {
                epsilon: acc.epsilon.max(m.epsilon),
                tolerance: acc.tolerance.max(m.tolerance),
            },
        )
    }

    /// Raise `imp`'s epsilon and tolerance to the carried values.
    pub fn apply(&self, imp: &mut ManifoldImpl) {
        imp.epsilon = imp.epsilon.max(self.epsilon);
        imp.tolerance = imp.tolerance.max(self.tolerance);
    }
}

/// Exact barycentric coordinates of `p` on triangle `tri` (p must lie on the
/// triangle's plane), computed in the dominant-axis projection. The three
/// weights sum to exactly 1.
fn barycentric_r(p: &R3, tri: &[R3; 3]) -> [Rational; 3] {
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
    let w2 = rat_one() - &w0 - &w1;
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
/// byte-identical to the pre-property behavior. The output is one fresh run
/// with no face IDs; [`assemble_with_refs`] keeps the source faces.
pub fn assemble<F: Fn(usize) -> bool>(
    pieces: &[Piece],
    verts: &[R3],
    verts_f64: &[Vec3],
    select: F,
    props: Option<&PropCtx>,
) -> Manifold {
    assemble_with_refs(pieces, verts, verts_f64, select, props, None)
}

/// [`assemble`] where, given a `RefCtx`, each output triangle carries the
/// `TriRef` of the input triangle its piece lies on (robust/provenance.rs).
pub fn assemble_with_refs<F: Fn(usize) -> bool>(
    pieces: &[Piece],
    verts: &[R3],
    verts_f64: &[Vec3],
    select: F,
    props: Option<&PropCtx>,
    refs: Option<&RefCtx>,
) -> Manifold {
    assemble_with_tolerance(pieces, verts, verts_f64, select, props, refs, None)
}

/// [`assemble_with_refs`] whose result keeps the operands' epsilon and
/// tolerance (`carried`) instead of the automatic ones its import computes.
/// Every robust boolean, rebuild and compose assembles through here with its
/// operands' values; `None` leaves the automatic ones.
pub fn assemble_with_tolerance<F: Fn(usize) -> bool>(
    pieces: &[Piece],
    verts: &[R3],
    verts_f64: &[Vec3],
    select: F,
    props: Option<&PropCtx>,
    refs: Option<&RefCtx>,
    carried: Option<CarriedTolerance>,
) -> Manifold {
    let out_prop = props.map_or(0, |p| p.out_num_prop());

    let selected: Vec<&Piece> = pieces
        .iter()
        .enumerate()
        .filter(|(pi, _)| select(*pi))
        .map(|(_, piece)| piece)
        .collect();
    if selected.is_empty() {
        return Manifold::empty();
    }

    // A boundary that touches itself along an edge carries more than two
    // half-edges on that vertex-id edge, which the import's id-based pairing
    // can only guess at. Splitting the pinched vertices into one copy per
    // geometric fan makes that pairing reproduce the geometry. The plan is
    // `None` — and everything below unchanged — for every mesh without such
    // an edge.
    let tris: Vec<[u32; 3]> = selected.iter().map(|piece| piece.vi).collect();
    let plan = super::pairing::plan_vertex_splits(&tris, VertTables { verts, verts_f64 });

    // Property-vertex identity: interned position id + fan copy + property
    // bit pattern (id equality is exact geometric identity — see
    // VertInterner).
    type Key = (u32, u32, Vec<u64>);
    // Fx hashing (unseeded): probe-only map — output vertex ids come from
    // `vert_order.len()` at first sight, i.e. from triangle/corner order.
    let mut vert_index: rustc_hash::FxHashMap<Key, u64> = rustc_hash::FxHashMap::default();
    let mut vert_order: Vec<(u32, u32, Vec<f64>)> = Vec::new();
    let mut tri_verts: Vec<u64> = Vec::new();

    for (t, piece) in selected.iter().enumerate() {
        for (c, &vid) in piece.vi.iter().enumerate() {
            let split = plan.as_ref().map_or(0, |p| p[3 * t + c]);
            let pvals = match props {
                Some(ctx) if out_prop > 0 => {
                    interpolate_props(ctx, piece, &verts[vid as usize], out_prop)
                }
                _ => Vec::new(),
            };
            let key = (vid, split, pvals.iter().map(|x| x.to_bits()).collect());
            let next = vert_order.len() as u64;
            let id = *vert_index.entry(key).or_insert_with(|| {
                vert_order.push((vid, split, pvals));
                next
            });
            tri_verts.push(id);
        }
    }

    let stride = 3 + out_prop;
    let mut mesh = MeshGL64::default();
    mesh.num_prop = stride as u64;
    mesh.vert_properties = Vec::with_capacity(stride * vert_order.len());
    for (vid, _, pvals) in &vert_order {
        let p = verts_f64[*vid as usize];
        mesh.vert_properties.extend([p.x, p.y, p.z]);
        mesh.vert_properties.extend(pvals.iter());
    }
    mesh.tri_verts = tri_verts;

    // Coincident positions with different properties are distinct property
    // vertices; merge vectors tell the import they are topologically one.
    // Keyed on the fan copy too, so split copies of a pinched vertex stay
    // separate geometric vertices.
    if out_prop > 0 {
        // Probe-only; merge pairs are emitted in `vert_order` index order.
        let mut by_pos: rustc_hash::FxHashMap<(u32, u32), u64> = rustc_hash::FxHashMap::default();
        for (i, (vid, split, _)) in vert_order.iter().enumerate() {
            match by_pos.get(&(*vid, *split)) {
                Some(&first) => {
                    mesh.merge_from_vert.push(i as u64);
                    mesh.merge_to_vert.push(first);
                }
                None => {
                    by_pos.insert((*vid, *split), i as u64);
                }
            }
        }
    }

    // The robust import handles everything rounding can produce: verts that
    // collapsed to identical f64 positions, exactly-degenerate triangles,
    // and non-manifold connectivity (kept as a soup impl).
    let tri_refs = refs.map(|ctx| {
        let (tri_ref, relations) = ctx.output_refs(&selected);
        crate::manifold::AssembledRefs { tri_ref, relations }
    });
    let out = Manifold::from_mesh_gl64_robust_assembled(&mesh, tri_refs);
    if out.status() != crate::types::Error::NoError || out.is_empty() {
        return out;
    }
    let mut imp = out.into_impl();

    // Manifold results get the same topology simplification the exact
    // engine's boolean_result applies: without it the CDT's coplanar
    // subdivision vertices survive and the output carries more (redundant)
    // vertices than the exact engine produces for the same inputs.
    //
    // The one stage held back is `swap_degenerates`, a deliberate difference
    // from the C++ (`Impl::SimplifyTopology`): the pieces of
    // `simplify_topology` are composed here without it, matching the import
    // above (`Manifold::from_mesh_gl64_robust_assembled`). `edge_op`'s own
    // functions are untouched and every other import path, the exact
    // engine's included, still runs the swap. The reason is an accuracy fix:
    // `face_op::set_normals_and_coplanar` (a faithful port of C++
    // `Impl::SetNormalsAndCoplanar`) flood-fills a seed triangle's normal onto
    // every coplanar neighbour with no orientation check, and a boolean result
    // legitimately contains coplanar *antiparallel* adjacencies (material
    // below the plane on one side of a shared edge, above it on the other).
    // A few large triangles get sign-flipped normals, the swap reads them as
    // degenerate (`ccw <= 0`) and swaps large non-planar quads, physically
    // moving the surface: −2.5e-3 of the robust union volume on Thingi10K
    // #301921 ∪ rotated-self (the extraction's 0.438010 survives cleanup and
    // is lost only in the swap). Skipping it keeps robust volumes
    // extraction-exact, at the cost of slightly larger triangle and vertex
    // counts than the exact engine's, which the sweep gate treats as
    // advisory. Pinned by
    // `robust::thingi_tests::thingi_301921_union_rotated_self_matches_exact`
    // and checked by the Monte-Carlo referee (`examples/volume_referee.rs`).
    // A twin of this assembly (manifold-sharp) must skip the same stage to
    // stay bit-equal.
    //
    // It runs at the automatic epsilon and tolerance the import set, not the
    // operands' carried ones, which are applied only afterwards. Collapsing
    // at a user's coarser tolerance (as the exact engine's boolean
    // simplification does for the vertices it creates) trades the robust
    // engine's volume exactness for fewer triangles, and without the
    // `swap_degenerates` clean-up that follows it there it measurably broke
    // results: an r=1 48-segment sphere and a copy rotated (13°, 7°, 3°)
    // and shifted 5e-4, unioned at tolerance 1e-2, came out self-intersecting
    // with genus -2, and their difference at an imported tolerance 1e-3 —
    // which raises epsilon — self-intersecting with genus 14, where the
    // automatic values give clean genus-0 / genus-20 results. A coarser
    // simplification stays one explicit `simplify()` away, at the carried
    // tolerance.
    if !imp.is_soup {
        crate::edge_op::cleanup_topology(&mut imp);
        crate::edge_op::collapse_short_edges(&mut imp, 0);
        crate::edge_op::collapse_colinear_edges(&mut imp, 0);
        crate::face_op::calculate_vert_normals(&mut imp);
        imp.remove_unreferenced_verts();
        imp.calculate_bbox();
        imp.sort_geometry();
    }
    if let Some(c) = carried {
        c.apply(&mut imp);
    }
    Manifold::from_impl(imp)
}
