// robust/mod.rs — Robust boolean engine for general (possibly non-manifold)
// closed, orientable triangle meshes.
//
// Implements Barki, Guennebaud, Foufou 2015, "Exact, robust, and efficient
// regularized Booleans on general 3D meshes" (docs/Exact, robust, and
// efficient booleans.pdf). This engine is a parallel alternative to the
// ported exact pipeline in src/boolean3.rs: it requires inputs only to be
// geometrically closed and orientable (triangle soup is fine — connectivity
// is never trusted), at the cost of exact rational arithmetic on the hard
// predicate/construction cases.
//
// Selection between the two engines is via `types::BooleanEngine`
// (Exact | Robust | Auto); the exact engine remains the default and its
// behavior is byte-identical to before this module existed.
//
// Submodules (pipeline order):
//   exact              — rational points, filtered predicates, constructions
//   tri_tri            — exact triangle-triangle intersection (narrow phase)
//   arrangement        — per-triangle 2D arrangement of intersection prims
//   cdt                — exact constrained Delaunay triangulation
//   intersection_graph — broad phase, prim distribution, piece emission
//   cells              — arrangement cell complex + winding propagation
//   ray_shoot          — exact winding numbers (residual component seeds)
//   soup               — triangle-soup import (closed/orientable validation)
//
// Classification follows the mesh-arrangement formulation (Zhou, Grinspun,
// Zorin, Jacobson 2016 — what libigl's mesh_boolean uses), which subsumes
// both the paper's local Prop 2/3 ring walk and the per-component winding
// queries this engine used before. The arrangement's cells carry a winding
// number per operand, propagated combinatorially from the unbounded cell, so
// adjacent regions cannot disagree; each operand's solid is {w ≥ 1}, so a
// negative winding (an inverted region of a self-intersecting scan) is never
// material. The output keeps a wall exactly where the operation's predicate
// differs across it, wound from the cell labels rather than from the input
// face — which is what makes the result closed and consistently oriented no
// matter how the input was wound.
//
// `classify` and `propagate` are retained for the ring regularization and
// flood-fill they still provide to other callers; the boolean path no longer
// needs either, since thin material now cancels arithmetically in the
// winding sum instead of being discarded up front.

pub mod arrangement;
pub mod assemble;
pub mod cdt;
pub mod cells;
pub mod classify;
pub mod exact;
pub mod intersection_graph;
pub mod propagate;
pub mod ray_shoot;
pub mod soup;
pub mod tri_tri;

use crate::cancel::CancelToken;
use crate::impl_mesh::ManifoldImpl;
use crate::linalg::Vec3;
use crate::types::{Error, OpType};





fn is_cancelled(token: Option<&CancelToken>) -> bool {
    token.is_some_and(|t| t.is_cancelled())
}

fn cancelled_impl() -> ManifoldImpl {
    let mut out = ManifoldImpl::new();
    out.make_empty(Error::Cancelled);
    out
}

/// Robust boolean of two impls (manifold or soup). Same observable contract
/// as `boolean3::boolean_with_token`: intersect exactly, arrange +
/// retriangulate, build the arrangement's cell complex, propagate winding
/// numbers, and keep the walls the operation's predicate separates.
///
/// Every operation is one predicate on a cell's winding vector, so `Subtract`
/// needs no operand flip — it is simply "inside P and not inside Q".
pub fn boolean(
    a: &ManifoldImpl,
    b: &ManifoldImpl,
    op: OpType,
    token: Option<&CancelToken>,
) -> ManifoldImpl {
    if is_cancelled(token) {
        return cancelled_impl();
    }
    // Fast paths mirror the exact engine's observable behavior.
    if a.is_empty() {
        return match op {
            OpType::Add => b.clone(),
            OpType::Intersect | OpType::Subtract => ManifoldImpl::new(),
        };
    }
    if b.is_empty() {
        return match op {
            OpType::Add | OpType::Subtract => a.clone(),
            OpType::Intersect => ManifoldImpl::new(),
        };
    }
    let p_tris = soup::impl_to_tris(a);
    let q_tris = soup::impl_to_tris(b);
    let p_props = soup::impl_to_corner_props(a);
    let q_props = soup::impl_to_corner_props(b);

    if !a.bbox.does_overlap_box(&b.bbox) {
        match op {
            OpType::Add => {
                // Disjoint union: concatenate the soups and re-import. The
                // property context tags the two halves so each keeps its own
                // interpolated properties.
                let mut tris = p_tris.clone();
                tris.extend(q_tris.iter().cloned());
                let mut interner = intersection_graph::VertInterner::default();
                let pieces: Vec<intersection_graph::Piece> = tris
                    .iter()
                    .enumerate()
                    .map(|(i, t)| intersection_graph::Piece {
                        mesh: if i < p_tris.len() { 0 } else { 1 },
                        tri: if i < p_tris.len() { i } else { i - p_tris.len() },
                        vi: [
                            interner.intern_f64(t[0]),
                            interner.intern_f64(t[1]),
                            interner.intern_f64(t[2]),
                        ],
                    })
                    .collect();
                let ctx = assemble::PropCtx {
                    num_prop: [a.num_prop, b.num_prop],
                    tris: [&p_tris, &q_tris],
                    props: [&p_props, &q_props],
                };
                let props = (ctx.out_num_prop() > 0).then_some(&ctx);
                return assemble::assemble(
                    &pieces,
                    &interner.verts,
                    &interner.verts_f64,
                    |_| true,
                    props,
                )
                .into_impl();
            }
            OpType::Intersect => return ManifoldImpl::new(),
            OpType::Subtract => return a.clone(),
        }
    }

    // Subtraction needs no operand flip: the cell predicate expresses it
    // directly as "inside P and not inside Q", so both operands keep their
    // own winding and their corner properties stay in their original order.
    let graph = intersection_graph::build_graph(&p_tris, &q_tris);
    if is_cancelled(token) {
        return cancelled_impl();
    }
    let t_cells = crate::timing::start();
    let complex = cells::build_cells(&graph);
    crate::timing::print("robust: cell complex", t_cells);
    if is_cancelled(token) {
        return cancelled_impl();
    }

    // One exact query anchors each connected component; the rest of its
    // cells follow combinatorially.
    let t_winding = crate::timing::start();
    let wind = cells::windings(&graph, &complex, [&p_tris, &q_tris]);
    crate::timing::print("robust: winding propagation", t_winding);
    if is_cancelled(token) {
        return cancelled_impl();
    }

    // Boundary of the result, wound from the cell labels.
    let pieces = cells::extract(&graph, &complex, &wind, op);
    let ctx = assemble::PropCtx {
        num_prop: [a.num_prop, b.num_prop],
        tris: [&p_tris, &q_tris],
        props: [&p_props, &q_props],
    };
    let props = (ctx.out_num_prop() > 0).then_some(&ctx);
    let t_asm = crate::timing::start();
    let out = assemble::assemble(&pieces, &graph.verts, &graph.verts_f64, |_| true, props);
    crate::timing::print("robust: assemble+import", t_asm);
    out.into_impl()
}

/// Import a raw triangle list as a boolean result (used by
/// `boolean3::compose_meshes` when any input is a soup; positions only —
/// the property-aware disjoint-union path in `boolean` builds its own
/// tagged pieces).
pub(crate) fn assemble_all(tris: &[[Vec3; 3]]) -> ManifoldImpl {
    let mut interner = intersection_graph::VertInterner::default();
    let pieces: Vec<intersection_graph::Piece> = tris
        .iter()
        .enumerate()
        .map(|(i, t)| intersection_graph::Piece {
            mesh: 0,
            tri: i,
            vi: [
                interner.intern_f64(t[0]),
                interner.intern_f64(t[1]),
                interner.intern_f64(t[2]),
            ],
        })
        .collect();
    assemble::assemble(&pieces, &interner.verts, &interner.verts_f64, |_| true, None).into_impl()
}

#[cfg(test)]
#[path = "engine_tests.rs"]
mod engine_tests;

#[cfg(test)]
#[path = "cross_validation_tests.rs"]
mod cross_validation_tests;

#[cfg(test)]
#[path = "nonmanifold_tests.rs"]
mod nonmanifold_tests;

#[cfg(test)]
#[path = "property_tests.rs"]
mod property_tests;

#[cfg(test)]
#[path = "thingi_tests.rs"]
mod thingi_tests;
