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
//   classify           — radial rings, Prop 2/3 union/intersection tagging
//   propagate          — per-mesh tag flood fill between intersection cuts
//   ray_shoot          — exact winding numbers for untouched components
//   soup               — triangle-soup import (closed/orientable validation)

pub mod arrangement;
pub mod assemble;
pub mod cdt;
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

use classify::Tag;
use ray_shoot::{piece_centroid, winding_number};

fn is_cancelled(token: Option<&CancelToken>) -> bool {
    token.is_some_and(|t| t.is_cancelled())
}

fn cancelled_impl() -> ManifoldImpl {
    let mut out = ManifoldImpl::new();
    out.make_empty(Error::Cancelled);
    out
}

/// Robust boolean of two impls (manifold or soup). Same observable contract
/// as `boolean3::boolean_with_token`, computed by the Barki 2015 pipeline:
/// intersect exactly, arrange + retriangulate, classify pieces into
/// union/intersection sets, assemble the requested one.
///
/// `Subtract` uses the identity P − Q = P ∩ Q^c: Q's winding is flipped on
/// the working copy and the intersection set is assembled; the ray-shooting
/// fallback interprets the flipped operand as its (unbounded) complement.
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
    if !a.bbox.does_overlap_box(&b.bbox) {
        match op {
            OpType::Add => {
                // Disjoint union: concatenate the soups and re-import.
                let mut tris = soup::impl_to_tris(a);
                tris.extend(soup::impl_to_tris(b));
                return assemble_all(&tris);
            }
            OpType::Intersect => return ManifoldImpl::new(),
            OpType::Subtract => return a.clone(),
        }
    }

    let p_tris = soup::impl_to_tris(a);
    let complement = op == OpType::Subtract;
    let mut q_tris = soup::impl_to_tris(b);
    if complement {
        for t in &mut q_tris {
            t.swap(1, 2);
        }
    }

    let graph = intersection_graph::build_graph(&p_tris, &q_tris);
    if is_cancelled(token) {
        return cancelled_impl();
    }
    let cls = classify::classify_rings(&graph);
    if is_cancelled(token) {
        return cancelled_impl();
    }
    let prop = propagate::propagate(&graph, &cls);
    let mut tags = prop.tags;
    for &(root, rep) in &prop.untagged {
        if is_cancelled(token) {
            return cancelled_impl();
        }
        let piece = &graph.pieces[rep];
        let (other, other_is_complement): (&[[Vec3; 3]], bool) = if piece.mesh == 0 {
            (&q_tris, complement)
        } else {
            (&p_tris, false)
        };
        let w = winding_number(&piece_centroid(&piece.v), other);
        let inside = if other_is_complement { w == 0 } else { w != 0 };
        let tag = if inside { Tag::Inter } else { Tag::Union };
        for pi in 0..graph.pieces.len() {
            if !cls.discarded[pi] && prop.component[pi] == root {
                tags[pi] = Some(tag);
            }
        }
    }

    let want = match op {
        OpType::Add => Tag::Union,
        OpType::Subtract | OpType::Intersect => Tag::Inter,
    };
    let out = assemble::assemble(&graph.pieces, |pi| {
        !cls.discarded[pi] && tags[pi] == Some(want)
    });
    out.into_impl()
}

/// Import a raw triangle list as a boolean result (used by the disjoint
/// union fast path and by `boolean3::compose_meshes` when any input is a
/// soup).
pub(crate) fn assemble_all(tris: &[[Vec3; 3]]) -> ManifoldImpl {
    use exact::rational::R3;
    let pieces: Vec<intersection_graph::Piece> = tris
        .iter()
        .enumerate()
        .map(|(i, t)| intersection_graph::Piece {
            mesh: 0,
            tri: i,
            v: [
                R3::from_vec3(t[0]),
                R3::from_vec3(t[1]),
                R3::from_vec3(t[2]),
            ],
        })
        .collect();
    assemble::assemble(&pieces, |_| true).into_impl()
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
