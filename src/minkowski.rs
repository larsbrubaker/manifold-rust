// Copyright 2026 Lars Brubaker
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// Phase 17: Minkowski Sum/Difference — ported from C++ minkowski.cpp (175 lines)
//
// Implements the Minkowski sum/difference using:
// - Convex+Convex: pairwise vertex sums → Hull
// - NonConvex+Convex: per-triangle vertex sums → Hull → BatchBoolean (in batches)
// - NonConvex+NonConvex: per-face-pair sums with coplanarity filtering → BatchBoolean
//
// Progress and cancellation (`minkowski_with_progress`): both arrive the way
// the boolean pipeline takes them, as `Option`s whose `None` path is
// byte-for-byte the code that ran before they existed —
// `maybe_par_map_ct_progress` with no reporter and no token IS
// `maybe_par_map`, and `evaluate_with_token(None)` IS `evaluate()`. Nothing
// about the geometry moves.
//
//   Cancellation answers with an empty `ManifoldImpl` carrying
//   `Error::Cancelled`, as the boolean and the CSG tree do. Checks sit at the
//   entry, at every hull-map element, between batches / faces, and — the
//   load-bearing one — after the final batch boolean, so `cancel.rs`'s
//   invariant holds here too: a cancelled token can never produce a NoError
//   result. Nothing is checked per vertex; the hull loops stay as they were.
//
//   Progress reports one `Phase::Minkowski` unit per hull, one per batch
//   reduction and one for the closing merge, and finishes with
//   `complete_phase`, so a watched run ends on exactly 1.0.
//
// Shared with manifold-sharp's Minkowski.cs (its commit cf2c170).

use crate::cancel::{is_cancelled, CancelToken};
use crate::csg_tree::{CsgLeafNode, CsgNode};
use crate::impl_mesh::ManifoldImpl;
use crate::linalg::{dot, Vec3};
use crate::progress::{
    begin_phase, complete_phase, maybe_par_map_ct_progress, Phase, ProgressReporter,
};
use crate::quickhull;
use crate::types::OpType;

const BATCH_SIZE: usize = 1000;
const REDUCE_THRESHOLD: usize = 200;
const K_COPLANAR_TOL: f64 = 1e-12;

/// Compute the Minkowski sum or difference of two meshes.
/// Port of C++ Manifold::Impl::Minkowski()
///
/// `inset`: if true, computes the Minkowski difference (erosion);
///          if false, computes the Minkowski sum (dilation).
pub fn minkowski(a: &ManifoldImpl, b: &ManifoldImpl, inset: bool) -> ManifoldImpl {
    minkowski_with_progress(a, b, inset, None, None)
}

/// [`minkowski`] with cooperative cancellation and progress reporting.
///
/// A cancelled run answers with an empty mesh whose status is
/// [`crate::types::Error::Cancelled`], as the boolean pipeline does. The
/// reporter sees [`Phase::Minkowski`] only: one unit per hull, per batch
/// reduction and for the closing merge. `None` for both is exactly
/// [`minkowski`].
pub fn minkowski_with_progress(
    a: &ManifoldImpl,
    b: &ManifoldImpl,
    inset: bool,
    token: Option<&CancelToken>,
    progress: Option<&ProgressReporter>,
) -> ManifoldImpl {
    // Entry gate before any work, the shape every cancellable entry point
    // opens with: a pre-cancelled token must not buy a single hull.
    if is_cancelled(token) {
        return crate::boolean3::cancelled_impl();
    }

    let mut a_impl = a;
    let mut b_impl = b;

    let mut a_convex = a_impl.is_convex();
    let mut b_convex = b_impl.is_convex();

    // If the convex manifold was supplied first, swap them
    let (a_ref, b_ref);
    if a_convex && !b_convex {
        a_ref = b;
        b_ref = a;
        std::mem::swap(&mut a_convex, &mut b_convex);
        a_impl = a_ref;
        b_impl = b_ref;
    }

    // Early-exit if either input is empty
    if b_impl.is_empty() {
        report_trivial_completion(progress);
        return a_impl.clone();
    }
    if a_impl.is_empty() {
        report_trivial_completion(progress);
        return b_impl.clone();
    }

    // Costed from the branch about to run, because the three strategies are
    // orders of magnitude apart: a bar driven by a branch-blind total would
    // move at three different speeds. `work_units` repeats the branch
    // conditions below and must stay in step with them.
    begin_phase(
        progress,
        Phase::Minkowski,
        work_units(a_impl, b_impl, a_convex, b_convex, inset),
    );

    let mut composed_hulls: Vec<ManifoldImpl> = Vec::new();
    composed_hulls.push(a_impl.clone());

    // Convex-Convex Minkowski: Very Fast
    if !inset && a_convex && b_convex {
        let mut simple_hull: Vec<Vec3> =
            Vec::with_capacity(b_impl.vert_pos.len() * a_impl.vert_pos.len());
        for &a_vert in &a_impl.vert_pos {
            for &b_vert in &b_impl.vert_pos {
                simple_hull.push(a_vert + b_vert);
            }
        }
        composed_hulls.push(quickhull::convex_hull(&simple_hull));
        if let Some(p) = progress {
            p.advance(1);
        }

    // Convex + Non-Convex (or inset): Slower
    } else if (inset || !a_convex) && b_convex {
        let num_tri = a_impl.num_tri();

        // Process in batches. Each per-triangle hull is independent (C++ runs
        // this loop via for_each_n); results are collected in index order so
        // the batch content matches sequential. C++ pushes every hull
        // unconditionally — no empty filter here (unlike the
        // non-convex×non-convex branch); filtering would shift BatchBoolean
        // serials and change the reduction order.
        let mut offset = 0;
        while offset < num_tri {
            // Per-batch gate, the granularity the CSG tree's batch boolean
            // checks at: the hull map polls per element on its own, so this
            // catches a cancel that landed inside the previous batch's
            // boolean.
            if is_cancelled(token) {
                return crate::boolean3::cancelled_impl();
            }
            let num_iter = (num_tri - offset).min(BATCH_SIZE);
            let new_hulls = maybe_par_map_ct_progress(num_iter, 8, token, progress, |iter| {
                let tri = offset + iter;
                let mut simple_hull: Vec<Vec3> = Vec::with_capacity(3 * b_impl.vert_pos.len());
                for i in 0..3 {
                    let a_vert = a_impl.vert_pos[a_impl.halfedge[tri * 3 + i].start_vert as usize];
                    for &b_vert in &b_impl.vert_pos {
                        simple_hull.push(a_vert + b_vert);
                    }
                }
                quickhull::convex_hull(&simple_hull)
            });
            // `None` is the map's "a worker saw the flag"; the partial results
            // go with it.
            let Some(new_hulls) = new_hulls else {
                return crate::boolean3::cancelled_impl();
            };

            composed_hulls.push(batch_boolean_impls(&new_hulls, OpType::Add, token));
            if let Some(p) = progress {
                p.advance(1);
            }
            offset += BATCH_SIZE;
        }

    // Non-Convex + Non-Convex: Very Slow
    } else if !a_convex && !b_convex {
        let num_tri_a = a_impl.num_tri();
        let num_tri_b = b_impl.num_tri();

        let mut accumulated: Vec<ManifoldImpl> = Vec::new();

        for a_face in 0..num_tri_a {
            // Per-face gate: one A-face is a whole map plus a batch boolean,
            // so this is the coarsest granularity that still bounds cancel
            // latency by one face.
            if is_cancelled(token) {
                return crate::boolean3::cancelled_impl();
            }
            let a1 = a_impl.vert_pos[a_impl.halfedge[a_face * 3].start_vert as usize];
            let a2 = a_impl.vert_pos[a_impl.halfedge[a_face * 3 + 1].start_vert as usize];
            let a3 = a_impl.vert_pos[a_impl.halfedge[a_face * 3 + 2].start_vert as usize];
            let n_a = a_impl.face_normal[a_face];

            // Per-B-face hulls are independent (C++ parallel for_each_n over
            // bFace); collect in index order, then filter like C++'s
            // validFaceHulls pass so batch content and order match sequential.
            let hulls: Option<Vec<Option<ManifoldImpl>>> =
                maybe_par_map_ct_progress(num_tri_b, 8, token, progress, |b_face| {
                    let n_b = b_impl.face_normal[b_face];
                    let dot_same = dot(n_a, n_b);
                    let dot_opp = dot(n_a, Vec3::new(-n_b.x, -n_b.y, -n_b.z));
                    let coplanar = (dot_same - 1.0).abs() < K_COPLANAR_TOL
                        || (dot_opp - 1.0).abs() < K_COPLANAR_TOL;
                    if coplanar {
                        return None;
                    }

                    let b1 = b_impl.vert_pos[b_impl.halfedge[b_face * 3].start_vert as usize];
                    let b2 = b_impl.vert_pos[b_impl.halfedge[b_face * 3 + 1].start_vert as usize];
                    let b3 = b_impl.vert_pos[b_impl.halfedge[b_face * 3 + 2].start_vert as usize];

                    Some(quickhull::convex_hull(&[
                        a1 + b1,
                        a1 + b2,
                        a1 + b3,
                        a2 + b1,
                        a2 + b2,
                        a2 + b3,
                        a3 + b1,
                        a3 + b2,
                        a3 + b3,
                    ]))
                });
            let Some(hulls) = hulls else {
                return crate::boolean3::cancelled_impl();
            };
            let mut face_hulls: Vec<ManifoldImpl> = Vec::new();
            for hull in hulls.into_iter().flatten() {
                if !hull.is_empty() {
                    face_hulls.push(hull);
                }
            }

            if !face_hulls.is_empty() {
                accumulated.push(batch_boolean_impls(&face_hulls, OpType::Add, token));
            }

            // Periodically reduce to limit memory
            if accumulated.len() >= REDUCE_THRESHOLD {
                let reduced = batch_boolean_impls(&accumulated, OpType::Add, token);
                accumulated.clear();
                accumulated.push(reduced);
            }

            // One unit per A-face, whether or not it contributed a hull: the
            // face is the work item, and a coplanar-filtered face still cost a
            // whole map.
            if let Some(p) = progress {
                p.advance(1);
            }
        }

        if !accumulated.is_empty() {
            composed_hulls.push(batch_boolean_impls(&accumulated, OpType::Add, token));
        }
    }

    // Final merge; C++ finishes with AsOriginal() = InitializeOriginal +
    // SetNormalsAndCoplanar.
    let op = if inset { OpType::Subtract } else { OpType::Add };
    let mut out = batch_boolean_impls(&composed_hulls, op, token);

    // The closing check `cancel.rs`'s invariant requires — "a cancelled token
    // can never produce a NoError result". Without it a cancel that landed
    // inside this last reduction would come back as an empty Cancelled leaf,
    // and the two calls below would dress it up into an empty NoError mesh.
    if is_cancelled(token) {
        return crate::boolean3::cancelled_impl();
    }
    out.initialize_original();
    out.set_normals_and_coplanar();

    // `complete_phase`, not `advance(1)`, for the closing merge's unit: the
    // throttle only emits on a step boundary, so on a run of more than 100
    // units the trailing advances are swallowed (510/514 on a 512-triangle
    // erosion). Reached only on success: the cancelled returns above all skip
    // it, because a full bar is a claim that the work was done.
    complete_phase(progress);
    out
}

/// How many progress units the branch that is about to run will report: one
/// per hull, one per batch reduction and one for the closing merge, counted
/// exactly, so the fractions are scaled against the work that actually runs.
/// The closing merge's unit is reported by `complete_phase` rather than an
/// `advance`. The conditions repeat the branch chain in
/// [`minkowski_with_progress`] and have to be changed with it.
fn work_units(
    a_impl: &ManifoldImpl,
    b_impl: &ManifoldImpl,
    a_convex: bool,
    b_convex: bool,
    inset: bool,
) -> u64 {
    // Every branch finishes with one batch boolean over `composed_hulls`.
    const FINAL_MERGE: u64 = 1;
    if !inset && a_convex && b_convex {
        // One hull over the pairwise vertex sums.
        return 1 + FINAL_MERGE;
    }
    if (inset || !a_convex) && b_convex {
        let num_tri = a_impl.num_tri() as u64;
        let batches = num_tri.div_ceil(BATCH_SIZE as u64);
        return num_tri + batches + FINAL_MERGE;
    }
    if !a_convex && !b_convex {
        // The periodic REDUCE_THRESHOLD merges are folded into the per-A-face
        // unit: how many run depends on how many faces survive the
        // coplanarity filter, which is not knowable up front.
        let num_tri_a = a_impl.num_tri() as u64;
        let num_tri_b = b_impl.num_tri() as u64;
        return num_tri_a * num_tri_b + num_tri_a + FINAL_MERGE;
    }
    // Unreachable: the swap above puts the convex operand second, so a convex
    // `a` with a non-convex `b` cannot arrive here. Costed as the merge alone
    // so that a future fourth branch reports something.
    FINAL_MERGE
}

/// A whole bar for a call that returns without computing anything — the
/// empty-operand early exits — so a watching caller sees a finished
/// operation instead of a bar that never started.
fn report_trivial_completion(progress: Option<&ProgressReporter>) {
    begin_phase(progress, Phase::Minkowski, 1);
    complete_phase(progress);
}

/// Helper: BatchBoolean on ManifoldImpl directly via the CSG tree.
/// `evaluate_with_token(None)` is what `evaluate()` calls, so the untokened
/// path is unchanged.
fn batch_boolean_impls(
    meshes: &[ManifoldImpl],
    op: OpType,
    token: Option<&CancelToken>,
) -> ManifoldImpl {
    if meshes.is_empty() {
        return ManifoldImpl::new();
    }
    if meshes.len() == 1 {
        return meshes[0].clone();
    }

    let children: Vec<CsgNode> = meshes
        .iter()
        .map(|m| CsgNode::leaf_node(CsgLeafNode::new(m.clone())))
        .collect();
    let tree = CsgNode::op_n(op, children);
    tree.evaluate_with_token(token)
}

/// Convenience wrapper: Minkowski sum (dilation).
pub fn minkowski_sum(a: &ManifoldImpl, b: &ManifoldImpl) -> ManifoldImpl {
    minkowski(a, b, false)
}

/// Convenience wrapper: Minkowski difference (erosion).
pub fn minkowski_difference(a: &ManifoldImpl, b: &ManifoldImpl) -> ManifoldImpl {
    minkowski(a, b, true)
}

/// [`minkowski_sum`] with cancellation and progress; see
/// [`minkowski_with_progress`].
pub fn minkowski_sum_with_progress(
    a: &ManifoldImpl,
    b: &ManifoldImpl,
    token: Option<&CancelToken>,
    progress: Option<&ProgressReporter>,
) -> ManifoldImpl {
    minkowski_with_progress(a, b, false, token, progress)
}

/// [`minkowski_difference`] with cancellation and progress; see
/// [`minkowski_with_progress`].
pub fn minkowski_difference_with_progress(
    a: &ManifoldImpl,
    b: &ManifoldImpl,
    token: Option<&CancelToken>,
    progress: Option<&ProgressReporter>,
) -> ManifoldImpl {
    minkowski_with_progress(a, b, true, token, progress)
}

#[cfg(test)]
#[path = "minkowski_union_regression_tests.rs"]
mod union_regression_tests;

#[cfg(test)]
#[path = "minkowski_progress_tests.rs"]
mod progress_tests;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linalg::{mat4_to_mat3x4, scaling_matrix, translation_matrix, Mat3x4};

    #[test]
    fn test_convex_convex_minkowski_sum() {
        let a = ManifoldImpl::cube(&Mat3x4::identity());
        let b = ManifoldImpl::cube(&Mat3x4::identity());
        let sum = minkowski_sum(&a, &b);
        assert!(
            sum.num_tri() > 0,
            "Minkowski sum should produce non-empty mesh"
        );
        // Two unit cubes: Minkowski sum should be a 2×2×2 cube
        let vol = sum.get_property(crate::properties::Property::Volume).abs();
        assert!(
            (vol - 8.0).abs() < 0.5,
            "Minkowski sum of two unit cubes should have volume ~8, got {}",
            vol
        );
    }

    #[test]
    fn test_convex_convex_minkowski_difference() {
        let a = ManifoldImpl::cube(&mat4_to_mat3x4(scaling_matrix(Vec3::splat(2.0))));
        let b = ManifoldImpl::cube(&mat4_to_mat3x4(
            translation_matrix(Vec3::splat(-0.25)) * scaling_matrix(Vec3::splat(0.5)),
        ));
        let diff = minkowski_difference(&a, &b);
        assert!(
            diff.num_tri() > 0,
            "Minkowski difference should produce non-empty mesh"
        );
    }

    #[test]
    fn test_empty_minkowski() {
        let a = ManifoldImpl::cube(&Mat3x4::identity());
        let b = ManifoldImpl::new();
        let sum = minkowski_sum(&a, &b);
        // If b is empty, result should be a
        assert_eq!(sum.num_tri(), a.num_tri());
    }
}
