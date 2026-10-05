// convex_dilation.rs — NOT A C++ PORT; a 1:1 mirror of manifold-sharp's
// ConvexDilation.cs (its RUST_DIVERGENCES.md entry 6; here
// docs/CPP_DIVERGENCES.md entry 13). minkowski.rs dilates a non-convex solid by
// a convex tool with one hull per triangle, then `batch_union` 1000 hulls at a
// time. This answers the same question with the same hulls but a different
// reduction — an ADDED entry point (`Manifold::try_dilate_by_convex` /
// `try_erode_by_convex`); `minkowski::minkowski` is not routed through it.
//
// ── Why a different reduction ───────────────────────────────────────────────
// batch_boolean pops the largest-vertex mesh first, so within a batch the
// growing union is a serial chain and most of a dilation runs on one core. A
// balanced tree exposes the work to the `parallel` feature (a single-threaded
// balanced tree was measured and was not faster). Measured in manifold-sharp:
// a 2642-triangle part 13.3 s → 2.96 s, identical volume.
//
// ── The tree ────────────────────────────────────────────────────────────────
//   1. Leaves: the solid on its own (dilation only), then runs of LEAF_SIZE
//      hull units in seed order — a unit being one triangle or, when dilating,
//      a convex patch (convex_patches.rs proves its hull lies in the dilation)
//      — each leaf building its units' hulls in minkowski.rs's vertex-sum order
//      and unioning them through `csg_tree::batch_union`.
//   2. Pairwise levels: node k is the union of nodes 2k and 2k+1 of the level
//      below; an odd last node is carried up unchanged. Each level is one map.
//
// ── Determinism ─────────────────────────────────────────────────────────────
// The leaf and level maps go through `maybe_par_map_ct(_progress)`, so the
// `parallel` feature governs them (single-threaded on wasm and without it).
// Each worker writes its own slot and reads only the input meshes or the
// previous, complete level. The engine is read once, before the tree. The
// tree's shape is a function of the input alone, so a parallel run performs
// the same booleans on the same operands as a sequential one — the same bits,
// up to the process-global mesh IDs, which the closing initialize_original
// replaces anyway.
//
// ── Erosion: the same tree, minus the solid ─────────────────────────────────
// minkowski.rs's inset branch computes A \ (boundary(A) ⊕ B) from the same
// per-triangle hulls. try_compute_erosion is that with this reduction: hull
// leaves WITHOUT the solid leaf, then one solid − union on the same engine.
// One routine (try_reduce) serves both. Erosion keeps one hull per triangle (a
// patch hull would carve kept material) and takes a convex solid.
//
// ── Nested and crossing shells ──────────────────────────────────────────────
// Shells that nest or cross wind 2 where they overlap, which the exact
// engine's unions are not defined for. That needs overlapping component boxes,
// so when any pair overlaps (has_overlapping_components) the robust
// `rebuild_with_rule(Positive)` first turns the solid into the union of its
// shells; a rebuild that is not a clean manifold declines.
//
// ── Progress from inside the top unions ─────────────────────────────────────
// Levels of at most SUB_PROGRESS_MAX_PAIRS unions (and erosion's closing
// subtraction) hand each exact boolean a stage sink (boolean_stage_progress.rs)
// feeding one NodeProgress per level, which reports finished nodes plus every
// running node's fraction through `report_units`. Monotone under the tracker's
// lock; those levels advance the reporter only after their map returns.
//
// ── What it is not ──────────────────────────────────────────────────────────
// Not bit-identical to minkowski_sum: the same hulls unioned in a different
// order round intersection vertices differently. Volume and genus agree.

use std::cell::Cell;
use std::sync::Mutex;

use crate::boolean3::{boolean_dispatch, boolean_with_token_and_stage, cancelled_impl};
use crate::cancel::{is_cancelled, CancelToken};
use crate::convex_patches;
use crate::csg_tree::{batch_union, CsgLeafNode};
use crate::disjoint_sets::DisjointSets;
use crate::impl_mesh::ManifoldImpl;
use crate::linalg::Vec3;
use crate::par::maybe_par_map_ct;
use crate::progress::{
    begin_phase, complete_phase, maybe_par_map_ct_progress, Phase, ProgressReporter,
};
use crate::quickhull::convex_hull;
use crate::types::{BooleanConfig, BooleanEngine, Error, OpType, WindingRule};

/// Hull units per leaf.
const LEAF_SIZE: usize = 16;

/// Parallel threshold for the leaf and level maps.
const UNION_PAR_THRESHOLD: usize = 2;

/// Levels with at most this many unions report from inside their booleans.
const SUB_PROGRESS_MAX_PAIRS: usize = 8;

thread_local! {
    /// Test hook: how many nesting rebuilds ran on this thread.
    pub(crate) static REBUILDS_RUN: Cell<usize> = const { Cell::new(0) };
    /// Test hook: overrides the dilation patch cap on this thread.
    pub(crate) static PATCH_SIZE_OVERRIDE: Cell<Option<usize>> = const { Cell::new(None) };
    /// Test hook: the unit (hull) count of the last run on this thread.
    pub(crate) static LAST_HULL_COUNT: Cell<usize> = const { Cell::new(0) };
}

/// The dilation `solid ⊕ tool` for a non-convex solid and a convex tool,
/// through the balanced union tree. `None` declines (the caller runs
/// `minkowski_sum`); a cancelled run is `Some` of an `Error::Cancelled` impl.
pub fn try_compute(
    solid: &ManifoldImpl,
    tool: &ManifoldImpl,
    token: Option<&CancelToken>,
    progress: Option<&ProgressReporter>,
) -> Option<ManifoldImpl> {
    try_reduce(solid, tool, false, token, progress)
}

/// The erosion `solid ⊖ tool` for a convex tool: the per-triangle hulls
/// reduced by the tree and subtracted from the solid once. `None` declines.
pub fn try_compute_erosion(
    solid: &ManifoldImpl,
    tool: &ManifoldImpl,
    token: Option<&CancelToken>,
    progress: Option<&ProgressReporter>,
) -> Option<ManifoldImpl> {
    try_reduce(solid, tool, true, token, progress)
}

fn try_reduce(
    solid_in: &ManifoldImpl,
    tool: &ManifoldImpl,
    inset: bool,
    token: Option<&CancelToken>,
    progress: Option<&ProgressReporter>,
) -> Option<ManifoldImpl> {
    if is_cancelled(token) {
        return Some(cancelled_impl());
    }

    if solid_in.is_empty()
        || tool.is_empty()
        || solid_in.is_soup
        || tool.is_soup
        || solid_in.status != Error::NoError
        || tool.status != Error::NoError
    {
        return None;
    }

    // minkowski.rs's middle branch without the operand swap.
    if (!inset && solid_in.is_convex()) || !tool.is_convex() {
        return None;
    }

    // Nested or crossing shells: rebuild as the union of shells first, with no
    // reporter (its phases each end on 1.0, which would read as finished).
    let rebuilt;
    let solid: &ManifoldImpl = if has_overlapping_components(solid_in) {
        REBUILDS_RUN.with(|c| c.set(c.get() + 1));
        rebuilt = crate::robust::rebuild_with_rule(solid_in, WindingRule::Positive, token, None);
        if is_cancelled(token) {
            return Some(cancelled_impl());
        }

        if rebuilt.is_empty() || rebuilt.is_soup || rebuilt.status != Error::NoError {
            return None;
        }

        &rebuilt
    } else {
        solid_in
    };

    let max_size = if inset {
        1
    } else {
        PATCH_SIZE_OVERRIDE
            .with(|c| c.get())
            .unwrap_or(convex_patches::MAX_PATCH_SIZE)
    };
    let Some(units) = convex_patches::build(solid, max_size, token) else {
        return Some(cancelled_impl());
    };

    let num_units = units.len();
    LAST_HULL_COUNT.with(|c| c.set(num_units));
    let num_hull_leaves = num_units.div_ceil(LEAF_SIZE);
    let solid_leaves = if inset { 0 } else { 1 };
    let num_leaves = num_hull_leaves + solid_leaves;

    // A binary reduction of L leaves performs exactly L - 1 unions. Erosion
    // adds one unit for its closing subtraction. (C# casts `numLeaves - 1` to
    // ulong; for L = 0 that wraps, and so does this.)
    let total: u64 = (num_units as u64)
        .wrapping_add(num_leaves as u64)
        .wrapping_add((num_leaves as u64).wrapping_sub(1))
        .wrapping_add(1)
        .wrapping_add(if inset { 1 } else { 0 });
    begin_phase(progress, Phase::Minkowski, total);

    // Read once so every node runs on the same engine.
    let engine = BooleanConfig::default_engine();

    let mut units_done: u64 = num_units as u64 + num_leaves as u64;
    let mut level: Option<Vec<ManifoldImpl>> =
        maybe_par_map_ct_progress(num_leaves, UNION_PAR_THRESHOLD, token, progress, |leaf| {
            if leaf < solid_leaves {
                return solid.clone();
            }

            let start = (leaf - solid_leaves) * LEAF_SIZE;
            let count = LEAF_SIZE.min(num_units - start);
            leaf_union(solid, tool, &units, start, count, engine, token, progress)
        });

    while let Some(below) = level.take() {
        if below.len() <= 1 {
            level = Some(below);
            break;
        }

        // Per-level gate: Cancelled leaves must not be unioned into a result.
        if is_cancelled(token) {
            level = Some(below);
            break;
        }

        let pairs = below.len() / 2;
        let merged: Option<Vec<ManifoldImpl>> = match progress {
            Some(p) if pairs <= SUB_PROGRESS_MAX_PAIRS => {
                let tracker = NodeProgress::new(p, units_done, pairs);
                let merged = maybe_par_map_ct(pairs, UNION_PAR_THRESHOLD, token, |node| {
                    let sink = |fraction: f64| tracker.stage(node, fraction);
                    let union = union(
                        &below[2 * node],
                        &below[2 * node + 1],
                        OpType::Add,
                        engine,
                        token,
                        Some(&sink),
                    );
                    tracker.complete(node);
                    union
                });
                p.advance(pairs as u64);
                merged
            }
            _ => maybe_par_map_ct_progress(pairs, UNION_PAR_THRESHOLD, token, progress, |node| {
                boolean_dispatch(
                    &below[2 * node],
                    &below[2 * node + 1],
                    OpType::Add,
                    engine,
                    token,
                )
            }),
        };

        units_done += pairs as u64;

        let Some(mut merged) = merged else {
            level = None;
            break;
        };

        if (below.len() & 1) == 1 {
            // The odd node rides up unchanged; no boolean, no unit.
            let mut below = below;
            if let Some(last) = below.pop() {
                merged.push(last);
            }
        }

        level = Some(merged);
    }

    // A cancelled token can never produce a NoError result.
    let Some(mut level) = level else {
        return Some(cancelled_impl());
    };
    if is_cancelled(token) || level.is_empty() {
        return Some(cancelled_impl());
    }

    let mut out_r = level.swap_remove(0);
    if inset {
        // minkowski.rs's closing merge with two operands, on the tree's engine.
        let tracker = progress.map(|p| NodeProgress::new(p, units_done, 1));
        let sink = |fraction: f64| {
            if let Some(t) = &tracker {
                t.stage(0, fraction);
            }
        };
        out_r = union(
            solid,
            &out_r,
            OpType::Subtract,
            engine,
            token,
            if tracker.is_some() { Some(&sink) } else { None },
        );
        if let Some(t) = &tracker {
            t.complete(0);
        }
        if let Some(p) = progress {
            p.advance(1);
        }
        if is_cancelled(token) {
            return Some(cancelled_impl());
        }
    }

    // minkowski.rs's closing as_original: one fresh original ID, normals and
    // coplanar faces set.
    out_r.initialize_original();
    out_r.set_normals_and_coplanar();

    complete_phase(progress);
    Some(out_r)
}

/// One tree node's boolean: the exact boolean with a stage sink when one is
/// given and the engine is exact, otherwise the plain dispatch.
fn union(
    a: &ManifoldImpl,
    b: &ManifoldImpl,
    op: OpType,
    engine: BooleanEngine,
    token: Option<&CancelToken>,
    stage: Option<&(dyn Fn(f64) + Sync)>,
) -> ManifoldImpl {
    if stage.is_some() && engine == BooleanEngine::Exact {
        return boolean_with_token_and_stage(a, b, op, token, stage);
    }

    boolean_dispatch(a, b, op, engine, token)
}

/// Whether any two connected components' bounding boxes overlap (closed), by
/// a sort and sweep on X — the only way shells can nest or cross.
fn has_overlapping_components(solid: &ManifoldImpl) -> bool {
    let num_vert = solid.num_vert();
    let sets = DisjointSets::new(num_vert as u32);
    for halfedge in solid.halfedge.iter() {
        if halfedge.is_forward() {
            sets.unite(halfedge.start_vert as u32, halfedge.end_vert as u32);
        }
    }

    let mut component: Vec<i32> = vec![0; num_vert];
    let num_components = sets.connected_components(&mut component) as usize;
    if num_components <= 1 {
        return false;
    }

    let mut min = vec![Vec3::new(0.0, 0.0, 0.0); num_components];
    let mut max = vec![Vec3::new(0.0, 0.0, 0.0); num_components];
    let mut seen = vec![false; num_components];
    for v in 0..num_vert {
        let c = component[v] as usize;
        let p = solid.vert_pos[v];
        if !seen[c] {
            min[c] = p;
            max[c] = p;
            seen[c] = true;
            continue;
        }

        // C# Math.Min/Max; NaN-free positions make f64::min/max identical.
        min[c] = Vec3::new(min[c].x.min(p.x), min[c].y.min(p.y), min[c].z.min(p.z));
        max[c] = Vec3::new(max[c].x.max(p.x), max[c].y.max(p.y), max[c].z.max(p.z));
    }

    // Stable sort, as LINQ OrderBy is; -0.0 and 0.0 compare equal, as there.
    let mut order: Vec<usize> = (0..num_components).filter(|&c| seen[c]).collect();
    order.sort_by(|&a, &b| {
        min[a]
            .x
            .partial_cmp(&min[b].x)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    for i in 0..order.len() {
        let a = order[i];
        let mut j = i + 1;
        while j < order.len() && min[order[j]].x <= max[a].x {
            let b = order[j];
            if min[b].y <= max[a].y
                && min[a].y <= max[b].y
                && min[b].z <= max[a].z
                && min[a].z <= max[b].z
            {
                return true;
            }
            j += 1;
        }
    }

    false
}

/// One leaf: the hulls of units `[start, start + count)`, built in
/// minkowski.rs's vertex-sum order (a patch lists its distinct vertices in
/// joining order), unioned through `batch_union` on `engine`.
#[allow(clippy::too_many_arguments)]
fn leaf_union(
    solid: &ManifoldImpl,
    tool: &ManifoldImpl,
    units: &[Vec<usize>],
    start: usize,
    count: usize,
    engine: BooleanEngine,
    token: Option<&CancelToken>,
    progress: Option<&ProgressReporter>,
) -> ManifoldImpl {
    let mut children: Vec<CsgLeafNode> = Vec::with_capacity(count);
    for i in 0..count {
        if is_cancelled(token) {
            return cancelled_impl();
        }

        let unit = &units[start + i];
        let mut simple_hull: Vec<Vec3> = Vec::with_capacity(3 * unit.len() * tool.vert_pos.len());
        let mut seen: Option<Vec<i32>> = if unit.len() > 1 {
            Some(Vec::new())
        } else {
            None
        };
        for &tri in unit {
            for k in 0..3 {
                let v = solid.halfedge[tri * 3 + k].start_vert;
                if let Some(s) = seen.as_mut() {
                    if s.contains(&v) {
                        continue;
                    }
                    s.push(v);
                }

                let a_vert = solid.vert_pos[v as usize];
                for &b_vert in tool.vert_pos.iter() {
                    simple_hull.push(a_vert + b_vert);
                }
            }
        }

        let hull = convex_hull(&simple_hull);
        if let Some(p) = progress {
            p.advance(1);
        }
        if count == 1 {
            return hull;
        }

        children.push(CsgLeafNode::new(hull));
    }

    batch_union(&mut children, token, Some(engine), None).get_impl()
}

/// Sub-unit progress for one level of few, long unions: finished nodes plus
/// every running node's fraction, reported as units through `report_units`,
/// monotone under one lock.
struct NodeProgress<'a> {
    reporter: &'a ProgressReporter,
    base_units: f64,
    state: Mutex<NodeState>,
}

struct NodeState {
    running: Vec<f64>,
    completed: usize,
    last_reported: f64,
}

impl<'a> NodeProgress<'a> {
    fn new(reporter: &'a ProgressReporter, base_units: u64, nodes: usize) -> Self {
        Self {
            reporter,
            base_units: base_units as f64,
            state: Mutex::new(NodeState {
                running: vec![0.0; nodes],
                completed: 0,
                last_reported: 0.0,
            }),
        }
    }

    fn stage(&self, node: usize, fraction: f64) {
        let Ok(mut s) = self.state.lock() else {
            return;
        };
        s.running[node] = s.running[node].max(fraction.clamp(0.0, 1.0));
        self.report_locked(&mut s);
    }

    fn complete(&self, node: usize) {
        let Ok(mut s) = self.state.lock() else {
            return;
        };
        s.running[node] = 0.0;
        s.completed += 1;
        self.report_locked(&mut s);
    }

    fn report_locked(&self, s: &mut NodeState) {
        let mut value = s.completed as f64;
        for &fraction in s.running.iter() {
            value += fraction;
        }

        if value > s.last_reported {
            s.last_reported = value;
            self.reporter.report_units(self.base_units + value);
        }
    }
}

#[cfg(test)]
#[path = "convex_dilation_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "convex_dilation_erosion_tests.rs"]
mod erosion_tests;

#[cfg(test)]
#[path = "convex_dilation_patch_tests.rs"]
mod patch_tests;
