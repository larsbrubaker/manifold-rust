// par.rs — determinism-preserving parallel execution helpers
//
// Mirrors the C++ `autoPolicy` pattern (parallel.h): each site opts into
// parallelism above a size threshold, staying sequential for small inputs.
// Unlike upstream MANIFOLD_PAR (which allows nondeterministic vertex order in
// some phases), only sites whose output is provably identical to the
// sequential build are parallelized — per-index maps with indexed writes, and
// collect-then-sort pipelines whose final sort is a total order. This keeps
// the `parallel` feature bit-exact with the sequential reference.
//
// The sites, with the size at which each goes parallel (manifold-sharp's
// Par.cs mirrors this list; keep it current):
//   face_op.rs            `calculate_vert_normals`, per-vertex map    10,000 verts
//   face_op_triangulate.rs `face2tri_ct` triangulation, map_ct       512 faces
//                         `face2tri` writes, runs of 4,096 faces into
//                         disjoint `split_at_mut` slices             2 runs
//   boolean3_kernels.rs   `intersect12` query map_ct; result stable
//                         sort and gathers                           10,000
//                         `winding03` unbroken-edge search, chunks of
//                         1,024 halfedges (unions stay sequential)   10 chunks
//                         `winding03` per-vert winding query map_ct  1,000 verts
//   boolean_result.rs     `EdgeGroups::new` stable sort of the
//                         AddNewEdgeVerts edge lists                 10,000 entries
//   sort.rs               `sort_geometry`: Morton codes, stable sorts,
//                         vertex/face/halfedge/tangent gathers        100,000
//                         (PAR_THRESHOLD; per-call n is verts, tris
//                         or 3 * tris)
//   edge_op.rs            edge-flag scans of `collapse_short_edges`,
//                         `collapse_colinear_edges`, `swap_degenerates`
//                         (filter; collapses/swaps stay sequential)  100,001 halfedges
//                                                                    (FLAG_PAR_THRESHOLD)
//   edge_op_orbits.rs     `orbit_owners` for `split_pinched_verts` and
//                         `dedupe_edges` (walks capped at 64 steps,
//                         then one sequential pass)                  100,001 halfedges
//                                                                    (ORBIT_PAR_THRESHOLD)
//                         `dedupe_edges` per-owner orbit map         12,500 owners
//   sdf.rs                level-set voxel evaluation                 10,000 voxels
//   minkowski.rs          per-face hull maps                         8
//   robust/intersection_graph.rs (via progress.rs) map_ct stages     64 / 16
//   robust/soup.rs        `compute_self_intersections` rows, any_ct 1,000 triangles
//
// One site is an existential "any" rather than a map: `maybe_par_any_ct`, the
// Auto engine's self-intersection pre-check
// (robust/soup.rs `compute_self_intersections`, rows of 1,000+ triangles).
// It returns no array but a boolean identical to the sequential loop's,
// because "some index hits" does not depend on which index is found first or
// how many ran (see its doc). Shared with manifold-sharp's
// `Par.MaybeParAnyCt` (its commit 1d9162c).

/// Map `f` over `0..n`, in parallel when the `parallel` feature is enabled and
/// `n >= threshold`. Results are returned in index order either way.
#[cfg(feature = "parallel")]
pub fn maybe_par_map<T, F>(n: usize, threshold: usize, f: F) -> Vec<T>
where
    T: Send,
    F: Fn(usize) -> T + Sync + Send,
{
    use rayon::prelude::*;
    if n >= threshold {
        (0..n).into_par_iter().map(f).collect()
    } else {
        (0..n).map(f).collect()
    }
}

/// Sequential fallback: identical output to the parallel version.
#[cfg(not(feature = "parallel"))]
pub fn maybe_par_map<T, F>(n: usize, _threshold: usize, f: F) -> Vec<T>
where
    F: Fn(usize) -> T,
{
    (0..n).map(f).collect()
}

/// The indices in `0..n` where `pred` holds, ascending, testing in parallel
/// when `n >= threshold`: the flag half of C++ `FlagStore::run_par`
/// (edge_op.cpp:54-84). rayon's `collect` keeps the order, so unlike C++ no
/// sort is needed.
#[cfg(feature = "parallel")]
pub fn maybe_par_filter<F>(n: usize, threshold: usize, pred: F) -> Vec<usize>
where
    F: Fn(usize) -> bool + Sync + Send,
{
    use rayon::prelude::*;
    if n >= threshold {
        (0..n).into_par_iter().filter(|&i| pred(i)).collect()
    } else {
        (0..n).filter(|&i| pred(i)).collect()
    }
}

/// Sequential fallback: identical output to the parallel version.
#[cfg(not(feature = "parallel"))]
pub fn maybe_par_filter<F>(n: usize, _threshold: usize, pred: F) -> Vec<usize>
where
    F: Fn(usize) -> bool,
{
    (0..n).filter(|&i| pred(i)).collect()
}

/// A stable sort by key, in parallel when `v.len() >= threshold`. A stable
/// sort is determined by the keys and input order, so rayon's gives the same
/// slice, provided `key` is a total order.
#[cfg(feature = "parallel")]
pub fn maybe_par_sort_by_key<T, K, F>(v: &mut [T], threshold: usize, key: F)
where
    T: Send,
    K: Ord,
    F: Fn(&T) -> K + Sync,
{
    use rayon::prelude::*;
    if v.len() >= threshold {
        v.par_sort_by_key(key);
    } else {
        v.sort_by_key(key);
    }
}

/// Sequential fallback: identical output to the parallel version.
#[cfg(not(feature = "parallel"))]
pub fn maybe_par_sort_by_key<T, K, F>(v: &mut [T], _threshold: usize, key: F)
where
    K: Ord,
    F: Fn(&T) -> K,
{
    v.sort_by_key(key);
}

/// Run `f` on every item, in parallel when `items.len() >= threshold`. Items
/// must own disjoint output, so the run order cannot change what is written.
#[cfg(feature = "parallel")]
pub fn maybe_par_for_each<T, F>(items: Vec<T>, threshold: usize, f: F)
where
    T: Send,
    F: Fn(T) + Sync + Send,
{
    use rayon::prelude::*;
    if items.len() >= threshold {
        items.into_par_iter().for_each(f);
    } else {
        items.into_iter().for_each(f);
    }
}

/// Sequential fallback: identical output to the parallel version.
#[cfg(not(feature = "parallel"))]
pub fn maybe_par_for_each<T, F>(items: Vec<T>, _threshold: usize, f: F)
where
    F: Fn(T),
{
    items.into_iter().for_each(f);
}

/// Apply `f` to every element in place, in parallel when `v.len() >=
/// threshold`. `f` sees only its own element, so the result does not depend on
/// the order.
#[cfg(feature = "parallel")]
pub fn maybe_par_for_each_mut<T, F>(v: &mut [T], threshold: usize, f: F)
where
    T: Send,
    F: Fn(&mut T) + Sync + Send,
{
    use rayon::prelude::*;
    if v.len() >= threshold {
        v.par_iter_mut().for_each(f);
    } else {
        v.iter_mut().for_each(f);
    }
}

/// Sequential fallback: identical output to the parallel version.
#[cfg(not(feature = "parallel"))]
pub fn maybe_par_for_each_mut<T, F>(v: &mut [T], _threshold: usize, f: F)
where
    F: Fn(&mut T),
{
    v.iter_mut().for_each(f);
}

/// [`maybe_par_map`] with cooperative cancellation: `None` means the token was
/// cancelled and the (necessarily incomplete) results were discarded.
///
/// Mirrors the ctx-aware `for_each` overload in C++ `parallel.h:400-437`, whose
/// contract is the same — "only safe when *skip the rest of the range* produces
/// a result the caller will discard via a post-loop `IsCancelled` check".
///
/// Two deliberate differences from the C++:
/// - A `None` token dispatches straight to [`maybe_par_map`], so the
///   uncancellable path is byte-for-byte the code that ran before cancellation
///   existed — stricter than C++, which still branches on `ctx != nullptr`
///   inside the loop.
/// - With a token we check the flag on *every* element rather than once per
///   `kSeqCancelChunk` (1024) elements. The flag is written at most once, so its
///   cache line stays shared and the relaxed load is an L1 hit; paying it per
///   element buys strictly better cancellation latency than the C++ chunking.
#[cfg(feature = "parallel")]
pub fn maybe_par_map_ct<T, F>(
    n: usize,
    threshold: usize,
    token: Option<&crate::cancel::CancelToken>,
    f: F,
) -> Option<Vec<T>>
where
    T: Send,
    F: Fn(usize) -> T + Sync + Send,
{
    use rayon::prelude::*;
    let Some(token) = token else {
        return Some(maybe_par_map(n, threshold, f));
    };
    // Collecting `Option<T>` into `Option<Vec<T>>` short-circuits on the first
    // `None` in both rayon and std, so a cancel stops the remaining work
    // instead of merely skipping each element's body.
    if n >= threshold {
        (0..n)
            .into_par_iter()
            .map(|i| {
                if token.is_cancelled() {
                    None
                } else {
                    Some(f(i))
                }
            })
            .collect()
    } else {
        (0..n)
            .map(|i| {
                if token.is_cancelled() {
                    None
                } else {
                    Some(f(i))
                }
            })
            .collect()
    }
}

/// Sequential fallback: identical output to the parallel version.
#[cfg(not(feature = "parallel"))]
pub fn maybe_par_map_ct<T, F>(
    n: usize,
    threshold: usize,
    token: Option<&crate::cancel::CancelToken>,
    f: F,
) -> Option<Vec<T>>
where
    F: Fn(usize) -> T,
{
    let Some(token) = token else {
        return Some(maybe_par_map(n, threshold, f));
    };
    (0..n)
        .map(|i| {
            if token.is_cancelled() {
                None
            } else {
                Some(f(i))
            }
        })
        .collect()
}

/// Does `predicate` hold for some index in `0..n`? Stops at the first hit, in
/// parallel when the `parallel` feature is enabled and `n >= threshold`.
/// `None` means the token was cancelled before a verdict.
///
/// The one existential site (the Auto engine's self-intersection pre-check,
/// `robust::soup::compute_self_intersections`). Its answer is a single
/// boolean rather than an array, so the index-order argument the maps above
/// rest on becomes this one: the answer is "some index's predicate is true",
/// and an existential does not depend on which index is found first or how
/// many ran. A hit, from any worker, answers `Some(true)`; a run in which
/// every index ran without a hit answers `Some(false)`. Both are exactly what
/// the sequential loop answers, so the verdict is the same in both builds;
/// only *which* hit stops the loop differs, and the caller never learns which.
///
/// `predicate` must be pure apart from its per-job scratch `L`, made by
/// `local_init` (rayon's `*_init`: one per split job in parallel, one for the
/// whole sequential loop).
///
/// Cancellation: a cancel any worker observed answers `None`, unless some
/// worker had already found a hit, which answers `Some(true)` — a real
/// verdict (the hit is a genuine witness), not a cancelled partial result.
/// The sequential loop polls the token before every index.
#[cfg(feature = "parallel")]
pub fn maybe_par_any_ct<L, I, F>(
    n: usize,
    threshold: usize,
    token: Option<&crate::cancel::CancelToken>,
    local_init: I,
    predicate: F,
) -> Option<bool>
where
    I: Fn() -> L + Sync + Send,
    F: Fn(usize, &mut L) -> bool + Sync + Send,
{
    use rayon::prelude::*;
    use std::sync::atomic::{AtomicBool, Ordering};
    if n < threshold {
        return seq_any_ct(n, token, local_init, predicate);
    }
    // Written only ever from false to true, so two racing hits are benign.
    let found = AtomicBool::new(false);
    let stopped = (0..n)
        .into_par_iter()
        .try_for_each_init(local_init, |local, i| {
            if token.is_some_and(|t| t.is_cancelled()) {
                return Err(());
            }
            if predicate(i, local) {
                found.store(true, Ordering::Relaxed);
                return Err(());
            }
            Ok(())
        })
        .is_err();
    if found.load(Ordering::Relaxed) {
        Some(true)
    } else if stopped {
        // Stopped without a hit means some worker stopped on the cancel.
        None
    } else {
        Some(false)
    }
}

/// Sequential fallback: identical verdict to the parallel version.
#[cfg(not(feature = "parallel"))]
pub fn maybe_par_any_ct<L, I, F>(
    n: usize,
    _threshold: usize,
    token: Option<&crate::cancel::CancelToken>,
    local_init: I,
    predicate: F,
) -> Option<bool>
where
    I: Fn() -> L,
    F: Fn(usize, &mut L) -> bool,
{
    seq_any_ct(n, token, local_init, predicate)
}

/// The sequential loop both builds' [`maybe_par_any_ct`] reduce to.
fn seq_any_ct<L, I, F>(
    n: usize,
    token: Option<&crate::cancel::CancelToken>,
    local_init: I,
    predicate: F,
) -> Option<bool>
where
    I: Fn() -> L,
    F: Fn(usize, &mut L) -> bool,
{
    let mut local = local_init();
    for i in 0..n {
        if token.is_some_and(|t| t.is_cancelled()) {
            return None;
        }
        if predicate(i, &mut local) {
            return Some(true);
        }
    }
    Some(false)
}

#[cfg(test)]
#[path = "par_tests.rs"]
mod tests;
