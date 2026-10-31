// robust/soup.rs — Triangle-soup import support for the robust boolean
// engine.
//
// When the strict halfedge pairing of `Manifold::from_mesh_gl` fails
// (non-manifold connectivity), `from_mesh_gl_robust` falls through to
// `soupify`: the geometry is kept as an unpaired-halfedge triangle soup
// inside `ManifoldImpl` (is_soup = true), provided it passes the one check
// the robust engine genuinely needs — the soup must be geometrically
// **closed and orientable**: after welding vertices by exact position and
// dropping exactly-degenerate triangles (paper §5), every directed edge
// must be balanced by its reverse. That is precisely the condition for the
// soup to bound a solid via winding numbers.

use std::collections::BTreeMap;

use crate::impl_mesh::ManifoldImpl;
use crate::linalg::{IVec3, Vec3};
use crate::types::{Error, Halfedge};

use super::exact::rational::R3;

/// Weld key: exact position identity (with -0.0 normalized so it equals 0.0).
fn pos_key(v: Vec3) -> (u64, u64, u64) {
    let norm = |x: f64| if x == 0.0 { 0.0f64 } else { x }.to_bits();
    (norm(v.x), norm(v.y), norm(v.z))
}

/// Exact zero-area test on f64 positions.
fn is_degenerate(a: Vec3, b: Vec3, c: Vec3) -> bool {
    R3::from_vec3(b)
        .sub(&R3::from_vec3(a))
        .cross(&R3::from_vec3(c).sub(&R3::from_vec3(a)))
        .is_zero()
}

/// Convert `imp` (with `vert_pos` and `mesh_relation.tri_ref` already
/// populated, and `halfedge` full of whatever strict pairing produced) into
/// a validated triangle soup:
///  - drop exactly-degenerate triangles (and their tri_refs),
///  - verify the remainder is closed and orientable (else `Error::NotClosed`),
///  - rebuild halfedges with best-effort pairing (`paired_halfedge == -1`
///    where no partner exists),
///  - recompute per-face normals (no pairing required),
///  - set `is_soup`.
///
/// `tri_prop` / `tri_vert` mirror the `create_halfedges` inputs: `tri_vert`
/// holds the position indices when properties are mapped separately, else
/// `tri_prop` does double duty.
pub fn soupify(
    imp: &mut ManifoldImpl,
    tri_prop: &[IVec3],
    tri_vert: &[IVec3],
) -> Result<(), Error> {
    let position_tris: &[IVec3] = if tri_vert.is_empty() {
        tri_prop
    } else {
        tri_vert
    };
    debug_assert_eq!(position_tris.len(), tri_prop.len());

    // Weld vertex ids by exact position for the closedness bookkeeping only
    // (vert_pos itself is left as imported).
    let mut weld: BTreeMap<(u64, u64, u64), i32> = BTreeMap::new();
    let mut welded_id = vec![0i32; imp.vert_pos.len()];
    for (i, &p) in imp.vert_pos.iter().enumerate() {
        let id = *weld.entry(pos_key(p)).or_insert(i as i32);
        welded_id[i] = id;
    }

    // Keep non-degenerate triangles; balance directed edges on welded ids.
    let mut keep: Vec<usize> = Vec::with_capacity(position_tris.len());
    let mut balance: BTreeMap<(i32, i32), i64> = BTreeMap::new();
    for (t, tv) in position_tris.iter().enumerate() {
        let (a, b, c) = (tv.x as usize, tv.y as usize, tv.z as usize);
        let (wa, wb, wc) = (welded_id[a], welded_id[b], welded_id[c]);
        if wa == wb
            || wb == wc
            || wc == wa
            || is_degenerate(imp.vert_pos[a], imp.vert_pos[b], imp.vert_pos[c])
        {
            continue; // paper §5: degenerate triangles are safe to drop
        }
        keep.push(t);
        for (u, v) in [(wa, wb), (wb, wc), (wc, wa)] {
            let key = (u.min(v), u.max(v));
            *balance.entry(key).or_insert(0) += if u < v { 1 } else { -1 };
        }
    }
    if keep.len() < 4 {
        return Err(Error::NotClosed);
    }
    if balance.values().any(|&n| n != 0) {
        return Err(Error::NotClosed);
    }

    // Rebuild halfedges for the kept triangles with best-effort pairing:
    // LIFO multimap on welded undirected edges, forward (start < end) pairs
    // with reverse.
    let has_props = !tri_vert.is_empty();
    let mut halfedges: Vec<Halfedge> = Vec::with_capacity(3 * keep.len());
    let mut tri_ref = Vec::with_capacity(keep.len());
    for (new_t, &old_t) in keep.iter().enumerate() {
        let tv = position_tris[old_t];
        let tp = tri_prop[old_t];
        for i in 0..3 {
            let j = (i + 1) % 3;
            halfedges.push(Halfedge {
                start_vert: tv[i],
                end_vert: tv[j],
                paired_halfedge: -1,
                prop_vert: if has_props { tp[i] } else { tv[i] },
            });
        }
        if old_t < imp.mesh_relation.tri_ref.len() {
            tri_ref.push(imp.mesh_relation.tri_ref[old_t]);
        }
        let _ = new_t;
    }
    let mut open: BTreeMap<(i32, i32), Vec<usize>> = BTreeMap::new();
    for (idx, he) in halfedges.iter().enumerate() {
        let (u, v) = (
            welded_id[he.start_vert as usize],
            welded_id[he.end_vert as usize],
        );
        open.entry((u.min(v), u.max(v))).or_default().push(idx);
    }
    for (_key, mut idxs) in open {
        // Pair forwards with reverses greedily; leftovers stay -1.
        let mut fwd: Vec<usize> = Vec::new();
        let mut bwd: Vec<usize> = Vec::new();
        for idx in idxs.drain(..) {
            let he = &halfedges[idx];
            if welded_id[he.start_vert as usize] < welded_id[he.end_vert as usize] {
                fwd.push(idx);
            } else {
                bwd.push(idx);
            }
        }
        while let (Some(f), Some(b)) = (fwd.pop(), bwd.pop()) {
            halfedges[f].paired_halfedge = b as i32;
            halfedges[b].paired_halfedge = f as i32;
        }
    }
    imp.halfedge = halfedges;
    if tri_ref.len() == keep.len() {
        imp.mesh_relation.tri_ref = tri_ref;
    } else {
        imp.mesh_relation.tri_ref.clear();
    }

    // Per-face normals need only the triangle itself.
    imp.face_normal = (0..imp.num_tri())
        .map(|t| {
            let a = imp.vert_pos[imp.halfedge[3 * t].start_vert as usize];
            let b = imp.vert_pos[imp.halfedge[3 * t + 1].start_vert as usize];
            let c = imp.vert_pos[imp.halfedge[3 * t + 2].start_vert as usize];
            let n = crate::linalg::cross(b - a, c - a);
            let len = crate::linalg::length(n);
            if len > 0.0 {
                n / len
            } else {
                Vec3::new(0.0, 0.0, 0.0)
            }
        })
        .collect();
    imp.vert_normal.clear();
    imp.is_soup = true;
    Ok(())
}

// ---------------------------------------------------------------------------
// Geometric self-intersection detection
// ---------------------------------------------------------------------------

/// Lazily-resolved "does this impl self-intersect" verdict, stored on
/// [`ManifoldImpl`]. `OnceLock` because the robust engine queries impls from
/// rayon workers under the `parallel` feature, so the cell must be
/// thread-safe. `Clone` carries the settled value across, which is why every
/// operation that clones an impl and then edits its geometry must call
/// `ManifoldImpl::invalidate_self_intersects`.
#[derive(Debug, Default)]
pub struct SelfIntersectCache(std::sync::OnceLock<bool>);

impl Clone for SelfIntersectCache {
    fn clone(&self) -> Self {
        let out = std::sync::OnceLock::new();
        if let Some(&v) = self.0.get() {
            let _ = out.set(v);
        }
        SelfIntersectCache(out)
    }
}

impl SelfIntersectCache {
    /// The settled verdict, if the detector has already run.
    pub fn get(&self) -> Option<bool> {
        self.0.get().copied()
    }

    /// Seed an already-known verdict (used when a transform carries the
    /// answer forward). No-op once the cell is settled.
    pub fn set(&self, value: bool) {
        let _ = self.0.set(value);
    }
}

/// True when two of `imp`'s own triangles genuinely intersect — they cross,
/// they overlap, or they are coincident surface — as opposed to merely
/// sharing an edge or a vertex as every closed mesh does.
///
/// Answers from the cache after the first call. The narrow phase is
/// `intersection_graph::real_self_contact`, the same predicate the robust
/// engine uses to decide which self-cuts it must make, plus the exact
/// duplicate-triangle case that predicate deliberately passes over (see
/// [`genuine_contact`]). Unlike the engine's phase 2b this stops at the
/// first genuine contact and builds no graph.
pub fn has_self_intersections(imp: &ManifoldImpl) -> bool {
    has_self_intersections_with_token(imp, None)
}

/// [`has_self_intersections`] with cancellation, for the boolean dispatcher.
///
/// A cancelled scan answers `true` (route to the robust engine, which then
/// reports `Error::Cancelled` itself) and caches nothing, so the verdict is
/// recomputed properly if the impl is used again.
pub fn has_self_intersections_with_token(
    imp: &ManifoldImpl,
    token: Option<&crate::cancel::CancelToken>,
) -> bool {
    if let Some(v) = imp.self_intersects.get() {
        return v;
    }
    match compute_self_intersections(imp, token) {
        Some(verdict) => {
            imp.self_intersects.set(verdict);
            verdict
        }
        None => true,
    }
}

/// Whether two triangles that share no vertex are proven apart before the
/// full test: `t2` strictly on one side of `t1`'s plane by the float filter
/// alone (`filtered::orient3d_filter`, the very filter
/// `real_self_contact`'s plane signs start with, so the two always agree
/// wherever it certifies a sign, including at magnitudes where its own
/// bound has known gaps (see filtered.rs); no exact fallback), or else, when the filter cannot settle that
/// or finds a corner on the plane (the pair lies in or near one plane),
/// their projections strictly
/// separated by an edge line (tri_tri.rs `coplanar_separated_2d`, the same
/// predicate and axis `real_self_contact`'s coplanar path uses).
///
/// Sound because two triangles that meet share a point, and so does any
/// projection of them; either proof leaves no contact, so
/// `real_self_contact` would answer `None` and no verdict moves. The order
/// keeps ordinary pairs cheap (three filtered signs, the test
/// `real_self_contact` would start with anyway) and spends the projection
/// test only where it pays: on exactly coplanar pairs every plane sign is
/// zero, and in `real_self_contact` each would escalate to the exact
/// `orient3d_i` before the coplanar path rejected the pair. A finely
/// revolved part's flat annuli are nearly all such pairs.
fn vertex_disjoint_apart(t1: [Vec3; 3], t2: [Vec3; 3]) -> bool {
    use super::exact::filtered::orient3d_filter;
    use super::exact::Sign;
    let side = |v: Vec3| orient3d_filter(t1[0], t1[1], t1[2], v);
    let signs = [side(t2[0]), side(t2[1]), side(t2[2])];
    if let [Some(a), Some(b), Some(c)] = signs {
        if a != Sign::Zero && a == b && b == c {
            return true;
        }
    }
    if signs.iter().all(|s| s.is_some_and(|s| s != Sign::Zero)) {
        // Decided, off the plane and not one-sided: `t2` crosses the
        // plane, so the pair is not near-coplanar and the projection test
        // is unlikely to help.
        return false;
    }
    super::tri_tri::coplanar_separated_2d(t1, t2, super::tri_tri::f64_drop_axis(t1))
}

/// Do these two triangles of one mesh meet in anything beyond ordinary
/// adjacency?
///
/// `real_self_contact` answers that for every case but one: it reports
/// exactly duplicated triangles (all three vertices coincide, either
/// winding) as benign, because the robust arrangement needs no cut there —
/// both copies emit identical pieces and the winding arithmetic resolves
/// them. They are still coincident surface, which is precisely what the
/// exact engine cannot integrate (Thingi10K #92068's shells are triple-wound
/// duplicates and nothing else), so the dispatch detector counts them.
///
/// Before that, a pair sharing exactly one vertex is benign when `t1`'s two
/// other corners lie strictly on one side of `t2`'s plane: then every point
/// of t1 but the shared vertex is strictly off that plane, so the pair meets
/// in that vertex alone, which `real_self_contact` also calls benign (by its
/// own shortcut or, after a full tri-tri test, as a point contact). Its
/// shortcut only tests t2's corners against t1's plane; this is the mirror,
/// kept here so `real_self_contact` itself is untouched, and it changes no
/// verdict. (manifold-sharp, its commit 1d9162c: on a 21.6k-triangle part it
/// settles about 171k of the 209k vertex-neighbour pairs that would
/// otherwise pay for the full test.)
///
/// First of all, a pair sharing no vertex that `vertex_disjoint_apart`
/// proves apart is benign: no point in common means no contact of any kind.
pub(super) fn genuine_contact(
    t1: [Vec3; 3],
    t2: [Vec3; 3],
    stats: &mut super::intersection_graph::SelfCutStats,
) -> bool {
    use super::exact::Sign;
    // Shared vertices counted as real_self_contact counts them: t1's corners
    // found in t2 by IEEE equality.
    let n_shared = t1.iter().filter(|v| t2.contains(v)).count();
    if n_shared == 3 {
        return true;
    }
    if n_shared == 0 && vertex_disjoint_apart(t1, t2) {
        stats.apart += 1;
        return false;
    }
    if n_shared == 1 {
        let mut first = Sign::Zero;
        let mut one_sided = true;
        for &v in t1.iter().filter(|v| !t2.contains(v)) {
            let s = super::intersection_graph::orient3d_plane(&t2, v);
            if s == Sign::Zero || (first != Sign::Zero && s != first) {
                one_sided = false;
                break;
            }
            first = s;
        }
        if one_sided {
            stats.vert_benign += 1;
            return false;
        }
    }
    super::intersection_graph::real_self_contact(t1, t2, stats).is_some()
}

/// Triangle count at or above which the parallel build scans in parallel.
pub(crate) const SELF_INTERSECT_PAR_THRESHOLD: usize = 1_000;

/// Uncached detector: BVH broad phase over the impl's own triangles, exact
/// narrow phase, early exit on the first genuine contact. `None` means the
/// scan was cancelled before it could reach a verdict.
///
/// The broad phase reuses `imp.collider` — the face BVH `sort_geometry`
/// already built, whose leaves are in face order — and only builds a private
/// morton-ordered tree (as `intersection_graph::build_graph`'s self-cut
/// phase does) when the impl carries no matching collider, which is the case
/// for soup impls.
fn compute_self_intersections(
    imp: &ManifoldImpl,
    token: Option<&crate::cancel::CancelToken>,
) -> Option<bool> {
    scan_triangles(&impl_to_tris(imp), Some(&imp.collider), token)
}

/// [`compute_self_intersections`] over plain triangles. `face_collider` is a
/// tree whose leaves are the triangles in order, reused when its leaf count
/// matches; any other tree (or `None`) builds a private one.
pub(super) fn scan_triangles(
    tris: &[[Vec3; 3]],
    face_collider: Option<&crate::collider::Collider>,
    token: Option<&crate::cancel::CancelToken>,
) -> Option<bool> {
    use super::intersection_graph::{is_degenerate as is_degenerate_tri, tri_box, SelfCutStats};
    use crate::types::Box;

    if tris.len() < 2 {
        return Some(false);
    }
    // Non-finite positions (a warp to NaN/inf, which no import check
    // rejects) have no exact rational form, so the narrow phase cannot run
    // on them. "Self-intersecting" is the safe verdict: it routes the
    // operand to the robust engine rather than letting the exact kernels
    // integrate garbage.
    if tris
        .iter()
        .flatten()
        .any(|v| !v.x.is_finite() || !v.y.is_finite() || !v.z.is_finite())
    {
        return Some(true);
    }

    let boxes: Vec<Box> = tris.iter().map(tri_box).collect();
    let live: Vec<bool> = tris.iter().map(|t| !is_degenerate_tri(t)).collect();

    // Flat fans settle their own pairs in one sort (robust/fan_star.rs); no
    // verdict moves, only which pairs are listed.
    let stars = super::fan_star::flat_stars(tris, &live);
    if !stars.owned.is_empty() {
        if token.is_some_and(|t| t.is_cancelled()) {
            return None;
        }
        return super::fan_star_scan::scan_with_flat_stars(tris, &boxes, &live, &stars, token);
    }

    // Leaf index -> triangle index; empty when the cached face collider is
    // used, whose leaves already are triangle indices.
    let mut leaf_tri: Vec<usize> = Vec::new();
    let owned;
    let collider = if let Some(face) = face_collider.filter(|c| c.num_leaves() == tris.len()) {
        face
    } else {
        let mut order: Vec<usize> = (0..tris.len()).filter(|&i| live[i]).collect();
        if order.len() < 2 {
            return Some(false);
        }
        let scene = boxes
            .iter()
            .enumerate()
            .filter(|(i, _)| live[*i])
            .fold(Box::new(), |acc, (_, b)| acc.union_box(b));
        order.sort_by_key(|&i| crate::sort::morton_code(boxes[i].center(), &scene));
        owned = crate::collider::Collider::new(
            order.iter().map(|&i| boxes[i]).collect(),
            order
                .iter()
                .map(|&i| crate::sort::morton_code(boxes[i].center(), &scene))
                .collect(),
        );
        leaf_tri = order;
        &owned
    };
    let mapped = !leaf_tri.is_empty();

    // One row per triangle i: does it genuinely touch some later triangle j?
    // The verdict is "some row hits", which `maybe_par_any_ct` answers the
    // same in both builds (par.rs). Each job gets its own candidate list and
    // stats; the stats are discarded, and the rows read only the shared,
    // read-only tris/boxes/live and the collider, which the robust engine's
    // self-cut map already queries concurrently. Cancellation keeps its
    // contract: `None` when the scan stopped on the token before finding a
    // contact.
    crate::par::maybe_par_any_ct(
        tris.len(),
        SELF_INTERSECT_PAR_THRESHOLD,
        token,
        || (Vec::<usize>::new(), SelfCutStats::default()),
        |i, (cands, stats)| {
            if !live[i] {
                return false;
            }
            cands.clear();
            collider.collisions_one(&boxes[i], i, |_, leaf| {
                cands.push(if mapped { leaf_tri[leaf] } else { leaf });
            });
            cands.sort_unstable();
            cands.iter().any(|&j| {
                j > i
                    && live[j]
                    && boxes[i].does_overlap_box(&boxes[j])
                    && genuine_contact(tris[i], tris[j], stats)
            })
        },
    )
}

#[cfg(test)]
#[path = "soup_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "self_intersection_scan_tests.rs"]
mod scan_tests;

#[cfg(test)]
#[path = "fan_star_scan_tests.rs"]
mod fan_star_tests;

#[cfg(test)]
#[path = "revolved_ring_scan_tests.rs"]
mod revolved_ring_tests;

/// Per-corner vertex properties of an impl, flattened as
/// `props[(3*tri + corner) * num_prop + channel]`, aligned with
/// `impl_to_tris` ordering. Empty when the impl carries no properties.
pub fn impl_to_corner_props(imp: &ManifoldImpl) -> Vec<f64> {
    let np = imp.num_prop;
    if np == 0 {
        return Vec::new();
    }
    let mut out = Vec::with_capacity(3 * imp.num_tri() * np);
    for t in 0..imp.num_tri() {
        for i in 0..3 {
            let pv = imp.halfedge[3 * t + i].prop_vert as usize;
            out.extend_from_slice(&imp.properties[pv * np..(pv + 1) * np]);
        }
    }
    out
}

/// The triangle list of any impl (soup or manifold) as position triples —
/// the robust engine's working form.
pub fn impl_to_tris(imp: &ManifoldImpl) -> Vec<[Vec3; 3]> {
    (0..imp.num_tri())
        .map(|t| {
            [
                imp.vert_pos[imp.halfedge[3 * t].start_vert as usize],
                imp.vert_pos[imp.halfedge[3 * t + 1].start_vert as usize],
                imp.vert_pos[imp.halfedge[3 * t + 2].start_vert as usize],
            ]
        })
        .collect()
}
