// convex_patches.rs — NOT A C++ PORT; a 1:1 mirror of manifold-sharp's
// ConvexPatches.cs (its RUST_DIVERGENCES.md entry 6; here
// docs/CPP_DIVERGENCES.md entry 13). Stage B of convex_dilation.rs's tree:
// adjacent triangles are grouped into patches P so one hull,
// hull(P ⊕ B) = hull(P) ⊕ B, replaces |P| per-triangle hulls. DILATION ONLY.
//
// ── Why a patch hull is sound ───────────────────────────────────────────────
// hull(P) ⊕ B ⊆ solid ⊕ B whenever H = hull(P) ⊆ solid, and it covers every
// t ⊕ B for t in P, so swapping the |P| hulls for one leaves the union the same
// set — up to QuickHull's epsilon (a hull can only shrink, never add material),
// which is why patched and per-triangle trees are compared at 1e-9 in volume.
// A patch is accepted only when exact checks prove H ⊆ solid:
//   (1) Supporting faces: for every triangle t of P, every vertex of P is on or
//       below t's plane (`orient3d`, never Pos) and at least one is strictly
//       below. A flat patch is refused.
//   (2) Guard: every surface triangle T outside P misses int H, shown by one
//       face plane of the computed hull with all of T's corners on or above it
//       (or T's own plane with every P vertex on one closed side). The hull's
//       faces are first checked to be real supporting planes of the patch, so a
//       zero-area face cannot wave a T through. Candidates come from the
//       solid's collider, queried with the patch's box.
//   Then S ∩ int H = ∅, int H is connected and lies just below each t in P, so
//   int H ⊆ solid and H ⊆ solid. Any failure regrows the patch at half size,
//   down to single triangles, so doubt only costs speed.
//
// Not used for erosion: hull(P) ⊕ B also contains int H ⊕ B, material the
// eroded answer keeps.
//
// ── Determinism ─────────────────────────────────────────────────────────────
// Sequential and a function of the mesh alone: seeds in triangle-index order,
// growth breadth-first across paired halfedges in halfedge order, each
// candidate decided by exact predicates. The unit list is in seed order.

use std::cell::Cell;
use std::collections::VecDeque;

use crate::cancel::{is_cancelled, CancelToken};
use crate::impl_mesh::ManifoldImpl;
use crate::linalg::Vec3;
use crate::quickhull::convex_hull;
use crate::robust::exact::filtered::orient3d;
use crate::robust::exact::Sign;
use crate::types::Box as BBox;

/// The most triangles one patch may hold.
pub(crate) const MAX_PATCH_SIZE: usize = 16;

thread_local! {
    /// Test-only, read on the calling thread: skip the guard (2), so a test can
    /// show a fixture the guard is what keeps open.
    pub(crate) static SKIP_GUARD_FOR_TESTS: Cell<bool> = const { Cell::new(false) };
}

/// The solid's triangles as hull units: each entry is a patch (triangles in
/// the order they joined) or a single triangle, in seed-index order. `None`
/// when cancelled. `max_size <= 1` gives one unit per triangle.
pub(crate) fn build(
    solid: &ManifoldImpl,
    max_size: usize,
    token: Option<&CancelToken>,
) -> Option<Vec<Vec<usize>>> {
    let num_tri = solid.num_tri();
    let mut units: Vec<Vec<usize>> = Vec::with_capacity(num_tri);
    if max_size <= 1 || solid.collider.num_leaves() != num_tri {
        for t in 0..num_tri {
            units.push(vec![t]);
        }
        return Some(units);
    }

    let mut taken = vec![false; num_tri];
    let mut patch: Vec<usize> = Vec::with_capacity(max_size);
    let mut patch_verts: Vec<usize> = Vec::with_capacity(3 * max_size);
    let mut frontier: VecDeque<usize> = VecDeque::new();
    for seed in 0..num_tri {
        if taken[seed] {
            continue;
        }

        if (seed & 63) == 0 && is_cancelled(token) {
            return None;
        }

        // A refused patch is released and regrown from the same seed at half
        // the cap, always below the refused patch's size, down to a single
        // triangle. Growth is deterministic, so the smaller patch is a prefix.
        let mut accepted = false;
        let mut cap = max_size;
        while cap > 1 && !accepted {
            grow(
                solid,
                seed,
                cap,
                &mut taken,
                &mut patch,
                &mut patch_verts,
                &mut frontier,
            );
            accepted = patch.len() > 1 && accept(solid, &patch, &patch_verts);
            if !accepted {
                for &t in patch.iter().skip(1) {
                    taken[t] = false;
                }

                if patch.len() <= 1 {
                    break;
                }
            }
            // C#: cap = Math.Min(cap / 2, patch.Count - 1); patch.len() >= 1 here.
            cap = (cap / 2).min(patch.len() - 1);
        }

        units.push(if accepted { patch.clone() } else { vec![seed] });
    }

    Some(units)
}

fn grow(
    solid: &ManifoldImpl,
    seed: usize,
    max_size: usize,
    taken: &mut [bool],
    patch: &mut Vec<usize>,
    patch_verts: &mut Vec<usize>,
    frontier: &mut VecDeque<usize>,
) {
    patch.clear();
    patch_verts.clear();
    frontier.clear();
    patch.push(seed);
    taken[seed] = true;
    add_verts(solid, seed, patch_verts);
    frontier.push_back(seed);
    while patch.len() < max_size {
        let Some(tri) = frontier.pop_front() else {
            break;
        };
        for k in 0..3 {
            if patch.len() >= max_size {
                break;
            }
            let paired = solid.halfedge[3 * tri + k].paired_halfedge;
            if paired < 0 {
                continue;
            }
            let neighbor = paired as usize / 3;
            if taken[neighbor] || !compatible(solid, neighbor, patch, patch_verts) {
                continue;
            }

            patch.push(neighbor);
            taken[neighbor] = true;
            add_verts(solid, neighbor, patch_verts);
            frontier.push_back(neighbor);
        }
    }
}

fn add_verts(solid: &ManifoldImpl, tri: usize, verts: &mut Vec<usize>) {
    for k in 0..3 {
        let v = solid.halfedge[3 * tri + k].start_vert as usize;
        if !verts.contains(&v) {
            verts.push(v);
        }
    }
}

/// Condition (1) kept incrementally: the candidate's corners are on or below
/// every patch plane, and every patch vertex is on or below the candidate's.
fn compatible(
    solid: &ManifoldImpl,
    candidate: usize,
    patch: &[usize],
    patch_verts: &[usize],
) -> bool {
    for &tri in patch {
        for k in 0..3 {
            let p = solid.vert_pos[solid.halfedge[3 * candidate + k].start_vert as usize];
            if side(solid, tri, p) == Sign::Pos {
                return false;
            }
        }
    }

    for &v in patch_verts {
        if side(solid, candidate, solid.vert_pos[v]) == Sign::Pos {
            return false;
        }
    }

    true
}

/// Which side of triangle `tri`'s outward plane a point is on, exactly.
fn side(solid: &ManifoldImpl, tri: usize, p: Vec3) -> Sign {
    orient3d(
        solid.vert_pos[solid.halfedge[3 * tri].start_vert as usize],
        solid.vert_pos[solid.halfedge[3 * tri + 1].start_vert as usize],
        solid.vert_pos[solid.halfedge[3 * tri + 2].start_vert as usize],
        p,
    )
}

/// The strict half of condition (1), then the guard (2).
fn accept(solid: &ManifoldImpl, patch: &[usize], patch_verts: &[usize]) -> bool {
    for &tri in patch {
        let below = patch_verts
            .iter()
            .any(|&v| side(solid, tri, solid.vert_pos[v]) == Sign::Neg);
        if !below {
            return false;
        }
    }

    let mut points: Vec<Vec3> = Vec::with_capacity(patch_verts.len());
    let mut bx = BBox::new();
    for &v in patch_verts {
        points.push(solid.vert_pos[v]);
        bx.union_point(solid.vert_pos[v]);
    }

    let hull = convex_hull(&points);
    let num_faces = hull.num_tri();
    if hull.is_empty() || num_faces < 4 {
        return false;
    }

    let mut face_corners: Vec<Vec3> = Vec::with_capacity(3 * num_faces);
    for f in 0..num_faces {
        for k in 0..3 {
            face_corners.push(hull.vert_pos[hull.halfedge[3 * f + k].start_vert as usize]);
        }
    }

    if !faces_support(&face_corners, num_faces, &points) {
        return false;
    }

    if SKIP_GUARD_FOR_TESTS.with(|c| c.get()) {
        return true;
    }

    let mut clear = true;
    solid.collider.collisions_one(&bx, 0, |_, tri| {
        if !clear || patch.contains(&tri) {
            return;
        }

        if !separated(solid, tri, &face_corners, num_faces, &points) {
            clear = false;
        }
    });

    clear
}

/// Whether every hull face is an exact supporting plane of the patch: every
/// patch vertex on or below it and at least one strictly below (refusing a
/// zero-area face, whose orient3d is zero for every point).
pub(crate) fn faces_support(face_corners: &[Vec3], num_faces: usize, points: &[Vec3]) -> bool {
    for f in 0..num_faces {
        let mut below = false;
        for &p in points {
            let s = orient3d(
                face_corners[3 * f],
                face_corners[3 * f + 1],
                face_corners[3 * f + 2],
                p,
            );
            if s == Sign::Pos {
                return false;
            }

            below |= s == Sign::Neg;
        }

        if !below {
            return false;
        }
    }

    true
}

/// Whether a plane shows `tri` missing the hull's interior: the triangle's own
/// plane with every patch vertex on one closed side, or one hull face with all
/// of the triangle's corners on or above it. Exact and sufficient only.
fn separated(
    solid: &ManifoldImpl,
    tri: usize,
    face_corners: &[Vec3],
    num_faces: usize,
    points: &[Vec3],
) -> bool {
    let a = solid.vert_pos[solid.halfedge[3 * tri].start_vert as usize];
    let b = solid.vert_pos[solid.halfedge[3 * tri + 1].start_vert as usize];
    let c = solid.vert_pos[solid.halfedge[3 * tri + 2].start_vert as usize];

    // A degenerate triangle puts every point at zero and falls through.
    let mut any_pos = false;
    let mut any_neg = false;
    for &p in points {
        let s = orient3d(a, b, c, p);
        any_pos |= s == Sign::Pos;
        any_neg |= s == Sign::Neg;
    }

    if any_pos != any_neg {
        return true;
    }
    for f in 0..num_faces {
        let f0 = face_corners[3 * f];
        let f1 = face_corners[3 * f + 1];
        let f2 = face_corners[3 * f + 2];
        if orient3d(f0, f1, f2, a) != Sign::Neg
            && orient3d(f0, f1, f2, b) != Sign::Neg
            && orient3d(f0, f1, f2, c) != Sign::Neg
        {
            return true;
        }
    }

    false
}
