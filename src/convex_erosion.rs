// convex_erosion.rs — the closed-form erosion of a *convex* solid. NOT A C++ PORT:
// the C++ (and `minkowski.rs`) has one erosion algorithm, the per-triangle sweep.
// This is a 1:1 mirror of manifold-sharp's `ConvexErosion.cs` (its
// RUST_DIVERGENCES.md entry 5, now taken here; recorded in this repo's
// docs/CPP_DIVERGENCES.md entry 13). It is an ADDED entry point, reached only
// through `Manifold::try_convex_erosion`; `minkowski::minkowski` is not routed
// through it, so every ported path still produces the C++'s bits. The two ports
// must stay bit-identical to each other: every operation below is written in
// the same order as the C#.
//
// ── The closed form ─────────────────────────────────────────────────────────
// A convex solid IS the intersection of its face halfspaces,
//
//     A = { x : n_i . x <= d_i }   over A's faces i, n_i outward and unit,
//
// and the erosion of a convex A by any B is that intersection with each plane
// pushed inward by B's support in that direction:
//
//     A (-) B = { x : n_i . x <= d_i - h_B(-n_i) },   h_B(u) = max_b u.b
//
// The sign is `h_B(-n_i)` because minkowski.rs's erosion is
// A \ (boundary(A) (+) B) — it sweeps B, not -B — so it computes
// { x : x - B subset A }. For the centred ball every caller uses, B = -B.
//
// ── Why a dual hull rather than a plane-by-plane clip ───────────────────────
// With p strictly inside the eroded body, the dual point q_i = n_i / c_i (c_i
// the plane's slack at p) has the property that a FACET of hull{q_i} names the
// three planes meeting at a VERTEX of the result. The dual hull only enumerates
// those triples; each vertex is solved in the primal from the original n and d,
// which keeps a box's corners exactly on their planes.
//
// ── Declines ────────────────────────────────────────────────────────────────
// A fast path: every case it is not sure of returns `None` and the caller runs
// the general sweep. It declines a non-convex solid or tool, a tool that does
// not contain the origin, a centroid not strictly inside the eroded body, a
// degenerate plane triple, and — the backstop — any output vertex that does not
// satisfy every constraint it was built from.
//
// ── Accuracy ────────────────────────────────────────────────────────────────
// A 20-cube eroded by a unit ball comes out at exactly 5832.0. Agreement with
// the sweep is ~1e-15 relative; on a 2048-triangle sphere the volume agrees to
// 8e-15 but the triangulation differs (4016 vs 4020 triangles), since QuickHull
// drops dual points within its epsilon of a facet (redundant halfspaces).

// `!(x > 0.0)` is deliberate throughout: a NaN length or slack must decline,
// exactly as the C# `!(x > 0.0)` it mirrors does.
#![allow(clippy::neg_cmp_op_on_partial_ord)]

use std::collections::HashMap;

use crate::boolean3::cancelled_impl;
use crate::cancel::{is_cancelled, CancelToken};
use crate::impl_mesh::ManifoldImpl;
use crate::linalg::{cross, dot, length, Vec3};
use crate::progress::{begin_phase, complete_phase, Phase, ProgressReporter};
use crate::quickhull::convex_hull;
use crate::types::Error;

/// How often the support pass polls the cancel flag, as a power-of-two mask on
/// the face index. One dot product per face per tool vertex is cheap enough
/// that a poll per face would cost more than the work it guards.
const CANCEL_POLL_MASK: usize = 63;

/// Relative slack a candidate interior point must clear, and the relative slop
/// the output vertices are verified within, both scaled by the solid's
/// bounding-box diagonal.
const RELATIVE_TOLERANCE: f64 = 1e-9;

/// The erosion of `solid` by `tool` in closed form, when `solid` is convex.
///
/// Returns `None` when the caller must run the general erosion
/// ([`crate::minkowski::minkowski_difference`]). A cancelled run returns
/// `Some` of an empty impl carrying [`Error::Cancelled`] — `Some`, because a
/// cancelled token must never send the caller off to the path it was
/// cancelled out of.
///
/// Progress: one [`Phase::Minkowski`] unit per face of `solid` (the support
/// pass) plus one for everything after it, closed by `complete_phase` on
/// success only. The cheap declines (non-convex operand, a tool missing the
/// origin) happen before the phase opens and report nothing; the rare numeric
/// declines after it leave the phase open and short, deliberately.
pub fn try_compute(
    solid: &ManifoldImpl,
    tool: &ManifoldImpl,
    token: Option<&CancelToken>,
    progress: Option<&ProgressReporter>,
) -> Option<ManifoldImpl> {
    if is_cancelled(token) {
        return Some(cancelled_impl());
    }

    // Empty or errored operands are the general path's business.
    if solid.is_empty()
        || tool.is_empty()
        || solid.is_soup
        || tool.is_soup
        || solid.status != Error::NoError
        || tool.status != Error::NoError
    {
        return None;
    }

    // Tool convexity is not needed by the mathematics, but it is by the
    // promise to agree with the general erosion, which swaps its operands for
    // a convex solid and a non-convex tool.
    if !solid.is_convex() || !tool.is_convex() {
        return None;
    }

    if !tool_contains_origin(tool) {
        return None;
    }

    begin_phase(progress, Phase::Minkowski, solid.num_tri() as u64 + 1);

    let (normals, offsets) = match support_planes(solid, tool, token, progress) {
        SupportPlanes::Planes(normals, offsets) => (normals, offsets),
        SupportPlanes::Cancelled => return Some(cancelled_impl()),
        SupportPlanes::Unusable => return None,
    };

    let diagonal = solid.bbox.size();
    let tolerance = length(diagonal) * RELATIVE_TOLERANCE;

    let (_interior, slacks) = interior_point(solid, &normals, &offsets, tolerance)?;

    let vertices = solve_vertices(&normals, &offsets, &slacks, tolerance)?;

    if is_cancelled(token) {
        return Some(cancelled_impl());
    }

    let eroded = convex_hull(&vertices);

    // A cancelled token can never produce a NoError result: a cancel that
    // landed during the final hull must not come back as a finished erosion.
    if is_cancelled(token) {
        return Some(cancelled_impl());
    }

    if eroded.is_empty() {
        // Collapsed to a point, segment or sheet; the general path has a
        // considered answer for those.
        return None;
    }

    // Reached only on success, so a full bar is never a claim about work that
    // was abandoned.
    complete_phase(progress);
    Some(eroded)
}

/// Whether the origin lies inside or on `tool` — the condition under which the
/// closed form and the general sweep are the same function.
///
/// The sweep drops a point of A exactly when its swept copy `x - B` meets the
/// boundary; with the origin in B that is the erosion. Without it, a swept copy
/// can land wholly outside A and be kept (a 20-cube by a unit ball at
/// (0.8,0.8,0) sweeps to 5832.547 where the erosion is 5832). Exact because the
/// tool is convex: it contains the origin precisely when every outward face
/// plane has a non-negative offset, within the relative epsilon.
fn tool_contains_origin(tool: &ManifoldImpl) -> bool {
    let tolerance = length(tool.bbox.size()) * RELATIVE_TOLERANCE;

    for tri in 0..tool.num_tri() {
        let mut normal = tool.face_normal[tri];
        let len = length(normal);
        if !(len > 0.0) {
            // Degenerate face, no plane to test against.
            continue;
        }

        normal /= len;
        let plane_offset = dot(
            normal,
            tool.vert_pos[tool.halfedge[tri * 3].start_vert as usize],
        );
        if plane_offset < -tolerance {
            return false;
        }
    }

    true
}

/// Outcome of [`support_planes`].
enum SupportPlanes {
    /// Unit outward normals and their pushed-in offsets.
    Planes(Vec<Vec3>, Vec<f64>),
    /// The pass stopped on the cancel flag.
    Cancelled,
    /// Fewer than four planes survived; the faces no longer describe a solid.
    Unusable,
}

/// The solid's face planes, each pushed inward by the tool's support against
/// `-n_i` — the `d_i - h_B(-n_i)` of the file header. One progress unit per face.
fn support_planes(
    solid: &ManifoldImpl,
    tool: &ManifoldImpl,
    token: Option<&CancelToken>,
    progress: Option<&ProgressReporter>,
) -> SupportPlanes {
    let num_tri = solid.num_tri();
    let mut normals: Vec<Vec3> = Vec::with_capacity(num_tri);
    let mut offsets: Vec<f64> = Vec::with_capacity(num_tri);

    let tool_verts = &tool.vert_pos;

    for tri in 0..num_tri {
        if (tri & CANCEL_POLL_MASK) == 0 && is_cancelled(token) {
            return SupportPlanes::Cancelled;
        }

        if let Some(p) = progress {
            p.advance(1);
        }

        let mut normal = solid.face_normal[tri];
        let len = length(normal);
        if !(len > 0.0) {
            // A degenerate face bounds nothing; on a convex solid every plane
            // that matters is also carried by a face with area.
            continue;
        }

        // Exact when the normal is already unit: dividing by 1.0 moves no bit.
        normal /= len;

        let plane_offset = dot(
            normal,
            solid.vert_pos[solid.halfedge[tri * 3].start_vert as usize],
        );

        let mut support = f64::NEG_INFINITY;
        for &b in tool_verts.iter() {
            let reach = -dot(normal, b);
            if reach > support {
                support = reach;
            }
        }

        normals.push(normal);
        offsets.push(plane_offset - support);
    }

    // A bounded solid needs four planes.
    if normals.len() >= 4 {
        SupportPlanes::Planes(normals, offsets)
    } else {
        SupportPlanes::Unusable
    }
}

/// A point strictly inside the eroded body (the solid's vertex centroid) and
/// each plane's slack there; `None` when the centroid is not strictly inside
/// every pushed-in plane (near-total erosion of a skew solid).
fn interior_point(
    solid: &ManifoldImpl,
    normals: &[Vec3],
    offsets: &[f64],
    tolerance: f64,
) -> Option<(Vec3, Vec<f64>)> {
    let verts = &solid.vert_pos;
    let mut sum = Vec3::new(0.0, 0.0, 0.0);
    for &v in verts.iter() {
        sum += v;
    }

    let interior = sum / verts.len() as f64;
    let mut slacks: Vec<f64> = Vec::with_capacity(normals.len());

    for i in 0..normals.len() {
        let slack = offsets[i] - dot(normals[i], interior);
        if !(slack > tolerance) {
            return None;
        }
        slacks.push(slack);
    }

    Some((interior, slacks))
}

/// Exact-bits key for a dual point, the equality manifold-sharp's
/// `Vec3.Equals` uses (identically encoded components; `0.0 != -0.0`).
fn bits_key(v: Vec3) -> [u64; 3] {
    [v.x.to_bits(), v.y.to_bits(), v.z.to_bits()]
}

/// The vertices of the halfspace intersection, via the polar dual.
///
/// The dual hull only names WHICH three planes meet at each vertex; the vertex
/// is solved from the original normals and offsets. Hull vertices are mapped
/// back to their planes by the exact bits of the dual point, which works
/// because `convex_hull` selects and compacts its input points rather than
/// moving them, and a dual point determines its halfspace uniquely. The map is
/// probe-only (a later duplicate overwrites, as the C# indexer does); its
/// iteration order never reaches an output.
fn solve_vertices(
    normals: &[Vec3],
    offsets: &[f64],
    slacks: &[f64],
    tolerance: f64,
) -> Option<Vec<Vec3>> {
    let mut vertices: Vec<Vec3> = Vec::new();

    let mut dual_points: Vec<Vec3> = Vec::with_capacity(normals.len());
    let mut plane_of_dual_point: HashMap<[u64; 3], usize> = HashMap::with_capacity(normals.len());
    for i in 0..normals.len() {
        let dual = normals[i] / slacks[i];
        dual_points.push(dual);
        plane_of_dual_point.insert(bits_key(dual), i);
    }

    let dual_hull = convex_hull(&dual_points);
    if dual_hull.is_empty() {
        return None;
    }

    let num_tri = dual_hull.num_tri();
    for tri in 0..num_tri {
        let corner = |k: usize| {
            let v = dual_hull.vert_pos[dual_hull.halfedge[tri * 3 + k].start_vert as usize];
            plane_of_dual_point.get(&bits_key(v)).copied()
        };
        // A hull vertex that is not one of the input points cannot name its
        // triple, so the whole answer is handed back.
        let a = corner(0)?;
        let b = corner(1)?;
        let c = corner(2)?;

        let vertex = solve_plane_triple(
            normals[a], normals[b], normals[c], offsets[a], offsets[b], offsets[c],
        )?;
        vertices.push(vertex);
    }

    // The backstop: check the answer itself against every constraint.
    for &v in vertices.iter() {
        for i in 0..normals.len() {
            if dot(normals[i], v) > offsets[i] + tolerance {
                return None;
            }
        }
    }

    if vertices.is_empty() {
        None
    } else {
        Some(vertices)
    }
}

/// The point where three planes meet, by Cramer's rule; `None` when the unit
/// normals are close enough to coplanar (|triple product| < 1e-9, a pure angle
/// measure) that the point is noise.
fn solve_plane_triple(na: Vec3, nb: Vec3, nc: Vec3, da: f64, db: f64, dc: f64) -> Option<Vec3> {
    let b_cross_c = cross(nb, nc);
    let det = dot(na, b_cross_c);
    if det.abs() < 1e-9 {
        return None;
    }

    let c_cross_a = cross(nc, na);
    let a_cross_b = cross(na, nb);
    Some(((b_cross_c * da) + (c_cross_a * db) + (a_cross_b * dc)) / det)
}

#[cfg(test)]
#[path = "convex_erosion_tests.rs"]
mod tests;
