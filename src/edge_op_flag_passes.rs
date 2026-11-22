// edge_op_flag_passes.rs — the three flag-then-act passes of `simplify_topology`
// and `remove_degenerates` (edge_op.rs): short-edge collapse, colinear-edge
// collapse and degenerate-edge swaps. Each flags its candidate halfedges (in
// parallel on large meshes) and then acts on them in ascending order through
// edge_op.rs's `collapse_edge` and `recursive_edge_swap`.

use super::{collapse_edge, is_01_longest, recursive_edge_swap, tri_of};
use crate::face_op::get_axis_aligned_projection;
use crate::impl_mesh::ManifoldImpl;
use crate::linalg::{dot, Vec2};
use crate::polygon::ccw;
use crate::types::next_halfedge;

/// Halfedge count from which the edge-flag scans run in parallel, as C++
/// `FlagStore::run`'s `n > 1e5` (edge_op.cpp:90).
const FLAG_PAR_THRESHOLD: usize = 100_001;

/// Collapses edges shorter than epsilon_ when at least one endpoint is new,
/// repeating the pass until it collapses nothing.
///
/// Output change (manifold-sharp must mirror it in `EdgeOp.cs`
/// `CollapseShortEdges`, so that every caller repeats, not only its Boolean
/// path): this used to run one pass over the edges flagged before it started,
/// as the C++ does. The repeat applies to every caller: Booleans
/// (`first_new_vert > 0`) and the epsilon-only mode (`first_new_vert == 0`)
/// that runs in every import through `remove_degenerates`, in compose
/// (boolean3.rs), in the robust engine's clean-up (robust/assemble.rs) and in
/// the robust-assembled import (manifold_meshgl.rs).
///
/// Why: a Boolean may collapse a new vert into an old one across up to the
/// tolerance, but a new-to-new edge only across epsilon. Where several new
/// verts land within tolerance of one old vert, collapsing the first moves the
/// old vert onto the end of another new vert's edge, which was new-to-new when
/// the pass flagged and so was never flagged. One pass left that edge as a
/// needle triangle a billionth of the scale wide, whose corners coincide once
/// the result is exported as a single-precision `MeshGL`. Later passes flag
/// and collapse what earlier ones exposed; some inputs take several
/// (CraycloudBool 7 passes, `offset_triangulation_failure` 5,
/// `openscad_crash` 4).
///
/// The loop ends: a pass repeats only after a successful `collapse_edge`, each
/// of which removes the live triangle of its edge through `collapse_tri`
/// (`form_loop` may add verts but never triangles), and nothing in the pass
/// adds triangles, so there are at most as many successful collapses as
/// triangles.
///
/// Evidence: `edge_op::tests::test_union_leaves_no_edge_shorter_than_tolerance`,
/// MatterCAD's union of a thread and its core cylinder
/// (`testdata/thread-core-union-*.txt`): one pass leaves two halfedges of
/// length 1.03e-9 under the union's 5.4e-7 tolerance; without the f32
/// tolerance floor at import (manifold_meshgl.rs) the tolerance is 4e-12 and
/// needles export as two f32 triangles with coincident corners. The
/// loop and the floor together leave none, and the union's volume still
/// matches inclusion-exclusion to 1e-9. Expected values that moved:
/// `TEST(BooleanComplex, CraycloudBool)`'s difference of complements is empty
/// straight out of the Boolean, where one pass left needles only `simplify`
/// removed (`test_cpp_complex_craycloud`, `test_cpp_craycloud_bool`).
/// Compose output can also move with no pinned value changing: extra passes
/// collapse edges at epsilon three times in
/// `compose_tolerance::compose_with_a_part_finer_than_the_combined_epsilon_does_not_panic`.
pub fn collapse_short_edges(mesh: &mut ManifoldImpl, first_new_vert: i32) {
    while collapse_short_edges_once(mesh, first_new_vert) {}
}

/// One pass of [`collapse_short_edges`]; true when it collapsed an edge.
fn collapse_short_edges_once(mesh: &mut ManifoldImpl, first_new_vert: i32) -> bool {
    let mut scratch = Vec::with_capacity(10);
    let n = mesh.halfedge.len();

    // Per #1671: outside a Boolean (first_new_vert == 0) collapse only up to
    // epsilon to avoid error stacking; in a Boolean we only touch new verts, so
    // we may collapse up to tolerance. The per-edge max length below still
    // restricts new->old edges (old verts only move by epsilon).
    let tol = if first_new_vert == 0 {
        mesh.epsilon
    } else {
        mesh.tolerance
    };

    // Flag in parallel (see `par::maybe_par_filter`).
    let m: &ManifoldImpl = mesh;
    let flagged = crate::par::maybe_par_filter(n, FLAG_PAR_THRESHOLD, |i| {
        let h = &m.halfedge[i];
        if h.paired_halfedge < 0 {
            return false;
        }
        if h.start_vert < first_new_vert && h.end_vert < first_new_vert {
            return false;
        }
        if h.start_vert < 0 || h.end_vert < 0 {
            return false;
        }
        let delta = m.vert_pos[h.end_vert as usize] - m.vert_pos[h.start_vert as usize];
        let len_sq = dot(delta, delta);
        let max_len = if h.end_vert < first_new_vert {
            tol * tol
        } else {
            m.epsilon * m.epsilon
        };
        len_sq < max_len
    });

    let mut collapsed = false;
    for &i in &flagged {
        scratch.clear();
        collapsed |= collapse_edge(mesh, i, &mut scratch, tol, first_new_vert);
    }
    collapsed
}

/// Collapses redundant colinear edges (edges where the start vertex touches only
/// two distinct face groups).
pub fn collapse_colinear_edges(mesh: &mut ManifoldImpl, first_new_vert: i32) {
    let mut scratch = Vec::with_capacity(10);
    loop {
        let n = mesh.halfedge.len();

        // Flag in parallel (see `par::maybe_par_filter`).
        let m: &ManifoldImpl = mesh;
        let flagged = crate::par::maybe_par_filter(n, FLAG_PAR_THRESHOLD, |i| {
            let h = &m.halfedge[i];
            if h.paired_halfedge < 0 || h.start_vert < first_new_vert {
                return false;
            }
            if h.start_vert < 0 {
                return false;
            }
            if m.mesh_relation.tri_ref.is_empty() {
                return false;
            }
            if i / 3 >= m.mesh_relation.tri_ref.len() {
                return false;
            }

            let ref0 = m.mesh_relation.tri_ref[i / 3];
            let mut current = next_halfedge(m.halfedge[i].paired_halfedge) as usize;
            if current >= m.halfedge.len() {
                return false;
            }
            let mut ref1 = if current / 3 < m.mesh_relation.tri_ref.len() {
                m.mesh_relation.tri_ref[current / 3]
            } else {
                return false;
            };
            let mut ref1_updated = !ref0.same_face(&ref1);

            let max_orbit = m.halfedge.len(); // orbit can't exceed total halfedges
            let mut orbit_count = 0;
            while current != i {
                orbit_count += 1;
                if orbit_count > max_orbit {
                    return false;
                }
                // FlagEdge traversal: current = NextHalfedge(halfedge[current].pairedHalfedge)
                let pair = m.halfedge[current].paired_halfedge;
                if pair < 0 {
                    return false;
                }
                current = next_halfedge(pair) as usize;
                if current >= m.halfedge.len() {
                    return false;
                }
                let tri = current / 3;
                if tri >= m.mesh_relation.tri_ref.len() {
                    return false;
                }
                let ref_cur = m.mesh_relation.tri_ref[tri];
                if !ref_cur.same_face(&ref0) && !ref_cur.same_face(&ref1) {
                    if !ref1_updated {
                        ref1 = ref_cur;
                        ref1_updated = true;
                    } else {
                        return false;
                    }
                }
            }
            true
        });

        if flagged.is_empty() {
            break;
        }
        let mut num_collapsed = 0;
        for &i in &flagged {
            scratch.clear();
            if collapse_edge(mesh, i, &mut scratch, -1.0, 0) {
                num_collapsed += 1;
            }
        }
        if num_collapsed == 0 {
            break;
        }
    }
}

/// Swaps long edges of degenerate triangles.
pub fn swap_degenerates(mesh: &mut ManifoldImpl, first_new_vert: i32) {
    if mesh.face_normal.is_empty() {
        return;
    }

    let n = mesh.halfedge.len();

    // Flag in parallel (see `par::maybe_par_filter`).
    let m: &ManifoldImpl = mesh;
    let flagged = crate::par::maybe_par_filter(n, FLAG_PAR_THRESHOLD, |i| {
        let h = &m.halfedge[i];
        if h.paired_halfedge < 0 {
            return false;
        }
        // Skip edges where all 4 involved verts are old
        let tri0edge = tri_of(i);
        let pair = h.paired_halfedge as usize;
        let tri1edge = tri_of(pair);
        if m.halfedge[i].start_vert < first_new_vert
            && m.halfedge[i].end_vert < first_new_vert
            && m.halfedge[next_halfedge(i as i32) as usize].end_vert < first_new_vert
            && m.halfedge[next_halfedge(pair as i32) as usize].end_vert < first_new_vert
        {
            return false;
        }
        let tri = i / 3;
        if tri >= m.face_normal.len() {
            return false;
        }
        let proj = get_axis_aligned_projection(m.face_normal[tri]);
        let mut v = [Vec2::new(0.0, 0.0); 3];
        for j in 0..3 {
            let sv = m.halfedge[tri0edge[j]].start_vert;
            if sv >= 0 && (sv as usize) < m.vert_pos.len() {
                v[j] = proj.apply(m.vert_pos[sv as usize]);
            }
        }
        if ccw(v[0], v[1], v[2], m.tolerance) > 0 || !is_01_longest(v[0], v[1], v[2]) {
            return false;
        }
        // Switch to the neighbor's projection — C++ projects the PAIR
        // triangle's verts (pairTriEdge), not tri0's.
        let tri_p = pair / 3;
        if tri_p >= m.face_normal.len() {
            return false;
        }
        let proj_p = get_axis_aligned_projection(m.face_normal[tri_p]);
        for j in 0..3 {
            let sv = m.halfedge[tri1edge[j]].start_vert;
            if sv >= 0 && (sv as usize) < m.vert_pos.len() {
                v[j] = proj_p.apply(m.vert_pos[sv as usize]);
            }
        }
        ccw(v[0], v[1], v[2], m.tolerance) > 0 || is_01_longest(v[0], v[1], v[2])
    });

    let mut visited = vec![-1i32; n];
    let mut edge_swap_stack: Vec<i32> = Vec::new();
    let mut scratch: Vec<usize> = Vec::new();
    let mut tag = 0i32;

    for &i in &flagged {
        tag += 1;
        recursive_edge_swap(
            mesh,
            i as i32,
            &mut tag,
            &mut visited,
            &mut edge_swap_stack,
            &mut scratch,
        );
        while let Some(e) = edge_swap_stack.pop() {
            recursive_edge_swap(
                mesh,
                e,
                &mut tag,
                &mut visited,
                &mut edge_swap_stack,
                &mut scratch,
            );
        }
    }
}
