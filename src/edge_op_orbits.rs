// edge_op_orbits.rs — vertex-orbit scans for `split_pinched_verts` and
// `dedupe_edges` (edge_op.rs): the owner search that lets both run in
// parallel, and the per-orbit duplicate search they share.

use crate::impl_mesh::ManifoldImpl;
use crate::types::{next_halfedge, Halfedge};

// Both scans walk each vertex orbit (stepping `next_halfedge(paired_halfedge)`)
// once, from its smallest eligible halfedge, its owner, using visited flags in
// ascending order. C++ runs both in parallel (edge_op.cpp:722-796, 903-924).
//
// If every orbit is a closed cycle, a halfedge owns its orbit iff its walk
// returns to it without meeting a smaller eligible halfedge. If
// `paired_halfedge` is an involution the step is injective, so a walk either
// returns to its start or reaches a halfedge with no pair, which an open
// orbit's smallest eligible halfedge must do. Otherwise, or on an open orbit,
// the callers fall back to their sequential scans.
//
// Walking a whole orbit from every eligible halfedge costs the sum of the
// squared orbit lengths, Θ(n²) for one high-valence vertex. So the parallel
// walks stop after `ORBIT_WALK_CAP` steps, and the orbits they leave
// unresolved are settled by one sequential visited-flag pass, which walks
// each such orbit once. Why the owners are still exactly the sequential
// scan's, in its order, given closed orbits:
//
// - An orbit of length L <= cap: every walk from an eligible halfedge returns
//   to its start within L steps unless it first meets a smaller eligible one,
//   so the walks alone mark exactly its smallest eligible halfedge owner.
// - An orbit of length L > cap: no walk returns within the cap. A walk from a
//   halfedge that is not the smallest eligible one ends `NOT_OWNER` or
//   `UNRESOLVED`; the walk from the smallest one meets nothing smaller and
//   ends `UNRESOLVED`. So the orbit's `UNRESOLVED` halfedges are eligible and
//   include its smallest eligible halfedge, which is therefore the smallest
//   `UNRESOLVED` one. The ascending pass reaches it first, makes it the owner
//   and marks the whole orbit visited, so no other halfedge of the orbit
//   becomes an owner.
// - The pass walks only from `UNRESOLVED` halfedges, and a walk never leaves
//   its orbit, so short orbits keep the roles their walks gave them.
//
// Every orbit with an eligible halfedge thus has exactly one owner, its
// smallest eligible halfedge, as in the sequential scan, and the owners are
// collected ascending, the order in which that scan reaches them. An open
// orbit with an eligible halfedge is still caught: its smallest eligible
// halfedge's walk meets no smaller one, so it either reaches the open end
// (`OPEN`) or the cap, and then the pass reaches it first and its full walk
// reaches the open end.

/// Halfedge count from which the orbit scans run in parallel. Above C++'s 1e4:
/// every halfedge walks its own orbit, and below 100k that costs more than it saves.
pub(super) const ORBIT_PAR_THRESHOLD: usize = 100_001;

/// Steps after which a parallel orbit walk stops and leaves its orbit to the
/// sequential pass. Vertex valences rarely exceed it, so that pass is usually
/// empty, and it bounds the parallel work at 64 steps per halfedge.
pub(super) const ORBIT_WALK_CAP: usize = 64;

const NOT_OWNER: u8 = 0;
const OWNER: u8 = 1;
const OPEN: u8 = 2;
const UNRESOLVED: u8 = 3;

/// The owner of every orbit with an eligible halfedge, ascending, or `None`
/// below `threshold` or if an orbit is not a closed cycle.
pub(super) fn orbit_owners<F>(
    halfedge: &[Halfedge],
    threshold: usize,
    eligible: F,
) -> Option<Vec<usize>>
where
    F: Fn(&Halfedge) -> bool + Sync + Send,
{
    let n = halfedge.len();
    if n < threshold || !cfg!(feature = "parallel") {
        return None;
    }
    let paired_back = crate::par::maybe_par_filter(n, threshold, |i| {
        let p = halfedge[i].paired_halfedge;
        p >= 0 && (p as usize >= n || halfedge[p as usize].paired_halfedge != i as i32)
    });
    if !paired_back.is_empty() {
        return None;
    }
    // The step from `current`, or `None` at an open end (no pair, or a count
    // not a multiple of 3): those cases are left to the sequential scan.
    let step = |current: usize| -> Option<usize> {
        let p = halfedge[current].paired_halfedge;
        if p < 0 {
            return None;
        }
        let next = next_halfedge(p) as usize;
        (next < n).then_some(next)
    };
    let mut role: Vec<u8> = crate::par::maybe_par_map(n, threshold, |i| {
        if !eligible(&halfedge[i]) {
            return NOT_OWNER;
        }
        let mut current = i;
        for _ in 0..ORBIT_WALK_CAP {
            let Some(next) = step(current) else {
                return OPEN;
            };
            current = next;
            if current == i {
                return OWNER;
            }
            if current < i && eligible(&halfedge[current]) {
                return NOT_OWNER;
            }
        }
        UNRESOLVED
    });
    if role.contains(&OPEN) {
        return None;
    }
    // The sequential pass over the orbits longer than the cap (see above).
    let unresolved: Vec<usize> = (0..n).filter(|&i| role[i] == UNRESOLVED).collect();
    if !unresolved.is_empty() {
        let mut visited = vec![false; n];
        for i in unresolved {
            if visited[i] {
                role[i] = NOT_OWNER;
                continue;
            }
            role[i] = OWNER;
            let mut current = i;
            loop {
                visited[current] = true;
                current = step(current)?;
                if current == i {
                    break;
                }
            }
        }
    }
    Some((0..n).filter(|&i| role[i] == OWNER).collect())
}

/// Appends the duplicate edges of the orbit walked from `i`, in walk order:
/// halfedges whose end vertex a smaller halfedge of the orbit also ends at.
/// `mark` sees each halfedge the walk visits.
pub(super) fn dedupe_orbit(
    mesh: &ManifoldImpl,
    i: usize,
    n_edges: usize,
    mut mark: impl FnMut(usize),
    duplicates: &mut Vec<usize>,
) {
    // Track all endVerts seen in this vertex's orbit, keeping smallest edge idx.
    // Uses ForVert traversal: current = next_halfedge(halfedge[current].paired_halfedge)
    let mut end_verts: Vec<(i32, usize)> = Vec::new(); // (endVert, min_edge_idx)
                                                       // Process i itself first
    mark(i);
    let c_ev0 = mesh.halfedge[i].end_vert;
    if c_ev0 >= 0 {
        end_verts.push((c_ev0, i));
    }
    // Then orbit (with safety bound to prevent infinite loops)
    let mut current = i;
    let mut orbit_steps = 0;
    loop {
        let pair = mesh.halfedge[current].paired_halfedge;
        if pair < 0 {
            break;
        }
        current = next_halfedge(pair) as usize;
        if current == i {
            break;
        }
        orbit_steps += 1;
        if orbit_steps > n_edges {
            break;
        } // safety
        mark(current);
        let c_sv = mesh.halfedge[current].start_vert;
        let c_ev = mesh.halfedge[current].end_vert;
        if c_sv >= 0 && c_ev >= 0 {
            if let Some(entry) = end_verts.iter_mut().find(|(v, _)| *v == c_ev) {
                if current < entry.1 {
                    entry.1 = current;
                }
            } else {
                end_verts.push((c_ev, current));
            }
        }
    }

    // Second pass: find edges that aren't the minimum for their endVert
    let c_ev0 = mesh.halfedge[i].end_vert;
    if c_ev0 >= 0 {
        if let Some(&(_, min_edge)) = end_verts.iter().find(|(v, _)| *v == c_ev0) {
            if min_edge != i {
                duplicates.push(i);
            }
        }
    }
    current = i;
    orbit_steps = 0;
    loop {
        let pair = mesh.halfedge[current].paired_halfedge;
        if pair < 0 {
            break;
        }
        current = next_halfedge(pair) as usize;
        if current == i {
            break;
        }
        orbit_steps += 1;
        if orbit_steps > n_edges {
            break;
        } // safety
        let c_ev = mesh.halfedge[current].end_vert;
        if c_ev >= 0 {
            if let Some(&(_, min_edge)) = end_verts.iter().find(|(v, _)| *v == c_ev) {
                if min_edge != current {
                    duplicates.push(current);
                }
            }
        }
    }
}
