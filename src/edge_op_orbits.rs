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
// A walk from every halfedge costs the square of an orbit's length when the
// orbit meets its halfedges in ascending order, so a vertex of valence 1e5
// would take 1e10 steps against the sequential scan's 1e5. Walks therefore
// stop after `OWNER_WALK_CAP` steps, and the orbits they leave undecided are
// walked once each, sequentially, which keeps the search linear.

/// Halfedge count from which the orbit scans run in parallel. Above C++'s 1e4:
/// every halfedge walks part of its orbit, and below 100k that costs more than it saves.
pub(super) const ORBIT_PAR_THRESHOLD: usize = 100_001;

/// Steps a parallel owner walk takes before it leaves its orbit to the
/// sequential pass. Far above a typical valence, so that pass is rare, and it
/// bounds the parallel walks at this many steps per halfedge.
const OWNER_WALK_CAP: usize = 64;

/// What one halfedge's capped walk found.
#[derive(Clone, Copy, PartialEq)]
enum Role {
    NotOwner,
    Owner,
    /// The orbit is open, so owners are not the sequential scan's.
    Open,
    /// The walk reached `OWNER_WALK_CAP`; the orbit is decided in a later pass.
    Long,
}

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
    // The next halfedge of an orbit, or `None` where it is open. A count that
    // is not a multiple of 3 also ends up here, and goes to the sequential scan.
    let step = |current: usize| {
        let p = halfedge[current].paired_halfedge;
        if p < 0 {
            return None;
        }
        let next = next_halfedge(p) as usize;
        (next < n).then_some(next)
    };
    let mut role: Vec<Role> = crate::par::maybe_par_map(n, threshold, |i| {
        if !eligible(&halfedge[i]) {
            return Role::NotOwner;
        }
        let mut current = i;
        for _ in 0..OWNER_WALK_CAP {
            let Some(next) = step(current) else {
                return Role::Open;
            };
            current = next;
            if current == i {
                return Role::Owner;
            }
            if current < i && eligible(&halfedge[current]) {
                return Role::NotOwner;
            }
        }
        Role::Long
    });
    if role.contains(&Role::Open) {
        return None;
    }
    // An orbit with a `Long` halfedge has more than `OWNER_WALK_CAP`
    // halfedges, so none of its walks returned, and its owner, whose walk
    // meets no smaller eligible halfedge, is `Long` too. As `Long` halfedges
    // are eligible, the owner is the orbit's smallest `Long` halfedge: the
    // first this ascending scan reaches. Each such orbit is walked once, to
    // mark it visited and to find an open one.
    let long = crate::par::maybe_par_filter(n, threshold, |i| role[i] == Role::Long);
    let mut visited = vec![false; if long.is_empty() { 0 } else { n }];
    for &owner in &long {
        if visited[owner] {
            continue;
        }
        role[owner] = Role::Owner;
        let mut current = owner;
        loop {
            visited[current] = true;
            current = step(current)?;
            if current == owner {
                break;
            }
        }
    }
    Some((0..n).filter(|&i| role[i] == Role::Owner).collect())
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
