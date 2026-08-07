// robust/cells.rs — Arrangement cell complex and combinatorial winding
// propagation (Zhou, Grinspun, Zorin, Jacobson 2016, "Mesh Arrangements for
// Solid Geometry" — the formulation libigl's mesh_boolean uses).
//
// This replaces the per-component winding queries of robust/mod.rs with the
// structure that makes local inconsistency unrepresentable. The pieces of
// robust/intersection_graph.rs already form an intersection-free arrangement
// of both operands; what this module adds is the *dual*:
//
//   1. Every piece has two sides (normal / anti). Around each arrangement
//      edge the incident half-faces are radially sorted (the same exact
//      basis + angle_cmp the ring regularization uses), and the wedge
//      between consecutive radial positions unions the sides it touches.
//      The resulting equivalence classes are the 3D cells of the
//      arrangement.
//   2. Winding numbers are then propagated cell-to-cell by breadth-first
//      search: crossing a piece from its normal side to its anti side enters
//      the solid its operand bounds, so w[mesh] increases by one. Because
//      every cell's winding is *derived from one traversal*, adjacent
//      regions cannot disagree — the failure mode of independent per-region
//      queries (two pieces meeting along a segment both classified "outside",
//      leaving a surface that cannot close) does not exist here.
//
// Coincident duplicate faces fall out correctly: they share a radial
// position, so they form one "wall" whose winding step is the signed sum of
// the stack. A doubled sheet steps the winding by two and a fold cancels to
// zero, without any explicit regularization pass.

use std::cmp::Ordering;
use std::collections::HashMap;

use num_rational::BigRational;
use num_traits::{Signed, Zero};

use crate::disjoint_sets::DisjointSets;
use crate::linalg::Vec3;

use super::classify::angle_cmp;
use super::exact::rational::R3;
use super::intersection_graph::{edge_key, EdgeKey, IntersectionGraph, Piece};
use crate::types::OpType;

/// Side of a piece: `NORMAL` is the half-space its outward normal points
/// into, `ANTI` the one behind it (inside the solid its operand bounds).
pub const NORMAL: usize = 0;
pub const ANTI: usize = 1;

/// Node id in the side union-find: two per piece.
#[inline]
fn node(piece: usize, side: usize) -> u32 {
    (2 * piece + side) as u32
}

/// The arrangement's cell decomposition.
pub struct CellComplex {
    /// Compact cell id per (piece, side); index with [`node`].
    pub cell_of: Vec<u32>,
    pub num_cells: usize,
}

impl CellComplex {
    #[inline]
    pub fn cell(&self, piece: usize, side: usize) -> usize {
        self.cell_of[node(piece, side) as usize] as usize
    }
}

/// One half-face incident to an arrangement edge.
struct Inc {
    piece: usize,
    /// Traversal runs key.0 → key.1 (the edge's canonical direction).
    forward: bool,
    /// Radial direction of the opposite vertex, in the edge's (u, v) basis.
    du: BigRational,
    dv: BigRational,
}

impl Inc {
    /// The side facing counter-clockwise (increasing radial angle).
    ///
    /// A forward-traversing face has normal ∝ w × d, which sits 90° CCW of
    /// its apex direction, so its normal side is the CCW one; a backward
    /// face's normal is 90° CW, so the relationship inverts.
    #[inline]
    fn ccw_side(&self) -> usize {
        if self.forward {
            NORMAL
        } else {
            ANTI
        }
    }

    #[inline]
    fn cw_side(&self) -> usize {
        1 - self.ccw_side()
    }
}

/// Build the cell complex over every piece of the graph.
///
/// Discarded/regularized pieces are deliberately *not* excluded: thin
/// material cancels arithmetically in the winding sum, which is both simpler
/// and more robust than deciding up front which sheets are real.
pub fn build_cells(graph: &IntersectionGraph) -> CellComplex {
    let n = graph.pieces.len();
    let ds = DisjointSets::new((2 * n).max(1) as u32);

    // Edge → incident half-faces. Radial directions are filled per edge,
    // since the (u, v) basis depends on the edge axis.
    let mut incident: HashMap<EdgeKey, Vec<(usize, bool, u32)>> = HashMap::new();
    for (pi, piece) in graph.pieces.iter().enumerate() {
        let vi = piece.vi;
        for e in 0..3 {
            let (a, b) = (vi[e], vi[(e + 1) % 3]);
            let apex = vi[(e + 2) % 3];
            incident
                .entry(edge_key(a, b))
                .or_default()
                .push((pi, a < b, apex));
        }
    }

    for (key, raw) in &incident {
        if raw.len() < 2 {
            continue; // a boundary edge bounds no wedge
        }
        let k0 = &graph.verts[key.0 as usize];
        let k1 = &graph.verts[key.1 as usize];
        let mut incs = radial_directions(k0, k1, raw, graph);
        if incs.len() < 2 {
            continue;
        }
        incs.sort_by(|a, b| {
            angle_cmp((&a.du, &a.dv), (&b.du, &b.dv)).then_with(|| a.piece.cmp(&b.piece))
        });

        // Group coincident radial directions into single walls, then link
        // consecutive walls around the cyclic order.
        let groups = group_by_angle(&incs);
        for gi in 0..groups.len() {
            let (s, e) = groups[gi];
            // All faces of a wall share the cell on each of its two sides.
            for k in s..e {
                ds.unite(
                    node(incs[s].piece, incs[s].ccw_side()),
                    node(incs[k].piece, incs[k].ccw_side()),
                );
                ds.unite(
                    node(incs[s].piece, incs[s].cw_side()),
                    node(incs[k].piece, incs[k].cw_side()),
                );
            }
            // The wedge between this wall and the next: CCW side of this
            // wall meets the CW side of the next.
            let (ns, _) = groups[(gi + 1) % groups.len()];
            if groups.len() > 1 {
                ds.unite(
                    node(incs[s].piece, incs[s].ccw_side()),
                    node(incs[ns].piece, incs[ns].cw_side()),
                );
            }
        }
    }

    // Compact the union-find roots into dense cell ids.
    let mut cell_of = vec![0u32; 2 * n];
    let mut remap: HashMap<u32, u32> = HashMap::new();
    for i in 0..(2 * n) {
        let root = ds.find(i as u32);
        let next = remap.len() as u32;
        cell_of[i] = *remap.entry(root).or_insert(next);
    }
    CellComplex {
        num_cells: remap.len(),
        cell_of,
    }
}

/// Radial direction of each incident face's apex, in a right-handed (u, v, w)
/// basis with w along the edge — the same construction the ring
/// regularization uses, so both agree on angular order.
fn radial_directions(
    k0: &R3,
    k1: &R3,
    raw: &[(usize, bool, u32)],
    graph: &IntersectionGraph,
) -> Vec<Inc> {
    let w = k1.sub(k0);
    let (ax, ay, az) = (w.x.abs(), w.y.abs(), w.z.abs());
    let unit = |i: usize| {
        let z = BigRational::zero;
        let o = || BigRational::from_integer(1.into());
        match i {
            0 => R3::new(o(), z(), z()),
            1 => R3::new(z(), o(), z()),
            _ => R3::new(z(), z(), o()),
        }
    };
    let k = if ax <= ay && ax <= az {
        0
    } else if ay <= az {
        1
    } else {
        2
    };
    let u = w.cross(&unit(k));
    let v = w.cross(&u);

    let mut out = Vec::with_capacity(raw.len());
    for &(piece, forward, apex) in raw {
        let a = &graph.verts[apex as usize];
        let (du_n, du_d) = super::exact::predicates::dot_diff_raw(a, k0, &u);
        let (dv_n, dv_d) = super::exact::predicates::dot_diff_raw(a, k0, &v);
        let du = BigRational::new_raw(du_n, du_d);
        let dv = BigRational::new_raw(dv_n, dv_d);
        if du.is_zero() && dv.is_zero() {
            continue; // apex on the axis: a degenerate sliver, no wedge
        }
        out.push(Inc {
            piece,
            forward,
            du,
            dv,
        });
    }
    out
}

/// Half-open `[start, end)` ranges of equal radial direction.
fn group_by_angle(incs: &[Inc]) -> Vec<(usize, usize)> {
    let mut groups = Vec::new();
    let mut i = 0;
    while i < incs.len() {
        let mut j = i + 1;
        while j < incs.len()
            && angle_cmp((&incs[i].du, &incs[i].dv), (&incs[j].du, &incs[j].dv)) == Ordering::Equal
        {
            j += 1;
        }
        groups.push((i, j));
        i = j;
    }
    groups
}

/// Per-cell winding numbers, one entry per operand.
pub struct Windings {
    /// `w[cell] = [w_P, w_Q]`, valid only where `known[cell]`.
    pub w: Vec<[i32; 2]>,
    pub known: Vec<bool>,
}

/// Propagate winding numbers outward from `seed` (whose winding is
/// `seed_w`), crossing one piece at a time.
///
/// Crossing a piece from its normal side to its anti side enters the solid
/// its operand bounds, so that operand's winding gains one. Cells in other
/// connected components of the arrangement stay `known == false`; the caller
/// resolves those with one point query each.
pub fn propagate_windings(
    graph: &IntersectionGraph,
    complex: &CellComplex,
    seed: usize,
    seed_w: [i32; 2],
) -> Windings {
    let adj = cell_adjacency(graph, complex);
    let mut out = Windings {
        w: vec![[0i32; 2]; complex.num_cells],
        known: vec![false; complex.num_cells],
    };
    if seed < complex.num_cells {
        self::seed(&mut out, seed, seed_w);
        bfs(&adj, &mut out, seed);
    }
    out
}

/// Winding numbers for *every* cell.
///
/// The combinatorial BFS only reaches cells connected through shared
/// arrangement edges, so disjoint or nested components need a seed each:
/// exactly the residual ray-shooting libigl's `propagate_winding_numbers`
/// performs. One exact query pair per component, not per surface region.
pub fn propagate_all(
    graph: &IntersectionGraph,
    complex: &CellComplex,
    tris_f64: [&[[Vec3; 3]]; 2],
    tris_r: [&[[R3; 3]]; 2],
    boxes: [&[crate::types::Box]; 2],
) -> Windings {
    let adj = cell_adjacency(graph, complex);
    let mut out = Windings {
        w: vec![[0i32; 2]; complex.num_cells],
        known: vec![false; complex.num_cells],
    };

    // A representative (piece, side) per cell, for seeding by point query.
    let mut rep: Vec<Option<(usize, usize)>> = vec![None; complex.num_cells];
    for pi in 0..graph.pieces.len() {
        for side in [NORMAL, ANTI] {
            let c = complex.cell(pi, side);
            if rep[c].is_none() {
                rep[c] = Some((pi, side));
            }
        }
    }

    if let Some(outer) = outer_cell(graph, complex) {
        seed(&mut out, outer, [0, 0]);
        bfs(&adj, &mut out, outer);
    }
    // Any component the traversal could not reach gets its own exact query.
    for c in 0..complex.num_cells {
        if out.known[c] {
            continue;
        }
        let Some((pi, side)) = rep[c] else { continue };
        let pv = graph.piece_verts(pi);
        let point = super::ray_shoot::piece_centroid(pv);
        let n = pv[1].sub(pv[0]).cross(&pv[2].sub(pv[0]));
        let outward = if side == NORMAL {
            n
        } else {
            R3::new(-&n.x, -&n.y, -&n.z)
        };
        let mut w = [0i32; 2];
        for m in 0..2 {
            w[m] = super::ray_shoot::winding_off_surface(
                &point,
                &outward,
                tris_r[m],
                tris_f64[m],
                boxes[m],
            );
        }
        seed(&mut out, c, w);
        bfs(&adj, &mut out, c);
    }
    out
}

fn seed(out: &mut Windings, cell: usize, w: [i32; 2]) {
    out.w[cell] = w;
    out.known[cell] = true;
}

/// Winding step between adjacent cells, **summed** over every piece that
/// separates them.
///
/// Aggregation is what gives coincident stacks their multiplicity: a doubled
/// sheet contributes +2 and a fold cancels to 0, so self-overlapping input
/// classifies correctly with no separate regularization pass. Treating the
/// coincident pieces as independent adjacencies would apply only the first
/// and silently lose the rest.
fn cell_adjacency(
    graph: &IntersectionGraph,
    complex: &CellComplex,
) -> Vec<Vec<(usize, [i32; 2])>> {
    let mut agg: HashMap<(usize, usize), [i32; 2]> = HashMap::new();
    for Wall { rep, delta, .. } in walls(graph) {
        let (cn, ca) = (complex.cell(rep, NORMAL), complex.cell(rep, ANTI));
        if cn == ca {
            continue; // both sides in one cell: the sheet bounds nothing
        }
        let (lo, hi, s) = if cn < ca { (cn, ca, 1) } else { (ca, cn, -1) };
        agg.insert((lo, hi), [delta[0] * s, delta[1] * s]);
    }

    let mut adj: Vec<Vec<(usize, [i32; 2])>> = vec![Vec::new(); complex.num_cells];
    for ((lo, hi), d) in agg {
        adj[lo].push((hi, d));
        adj[hi].push((lo, [-d[0], -d[1]]));
    }
    adj
}

/// Does a cell with this winding vector lie inside the operation's result?
///
/// Each operand's solid is `{w ≥ 1}` — a negative winding is inverted
/// geometry, never material. Expressing every operation as one predicate on
/// the winding vector is what lets subtraction drop its operand-flip trick:
/// P − Q is just "inside P and not inside Q".
pub fn in_result(op: OpType, w: [i32; 2]) -> bool {
    let (a, b) = (w[0] >= 1, w[1] >= 1);
    match op {
        OpType::Add => a || b,
        OpType::Intersect => a && b,
        OpType::Subtract => a && !b,
    }
}

/// The boundary of the result: every wall whose two cells disagree about
/// containment, wound so its normal points out of the solid.
///
/// Orientation is *derived* from the cell labels rather than inherited from
/// the input face. That is what makes the output closed and consistently
/// oriented no matter how the input was wound — an inverted region of a
/// self-intersecting scan still lands in the right cells, and its emitted
/// orientation is corrected on the way out. One representative per wall is
/// emitted, so a coincident stack contributes a single face.
pub fn extract(
    graph: &IntersectionGraph,
    complex: &CellComplex,
    wind: &Windings,
    op: OpType,
) -> Vec<Piece> {
    let mut out = Vec::new();
    for Wall { rep, .. } in walls(graph) {
        let (cn, ca) = (complex.cell(rep, NORMAL), complex.cell(rep, ANTI));
        if cn == ca || !wind.known[cn] || !wind.known[ca] {
            continue;
        }
        let (in_n, in_a) = (in_result(op, wind.w[cn]), in_result(op, wind.w[ca]));
        if in_n == in_a {
            continue; // same material both sides: not a boundary
        }
        let piece = &graph.pieces[rep];
        // The representative's normal points from its anti side toward its
        // normal side. Material belongs behind the emitted normal, so the
        // winding reverses when the solid is on the normal side instead.
        let vi = if in_a {
            piece.vi
        } else {
            [piece.vi[0], piece.vi[2], piece.vi[1]]
        };
        out.push(Piece {
            mesh: piece.mesh,
            tri: piece.tri,
            vi,
        });
    }
    out
}

/// One distinct triangle of the arrangement, with the coincident stack that
/// occupies it collapsed into a single winding step.
pub struct Wall {
    /// Representative piece; its winding fixes the wall's normal side.
    pub rep: usize,
    /// Winding change per operand, crossing the representative's normal side
    /// to its anti side.
    pub delta: [i32; 2],
}

/// Group pieces into walls by exact triangle identity.
///
/// Pieces occupying the same triangle are a coincident stack whose
/// contributions add — a doubled sheet steps by two, a fold cancels to zero.
/// Two *different* triangles between the same pair of cells are alternative
/// crossings of one boundary and are never summed; that distinction is why
/// aggregation keys on the triangle rather than the cell pair.
fn walls(graph: &IntersectionGraph) -> Vec<Wall> {
    let mut by_tri: HashMap<[u32; 3], (usize, bool, [i32; 2])> = HashMap::new();
    for (pi, piece) in graph.pieces.iter().enumerate() {
        let (key, parity) = canonical(piece.vi);
        let m = piece.mesh as usize;
        let entry = by_tri.entry(key).or_insert((pi, parity, [0; 2]));
        // Opposite winding means this piece's normal side is the
        // representative's anti side, so it steps the other way.
        entry.2[m] += if parity == entry.1 { 1 } else { -1 };
    }
    by_tri
        .into_values()
        .map(|(rep, _, delta)| Wall { rep, delta })
        .collect()
}

/// Canonical key for a triangle (sorted vertex ids) plus the parity of the
/// piece's winding against that order — the same identity the coincident
/// binding uses, so both agree on what "the same triangle" means.
fn canonical(vi: [u32; 3]) -> ([u32; 3], bool) {
    let mut sorted = vi;
    sorted.sort_unstable();
    let i = (0..3).min_by_key(|&i| vi[i]).unwrap_or(0);
    let rotated = [vi[i], vi[(i + 1) % 3], vi[(i + 2) % 3]];
    (sorted, rotated == sorted)
}

fn bfs(adj: &[Vec<(usize, [i32; 2])>], out: &mut Windings, start: usize) {
    let mut queue = std::collections::VecDeque::from([start]);
    while let Some(c) = queue.pop_front() {
        let base = out.w[c];
        for &(next, d) in &adj[c] {
            if out.known[next] {
                continue;
            }
            seed(out, next, [base[0] + d[0], base[1] + d[1]]);
            queue.push_back(next);
        }
    }
}

/// The unbounded cell, found via an outer facet: at the lexicographically
/// largest vertex of the arrangement, the incident face whose normal is most
/// nearly parallel to +x has unbounded space on its +x side.
///
/// Returns `None` when every incident face is parallel to x (a knife edge at
/// the extreme vertex), leaving the caller to seed by point query instead.
pub fn outer_cell(graph: &IntersectionGraph, complex: &CellComplex) -> Option<usize> {
    let n = graph.pieces.len();
    if n == 0 {
        return None;
    }
    // Lexicographically largest vertex actually used by a piece.
    let mut best_v: Option<u32> = None;
    for piece in &graph.pieces {
        for &v in &piece.vi {
            let better = match best_v {
                None => true,
                Some(b) => lex_gt(&graph.verts[v as usize], &graph.verts[b as usize]),
            };
            if better {
                best_v = Some(v);
            }
        }
    }
    let vmax = best_v?;

    // Among faces at that vertex, maximize nx² / (n·n) — compared by
    // cross-multiplication so no roots or divisions are needed.
    let mut best: Option<(usize, BigRational, BigRational, BigRational)> = None;
    for (pi, piece) in graph.pieces.iter().enumerate() {
        if !piece.vi.contains(&vmax) {
            continue;
        }
        let pv = graph.piece_verts(pi);
        let nrm = pv[1].sub(pv[0]).cross(&pv[2].sub(pv[0]));
        let nx2 = &nrm.x * &nrm.x;
        let nn = &nrm.x * &nrm.x + &nrm.y * &nrm.y + &nrm.z * &nrm.z;
        if nn.is_zero() {
            continue;
        }
        let better = match &best {
            None => true,
            // nx2/nn > bx2/bn  ⇔  nx2·bn > bx2·nn   (nn, bn > 0)
            Some((_, bx2, bn, _)) => &nx2 * bn > bx2 * &nn,
        };
        if better {
            best = Some((pi, nx2, nn, nrm.x.clone()));
        }
    }
    let (pi, nx2, _, nx) = best?;
    if nx2.is_zero() {
        return None;
    }
    // Unbounded space lies on whichever side the +x-facing normal points to.
    let side = if nx.is_positive() { NORMAL } else { ANTI };
    Some(complex.cell(pi, side))
}

fn lex_gt(a: &R3, b: &R3) -> bool {
    match a.x.cmp(&b.x) {
        Ordering::Greater => return true,
        Ordering::Less => return false,
        Ordering::Equal => {}
    }
    match a.y.cmp(&b.y) {
        Ordering::Greater => return true,
        Ordering::Less => return false,
        Ordering::Equal => {}
    }
    a.z > b.z
}

#[cfg(test)]
#[path = "cells_tests.rs"]
mod tests;
