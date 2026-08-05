// robust/classify.rs — Radial regularization around intersection segments
// and coincident-piece binding (paper §7.1).
//
// For every exact intersection sub-edge (robust/intersection_graph.rs), the
// incident pieces from both meshes are sorted radially around the edge with
// pure sign arithmetic (quadrant + 2D cross, no trigonometry). Within each
// coincident-direction group, opposite-traversal pairs are infinitely thin
// material and are discarded (regularization, §7.1). Exactly coincident
// whole pieces (coplanar overlaps): opposite-orientation pairs are likewise
// discarded; same-orientation pairs get one ∪ and one ∩ binding so each
// shared region survives exactly once per output.
//
// The absolute ∪/∩ tags for everything else come from exact winding-number
// queries in robust/mod.rs (one per surface component, per piece for
// self-intersecting operands). The paper's local Prop 2/3 ring walk was
// replaced by those queries: it silently misclassifies pieces where an
// operand's winding exceeds 1 (self-overlapping sheets make a 1↔2 crossing
// locally indistinguishable from the 0↔1 crossing the walk assumes).

use std::cmp::Ordering;
use std::collections::BTreeMap;

use num_rational::BigRational;
use num_traits::{Signed, Zero};

use super::exact::rational::R3;
use super::exact::Sign;
use super::intersection_graph::{edge_key, EdgeKey, IntersectionGraph};

/// Which boolean output a piece belongs to.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Tag {
    Union,
    Inter,
}

/// Per-piece classification result. Tags are set only for coincident-bound
/// pieces; everything else is `None`, resolved by flood fill + winding
/// queries in robust/mod.rs.
pub struct Classification {
    pub tags: Vec<Option<Tag>>,
    pub discarded: Vec<bool>,
}

struct Incident {
    piece: usize,
    /// Piece winding traverses key.0 → key.1.
    forward: bool,
    /// The vertex opposite the ring edge; source of the radial direction.
    apex: R3,
    /// Radial direction of the apex: exact (d·u, d·v) coordinates, filled in
    /// once the ring's basis exists.
    du: BigRational,
    dv: BigRational,
}

/// CCW angular comparison of two nonzero direction vectors.
fn angle_cmp(a: (&BigRational, &BigRational), b: (&BigRational, &BigRational)) -> Ordering {
    fn quadrant(du: &BigRational, dv: &BigRational) -> u8 {
        let (su, sv) = (Sign::of_rat(du), Sign::of_rat(dv));
        debug_assert!(!(su == Sign::Zero && sv == Sign::Zero), "zero direction");
        match (su, sv) {
            (Sign::Pos, Sign::Pos) | (Sign::Pos, Sign::Zero) => 0,
            (Sign::Zero, Sign::Pos) | (Sign::Neg, Sign::Pos) => 1,
            (Sign::Neg, Sign::Zero) | (Sign::Neg, Sign::Neg) => 2,
            _ => 3,
        }
    }
    let (qa, qb) = (quadrant(a.0, a.1), quadrant(b.0, b.1));
    if qa != qb {
        return qa.cmp(&qb);
    }
    // Same quadrant: CCW order by cross-product sign.
    let cross = a.0 * b.1 - a.1 * b.0;
    match Sign::of_rat(&cross) {
        Sign::Pos => Ordering::Less,
        Sign::Neg => Ordering::Greater,
        Sign::Zero => Ordering::Equal, // coincident direction
    }
}

/// Coincident-piece regularization and binding (paper §7.1 done globally
/// rather than per ring, so every ring sees one consistent decision):
/// exactly-equal pieces from the two meshes — the copies coplanar overlap
/// regions produce on each side, identical thanks to the common-subdivision
/// guarantees — are paired up. Opposite-winding pairs are infinitely thin
/// material and get discarded outright; same-winding pairs are the shared
/// regions both outputs keep once, so the P copy is bound to Union and the
/// Q copy to Inter. (The paper notes the choice is arbitrary; making it
/// globally per piece, not per ring, is what keeps it conflict-free.)
fn bind_coincident(
    graph: &IntersectionGraph,
    tags: &mut [Option<Tag>],
    discarded: &mut [bool],
) {
    // Canonical key: sorted vertex triple; parity bit = winding orientation
    // relative to the sorted order.
    let mut by_key: BTreeMap<[R3; 3], Vec<(usize, bool)>> = BTreeMap::new();
    for (pi, piece) in graph.pieces.iter().enumerate() {
        let mut sorted = piece.v.clone();
        sorted.sort();
        // parity: does the winding cycle (v0,v1,v2), rotated to start at the
        // smallest vertex, match (sorted0, sorted1, sorted2)?
        let start = piece.v.iter().position(|v| *v == sorted[0]).unwrap();
        let same = piece.v[(start + 1) % 3] == sorted[1];
        by_key.entry(sorted).or_default().push((pi, same));
    }
    for owners in by_key.values() {
        if owners.len() < 2 {
            continue;
        }
        // Cancel opposite-parity cross-mesh pairs first (regularization)...
        let mut p_side: Vec<(usize, bool)> = Vec::new();
        let mut q_side: Vec<(usize, bool)> = Vec::new();
        for &(pi, parity) in owners {
            if graph.pieces[pi].mesh == 0 {
                p_side.push((pi, parity));
            } else {
                q_side.push((pi, parity));
            }
        }
        let mut used_q = vec![false; q_side.len()];
        for &(pi, parity) in &p_side {
            if let Some(k) = (0..q_side.len())
                .find(|&k| !used_q[k] && q_side[k].1 != parity)
            {
                used_q[k] = true;
                discarded[pi] = true;
                discarded[q_side[k].0] = true;
            }
        }
        // ...then bind same-parity cross-mesh pairs: P → Union, Q → Inter.
        for &(pi, parity) in &p_side {
            if discarded[pi] {
                continue;
            }
            if let Some(k) = (0..q_side.len())
                .find(|&k| !used_q[k] && q_side[k].1 == parity)
            {
                used_q[k] = true;
                tags[pi] = Some(Tag::Union);
                tags[q_side[k].0] = Some(Tag::Inter);
            }
        }
    }
}

/// Regularize every intersection ring (cancel coincident opposite-traversal
/// pieces) and bind exactly-coincident cross-mesh piece pairs.
pub fn classify_rings(graph: &IntersectionGraph) -> Classification {
    let n = graph.pieces.len();
    let mut tags: Vec<Option<Tag>> = vec![None; n];
    let mut discarded = vec![false; n];

    bind_coincident(graph, &mut tags, &mut discarded);

    // Ring construction: intersection edge key → incident pieces (globally
    // discarded pieces excluded up front).
    let mut rings: BTreeMap<EdgeKey, Vec<Incident>> = BTreeMap::new();
    for (pi, piece) in graph.pieces.iter().enumerate() {
        if discarded[pi] {
            continue;
        }
        for e in 0..3 {
            let a = &piece.v[e];
            let b = &piece.v[(e + 1) % 3];
            let key = edge_key(a, b);
            if !graph.isect_edges.contains(&key) {
                continue;
            }
            let forward = *a == key.0; // winding visits key.0 → key.1
            let apex = &piece.v[(e + 2) % 3];
            rings.entry(key).or_default().push(Incident {
                piece: pi,
                forward,
                apex: apex.clone(),
                du: BigRational::zero(),
                dv: BigRational::zero(),
            });
        }
    }

    for (key, incidents) in rings.iter_mut() {
        regularize_one_ring(key, incidents, &mut discarded);
    }

    Classification { tags, discarded }
}

fn regularize_one_ring(key: &EdgeKey, incidents: &mut [Incident], discarded: &mut [bool]) {
    // Radial basis: w along the edge, u ⊥ w via the axis of w's smallest
    // |component| (never parallel), v = w × u; (u, v, w) is right-handed.
    let w = key.1.sub(&key.0);
    let ax = w.x.abs();
    let ay = w.y.abs();
    let az = w.z.abs();
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
    debug_assert!(!u.is_zero());
    let v = w.cross(&u);

    for inc in incidents.iter_mut() {
        let d = inc.apex.sub(&key.0);
        inc.du = d.dot(&u);
        inc.dv = d.dot(&v);
        debug_assert!(
            !(inc.du.is_zero() && inc.dv.is_zero()),
            "apex on the ring axis"
        );
    }

    // CCW radial sort; coincident directions tie-break by piece id for
    // determinism.
    incidents.sort_by(|a, b| {
        angle_cmp((&a.du, &a.dv), (&b.du, &b.dv)).then_with(|| a.piece.cmp(&b.piece))
    });

    // Regularization: within each coincident-direction group, cancel
    // opposite-traversal pairs (discard both — infinitely thin material).
    let m = incidents.len();
    let mut cancelled = vec![false; m];
    let mut i = 0;
    while i < m {
        let mut j = i;
        while j + 1 < m
            && angle_cmp(
                (&incidents[i].du, &incidents[i].dv),
                (&incidents[j + 1].du, &incidents[j + 1].dv),
            ) == Ordering::Equal
        {
            j += 1;
        }
        // Group [i, j]: pair up +w with -w traversals.
        let mut fwd: Vec<usize> = (i..=j).filter(|&x| incidents[x].forward).collect();
        let mut bwd: Vec<usize> = (i..=j).filter(|&x| !incidents[x].forward).collect();
        while let (Some(f), Some(bk)) = (fwd.pop(), bwd.pop()) {
            cancelled[f] = true;
            cancelled[bk] = true;
            discarded[incidents[f].piece] = true;
            discarded[incidents[bk].piece] = true;
        }
        i = j + 1;
    }

}

#[cfg(test)]
#[path = "classify_tests.rs"]
mod tests;
