// robust/classify.rs — Radial sort and local ∪/∩ classification around
// intersection segments (paper §7.1–7.2).
//
// For every exact intersection sub-edge (robust/intersection_graph.rs), the
// incident pieces from both meshes are sorted radially around the edge with
// pure sign arithmetic (quadrant + 2D cross, no trigonometry). Successive
// pairs then classify by Propositions 2/3:
//   - opposite edge-traversal (orientable pair)  → same tag,
//   - same edge-traversal (non-orientable pair)  → opposite tags, and the
//     absolute assignment follows from the traversal direction: with the
//     radial order counterclockwise around w = key.1 - key.0, a pair both
//     traversing +w has its CCW-later member in the union; both traversing
//     -w, its CCW-earlier member. (Derivation in the module tests.)
// Exactly coincident pieces (coplanar overlaps): opposite orientation pairs
// are discarded outright (regularization, §7.1); same-orientation pairs get
// one ∪ and one ∩ so each shared region survives exactly once per output.

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

impl Tag {
    fn flip(self) -> Tag {
        match self {
            Tag::Union => Tag::Inter,
            Tag::Inter => Tag::Union,
        }
    }
}

/// Per-piece classification result. `None` tag = not adjacent to any ring
/// (resolved later by flood fill + ray shooting).
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

/// Classify every piece adjacent to an intersection ring.
pub fn classify_rings(graph: &IntersectionGraph) -> Classification {
    let n = graph.pieces.len();
    let mut tags: Vec<Option<Tag>> = vec![None; n];
    let mut discarded = vec![false; n];

    // Ring construction: intersection edge key → incident pieces.
    let mut rings: BTreeMap<EdgeKey, Vec<Incident>> = BTreeMap::new();
    for (pi, piece) in graph.pieces.iter().enumerate() {
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
        classify_one_ring(key, incidents, &mut tags, &mut discarded);
    }

    Classification { tags, discarded }
}

fn classify_one_ring(
    key: &EdgeKey,
    incidents: &mut [Incident],
    tags: &mut [Option<Tag>],
    discarded: &mut [bool],
) {
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

    // Survivors in CCW order.
    let ring: Vec<usize> = (0..m).filter(|&x| !cancelled[x]).collect();
    if ring.len() < 2 {
        return; // nothing to classify locally
    }

    // Successive-pair relation: coincident → non-orientable iff same
    // traversal (after cancellation, coincident survivors share traversal);
    // distinct angles → non-orientable iff same traversal. So the relation
    // is uniformly "same traversal = non-orientable".
    let coincident = |a: &Incident, b: &Incident| {
        angle_cmp((&a.du, &a.dv), (&b.du, &b.dv)) == Ordering::Equal
    };

    // Seed: first non-coincident successive pair with same traversal
    // (Prop 2). Both forward → CCW-later is Union; both backward →
    // CCW-earlier is Union.
    let mut seed: Option<(usize, Tag)> = None; // (ring position, its tag)
    for s in 0..ring.len() {
        let a = &incidents[ring[s]];
        let b = &incidents[ring[(s + 1) % ring.len()]];
        if a.forward == b.forward && !coincident(a, b) {
            let (pos, tag) = if a.forward {
                ((s + 1) % ring.len(), Tag::Union)
            } else {
                (s, Tag::Union)
            };
            seed = Some((pos, tag));
            break;
        }
    }

    let assign = |pos: usize, tag: Tag, tags: &mut [Option<Tag>]| {
        let piece = incidents[ring[pos]].piece;
        debug_assert!(
            tags[piece].is_none() || tags[piece] == Some(tag),
            "conflicting ring tags for piece {piece}"
        );
        tags[piece] = Some(tag);
    };

    match seed {
        Some((start, start_tag)) => {
            // Walk the full cycle propagating Prop 2/3.
            assign(start, start_tag, tags);
            let mut tag = start_tag;
            for step in 1..ring.len() {
                let prev = &incidents[ring[(start + step - 1) % ring.len()]];
                let cur = &incidents[ring[(start + step) % ring.len()]];
                tag = if prev.forward == cur.forward {
                    tag.flip() // non-orientable successive pair
                } else {
                    tag // orientable pair continues the same solid
                };
                assign((start + step) % ring.len(), tag, tags);
            }
        }
        None => {
            // Undefined configuration (§7.2.3): only orientable pairs and/or
            // coincident same-traversal groups. Tag coincident groups one
            // ∪ / one ∩ (deterministic by radial+id order); leave the rest
            // for flood fill / ray shooting.
            let mut s = 0;
            while s < ring.len() {
                let mut e = s;
                while e + 1 < ring.len()
                    && coincident(&incidents[ring[s]], &incidents[ring[e + 1]])
                {
                    e += 1;
                }
                if e > s {
                    debug_assert!(
                        e - s == 1,
                        "coincident multiplicity > 2 in undefined configuration"
                    );
                    assign(s, Tag::Union, tags);
                    assign(e, Tag::Inter, tags);
                }
                s = e + 1;
            }
        }
    }
}

#[cfg(test)]
#[path = "classify_tests.rs"]
mod tests;
