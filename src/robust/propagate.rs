// robust/propagate.rs — Component flood fill across each mesh's surface
// (paper §7.3, "union/intersection mesh browse").
//
// Pieces of the same mesh that share a non-constraint edge belong to the
// same boolean output: the surface is only ever "cut" along intersection
// segments (cross-mesh and self-intersection alike). Union-find
// (crate::disjoint_sets) groups pieces into components bounded by those
// cuts; a component takes the tag of any coincident-bound member
// (robust/classify.rs), and every other component is reported untagged for
// the per-component winding classification in robust/mod.rs.
//
// Edge identity is exact geometry — soups have no trustworthy connectivity;
// robust/intersection_graph.rs interns every distinct rational vertex, so
// integer id pairs are exact geometric edge identity here.

use std::collections::{BTreeMap, HashMap};

use crate::disjoint_sets::DisjointSets;

use super::classify::{Classification, Tag};
use super::intersection_graph::{edge_key, EdgeKey, IntersectionGraph};

/// Result of propagation: per-piece tags (None only for pieces of untagged
/// components), plus one representative piece per still-untagged component
/// for the ray-shooting fallback.
pub struct Propagation {
    pub tags: Vec<Option<Tag>>,
    /// (component root, representative piece) for every untagged component,
    /// discarded pieces excluded.
    pub untagged: Vec<(usize, usize)>,
    /// Component root per piece (valid where not discarded).
    pub component: Vec<usize>,
}

pub fn propagate(graph: &IntersectionGraph, cls: &Classification) -> Propagation {
    let n = graph.pieces.len();
    let ds = DisjointSets::new(n as u32);

    // Group pieces by non-intersection shared edges, per mesh.
    let mut edge_owner: HashMap<(u8, EdgeKey), Vec<usize>> = HashMap::new();
    for (pi, piece) in graph.pieces.iter().enumerate() {
        if cls.discarded[pi] {
            continue;
        }
        for e in 0..3 {
            let key = edge_key(piece.vi[e], piece.vi[(e + 1) % 3]);
            if graph.isect_edges.contains(&key) {
                continue; // never flood across a cut
            }
            edge_owner.entry((piece.mesh, key)).or_default().push(pi);
        }
    }
    for owners in edge_owner.values() {
        for w in owners.windows(2) {
            ds.unite(w[0] as u32, w[1] as u32);
        }
    }

    // Merge ring tags into components.
    let mut comp_tag: BTreeMap<usize, Tag> = BTreeMap::new();
    for pi in 0..n {
        if cls.discarded[pi] {
            continue;
        }
        if let Some(tag) = cls.tags[pi] {
            let root = ds.find(pi as u32) as usize;
            let prev = comp_tag.insert(root, tag);
            debug_assert!(
                prev.is_none() || prev == Some(tag),
                "conflicting tags within one surface component"
            );
        }
    }

    let mut tags: Vec<Option<Tag>> = vec![None; n];
    let mut component = vec![0usize; n];
    let mut untagged_roots: BTreeMap<usize, usize> = BTreeMap::new();
    for pi in 0..n {
        if cls.discarded[pi] {
            continue;
        }
        let root = ds.find(pi as u32) as usize;
        component[pi] = root;
        match comp_tag.get(&root) {
            Some(tag) => tags[pi] = Some(*tag),
            None => {
                untagged_roots.entry(root).or_insert(pi);
            }
        }
    }

    Propagation {
        tags,
        untagged: untagged_roots.into_iter().collect(),
        component,
    }
}
