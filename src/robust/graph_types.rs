// robust/graph_types.rs — Shared vocabulary of the intersection graph: the
// edge key spaces, the exact-point interner, and the `Piece` /
// `IntersectionGraph` output types.
//
// Split out of robust/intersection_graph.rs, which builds these values;
// robust/cells.rs, robust/pairing.rs, robust/propagate-style flood fills and
// robust/assemble.rs consume them (all through the `intersection_graph`
// re-exports, so the public paths are unchanged). The exact rational point
// type lives in robust/exact/rational.rs.

// Fx hashing instead of SipHash. Every map/set here is probe-only (documented
// per site); the hasher is unseeded, so even iteration order is stable across
// runs — output cannot depend on it.
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};

use crate::linalg::Vec3;

use super::exact::rational::{r3_eq, R3, R3Key};

/// Canonical (sorted) edge between two interned vertex ids. Downstream
/// stages (classify rings, propagate flood fill) key their maps on these
/// integers instead of exact rational point pairs — vertex interning at
/// piece-emission time makes id equality coincide with exact geometric
/// identity.
pub type EdgeKey = (u32, u32);

pub fn edge_key(a: u32, b: u32) -> EdgeKey {
    if a <= b {
        (a, b)
    } else {
        (b, a)
    }
}

/// Canonical (lexicographically sorted) exact edge between two points —
/// local key for the split-point registries built before interning exists.
/// Wrapped in [`R3Key`] so the registry hashes structurally instead of
/// comparing rationals per probe.
pub(super) type GeoEdgeKey = (R3Key, R3Key);

pub(super) fn geo_edge_key(a: &R3, b: &R3) -> GeoEdgeKey {
    if a <= b {
        (R3Key(a.clone()), R3Key(b.clone()))
    } else {
        (R3Key(b.clone()), R3Key(a.clone()))
    }
}

/// Canonical original-mesh edge keyed by raw coordinate bits — original
/// edges always join exact f64 vertices, so the boundary-split registry
/// never needs rational keys (and untouched triangles probe it for free).
pub(super) type BitEdgeKey = ([u64; 3], [u64; 3]);

pub(super) fn bit_edge_key(a: Vec3, b: Vec3) -> BitEdgeKey {
    let (ka, kb) = (f64_key(a), f64_key(b));
    if ka <= kb {
        (ka, kb)
    } else {
        (kb, ka)
    }
}

/// One output fragment: a sub-triangle of an arranged input triangle, or an
/// untouched whole triangle. `v` is wound to match the input mesh's outward
/// orientation; `vi` are the interned ids of the same three vertices.
#[derive(Clone, Copy, Debug)]
pub struct Piece {
    /// 0 = first operand (P), 1 = second operand (Q).
    pub mesh: u8,
    /// Index of the originating triangle in its soup.
    pub tri: usize,
    /// Interned vertex ids (indices into `IntersectionGraph::verts`), wound
    /// to the input mesh's outward orientation. Pieces carry no coordinates
    /// of their own — the shared tables keep untouched triangles free of
    /// rational clones entirely.
    pub vi: [u32; 3],
}

/// Everything classification and assembly need.
pub struct IntersectionGraph {
    pub pieces: Vec<Piece>,
    /// Interned unique vertices; `Piece::vi` and `EdgeKey` index into this.
    pub verts: Vec<R3>,
    /// Correctly rounded f64 approximation per interned vertex (exact for
    /// input vertices) — float filters and output assembly read these
    /// instead of re-rounding rationals.
    pub verts_f64: Vec<Vec3>,
    /// Canonical keys of every arrangement constraint edge — the exact
    /// intersection sub-segments the classification rings live on.
    pub isect_edges: HashSet<EdgeKey>,
    /// True when any P×Q pair intersected at all.
    pub any_intersections: bool,
}

impl IntersectionGraph {
    /// The three exact vertices of a piece.
    pub fn piece_verts(&self, pi: usize) -> [&R3; 3] {
        let vi = self.pieces[pi].vi;
        [
            &self.verts[vi[0] as usize],
            &self.verts[vi[1] as usize],
            &self.verts[vi[2] as usize],
        ]
    }
}

/// Exact-point interner: one id per distinct point, with two disjoint key
/// spaces. f64-representable points (all input vertices, and any constructed
/// point that rounds exactly) key on their coordinate bits — no rational
/// hashing, so untouched input triangles intern for the cost of a HashMap
/// probe. Only genuinely non-representable constructed points use the
/// rational map. `verts_f64` caches the correctly rounded approximation of
/// every id (exact for bit-keyed points), which downstream float filters
/// and output assembly reuse instead of re-rounding.
/// Order invariance: both maps are probe-only (`get`/`entry`, never
/// iterated); ids come from `verts.len()` at insertion time, so they depend
/// only on the sequential call order, not on the hasher.
#[derive(Default)]
pub struct VertInterner {
    map: HashMap<R3Key, u32>,
    fmap: HashMap<[u64; 3], u32>,
    pub verts: Vec<R3>,
    pub verts_f64: Vec<Vec3>,
}

pub(super) fn f64_key(v: Vec3) -> [u64; 3] {
    // Normalize -0.0 so it shares an id with +0.0 (they are the same
    // rational point).
    let norm = |x: f64| if x == 0.0 { 0.0f64 } else { x }.to_bits();
    [norm(v.x), norm(v.y), norm(v.z)]
}

impl VertInterner {
    /// Intern an exact-f64 point (input mesh vertices): zero rational work
    /// on hits; one `R3::from_vec3` on first sight, for the exact table.
    pub fn intern_f64(&mut self, v: Vec3) -> u32 {
        let key = f64_key(v);
        if let Some(&id) = self.fmap.get(&key) {
            return id;
        }
        let id = self.verts.len() as u32;
        self.fmap.insert(key, id);
        self.verts.push(R3::from_vec3(v));
        self.verts_f64.push(v);
        id
    }

    /// Intern an exact rational point. Representable points route to the
    /// f64 key space so both paths agree on ids.
    pub fn intern(&mut self, p: &R3) -> u32 {
        let rounded = p.to_vec3_rounded();
        if r3_eq(&R3::from_vec3(rounded), p) {
            return self.intern_f64(rounded);
        }
        let next = self.verts.len() as u32;
        match self.map.entry(R3Key(p.clone())) {
            std::collections::hash_map::Entry::Occupied(e) => *e.get(),
            std::collections::hash_map::Entry::Vacant(e) => {
                e.insert(next);
                self.verts.push(p.clone());
                self.verts_f64.push(rounded);
                next
            }
        }
    }
}
