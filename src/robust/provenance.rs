// robust/provenance.rs — Which input face each robust-engine output triangle
// lies on, in the fast engine's terms (`TriRef` + `Relation`).
//
// Every output triangle of the robust pipeline is a piece of one input
// triangle (robust/graph_types.rs `Piece::mesh`/`Piece::tri`): the
// intersection graph splits input triangles, each arrangement sub-triangle
// keeps its parent's index, cell extraction keeps the representative piece's
// index for a coincident stack, and assembly (robust/assemble.rs) hands the
// resulting per-triangle `TriRef`s to the import, so simplification merges
// only triangles of one source face (`TriRef::same_face`) and the export
// writes each triangle's run original ID, face ID, run transform and run
// flags exactly as it does for a fast-engine result.
//
// The mesh IDs follow the fast engine's `update_reference`
// (boolean_result_assemble.rs): operand P keeps its mesh IDs, operand Q's are
// shifted past a freshly reserved block, and a subtracted Q's relations are
// marked backside. The robust engine winds its output from the cell labels
// rather than from the input faces, so a triangle can also come out reversed
// where the fast engine has no counterpart (an inside-out region kept under
// the nonzero rule, or a mis-wound shell rebuilt): such a triangle moves to a
// sibling mesh ID whose relation differs only in its backside flag (run_flags
// bit 0), reserved in ascending base-ID order, so the export still tells a
// consumer exactly which triangles are flipped.
//
// A Rust extension: the robust engine is itself Rust-only, so nothing in the
// C++ constrains these values; the rules follow the fast engine's
// `update_reference` where it has a counterpart. Where they differ from it:
//
// - **Which operand wins a coincident face.** A coincident stack (coplanar
//   faces overlapping, from either operand) emits one triangle with the
//   provenance of its lowest-index piece, and the graph emits all of P's
//   pieces before Q's, so where a P face and a Q face coincide and the wall
//   survives, P's face wins for every operation. The fast engine agrees for
//   intersection (and a surviving coincident face in a difference is P's in
//   both); for a union of same-facing coplanar faces its choice follows its
//   symbolic perturbation and varies with the face's direction and position,
//   even within one face, so there the engines can name different operands.
//   Both name a face the triangle lies on.
// - **Coplanar groups.** The assembled import
//   (`Manifold::from_mesh_gl64_robust_assembled`) keeps the supplied coplanar
//   IDs instead of its own flood fill, as the fast engine keeps its inputs',
//   so simplification never merges across a seam between source faces, and
//   robust results keep a few more triangles where coplanar faces of
//   different sources meet than a single fresh run would.
// - **Soup concatenation.** `compose_meshes`' soup branch
//   (`robust::assemble_all`) keeps every triangle's reference, renumbering mesh
//   IDs and coplanar groups as its strict branch does.
//
// Pinned by `manifold::tests::face_id_provenance_robust` (single, coplanar,
// batch and nested cases under `BooleanEngine::Robust`; `Auto` with a soup and
// a self-intersecting operand; backside runs for a subtracted operand; a
// fast/robust colour comparison by point sampling) and by the re-frozen hash in
// `robust::coplanar_cross_copy_tests`. manifold-sharp's robust results must
// carry the same references to produce that hash.

use std::collections::{BTreeMap, BTreeSet};

use crate::impl_mesh::{reserve_ids, ManifoldImpl};
use crate::types::{Relation, TriRef};

use super::intersection_graph::Piece;

/// Output-ready provenance of every input triangle, per operand slot.
pub struct RefCtx {
    /// `tri_ref[m][t]`: the `TriRef` an output triangle lying on triangle `t`
    /// of operand `m` carries, its mesh ID already moved into the output's
    /// ID space.
    pub tri_ref: [Vec<TriRef>; 2],
    /// The output's relation per mesh ID (before backside siblings).
    pub relations: BTreeMap<i32, Relation>,
    /// Whether an operand's triangles are expected to come out reversed (a
    /// subtracted Q), which its relations already record as backside.
    pub expect_flip: [bool; 2],
}

/// An operand's triangle references and relations. An impl built without
/// references (none of the library's constructors or imports do that, but a
/// hand-built impl can) gets one fresh mesh ID covering every triangle, with
/// each triangle its own coplanar group — what the import gives a mesh with
/// no runs.
fn operand_refs(imp: &ManifoldImpl) -> (Vec<TriRef>, BTreeMap<i32, Relation>) {
    if imp.mesh_relation.tri_ref.len() == imp.num_tri() {
        return (
            imp.mesh_relation.tri_ref.clone(),
            imp.mesh_relation.mesh_id_transform.clone(),
        );
    }
    let id = reserve_ids(1) as i32;
    let refs = (0..imp.num_tri())
        .map(|t| TriRef {
            mesh_id: id,
            original_id: id,
            face_id: -1,
            coplanar_id: t as i32,
        })
        .collect();
    let rel = Relation {
        original_id: id,
        ..Relation::default()
    };
    (refs, BTreeMap::from([(id, rel)]))
}

impl RefCtx {
    /// Provenance for a two-operand boolean; `invert_b` for a subtraction,
    /// whose second operand's faces come out reversed.
    pub fn for_operands(a: &ManifoldImpl, b: &ManifoldImpl, invert_b: bool) -> Self {
        let (ref_a, mut relations) = operand_refs(a);
        let (mut ref_b, rel_b) = operand_refs(b);
        let offset = reserve_ids(rel_b.len() as u32) as i32;
        for r in &mut ref_b {
            r.mesh_id += offset;
        }
        for (k, mut rel) in rel_b {
            rel.back_side ^= invert_b;
            relations.insert(k + offset, rel);
        }
        RefCtx {
            tri_ref: [ref_a, ref_b],
            relations,
            expect_flip: [false, invert_b],
        }
    }

    /// Provenance for the single-operand rebuild: the operand in slot 0,
    /// slot 1 empty.
    pub fn for_single(a: &ManifoldImpl) -> Self {
        let (ref_a, relations) = operand_refs(a);
        RefCtx {
            tri_ref: [ref_a, Vec::new()],
            relations,
            expect_flip: [false, false],
        }
    }

    /// Provenance for concatenating `meshes` into slot 0 in order (the soup
    /// branch of `boolean3::compose_meshes`). Like the strict branch there,
    /// each input's mesh IDs become local IDs 1, 2, 3, … in input order and
    /// its coplanar groups are offset by the triangles before it; the caller
    /// then hands out fresh IDs with `increment_mesh_ids`.
    pub fn for_compose(meshes: &[ManifoldImpl]) -> Self {
        let mut tri_ref = Vec::new();
        let mut relations = BTreeMap::new();
        let mut next_local = 1i32;
        let mut tri_offset = 0i32;
        for mesh in meshes {
            let (refs, rels) = operand_refs(mesh);
            let mut local = BTreeMap::new();
            for (id, rel) in rels {
                local.insert(id, next_local);
                relations.insert(next_local, rel);
                next_local += 1;
            }
            tri_ref.extend(refs.into_iter().map(|r| TriRef {
                mesh_id: local.get(&r.mesh_id).copied().unwrap_or(r.mesh_id),
                coplanar_id: r.coplanar_id + tri_offset,
                ..r
            }));
            tri_offset += mesh.num_tri() as i32;
        }
        RefCtx {
            tri_ref: [tri_ref, Vec::new()],
            relations,
            expect_flip: [false, false],
        }
    }

    /// The `TriRef` of each output piece (in order) and the output's
    /// relations. A piece whose winding against its source triangle differs
    /// from its operand's expectation goes to a backside sibling of its mesh
    /// ID; siblings are reserved in ascending base-ID order, so the result
    /// does not depend on piece order beyond which IDs are needed.
    pub fn output_refs(&self, pieces: &[&Piece]) -> (Vec<TriRef>, BTreeMap<i32, Relation>) {
        let deviates = |p: &Piece| p.flipped != self.expect_flip[p.mesh as usize];
        let source = |p: &Piece| self.tri_ref[p.mesh as usize][p.tri];
        let bases: BTreeSet<i32> = pieces
            .iter()
            .filter(|p| deviates(p))
            .map(|p| source(p).mesh_id)
            .collect();
        let mut relations = self.relations.clone();
        let mut sibling: BTreeMap<i32, i32> = BTreeMap::new();
        if !bases.is_empty() {
            let first = reserve_ids(bases.len() as u32) as i32;
            for (id, base) in (first..).zip(bases) {
                let mut rel = relations.get(&base).cloned().unwrap_or_default();
                rel.back_side = !rel.back_side;
                relations.insert(id, rel);
                sibling.insert(base, id);
            }
        }
        let refs = pieces
            .iter()
            .map(|p| {
                let mut r = source(p);
                if deviates(p) {
                    r.mesh_id = sibling[&r.mesh_id];
                }
                r
            })
            .collect();
        (refs, relations)
    }
}
