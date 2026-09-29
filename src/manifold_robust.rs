// Copyright 2026 Lars Brubaker
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// manifold_robust.rs — Manifold entry points backed by the robust engine
// (src/robust) rather than the exact C++ port: the self-intersection query
// that `BooleanEngine::Auto` consults, winding repair of inside-out shells,
// and the full solid rebuild of arbitrary closed triangle soup. None of
// these has a C++ counterpart. A child module of manifold.rs so it can reach
// the private `imp` field.

use super::Manifold;

impl Manifold {
    /// True when two of this mesh's own triangles genuinely intersect —
    /// they cross, they overlap, or they are coincident surface — rather
    /// than merely sharing edges and vertices as every closed mesh does.
    ///
    /// Topologically manifold meshes can still be self-intersecting; those
    /// inputs break the exact boolean engine's assumptions, so
    /// [`crate::types::BooleanEngine::Auto`] routes them to the robust
    /// engine. A mesh carrying non-finite positions (e.g. after a warp to
    /// NaN) answers `true`, that being the safe verdict for geometry no
    /// exact predicate can evaluate.
    ///
    /// The scan is a BVH self-query with an exact narrow phase; the verdict
    /// is cached on the impl, so repeat queries (and the booleans that
    /// consult it) are free until the geometry changes.
    pub fn has_self_intersections(&self) -> bool {
        crate::robust::soup::has_self_intersections(&self.imp)
    }

    /// Repair the winding of inside-out shells so every body reads as solid
    /// material under the robust engine's {winding >= 1} semantics.
    ///
    /// Connected shells whose exact winding shows them inverted relative to
    /// their nesting are rewound: outermost shells end up winding +1 and
    /// cavity shells stay (or become) correctly inward-wound — legitimate
    /// voids are preserved, unlike a blanket flip of negative-signed-volume
    /// shells. Coincident/doubled sheets are deliberately left untouched;
    /// the robust boolean's winding-stack arithmetic already handles them.
    ///
    /// Works standalone (no boolean required) on both manifold and
    /// soup-backed impls; positions, properties, and mesh relations are
    /// untouched, only triangle winding changes. Returns `self` unchanged
    /// when nothing needs flipping.
    pub fn repair_orientation(&self) -> Self {
        if self.is_empty() {
            return self.clone();
        }
        let tris = crate::robust::soup::impl_to_tris(&self.imp);
        let plan = crate::robust::repair::plan_repair(&tris);
        if plan.is_noop() {
            return self.clone();
        }
        let mut out = self.imp.clone();
        crate::robust::repair::apply_flips(&mut out, &plan.flip);
        // Winding-only edit, but it rewrites halfedges in place; re-deriving
        // the verdict keeps the invalidate-on-in-place-edit rule absolute.
        out.invalidate_self_intersects();
        Self::from_impl(out)
    }

    /// Rebuild this mesh into a fresh, properly paired 2-manifold enclosing
    /// the same solid region under `rule`.
    ///
    /// The full robust pipeline — exact intersection (including the mesh
    /// against itself), arrangement, cell complex, winding-number
    /// classification, reassembly — run on this one mesh. Arbitrary triangle
    /// soup is fair game: self-intersections, T-junctions, duplicated or
    /// coincident sheets, more than two faces on an edge, interior walls.
    /// Every wall the winding numbers say has material on both sides
    /// dissolves, every surviving wall is rewound from the cell labels, and
    /// the output is re-imported with real halfedge pairing.
    ///
    /// What is *not* fair game is a surface with a hole in it. Winding numbers
    /// are only defined for a closed surface, and the soup import enforces it:
    /// [`Manifold::from_mesh_gl_robust`] balances directed edges on
    /// position-welded vertices and rejects anything left over with
    /// [`Error::NotClosed`](crate::types::Error::NotClosed), so an open or
    /// non-orientable mesh never reaches this method — it is already an empty
    /// `Manifold` carrying that status,
    /// and the rebuild is a no-op on it. Closed and orientable is the
    /// admission requirement; everything past that the pipeline will fix.
    ///
    /// Choose between this and the cheaper repairs by what is actually wrong:
    ///
    ///  * [`Manifold::repair_orientation`] when only the *winding* is wrong —
    ///    inside-out shells on geometry that is otherwise a clean manifold.
    ///    It touches nothing but triangle orientation, so it is fast, exact,
    ///    and preserves triangle count, properties and relations verbatim.
    ///  * `rebuild_solid` when the *geometry* is wrong — anything that cannot
    ///    be fixed by flipping triangles. It re-triangulates, so vertex and
    ///    triangle counts change and properties are re-interpolated.
    ///
    /// [`crate::types::WindingRule::Positive`] keeps `{w >= 1}`: an inverted body is not
    /// material and disappears. [`crate::types::WindingRule::Nonzero`] keeps `{w != 0}`,
    /// which reads an inside-out body as solid and rewinds it — the right
    /// choice for scans and CAD exports whose shells are wound arbitrarily.
    ///
    /// Empty input returns empty. A cancelled run returns an empty mesh with
    /// [`Error::Cancelled`](crate::types::Error::Cancelled); other pipeline
    /// failures surface through [`Manifold::status`] as usual.
    pub fn rebuild_solid(&self, rule: crate::types::WindingRule) -> Self {
        self.rebuild_solid_with_token(rule, None)
    }

    /// [`Manifold::rebuild_solid`] with cooperative cancellation. Soup
    /// rebuilds are as expensive as a boolean against a partner, so anything
    /// interactive wants this form.
    pub fn rebuild_solid_with_token(
        &self,
        rule: crate::types::WindingRule,
        token: Option<&crate::cancel::CancelToken>,
    ) -> Self {
        if self.is_empty() {
            return self.clone();
        }
        Self::from_impl(crate::robust::rebuild_with_rule(
            &self.imp, rule, token, None,
        ))
    }
}
