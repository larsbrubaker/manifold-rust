# Deliberate divergences from the C++ reference

Per CLAUDE.md, byte-matching `cpp-reference/manifold` is the default; the
entries below are the deliberate exceptions, each with the evidence that
justified it. Trace-diff debugging against the C++ must expect these.

Two kinds of entry live here, and they are not the same claim:

- **Justified divergence** — an accuracy fix, a real bug fix in the C++, or a
  measured improvement. CLAUDE.md's three qualifying reasons. Entries 1, 3 and 4.
- **Inherited divergence, documented and scheduled for coordinated
  harmonization** — an output shape this port already shipped, which we would
  resolve toward the C++ on the merits but cannot change unilaterally, because a
  downstream consumer verifies against this tree bit-for-bit. These are
  *disclosures*, not justifications: the entry states what differs, why it
  cannot be fixed here alone, and what a coordinated fix would take. Entries 2,
  5 and 6.

The second kind is deliberately uncomfortable to write, which is the point — it
is a debt with a name attached, not a decision that ends the discussion. Nothing
belongs in either category for convenience.

## 1. Robust-engine outputs skip `swap_degenerates` (2026-08-08)

**What differs:** meshes assembled by the robust boolean engine
(`src/robust/assemble.rs`) do not run `edge_op::swap_degenerates` during
their import/simplification. The exact engine's pipeline is unchanged and
still byte-matches C++ (`Impl::SimplifyTopology`).

**Where:** `Manifold::from_mesh_gl64_robust_assembled`
(`src/manifold_meshgl.rs`, `pub(crate)`, used only by `robust::assemble`)
runs `cleanup_topology` + `collapse_short_edges` + `calculate_vert_normals`
in place of `edge_op::remove_degenerates`, and `assemble` composes the same
pieces plus `collapse_colinear_edges` in place of `simplify_topology`.
`edge_op`'s own functions are untouched; only the composition differs. Every
other import path, including user meshes arriving through
`from_mesh_gl_robust` / `from_mesh_gl64_robust`, is byte-identical to before.

**Why:** `face_op::set_normals_and_coplanar` (faithful port of C++
`Impl::SetNormalsAndCoplanar`, impl.cpp:214) flood-fills a seed triangle's
normal onto every coplanar neighbor with no orientation check. A boolean
result legitimately contains coplanar *antiparallel* adjacencies (material
below the plane on one side of a shared edge, above it on the other), so a
handful of large triangles receive sign-flipped normals. `swap_degenerates`
then misclassifies them as degenerate (`ccw <= 0`) and swaps large
non-planar quads — physically moving the surface. On Thingi10K #301921 ∪
rotated-self this moved the robust union volume by −2.5e-3 (and the exact
engine's own result by +5.6e-5; both engines' outputs are exposed, exact
merely got tessellation-lucky). Skipping the swap keeps robust volumes
extraction-exact. The cost is slightly larger tri/vert-count drift versus
the exact engine, which the sweep gate already treats as advisory.

**Evidence:** `robust::thingi_tests::thingi_301921_union_rotated_self_matches_exact`,
the staged volume traces in the session log (extraction 0.438010 preserved
through cleanup, lost only in swap), and the Monte-Carlo referee
(`examples/volume_referee.rs`).

## 2. The centered cylinder centers in place, not through `Translate().AsOriginal()` (2026-08-30)

**What differs:** `constructors::cylinder`'s `center` branch shifts `vert_pos.z`
in place and repairs the derived caches. The C++ re-centers by composing a
transform, which additionally assigns a fresh original ID, re-marks coplanar
faces, and recomputes every coordinate. Unlike entry 1 this is not an accuracy
fix — it is an **inherited** output shape being deliberately preserved, and the
one entry here that we would resolve toward the C++ if we could do it alone.

**Where:** `src/constructors.rs`, the `if center` block at the end of `cylinder`.
Reached by `Manifold::cylinder_centered(…, center: true)` and — via the recursive
`cylinder(height, radius_high, 0.0, n, true)` its cone branch starts from — every
cone, centered or not. `Manifold::cylinder` does *not* reach it:
`manifold_shape.rs:48` forwards to `cylinder_centered` with `center: false`.
C++ v3.5.2 is `src/constructors.cpp:155-157`:

```cpp
Manifold cylinder = Manifold::Extrude({circle}, height, 0, 0.0, vec2(scale));
if (center)
  cylinder = cylinder.Translate(vec3(0.0, 0.0, -height / 2.0)).AsOriginal();
return cylinder;
```

`AsOriginal` (`src/manifold.cpp:449-457`) copies the impl, then runs
`InitializeOriginal()` and `SetNormalsAndCoplanar()`.

**What is observably different.** Measured on `cylinder(4.0, 1.0, 1.0, 8, true)`
against the same mesh built uncentered and put through `transform`:

| | in place (ours) | through a transform (C++'s shape) |
|---|---|---|
| `mesh_relation.original_id` | `1` — extrude's | `-1` from the transform, then a fresh ID from `InitializeOriginal` |
| `epsilon` | `4e-12` | `3.999999999999999e-12` |
| `tolerance` | `4e-12` | `4e-12` |
| positions | 4 coordinates are `-0.0` | the same 4 are `+0.0`; **zero** value differences |

The signed zeros are on x and y, never z, and the mechanism is exactly that: the
in-place edit touches only z, so x and y keep the `-0.0` that `cosd`/`sind`
produced at the quarter angles, where a transform recomputes them as
`1.0*x + 0.0*y + 0.0*z + 0.0` and normalizes `-0.0` to `+0.0`. The epsilon gap is
one ULP, from `Impl::Transform` scaling epsilon by a spectral norm that comes back
`0.9999999999999998` rather than exactly `1.0` for a translation. Beyond the
table, `AsOriginal`'s `SetNormalsAndCoplanar()` re-marks coplanar faces, which we
do not re-run.

**Why we keep it.** This divergence **predates** the cache repair that sits at
the same lines: the in-place shim shipped without `AsOriginal`'s semantics, so
`cylinder_centered`'s originalID, epsilon and signed zeros have differed from the
C++ in every released 0.14.x. The repair deliberately fixed only the stale
caches — the face BVH was describing the pre-shift positions, and every boolean
against a centered cylinder tripped `pair_up`'s non-manifold assert — while
retaining the existing output semantics exactly. `sort_geometry` was the right
size of claim for that: it rebuilds the derived caches from the positions already
there, where re-centering through a transform would also have moved the positions.
`set_epsilon` is deliberately not re-run for the same reason, and that is what
keeps epsilon at `4e-12`.

Switching to the C++ form is therefore not a free correction but a breaking
behavioral change, and it has a downstream consumer that would break:
[manifold-sharp](https://github.com/larsbrubaker/manifold-sharp) is a pure C#
port of this crate whose test contract is bit-exactness against this tree. It
transcribes the same in-place centering (`ManifoldSharp/Constructors.cs`, the
`if (center)` block of its cylinder), and the "Reference and oracle" section of
its `CLAUDE.md` names `fa18cc5`, this tree's cache repair, as the floor for
bit-agreement on this path. Its ledger once carried both stale-cache defects (the centered cylinder and `SubdivideImpl`) as its
own entries 4 and 5; those were retired in manifold-sharp `13e0a87` when
`fa18cc5` landed here, and its entries 4 and 5 today are unrelated (a progress
`Phase`, a convex erosion). Its verification nets compare with no slack, so a
signed zero moved on one side alone is a failure there.

**The two halves of `cylinder` disagree today, and we say so plainly.** The cone
branch a few lines above *does* finish with `cone.initialize_original()` and
`set_normals_and_coplanar(&mut cone)` — it models `AsOriginal` faithfully,
because it was written from the C++'s Mirror/Translate/AsOriginal chain. So a
cone reports a fresh original ID and a centered cylinder reports `extrude`'s.
That inconsistency is real, it is not defensible on its merits, and resolving it
belongs to a harmonization pass coordinated with manifold-sharp so both trees
move together. Until then, do not "fix" one half in isolation.

**Evidence:** the table above, measured with a scratch probe comparing the two
constructions coordinate by coordinate on raw f64 bits (4 signed-zero
differences, 0 value differences); the C++ sources cited above;
`constructors::tests::centered_cylinder_collider_matches_its_vertex_positions`
and `centered_cylinder_is_usable_in_a_boolean`, which pin the cache repair that
this entry's semantics were preserved *through*; and manifold-sharp's
`ConstructorsTests` (`CenteredCylinderColliderMatchesItsVertexPositions`,
`CenteredCylinderIsUsableInABoolean`, `CenteredConeIsUsableInABoolean`), which
pin the same repair on its transcription. manifold-sharp's oracle lane (34/34 at
`13e0a87`) is *not* evidence for this entry: it consumes the published
`ManifoldRust` 0.5.0 NuGet (`ManifoldSharp.OracleTests.csproj:22`), whose natives
were built from `43f377f` — before `fa18cc5` — and, per that repo's `CLAUDE.md`,
no oracle row exercises the centered cylinder on the native side.

## 3. `dedupe_edges` skips duplicate entries an earlier repair already resolved (2026-09-26)

**What differs:** `edge_op::dedupe_edges` checks each collected duplicate is
still a duplicate (`is_still_duplicated`: another halfedge leaving its start
vertex still ends at its end vertex) before calling `dedupe_edge` on it. The C++
(`src/edge_op.cpp`, `DedupeEdges`) collects every duplicated edge once per pass
and repairs each entry with no re-check.

**Why:** a real bug in the C++. An earlier repair in the same pass rewires the
neighbourhood, so a later entry can be stale; repairing it copies the position
of a vertex that is not the orbit's own into the new vertex it relabels the
orbit to, so triangle corners jump and solid disappears — an exact union of two
`NoError` operands lost 2.3e-5 of one of them, reported `NoError`. The outer
loop re-collects each pass, so skipping a stale entry loses nothing a later pass
would not catch.

**What is observably different:** Thingi10K #1147177 through the demo import
goes from 3201 verts / 6418 tris / genus 5 / volume 0.047535 to 3206 / 6424 /
genus 4 / 0.047563, with no corner moved. Thingi10K #939888 goes from 860 / 1716
to 861 / 1718. manifold-sharp made the same fix in its commit `dd55f9f`
(`ManifoldSharp/EdgeOp.Dedupe.cs`, pinned by `DedupeEdgesRegressionTests`) and
pins the same counts, so the two ports agree. It is not a ledger entry there,
since it does not diverge from this tree; its `docs/RUST_DIVERGENCES.md` mentions
it only in a parenthetical under entry 6, as "not a divergence".

**Evidence:** `edge_op::tests::test_dedupe_edges_never_moves_a_triangle_corner`
on `src/testdata/dedupe-stale-duplicate.txt` (the 852-triangle fixture shared
with manifold-sharp's `DedupeEdgesRegressionTests`): 16 corners moved before the
fix, 0 after; `robust::thingi_tests::thingi_1147177_import_counts` and
`thingi_939888_import_counts` pin the new counts.

## 4. Mirroring keeps each property with its corner (2026-09-28)

**What differs:** when `ManifoldImpl::transform` (`src/impl_mesh.rs`) flips
triangle winding for a negative-determinant transform, each halfedge's
`prop_vert` is reassigned so the new corners take the props of old corners
(0, 2, 1) — the old corner whose start vertex becomes the new start vertex. The
pinned C++ v3.5.2 `FlipTris` (`src/mesh_fixes.h:51-67`) swaps halfedges `3t` and
`3t+2` whole, so the props travel with the halfedges and new corner `i` takes the
prop of old corner `2-i`: corners 0 and 2 exchange properties while their
vertices do not.

**Why:** a real bug in the C++, fixed upstream on master in 422ab6fc ("Fixes
#1781 - Manifold::Transform() does not correctly handle mirror transforms",
upstream issue #1781). We take upstream's change verbatim — the same
`{Prop(3t), Prop(3t+2), Prop(3t+1)}` ordering — ahead of the pinned submodule.
With the bug, a mirrored mesh with properties has normals/UVs on the wrong
corners, and `get_mesh_gl` splits vertices that should be shared. The C++
`CsgLeafNode::Compose` path (`src/csg_tree.cpp`) also calls `FlipTris`; ours
composes leaves through the same `ManifoldImpl::transform`, so the one fix
covers both. This entry retires once the submodule pin moves past 422ab6fc.

**What is observably different:** a mirrored mesh with properties exports
differently: `Manifold::sphere(1.0, 32).calculate_normals(0, 180.0)` has 258
`MeshGL` verts; mirrored over x it emitted 1536 before the fix (every triangle
corner split off with a misplaced normal) and emits 258 after. A mirrored mesh
*without* properties exports identically, because `get_mesh_gl` ignores
`prop_vert` when `num_prop == 0`. But the old flip still broke the
`prop_vert == start_vert` identity on such meshes, and three property-adding
readers index by `prop_vert` even when `num_prop == 0`:
`Manifold::set_properties` (`src/manifold.rs`), `calculate_curvature`
(`src/properties.rs`) and `calculate_normals` (`src/smoothing.rs`). So any of
those applied after a mirror now sees corrected corners, including meshes the
user never mirrored: `Manifold::cylinder(2.0, 0.0, 1.0, 16)` (apex-bottom cone)
mirrors internally (`src/constructors.rs`). Measured before -> after: that cone
then `set_properties(1, ..)` exports 90 -> 17 verts; the cone then
`calculate_normals(0, 60.0)` 90 -> 33; a mirrored `cube` then `set_properties`
36 -> 8; a mirrored `sphere(1.0, 32)` then `set_properties`,
`calculate_curvature` or `calculate_normals` 1536 -> 258. All of these outputs
now diverge from pinned v3.5.2, whose `FlipTris` has the old behavior.
Positive-determinant transforms are unchanged.

**Evidence:** `manifold::tests::normals::test_cpp_mirrored_normals`, a port of
upstream's `TEST(Manifold, MirroredNormals)`: MeshGL vertex count unchanged by
the mirror and every normal outward-pointing. Plus, in the same file,
`test_mirrored_cone_set_properties_shares_verts` (17 verts),
`test_mirrored_cube_set_properties_shares_verts` (8 verts) and
`test_mirrored_cone_calculate_normals_vert_count` (33 verts). All four fail
with the prop capture removed (1536, 90, 36, 90 verts) and pass with it.
manifold-sharp carries the same fix in its commit 674b6ce.

## 5. `CrossSection::decompose` groups holes by bounding box, not by a `PolyTree` (2026-09-28)

**What differs:** `CrossSection::decompose` (`src/cross_section.rs:286-348`)
normalizes through `union` with an empty section, calls every contour with
non-negative signed area an outline, and gives each hole to the outline with the
smallest bounding box containing the hole's *first vertex*. C++ v3.5.2
`CrossSection::Decompose` (`src/cross_section/cross_section.cpp:475-494`, with
`decompose_outline` / `decompose_hole` at 126-151) runs
`C2::BooleanOp(Union, FillRule::Positive, …)` into a `C2::PolyTreeD`, whose
parent/child links are Clipper's own containment result, and emits one section
per outline node with exactly that node's children as holes. This is not an
accuracy fix or a bug fix on our side — the bounding-box heuristic is a
simplification the port shipped, and it is the less correct of the two. It is an
**inherited** divergence: debt, not a decision.

**What is observably different:**

- *Hole ownership.* A bounding box is not containment. When a hole's first vertex
  also falls in a smaller outline's box, the hole goes to the wrong component.
  Measured with a scratch probe: a bar `[0,10]×[0,2]` with a hole `[8,9]×[0.5,1.5]`,
  unioned with a U-shaped outline (bbox `[7,11]×[-0.5,2.5]`, area 12) whose
  opening embraces the bar's right end. `decompose` returns the U *carrying the
  bar's hole* and the bar as a solid rectangle with no hole; the PolyTree puts the
  hole under the bar. (`compose` of the two components still restores the input,
  because the union re-derives winding from all contours together — the
  components themselves are wrong.) Islands nested inside holes are separate
  outlines in both implementations, so they are not the failure mode on their
  own; any hole whose first vertex a smaller, unrelated outline's box covers is.
- *Component order.* C++ emits the reversed post-order of the tree walk (islands
  inside a node's holes are pushed before the node, siblings after); ours follows
  the path order Clipper's flat union output happens to have. The two existing
  tests (`test_cpp_cross_section_decompose` in `manifold_tests/cross_section2.rs`
  and `manifold_tests/advanced.rs`) check only counts, which agree.
- *Short-circuit.* C++ returns a copy of `*this` unchanged when
  `NumContour() < 2` — so an empty section decomposes to one empty section, and a
  single contour is not re-normalized. Ours always normalizes, and returns an empty
  `Vec` for an empty section (probe: `CrossSection::default().decompose().len()`
  is `0`).

**Why it stays for now.** manifold-sharp transcribes the same heuristic
(`ManifoldSharp/CrossSection.cs`, `Decompose`, with the same comments) and
verifies bit-for-bit against this tree, so fixing it here alone breaks its
parity.

**Harmonization path:** port the PolyTree grouping in both trees together.
`clipper2-rust` 1.0.3, already a dependency, exposes `boolean_op_tree_d` and
`PolyTreeD`, so the C++ `decompose_outline` / `decompose_hole` recursion ports
directly — including the `NumContour() < 2` short-circuit and the reversed
emission order — with no new dependency. The shared regression should be the
U-shape case above, asserting the hole stays with the bar.

**Evidence:** source reading of both implementations at the lines above, and the
scratch probe described (not checked in; it is the U-shape construction from
`CrossSection::from_polygons_fill` rectangles and `difference`/`union`).

## 6. `MeshGL::merge` dedupes open edges and open vertices (2026-09-28)

**What differs:** `MeshGLP<f32, u32>::merge` (`src/types_meshgl.rs:160-267`)
collects open halfedges into a `BTreeSet<(usize, usize)>` (175-191) and then
dedupes their start vertices through a second `BTreeSet` (196-206). C++ v3.5.2
`MergeMeshGLP` (`src/sort.cpp:62-98`) uses a `std::multiset<std::pair<int,int>>`,
erases one matching copy per reverse halfedge, and builds `openVerts` with one
entry per remaining open edge, duplicates kept. Like entry 5 this is neither an
accuracy fix nor a bug fix; it is an **inherited** simplification.

**What is observably different.** Two cases, both requiring input that is not
already a clean manifold (which is exactly what `merge` exists to repair):

- *Duplicate same-direction halfedges* change which edges are open, so they can
  change which vertices merge. If a halfedge `s→e` occurs twice before its reverse
  `e→s` arrives, the multiset holds two copies and the reverse erases one, leaving
  `s→e` open; the set holds one and the reverse leaves nothing. Whether a
  divergence appears depends on triangle order: if the reverse arrives between the
  copies, both implementations leave one open. Scratch probe: a tetrahedron whose
  face `(0,2,1)` is listed twice first, plus a separate open triangle `(4,5,6)`
  with vertex 4 coincident with vertex 0. Ours returns `true` with no merges (only
  4, 5, 6 are open); by the C++ source, vertices 0, 1, 2 are also open and 4
  merges to 0 (`mergeFromVert = [4]`, `mergeToVert = [0]`). No C++ build was run
  for this; the C++ result is traced from source.
- *A vertex that starts two or more open halfedges* — a pinched boundary, such as
  two open fans touching at one vertex — appears once in our `open_verts` and
  once per edge in the C++'s, with no duplicate edge involved. The collider then
  holds a different number of leaves, so collision pairs arrive in a different
  order. The resulting partition (which vertices end up together) is the same,
  but `DisjointSets::unite` is union-by-rank, so the representative written to
  `merge_to_vert` can in principle differ. We have not constructed a case that
  shows it.

**Why it stays for now.** manifold-sharp transcribes the same two sets
(`ManifoldSharp/MeshGL.cs`, `SortedSet<(int, int)>` and `SortedSet<int>`, with a
comment citing the Rust `BTreeSet`), and verifies against this tree bit-for-bit.

**Harmonization path:** in both trees together, replace the edge set with a
multiset (a `BTreeMap<(usize, usize), usize>` count, erasing one per reverse
match), iterate it in `(start, end)` order to emit one `open_verts` entry per open
edge including repeats, and keep the stable Morton sort that follows. The probe
above is the shared regression.

**Evidence:** source reading of both implementations at the lines above; the
scratch probe for the Rust half of the first case.
