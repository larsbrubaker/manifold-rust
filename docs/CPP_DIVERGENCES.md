# Deliberate divergences from the C++ reference

Per CLAUDE.md, byte-matching `cpp-reference/manifold` is the default; the
entries below are the deliberate exceptions, each with the evidence that
justified it. Trace-diff debugging against the C++ must expect these.

Two kinds of entry live here, and they are not the same claim:

- **Justified divergence** — an accuracy fix, a real bug fix in the C++, or a
  measured improvement. CLAUDE.md's three qualifying reasons. Entries 1, 3, 4
  and 11.
- **Inherited divergence, documented and scheduled for coordinated
  harmonization** — an output shape this port already shipped, which we would
  resolve toward the C++ on the merits but cannot change unilaterally, because a
  downstream consumer verifies against this tree bit-for-bit. These are
  *disclosures*, not justifications: the entry states what differs, why it
  cannot be fixed here alone, and what a coordinated fix would take. Entries 2
  and 12 (entries 5 and 6 were of this kind and are retired).

A third, narrower kind records where the C++ output is itself not pinned:

- **Implementation-defined in the C++** — the reference's result depends on
  standard-library behavior the C++ standard leaves unspecified (for example
  `std::unordered_set` iteration order), so "the C++ output" differs between
  toolchains and no single bit pattern exists to match. The entry states which
  part of the output is affected and what is still compared exactly. Entries 7
  and 8.

A fourth kind makes no numerical claim at all:

- **API shape or extension** — a public signature that differs from the C++
  one, or a Rust-only method the C++ does not have. The entry states what the
  numbers match (or that nothing in the C++ constrains them). Entries 9, 10
  and 13.

The second kind is deliberately uncomfortable to write, which is the point — it
is a debt with a name attached, not a decision that ends the discussion. Nothing
belongs in any category for convenience.

**Remaining divergences at a glance** (retired entries 5 and 6 omitted):

| # | Kind | What |
|---|------|------|
| 1 | Justified | Robust-engine outputs skip `swap_degenerates` |
| 2 | Inherited | Centered cylinder centers in place, not via `Translate().AsOriginal()` |
| 3 | Justified | `dedupe_edges` skips duplicates an earlier repair resolved |
| 4 | Justified | Mirroring keeps each property with its corner |
| 7 | Implementation-defined | `Impl::slice` contour start vertex (and contour order) |
| 8 | Implementation-defined | Hull: order of `+0.0` / `-0.0` ties |
| 9 | API shape | `Manifold::slice` / `project` return a `CrossSection`, not `Polygons` |
| 10 | Extension | `CrossSection::minkowski_sum` (no C++ counterpart) |
| 11 | Justified | QuickHull decides above-a-face exactly |
| 12 | Inherited (latent) | Keyhole bridge searches skip a degenerate outer ring whole |
| 13 | Extension | Mirrored from manifold-sharp: `try_convex_erosion`, `try_dilate_by_convex` / `try_erode_by_convex` (union tree, convex patches), exact-boolean stage sink, `ProgressReporter::phase_total` / `report_units` |

Anything else that differs from the C++ is a bug, not an entry. The ones known
and not yet fixed are listed under *Known unresolved mismatches* at the end, so
a trace-diff session does not mistake them for intent.

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

**What differs:** when `ManifoldImpl::transform` (`src/impl_transform.rs`) flips
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

## 5. `CrossSection::decompose` groups holes by bounding box, not by a `PolyTree` (2026-09-28) — retired 2026-09-29

**Retired:** `CrossSection::decompose` (`src/cross_section_ops.rs`) now ports
C++ v3.5.2 `CrossSection::Decompose` (`src/cross_section/cross_section.cpp:475-494`,
with `decompose_outline` / `decompose_hole` at 126-151) exactly: the
`NumContour() < 2` short-circuit returning the section unchanged (so an empty
section decomposes to one empty section, and a single contour is not snapped
through Clipper2), a `FillRule::Positive` union into `clipper2_rust::PolyTreeD`,
and the outline/hole walk emitted in reversed push order. The bounding-box
heuristic this entry disclosed gave the hole of a bar `[0,10]×[0,2]` (hole
`[8,9]×[0.5,1.5]`) to a U-shaped outline (bbox `[7,11]×[-0.5,2.5]`) embracing
the bar's end; the PolyTree keeps it with the bar. The crate's tree construction
(`build_tree_d` / `recursive_check_owners`) appends children in the same order
as Clipper2 46f6391's `BuildTreeD` / `RecursiveCheckOwners`. Evidence: the C++
reference's `cross_section.cpp` compiled (MSVC, `MANIFOLD_PAR=-1`) against
Clipper2 46f6391 agrees contour-for-contour and in component order with
`cross_section::tests::test_decompose_keeps_hole_with_its_outline` (the U case),
`test_decompose_order_matches_cpp` (an island inside a hole plus a separate
square) and `test_decompose_short_circuits_below_two_contours`. manifold-sharp
must twin this change (its `CrossSection.cs` `Decompose` still carries the
heuristic) to keep bit-agreement on `Decompose`.

## 6. `MeshGL::merge` dedupes open edges and open vertices (2026-09-28) — retired 2026-09-29

**Retired:** `MeshGLP::merge` (`src/types_meshgl.rs`) now ports C++ v3.5.2
`MergeMeshGLP` (`src/sort.cpp:62-180`) exactly: open halfedges live in a counted
multiset (`BTreeMap<(start, end), count>`, erasing one copy per reverse match),
and `open_verts` gets one entry per remaining open halfedge, duplicates kept, in
`(start, end)` order. It is also generic now, so `MeshGL64::merge` exists, as the
C++ template provides. Evidence: the C++ `MergeMeshGLP` compiled standalone
against the reference headers (MSVC, `MANIFOLD_PAR=-1`) agrees with this port on
the doubled-face probe (`from = [4]`, `to = [0]`, both precisions), on a pinched
boundary where the old port chose a different representative (`from = [4, 5]`,
`to = [0, 0]`; the old port gave `[0, 5] -> [4, 4]`), and on all 3000 cases of a
randomized open-mesh sweep. Regressions: `src/types_meshgl_merge_tests.rs`.
manifold-sharp must twin this change (its `MeshGL.cs` still uses `SortedSet`s)
to keep bit-agreement on `Merge`.

## 7. `Impl::slice` contour start vertex (2026-09-29) — implementation-defined in the C++

**What differs:** the first vertex of each raw contour returned by
`ManifoldImpl::slice` (`src/face_op.rs`). Every contour is the same cyclic
vertex sequence as the C++, bit-for-bit, in the same orientation; only the
rotation of the cycle can differ.

**Why:** C++ v3.5.2 `Manifold::Impl::Slice` (`src/face_op.cpp:370-430`) collects
the straddling triangles in a `std::unordered_set<int>` and starts each contour
at `*tris.begin()`. That iteration order is unspecified by the standard (MSVC,
libstdc++ and libc++ each give their own), so the C++ start vertex is a property
of the toolchain, not of the algorithm. When several contours exist, the order
of the contours follows the same set iteration. This port runs the same loop
over a `std::collections::BTreeSet<usize>`, so `begin()` is the lowest-indexed
untraced triangle: the start vertex and the contour order are deterministic
(ascending triangle index) and reproducible across processes and platforms,
but may still differ from any particular C++ toolchain's order. (Until
2026-09-29 the port used a `HashSet` with the default `RandomState` hasher, so
both varied from process to process.) Everything downstream of the public API is unaffected in shape:
`Manifold::slice` wraps the contours in a `FillRule::Positive` union
(`CrossSection::new`), as every C++ caller does.

**Evidence:** `cross_section::ctor_tests::test_raw_slice_matches_cpp_lerp_bits`
compares `Manifold::sphere(1.0, 8).as_impl().slice(0.3)` against the C++ (MSVC,
`MANIFOLD_PAR=-1`, Clipper2 46f6391) as a cyclic sequence, after the crossing
interpolation was corrected to `la::lerp(below, above, a)` =
`below * (1 - a) + above * a` (it was `below + a * (above - below)`, which
differed by one ULP in 10 of 24 coordinates). The MSVC build starts that contour
at a different triangle than this port does; the cycles agree bit-for-bit.
`cross_section::ctor_tests::test_raw_slice_contour_order_is_deterministic`
pins the contour order and start vertices of a three-sphere raw slice; it
failed on every run against the `HashSet` version. manifold-sharp must twin
both the `lerp` form and the ascending-triangle-index start (a sorted set, or
the minimum remaining index) to keep its raw slice bit-equal, cycle rotation
and contour order included.

## 8. `CrossSection` hull: order of `+0.0` / `-0.0` ties (2026-09-29) — implementation-defined in the C++

**What differs:** potentially, the sign of a zero coordinate in a hull vertex
when the input holds two points that differ only in the sign of a zero
coordinate. Nothing else: `hull_points` / `hull_cross_sections`
(`src/cross_section_ops.rs`) port C++ v3.5.2 `HullImpl`
(`src/cross_section/cross_section.cpp:183-206`) step for step.

**Why:** C++ sorts the points with `std::sort` and `V2Lesser`, under which
`(0.0, y)` and `(-0.0, y)` are equivalent; `std::sort` is not stable and its
permutation of equivalent elements is left to the library. MSVC's `std::sort`
uses insertion sort below 32 elements, which is stable, and introsort above, so
against MSVC the tie order (and therefore this divergence) can only show for
inputs of 32 or more points; other libraries' small-input paths differ. This
port sorts with the stable `slice::sort_by` on the same comparator. The `CCW(.., 0.0)` backtrack then
keeps one of the tied points per chain (the later one in the lower chain, the
earlier one in the upper), so the kept zero's sign follows the library's
permutation. No other input is affected, because equivalent
points under `V2Lesser` are otherwise bit-identical.

**Evidence:** `cross_section::tests::test_hull_matches_cpp_hull_impl` pins the
C++ (MSVC) result bit-for-bit on degenerate (fewer than three, collinear,
coincident), near-duplicate, underflowing-`CCW` and multi-section inputs.

## 9. `Manifold::slice` / `project` return a `CrossSection`, not `Polygons` (API shape)

**What differs:** C++ v3.5.2 `Manifold::Slice(double)` and `Manifold::Project()`
return raw `Polygons` (`src/manifold.cpp`, delegating to `Impl::Slice` /
`Impl::Project` in `src/face_op.cpp`). This port's `Manifold::slice` /
`Manifold::project` (`src/manifold.rs`) return a `CrossSection` built with
`CrossSection::new`, i.e. C++ `CrossSection(m.Slice(height))` /
`CrossSection(m.Project())` with the default `FillRule::Positive` union.

**What still matches:** the numbers. The returned section's contours are
bit-for-bit those of the C++ `CrossSection(Polygons)` constructor applied to the
C++ raw polygons (`cross_section::ctor_tests::test_slice_and_project_wrap_like_cpp`).
The raw polygons are reachable in-crate as `as_impl().slice(h)` /
`as_impl().project()`, which the ported `TEST(Smooth, Fillet)` uses to extrude
the un-unioned slice exactly as the C++ test does; the raw slice is pinned by
`cross_section::ctor_tests::test_raw_slice_matches_cpp_lerp_bits` (modulo
entry 7) and the raw projection by `test_slice_and_project_wrap_like_cpp`.

**Why kept:** the signature shipped long before this audit, and manifold-sharp
mirrors it. Changing it is a public-API decision, not a numerical fix.

## 10. `CrossSection::minkowski_sum` has no C++ counterpart (extension)

**What differs:** C++ v3.5.2 `CrossSection` has no Minkowski operation. This
port's `CrossSection::minkowski_sum` (`src/cross_section_ops.rs`) runs Clipper2's
`MinkowskiSum` (closed, at `precision_`) for every pair of contours and
concatenates the per-pair results. Clipper2 unions each pair's output with
`FillRule::NonZero`, but the concatenation is not unioned again, so the result
can hold overlapping contours when either operand has more than one contour.

**Why kept:** it is an extension, so nothing in the reference constrains it;
callers that need a clean section can pass the result through
`CrossSection::new`. Should the C++ ever gain a Minkowski operation, this entry
becomes a porting task.

## 11. QuickHull decides above-a-face exactly (2026-09-30)

**What differs:** the QuickHull flood fill (`src/quickhull_algo.rs`) calls a face
visible from the apex only when the exact `orient3d` of the face's three corners
and the apex is strictly positive (`quickhull::is_above`); exactly coplanar is
hidden. `add_point_to_face` queues a point on a face only when it is also exactly
above it; the float epsilon test that decides whether a point is outside at all
is unchanged. The C++ (`src/quickhull.cpp:391`, and its `AddPointToFace`) reads
both off the float distance to the face's stored plane (`d > 0.0`).

**Why:** a real bug in the C++. When the apex is collinear with a hull edge it
lies in the plane of both faces on that edge, and rounding can put one face at
+5.6e-17 (visible) and the other at 0 (hidden). The edge becomes a horizon edge,
the new face coned to it has zero area and a noise normal, and every later
visibility test against it is noise: the hull ends non-convex. Minkowski sums
hull each triangle swept by the tool, so a folded hull loses solid from the sum
(Thingi10K 63451, 641145 and 287448 in manifold-sharp's dilation). The exact
orientation gives both faces the same answer, so no zero-area face is built.

**What is observably different:** hulls whose points hold exactly collinear or
coplanar sets can come out with different triangles. No existing expected value
in this tree moved. manifold-sharp made the same fix in its commit `1931e87` (its
`RUST_DIVERGENCES.md` entry 7), so the two ports agree. (This tree's commit
`a69e579` cites it as `b178563`, its hash before manifold-sharp's history was
rewritten; `b178563` is on no branch there.)

**Evidence:** `quickhull::tests::test_hull_of_a_flat_triangle_swept_by_sphere_is_convex`
(Thingi10K 63451 triangle 163 swept by `sphere(0.3, 8)`) and
`test_thingi641145_triangle109_swept_hull_is_convex` (641145 triangle 109 swept by
`sphere(0.05781898171099809, 12)`): an input point lay 0.581 and 0.2336 outside a
hull face before the fix, under 1e-12 after. Both tests are shared with
manifold-sharp's `QuickHullContainmentTests`.

## 12. The keyhole bridge searches skip a degenerate outer ring whole (2026-10-03)

**What differs:** `EarClip::cut_keyhole` and `find_closer_bridge`
(`src/polygon_earclip_keyhole.rs`) walk every outer ring with
`for_each_loop_vert`. When the walk reports the ring degenerate (a vert with
`right == left`), the search restores the connector it had before that ring, so
a degenerate ring contributes nothing. The C++ (`src/polygon.cpp:544`, `Loop`)
applies the function to each vert up to the degenerate one and returns
`polygon_.end()`; its two callers here (`CutKeyhole`, `:724`, and
`FindCloserBridge`, `:773`) ignore that return value, so whatever the lambda
did before the degenerate vert stands.

**Why we keep it.** It is inherited. The first Rust port collected each ring
with `loop_verts`, which returns `None` for a degenerate ring, and the searches
skipped it with `continue`; that shipped in every release so far, and
manifold-sharp transcribes it (`ManifoldSharp/PolygonEarclip.Algorithm.cs`,
which verifies against this tree bit-for-bit). External PR #6 replaced the
collection with a visitor for speed and kept the skip-whole behaviour exactly
by saving and restoring the connector around each ring, rather than moving
toward the C++ in a speed change.

**What is observably different: nothing, and here is why.** The two behaviours
can differ only if a walk applies the function to at least one vert and then
meets a degenerate one. That cannot happen. `Link` (`polygon.cpp:531`) and
`JoinPolygons` (`:789`-`:792`) are the only pointer writes, and both keep every
unclipped vert's `left->right` and `right->left` pointing back at it. So the
unclipped verts form closed rings, and a vert with `right == left` belongs to a
ring of one or two verts in which every vert has `right == left`. `Loop` checks
`right == left` on each unclipped vert *before* calling the function on it, and
once it reaches an unclipped vert it only follows `right` within that vert's
ring. A degenerate ring therefore returns at the first vert the function would
have seen, in both the C++ and the port, and no vert is ever visited. The
save-and-restore never changes the connector. It stays as a guard in case the
walk or the ring invariant changes.

The entry is filed as inherited, not as "no difference", because the code shape
does differ from the C++ and a trace-diff session reading the two side by side
should not have to re-derive the argument above. Harmonizing is free: dropping
the restore (taking the C++ shape) cannot change any output, so manifold-sharp
can follow at its own pace and neither tree's pins move.

**Evidence:** the invariant argument above, from `polygon.cpp:531-566` and
`:785-797`. A throwaway fuzz (not committed) instrumented `for_each_loop_vert`
to count degenerate walks and walks that had called the function before
reporting degenerate. It triangulated 2,000,000 random sets of 1-5 rings with
1-6 verts each on a 5x5 integer grid, a third of the rings being reversed or
rotated copies of earlier ones (coincident outer and hole rings collapse to
two-vert rings when joined), at epsilon 0, 1e-9, 0.3 and automatic. The bridge
searches met a degenerate outer ring 181,312 times; every walk anywhere in the
triangulator that reported degenerate (2,921,184) had visited zero verts first.
`polygon_earclip::tests::keyholing_many_holes_keeps_its_triangles` pins the
triangles that the visitor change left unchanged.

**Related, and not a divergence:** both searches also skip any outer ring whose
bounding box (`EarClip::outer_bbox`, not in C++) shows that no vert in it could
pass the search's tests. That is a speed difference, not an output difference:
a skipped ring is one the C++ walk could not have taken a connector from, so
the bridges and triangles are the same. A trace that instruments the walks will
see fewer rings visited. The determinant bound in `find_closer_bridge` applies
only inside a magnitude window (|connector - start| >= 1e-60, box distance and
epsilon <= 1e60), because outside it `ccw`'s squares can underflow or overflow
and call an outside vert collinear; the comment there gives the derivation, and
`keyhole_cull_keeps_a_bridge_whose_ccw_underflows` and `..._overflows` pin it.

## 13. Convex erosion, the convex dilation tree and their progress hooks, mirrored from manifold-sharp (extension, 2026-10-05)

**What differs:** three Rust-only APIs with no C++ counterpart, each a 1:1
mirror of a manifold-sharp addition that was until now a sharp-only entry in its
`docs/RUST_DIVERGENCES.md`:

- `Manifold::try_convex_erosion(&self, other, token, progress) -> Option<Manifold>`
  (`src/convex_erosion.rs`, sharp `ConvexErosion.cs` / `Manifold.TryConvexErosion`,
  sharp entry 5). The Minkowski difference of a *convex* solid in closed form:
  the intersection of the solid's face halfspaces, each pushed inward by the
  tool's support, enumerated through the polar dual hull and solved per vertex
  by Cramer's rule. `None` (decline) for a non-convex solid or tool, a tool not
  containing the origin, a centroid not strictly inside the eroded body, a
  degenerate triple, or a vertex failing the closing feasibility check; a
  cancelled run is `Some` of an empty `Error::Cancelled` impl.
- `ProgressReporter::phase_total()` and `ProgressReporter::report_units(f64)`
  (`src/progress.rs`, sharp `Progress.cs` `PhaseTotal` / `ReportUnits`, sharp
  entry 6's reporter hooks). Read-only / write-only from the kernel's side.
- The exact boolean's stage sink (`src/boolean_stage_progress.rs`, sharp
  `BooleanStageProgress.cs`): `boolean3::boolean_with_token_and_stage`,
  `Boolean3::new_with_token_and_stage` and
  `boolean_result::boolean_result_with_token_and_stage` take an optional
  `&dyn Fn(f64)` and invoke it with eight constant cumulative marks at the
  cancel gates closing each heavy stage. It reads nothing and is handed only
  constants; the variants without it pass `None`, which is the prior code.
  `stage_sink_hears_every_mark_and_changes_no_bit` pins both halves.

**What does not differ:** nothing is rerouted. `minkowski::minkowski`,
`Manifold::minkowski_difference` and every other ported path still run the
ported sweep and produce the C++'s bits, including for a convex solid where the
closed form would apply. The callback-side additions compute nothing.

**Numbers:** nothing in the C++ constrains them; the constraint is the twin
port. The Rust performs the same floating-point operations in the same order as
`ConvexErosion.cs`, and maps dual-hull vertices back to planes by exact bits, the
equality sharp's `Vec3.Equals` uses. `src/convex_erosion_tests.rs` ports
`ConvexErosionTests.cs` and `ConvexErosionTests.Contract.cs` 1:1 (a 20-cube by a
unit ball is exactly 5832.0; a 2048-triangle sphere agrees with the sweep to
< 1e-12 relative).

- `Manifold::try_dilate_by_convex` / `Manifold::try_erode_by_convex`
  (`src/convex_dilation.rs`, `src/convex_patches.rs`; sharp `ConvexDilation.cs`,
  `ConvexPatches.cs`, the rest of sharp entry 6). The same per-triangle (or,
  when dilating, exactly-proven convex-patch) hulls `minkowski.rs` builds,
  reduced through a balanced union tree (leaves of 16 units, pairwise levels)
  whose leaf and level maps run through `maybe_par_map_ct(_progress)`, so the
  `parallel` feature is the Rust's `ManifoldParallel.Enabled` (single-threaded
  on wasm32 and without the feature). Solids whose component boxes overlap are
  first rebuilt by `robust::rebuild_with_rule(Positive)`. Erosion subtracts the
  tree's union from the solid once. To let every node run on one engine,
  `csg_tree`'s private `batch_union` / `batch_boolean` / `simple_boolean` take
  an optional engine (`None` = the process default, the prior behavior).
  Not bit-identical to `minkowski_sum` / `minkowski_difference` (a different
  union order rounds differently): volume and genus agree at 1e-9. Sequential
  and parallel builds are bit-identical to each other; `par_tests.rs` pins the
  drilled-part dilation and erosion fingerprints (0xe4a5747b38ef2b89,
  0x565d3e8a3c515759) in both builds.

Tests ported 1:1 from sharp: `ConvexDilationTests.cs`, `.Erosion.cs`,
`.Patches.cs`, `.Nested.cs` (`src/convex_dilation_*tests.rs`) and
`ParallelismTests.ConvexDilation.cs` (`src/par_tests.rs`; the C# runtime
toggle becomes 1-vs-8 rayon threads plus a hash pinned across both builds).

Sharp entry 10 (CSG tree evaluation with a progress reporter) was taken by
e444fb9 (`CsgNode::evaluate_with_token_and_progress`), and `batch_boolean`
now also serializes a round's pairs when a reporter is attached, as sharp
does. With this entry nothing in sharp's entries 5, 6 or 10 remains
sharp-only.

## Known unresolved mismatches (bugs, not entries)

Found while porting `TEST(Smooth, Fillet)` exactly (2026-09-29), against the C++
compiled with MSVC (`MANIFOLD_PAR=-1`):

- **Chained booleans evaluate eagerly.** `Manifold` holds its `ManifoldImpl`
  directly, so `&(&a + &b) + &c` runs two pairwise booleans, and
  `Manifold::batch_boolean` is a left fold. C++ builds a lazy `CsgOpNode` and
  flattens `a + b + c` into one `BatchBoolean` / `BatchUnion`. On the Fillet
  inputs the C++ inline `cylinder + chamfer + base` refines to 3770 triangles;
  with the first union forced (as this port always does) the C++ gives 3822, as
  does this port. `src/csg_tree.rs` ports the C++ tree but only
  `minkowski.rs` uses it.
- **`smooth_by_normals` + `refine_to_tolerance` drift by a few ULPs.** With the
  union order matched (3822 triangles, 1913 vertices), every stage through
  `SmoothByNormals` agrees bit-for-bit in volume and surface area, and the
  refined surface area agrees too, but most refined vertex positions differ in
  the last few bits and the volume by about 1e-9 relative
  (`0x40be40dfe1b35e91` here vs `0x40be40dfe1b7737f` in C++). The test asserts
  the C++ test's own `EXPECT_NEAR` bounds, which both meet.

Found by manifold-sharp while fixing QuickHull (entry 11), 2026-09-30:

- **`swap_degenerates` cuts solid beside a zero-thickness fin.** On the
  operand pair in `src/testdata/minkowski-641145-union-{a,b}.txt` (the 19th
  pairwise union of Thingi10K 641145's Minkowski sum before entry 11, both
  `NoError`) the exact union comes out 3.9e-3 relative short of
  inclusion-exclusion, with 4.2e-4 of B outside it (manifold-sharp, which
  traced it, reports the same in C++ v3.5.2). The boolean itself is right; `simplify_topology`'s
  `swap_degenerates` loses the solid. `set_normals_and_coplanar`'s flood fill
  (C++ `SetNormalsAndCoplanar`, the same orientation-blind fill entry 1
  describes) hands a sound triangle 0.041 tall the reversed normal of an
  opposite-facing coplanar seed; projected through it the triangle reads as
  inverted, so `recursive_edge_swap` swaps its long edge into a neighbour in
  another plane and cuts out a wedge. Two narrow fixes were rejected (they
  break `test_cpp_simplify` and `test_cpp_nonconvex_convex_minkowski_sum`).
  Nothing in this crate produces these operands since entry 11, but the
  defect stands. Pinned, ignored, by
  `minkowski::union_regression_tests::exact_union_of_thingi641145_partial_unions_contains_both_operands`,
  shared with manifold-sharp's `MinkowskiUnionRegressionTests` (its commit
  `1931e87`).

Found while landing external PR #7 (parallel batch rounds), 2026-10-03:

- **`boolean_with_token` composes disjoint operands of an Add.** When the two
  operands' boxes do not overlap, `boolean3::boolean_with_token` returns
  `compose_meshes(&[a, b])` (`src/boolean3.rs`, the non-overlapping fast
  path). C++ `Boolean3` has no such shortcut: it early-outs only the
  intersection phase (`boolean3.cpp:509`) and still assembles the result
  through `Boolean3::Result`, which orders verts, triangles and face IDs
  differently from `Compose`. The run table is the same either way (since the
  compose fix below). A plain `a + b` of disjoint meshes composes in both,
  because the CSG tree's `BatchUnion` composes leaves whose boxes do not
  overlap before any boolean runs; the difference shows only when leaf boxes
  overlap (as transformed boxes of rotated leaves can) while the meshes
  themselves do not. Example: the four rotated spheres of
  `manifold::tests::compose::test_batch_union_of_rotated_instanced_leaves_is_independent_of_scheduling`
  agree with C++ v3.5.2 run for run, but `vert_properties`, `tri_verts` and
  `face_id` come out in a different order (C++ whole-mesh hash
  `0x5a5fb492e1484228`). Removing the fast path would move shipped output
  that manifold-sharp verifies bit-for-bit, so it stays until both ports
  change together.

## Fixed mismatches that change shipped output (bugs, not entries)

Bugs that were never ledger entries but whose fix changes output manifold-sharp
verifies bit-for-bit, so its twin knows what to port.

- **Composing instanced copies merged their runs** (fixed 2026-10-03).
  `boolean3::compose_meshes` (behind `Manifold::compose`, the disjoint-union
  fast path in `boolean3::boolean_with_token` and `csg_tree::batch_union`)
  merged every input's raw mesh IDs into one map. Two copies of one mesh share
  their mesh IDs, so they became a single MeshGL run carrying only the last
  copy's `run_transform`. C++ `CsgLeafNode::Compose` shifts node `i`'s meshIDs
  by `i * meshIDCounter` (`src/csg_tree.cpp:289, 388, 400`) before
  `IncrementMeshIDs`, so each copy keeps its own run and transform, node by
  node. The port now renumbers each input's mesh IDs to local IDs `1, 2, 3, ...`
  in that order (ascending within an input) before `increment_mesh_ids`, which
  gives the same ranks. Evidence: the v3.5.2 reference built as a static
  library gives `Compose({cube, cube.Translate({3, 0, 0})})`, `cube +
  cube.Translate({3, 0, 0})`, a three-copy `Compose` and a three-copy
  `BatchBoolean(Add)` one run per copy in input order;
  `manifold::tests::compose` asserts those run tables and failed on the old
  code (one run). No other expected value in this tree moved. manifold-sharp
  must twin this (`Boolean3.Functions.cs` `ComposeMeshes`).
