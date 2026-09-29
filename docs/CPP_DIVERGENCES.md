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
  cannot be fixed here alone, and what a coordinated fix would take. Entry 2
  (entries 5 and 6 were of this kind and are retired).

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
  numbers match (or that nothing in the C++ constrains them). Entries 9 and 10.

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
