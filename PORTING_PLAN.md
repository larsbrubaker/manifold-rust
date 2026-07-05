# Manifold C++ → Rust Porting Plan

This is a **roadmap of remaining work** to finish porting
[Manifold](https://github.com/elalish/manifold) to Rust. It tracks what is left to do,
not what has already been done (use `git log` for history). Every change must reproduce the
C++ reference with **exact numerical match** — identical results on identical inputs.

**Status:** 517 passing, 0 failing, 8 ignored.
**C++ reference target:** v3.5.0 (submodule at tag `v3.5.0`, commit `541c33bd`).
**Core engine:** all 18 phases (linalg → boolean → CSG → cross-section → SDF → minkowski →
WASM) are implemented. Remaining work is the v3.5.0 deltas below plus the ignored-test
backlog.

---

## Remaining v3.5.0 work

All substantive v3.5.0 deltas are now ported or confirmed N/A. What's left is
optional/peripheral:

### Likely N/A — confirm before skipping
- Import/Export fixes + topological 3MF sort (#1705, #1685, #1692) — no 3MF in Rust. But
  #1685 "improve import transformation" may fix the non-manifold OBJ load behind the
  `craycloud` tests; worth checking when tackling import.

### Out of scope (do not re-litigate)
- ExecutionContext (progress + mid-boolean cancellation, #1663/#1669/#1699/#1704) —
  non-numerical; only needed if cancellation is wanted.
- Halfedge Refactoring (#1709) — structural cleanup, no behavior change.
- CrossSection backend selector (#1710), concurrent-const-access safety (#1636),
  affinity_partitioner race (#1664) — threading concerns moot in the sequential port.

---

## Ignored tests (8) — grouped by the work needed to clear them

> **Audit note:** C++ has zero `DISABLED_` tests, so the fair bar is Rust-release ≈
> C++-release. The old "9 slow in debug" bucket was a mislabeled mix; broken out below.

### Slow in debug only — fast & passing in release (4)
`sdf_blobs` (release 5.6s after the HashMap→Vec voxel-grid fix, was 62.6s), `hull_sphere`
(~7s), `generic_twin_7081` + `complex_generic_twin_7081` (~18s each in release; pass after
the RecursiveEdgeSwap fix). Legitimately debug-slow; no further work unless we run a release
test lane.

### Heavy CSG/boolean — gap is deferred parallelism (2)
`hull_menger_sponge` (depth-4 CSG, >60s in release), `properties_mingap_stretchy_bracelet`.
C++ runs Boolean/CSG via TBB; Rust is sequential-first (rayon behind the `parallel` feature).
Closing the gap = enabling/validating rayon, not an identified algorithmic bug.

### Minkowski exactness — 2 of 3 FIXED (2026-07), 1 near-match left
The instrumented-trace session (Rust vs local C++ MSVC build, per-boolean diff of
BatchBoolean order → Boolean3 kernels → vert normals → PairUp inputs) found and fixed four
real divergences; `nonconvex_convex_minkowski_sum/difference` now pass and are un-ignored:
1. **BatchBoolean reduction order** (`csg_tree.rs`): C++ v3.5.0 heap entries are
   `(leaf, serial)` pairs — pops the MOST-vert node (ties → largest serial), processes up to
   4 pairs per round before re-pushing (even sequentially). Rust popped smallest with no
   serials, one pair at a time. Also `batch_union` now push-back + front/back-swaps like C++
   (matters when chunking >1000 children).
2. **Vertex-normal accumulation order** (`face_op.rs calculate_vert_normals`): C++
   `ForVert` STEPS FIRST then calls func, so `firstEdge` is processed LAST. The Rust port
   processed it first — same triangles, shifted float-sum order, ~1-ULP different vertex
   normals, which flip `Shadows` SOS ties in Boolean3's kernels on exactly-coincident
   geometry (verified: `kernel02` attributed a vert-on-edge crossing to the wrong coplanar
   face).
3. **`append_partial_edges` edgeVec** (`boolean_result.rs`): C++ projects along
   `inP.vertPos[vEnd] - inP.vertPos[vStart]` (INPUT positions); Rust used output positions
   via `vp2r`, which are garbage for non-retained (i03==0) endpoints → wildly wrong edgePos
   → wrong PairUp pairing on degenerate edges (verified with PairUp input dumps).
4. **`pair_up` partition** (`boolean_result.rs`): C++ `std::partition` is the unstable
   Hoare two-pointer form; its swap pattern decides the order of fully-tied entries (same
   edgePos AND collisionId, e.g. duplicated retained verts), which the stable sorts keep.
   Ported as `cpp_partition`.

**FIXED (2026-07, follow-up session): `nonconvex_nonconvex_minkowski_sum` un-ignored.**
The trace-diff loop (BIN/ROUT input-output mesh fingerprints per boolean → B3 kernel dumps →
per-phase simplify hashes → per-face triangulation dumps) found and fixed, in order:
1. **Ear-queue FIFO tie-break** (`polygon_earclip.rs`): C++'s `std::multiset` pops
   equal-cost ears in insertion order; the Rust `BinaryHeap` popped ties in arbitrary heap
   order. Costs tie constantly (`kBest = -inf` for short ears). Added a monotonic `seq` to
   `EarEntry`, tie-breaking older-first.
2. **`cut_keyhole` connector selection** (`polygon_earclip.rs`): the C++ CheckEdge tests
   `start->InsideEdge(edge)` per candidate; the Rust port tested `start` against its own
   left neighbor (loop-invariant!) and only applied the vertical-ordering tie-break when
   the CCW test returned 0 (C++: whenever it isn't 1). This mispicked keyhole bridges on
   degenerate faces with holes — the last divergence blocking the test.
3. **`swap_degenerates` neighbor check** (`edge_op.rs`): C++ `swappableEdge` projects the
   PAIR triangle's verts (`pairTriEdge`) under the neighbor's projection; Rust projected
   tri0's. (Note: inside `RecursiveEdgeSwap` C++ genuinely uses tri0edge — the two sites
   differ on purpose.)
4. **`recursive_edge_swap` tag by-reference** (`edge_op.rs`): C++ takes `int& tag` and
   increments it in the facing-degenerates collapse branch so previously-visited edges
   become re-processable; the Rust port took `tag` by value.
5. **`collapse_edge` prop update** is an else-if chain in C++ (tri0's prop wins when a
   triangle matches both face groups); `SwapEdge`'s prop writes reordered to match C++.

Fixing #3 exposed two **construction-exactness** bugs previously masked by compensating
behavior (`cylinder_zero_radius_low`, `refine_quads` regressed, then both root-caused):
- **Global `sind` now uses remquo reduction** (`types.rs`), matching C++ `common.h` — the
  floor-based reduction differed by ~1 ULP for reduced args in (45°, 90°), putting cylinder
  circle verts off by 1 ULP. The old "global remquo-sind regresses RefineQuads" note was a
  compensating-bug artifact; with the fixes below it holds exactly. The local
  `sind_remquo` in `manifold.rs` was folded back into `types::sind`.
- **`extrude` cap triangulation** (`constructors.rs`): C++ `TriangulateIdx` defaults to
  `allowConvex=true` (alternating-fan for convex caps); Rust passed `false` (ear-clip) →
  different cap triangle sets on every cylinder/extrude. Also the side-wall scale·rotation
  is now built as matrix entries first (C++ `mat2` multiply association), and the
  zero-radius cone branch now applies C++'s composed Mirror∘Translate via
  `Impl::Transform` + `AsOriginal` instead of negating z in place (which skipped the
  winding flip). Cylinder(5,0,3,256) and its boolean results are now **bit-identical** to
  the C++ reference.

### SDF thin-shell marching topology — FIXED (2026-07)
`sdf_sphere_shell` now passes with genus **13396, exactly matching the local C++ build**
(still `#[ignore]`d only for debug-suite speed, ~1min debug / ~5s release). Four findings:
1. **Test-porting error (the "genus bug"):** the C++ test passes `tolerance=0.0001` to
   `LevelSet` (enabling FindSurface vertex refinement) and uses `r - 0.995f` (a FLOAT
   literal). The Rust test omitted the tolerance → pure interpolation → genus ~9576. The
   marching implementation itself was never topologically wrong.
2. **Nondeterminism:** `level_set` iterated a std `HashMap` for ComputeVerts/BuildTris —
   random order per process; three runs gave three different genus values. Replaced with a
   faithful port of C++ `hashtable.h` (`GridHashTable` in sdf.rs: splitmix64 `hash64bit`,
   linear probing, power-of-2 sizing, used*2>size Full check, open-slot default-value
   lookups, and the LevelSet resize-and-rerun protocol) iterated in SLOT order.
3. **`la::lerp` association:** C++ lerp is `a*(1-t) + b*t`, not `a + (b-a)*t` — fixed in
   `find_surface`'s scalar and vector lerps.
4. **`ComputeGridPow`:** C++ is `CeilLog2(n+3)`; Rust computed `ceil_log2(n+2)+1` (one bit
   wider) — encoded indices are the hash keys, so this changed slot order. Fixing it made
   the genus match C++ exactly.

### Coplanar/coincident-face geometry — FIXED (2026-07)
The whole cluster landed in one session; root causes and fixes:
1. **`Manifold::rotate` ported to the C++ matrix construction** (`csg_tree.cpp`: degree-based
   `sind`/`cosd` axis matrices composed `rZ*rY*rX`) with a **local remquo-`sind`**
   (`sind_remquo` in `manifold.rs`). The old quaternion path differed by ~1 ULP, flipping an
   exact-coplanar SOS tie in Boolean3's `kernel02` (the spurious 21st vert in
   `almost_coplanar`). The remquo reduction is local to `rotate` because switching the global
   floor-`sind` regresses `RefineQuads`'s cylinder count.
2. **`face2tri` rewritten to the C++ v3.5.0 scheme** (`face_op.cpp`): pairing is face-local
   during triangulation (ported `HalfedgeTriangulation` + `triangulate_idx_halfedges` in
   `polygon.rs`, LIFO multimap robust to duplicate vertex pairs), boundary edges pair
   across faces via the *original* assembly `paired_halfedge` (`contour2tri`). The old port
   re-derived all pairs from a global `(start_vert, end_vert)` HashMap, which mispaired
   duplicate edges on exactly-coplanar faces (menger_sponge's exact 90° rotations) → unpaired
   halfedge → non-manifold at simplify entry. Also: boolean now passes `allowConvex=false`
   (matching `boolean_result.cpp:932`) and `prop_vert` propagates through triangulation.
3. **`update_vert` assert corrected** to C++ semantics (start==end is a graceful no-op; the
   infinite-loop check is `current != start_edge` *after* stepping).

Cleared and un-ignored: `almost_coplanar`, `convex_convex_minkowski_difference`,
`nonconvex_nonconvex_minkowski_difference`, `openscad_crash` (all pass in debug + release).

### Boolean produces a non-manifold intermediate — deep robustness — FIXED (3 of 4)
**Root cause (found):** `RecursiveEdgeSwap` in `edge_op.rs` was missing two pieces of the
C++ `edge_op.cpp`:
1. The `SwapEdge` lambda's tail block that, when a swap creates a *duplicate* edge, calls
   `FormLoop`/`RemoveIfFolded` to split the mesh and keep it manifold. Omitting it left a
   non-manifold edge that later sent `collapse_edge`/`update_vert` walking an unpaired
   (`-1`) halfedge → `usize::MAX` OOB (or the sort.rs:298 odd-halfedge assertion).
2. The C++ `if (edge < 0) return;` guard. C++'s `edgeSwapStack` holds `int` and the
   normal-swap path legitimately pushes pair values that can be `-1`; the Rust stack used
   `usize`, wrapping `-1` to `usize::MAX`. Stack + param are now `i32` with the guard.

This cleared `complex_sweep` (vol≈3757, exact C++ match), `complex_craycloud`,
`craycloud_bool` — all now un-ignored and passing in debug + release.

Still ignored: `openscad_crash` — genuinely needs **processOverlaps** (the C++
`TEST(Manifold, OpenscadCrash)` is itself gated behind `MANIFOLD_DEBUG` and sets
`ManifoldParams().processOverlaps = true`). Same blocker as the minkowski/genus group.

### Boolean hangs — FIXED (kept ignored for debug speed)
`complex_generic_twin_7081`, `generic_twin_7081`. These were the same non-manifold
infinite loop, fixed by the `RecursiveEdgeSwap` change above. They now pass in release
(~18s each; C++ ≈50s) and are kept `#[ignore]`d only for debug-suite speed, like
`sdf_blobs`/`hull_sphere`.

---

## Unported C++ tests (2)
| File | Test | Blocker |
|------|------|---------|
| boolean_complex_test.cpp | InterpolatedNormals | Large inline mesh property data |
| boolean_complex_test.cpp | Ring | Huge inline `mgl_0`/`mgl_1` mesh data (~600 lines) |

---

## After all tests pass: performance parity
- Benchmark identical operations (sphere booleans, complex OBJ booleans).
- Profile hotspots (collider queries, boolean result assembly).
- Add `rayon` parallelism behind the `parallel` feature.
- Target: within 2× of C++ for all operations.

---

## Guiding principles
1. **Exact match** — same floating-point results as C++; instrument both to compare.
2. **No stubs** — every function implemented; no `todo!()` / `unimplemented!()`.
3. **Tests first** — port the relevant C++ test before/with the implementation; when C++
   changes an expected value, update the Rust test to the new C++ value (it's not a
   regression).
4. **Dependency ordered** — never implement a function before its dependencies exist.

## Key architectural decisions
| Decision | Choice | Rationale |
|----------|--------|-----------|
| Linear algebra | Roll our own `linalg.rs` | Control exact f64 semantics; no `glam` |
| Trig | Own musl-derived `math.rs` | C++ uses deterministic trig (#1606); platform libm diverges |
| Float precision | f64 throughout; f32 only at MeshGL boundary | Matches C++ double/float usage |
| Parallelism | Sequential first; `rayon` behind `parallel` feature | Verify correctness before parallelizing |
| Cross-section | `clipper2-rust` crate | Our own Rust port; same API surface |
| BVH/Collider | Radix tree + Morton codes | Matches C++ algorithm exactly |
| CSG tree | `Arc<ManifoldImpl>` + automatic `Drop` | No manual refcount teardown (so #1673 is N/A) |
| Sorting | `sort_by_key` (stable) where C++ uses `stable_sort` | Tied keys must break on original index order |
