# Manifold C++ → Rust Porting Plan

This is a **roadmap of remaining work** to finish porting
[Manifold](https://github.com/elalish/manifold) to Rust. It tracks what is left to do,
not what has already been done (use `git log` for history). Every change must reproduce the
C++ reference with **exact numerical match** — identical results on identical inputs.

**Status:** 514 passing, 0 failing, 11 ignored.
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

## Ignored tests (11) — grouped by the work needed to clear them

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

### Minkowski area slightly high — near-match (3)
`nonconvex_convex_minkowski_sum/difference`, `nonconvex_nonconvex_minkowski_sum`.
**Fast in release — not slow.** After the 2026-07 matrix-`rotate` + v3.5.0-`face2tri` fixes
(below), these are no longer non-manifold/genus failures. **Volume matches C++** (sum even to
±1e-5) but **surface area is slightly high**: 34.073 vs 34.06, 16.742 vs 16.70, 31.3996 vs
31.17691 — extra (likely coplanar sliver) surface survives the pipeline. The remaining gap
is a small geometric divergence in the minkowski decompose→hull→batch-boolean chain; needs
instrumented comparison against a local C++ build (hull-by-hull, then per-boolean).

### SDF thin-shell marching topology (1)
`sdf_sphere_shell` — genus 9560 vs expected ~14235 (perf now fine, 12s in release after the
Vec fix). The half-voxel-thick shell exposes a marching-tetrahedra crossing/snap topology
difference vs C++. Real correctness bug; needs root-cause in the SDF crossing logic.

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
