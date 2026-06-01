# Manifold C++ → Rust Porting Plan

This is a **roadmap of remaining work** to finish porting
[Manifold](https://github.com/elalish/manifold) to Rust. It tracks what is left to do,
not what has already been done (use `git log` for history). Every change must reproduce the
C++ reference with **exact numerical match** — identical results on identical inputs.

**Status:** 510 passing, 0 failing, 15 ignored.
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

## Ignored tests (15) — grouped by the work needed to clear them

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

### Overlapping-polygon triangulation match (4)
`nonconvex_convex_minkowski_sum/difference`, `nonconvex_nonconvex_minkowski_sum/difference`.
**Fast in release (0.2–2.9s) — not slow.** Volume + area now match C++ exactly; only **genus**
differs.
> **Correction (2026-05):** this is *not* a missing "processOverlaps feature."
> `processOverlaps` (default **true** in C++ `common.h`) only gates a **debug-only assertion**
> — `polygon.cpp:950` runs `CheckGeometry` (which throws `geometryErr` on an overlapping
> triangulation) *only* when `processOverlaps == false` AND `MANIFOLD_DEBUG` AND
> `intermediateChecks`. It changes **no** triangulation output. The real gap: Rust's `EarClip`
> triangulator must, for non-ε-valid (self-overlapping) polygons, "always return a manifold
> result that matches the input edge directions" (per `TriangulateIdxHalfedges` docs) with the
> **same topology** as C++. The genus delta = a degenerate-input triangulation divergence in
> `polygon_earclip.rs`. Same blocker for `convex_convex_minkowski_difference` and
> `openscad_crash` (the latter's C++ test sets `processOverlaps=true` only to silence that
> debug assert on its known-overlapping input).

### SDF thin-shell marching topology (1)
`sdf_sphere_shell` — genus 9560 vs expected ~14235 (perf now fine, 12s in release after the
Vec fix). The half-voxel-thick shell exposes a marching-tetrahedra crossing/snap topology
difference vs C++. Real correctness bug; needs root-cause in the SDF crossing logic.

### Coplanar/coincident-face geometry (2)
`almost_coplanar`, `convex_convex_minkowski_difference` (collapse_edge exposes a
boolean-pipeline geometry difference).
> **Fully root-caused (2026-05) — `almost_coplanar` is a `rotate` (trig) bug. The Boolean3
> SOS kernel is CORRECT.** Traced sign-by-sign against the C++ reference (built locally):
> - The spurious 21st vert is `xv12[10]` in the 2nd union (M1 edge 25 × tet face 2). It comes
>   from `kernel02(vertA=16, face2)`: C++ returns `s=1` (so the two edge endpoints net
>   `x12=0`, no crossing) but Rust returns `s=0` (net `x12=1`, spurious crossing). vertA=16
>   lies *exactly on* the face plane (x−y+z=1), so the final `Shadows(vertA.z, z02, …)` is an
>   **exact tie** in C++ (`z02 == vertA.z`) but in Rust `vertA.z − z02 = 8.3e-17`, flipping the
>   sign. Inputs to that kernel call (positions, normals, face) were **bit-identical** to C++ —
>   interpolate/shadow01 are bit-exact.
> - The 8.3e-17 traces back to **M1's intersection verts differing from C++ by ~2e-14**, which
>   traces back to **`rotatedTet` differing from C++ by 1–2 ULP**. Root: Rust `Manifold::rotate`
>   uses **quaternions** (sin/cos of half-angles in radians); C++ `CsgNode::Rotate` builds
>   degree-based `sind`/`cosd` axis matrices composed `rZ*rY*rX` (`csg_tree.cpp`).
> - **Verified fix:** porting `rotate` to the C++ construction makes `rotatedTet` bit-exact and
>   `almost_coplanar` pass (20/36). Batch-vs-pairwise ordering was ruled out (the C++
>   `BatchBoolean` size-heap reduces to the same `(tet∪rotated)∪tet` as Rust's eager `+`).
>
> **Why it's not landed yet — one separate blocker (needs a dedicated session):**
> - The matrix `rotate` needs `sind` to use **remquo (round-to-nearest)** reduction, not
>   `floor` (so `cosd(small)=sind(90+small)` matches C++). The *global* `sind` change regresses
>   `RefineQuads` (cylinder count 17044→16892: floor-`sind` happens to refine to C++'s count,
>   global remquo-`sind` doesn't). **Solution found: a *local* remquo-`sind` inside `rotate`
>   only — VERIFIED to fix `almost_coplanar` with no cylinder regression.**
> - The real remaining blocker: `menger_sponge` uses **exact 90° rotations**
>   (`hole.rotate(90,0,0)`). With the bit-exact matrix `rotate`, `rotate(90)` yields exact
>   0/1/−1 entries (vs the quaternion's ~2e-16 error), so the hole becomes **exactly coplanar**
>   with the grid. The difference-boolean's **`face2tri`** (boolean-result face triangulation,
>   `face_op.rs`) then emits an **unpaired halfedge** on that degenerate input — so the boolean
>   output is **non-manifold at simplify entry** (instrumented: a dangling halfedge exists
>   *before any collapse*; `cleanup_topology` doesn't fix it; C++'s `CleanupTopology` asserts
>   `IsManifold()` at entry, i.e. C++'s output IS manifold). It is **NOT** `collapse_edge` —
>   the `-1` OOB at `edge_op.rs:247` is just where the pre-existing dangling is first walked.
>   The fix is in the boolean-result triangulation/pairing for degenerate-coplanar faces — the
>   **same root area as the overlapping-triangulation bucket above** (minkowski genus). Fixing
>   that likely clears `almost_coplanar` + the 4 minkowski + `convex_convex` together.

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
