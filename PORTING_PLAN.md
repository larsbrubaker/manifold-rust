# Manifold C++ → Rust Porting Plan

This is a **roadmap of remaining work** to finish porting
[Manifold](https://github.com/elalish/manifold) to Rust. It tracks what is left to do,
not what has already been done (use `git log` for history). Every change must reproduce the
C++ reference with **exact numerical match** — identical results on identical inputs.

**Status:** 507 passing, 0 failing, 18 ignored.
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

## Ignored tests (18) — grouped by the work needed to clear them

### Just slow in debug — pass in release; not bugs (9)
`nonconvex_convex_minkowski_sum/difference`, `nonconvex_nonconvex_minkowski_sum/difference`
(O(n²)), `sdf_blobs`, `sdf_sphere_shell`, `hull_sphere` (1500 seg), `hull_menger_sponge`
(depth-4), `properties_mingap_stretchy_bracelet`. → Speed up or run in release; no correctness
work.

### Coplanar merge / boolean-pipeline geometry (2)
`almost_coplanar` (21 vs 20 verts; C++ v3.5.0 still expects {20,36}, so this is a real
coplanar-merge difference in `SimplifyTopology`, *not* tolerance stacking — #1671 edge_op is
already ported and doesn't change it), `convex_convex_minkowski_difference` (collapse_edge
exposes a boolean-pipeline geometry difference).

### Boolean produces a non-manifold intermediate — deep robustness (4)
`complex_sweep`, `openscad_crash`, `complex_craycloud`, `craycloud_bool`. All four end in a
non-manifold intermediate during/after a boolean: either `update_vert` walks an unpaired
halfedge (`paired = -1` → `usize::MAX` OOB) or the sort.rs:298 odd-halfedge assertion fires.
**Diagnosed:** for complex_sweep, a probe at `simplify_topology` entry shows the boolean
output is fully manifold, so the `-1` is created *mid-simplify* during the collapse sequence;
`collapse_tri` / `form_loop` structurally match C++ (`PairUp`, same orbits), so it's a subtle
index/ordering divergence in the multi-collapse interaction. The craycloud pair is
FP-borderline (passes in release, fails in debug) — the OBJs are cleanly indexed (no
welding/import issue) and `processOverlaps` is unrelated (only a polygon-triangulation flag).
Next step: instrument Rust + the C++ reference to log each collapsed edge + resulting pairing
on these inputs and find the first divergence. (Avoid a `-1` guard band-aid — masks the bug,
risks incorrect geometry.)

### Boolean hangs (2)
`complex_generic_twin_7081`, `generic_twin_7081`. Loop-termination/convergence bug in boolean
(likely the same robustness area as above).

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
