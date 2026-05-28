# Manifold C++ → Rust Porting Plan

This is a **roadmap of remaining work** to finish porting
[Manifold](https://github.com/elalish/manifold) to Rust. It tracks what is left to do,
not what has already been done (use `git log` for history). Every change must reproduce the
C++ reference with **exact numerical match** — identical results on identical inputs.

**Status:** 495 passing, 0 failing, 22 ignored.
**C++ reference target:** v3.5.0 (submodule at tag `v3.5.0`, commit `541c33bd`).
**Core engine:** all 18 phases (linalg → boolean → CSG → cross-section → SDF → minkowski →
WASM) are implemented. Remaining work is the v3.5.0 deltas below plus the ignored-test
backlog.

---

## Remaining v3.5.0 work

### #1718 — Normals recorded on Manifold (largest remaining item)
Record normals on the Manifold and auto-substitute them on `get_mesh_gl`; round-trip via the
MeshGL `run_flags` hasNormals bit; negate normals for backside (subtractee) runs; apply the
#1602 consistent normal transformations. The `run_flags` plumbing already exists in
`types.rs` / `manifold_meshgl.rs`, but auto-substitution and backside negation do not.
- Also add `Error::InvalidTangents` (#1718) and `Error::Cancelled` (#1663) variants.
- **Unblocks** the ignored `boolean::test_cpp_normals` (its `RelatedGL` check runs with
  `checkNormals=true`). C++ ref: bd56861d (#1718), 46135097 (#1602).

### #1671 — edge_op / simplify part
The smoothing part of #1671 is done; the `edge_op.cpp` part is not: `CollapseEdge` gains
`tol`/`firstNewVert` params, reworked `CollapseShortEdges`/`CollapseColinearEdges`
("improved boolean edge collapse"), plus the properties.cpp simplify cleanup and "tolerance
stacking" fix.
- **May unblock** the tolerance/coplanar ignored tests (`almost_coplanar`,
  `properties_tolerance_sphere`). C++ ref: 58b328a9 (#1671).

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

## Ignored tests (22) — grouped by the work needed to clear them

### Just slow in debug — pass in release; not bugs (9)
`nonconvex_convex_minkowski_sum/difference`, `nonconvex_nonconvex_minkowski_sum/difference`
(O(n²)), `sdf_blobs`, `sdf_sphere_shell`, `hull_sphere` (1500 seg), `hull_menger_sponge`
(depth-4), `properties_mingap_stretchy_bracelet`. → Speed up or run in release; no correctness
work.

### Blocked on #1718 normals feature (1)
`boolean::test_cpp_normals`.

### N-way subdivision (2)
`sphere_tri_count_n25` (8192 vs 5000 tris), `properties_tolerance_sphere`. Current
subdivision is binary (always doubles edge count); C++ supports arbitrary n-way splits via
the `Partition` class. Non-power-of-2 segment counts need this.

### edge_op / simplify & tolerance (3) — overlaps with #1671 edge_op work
`almost_coplanar` (21 vs 20 verts — colinear-edge-collapse / epsilon difference in
`SimplifyTopology`), `properties_tolerance_sphere` (set_tolerance), `convex_convex_minkowski_difference`
(collapse_edge exposes a boolean-pipeline geometry difference).

### `edge_op::update_vert` robustness — needs processOverlaps (2)
`openscad_crash`, `complex_sweep`. `update_vert` walks `paired_halfedge` around a shared
vertex; an unpaired halfedge (paired = -1) indexes with `usize::MAX` and panics. C++ masks
this via `processOverlaps=true`; our boolean output produces a slightly non-manifold
intermediate for these inputs. Root-cause in boolean result assembly.

### Boolean hangs (2)
`complex_generic_twin_7081`, `generic_twin_7081`. Loop-termination/convergence bug in boolean.

### Non-manifold OBJ import (2)
`complex_craycloud`, `craycloud_bool` — OBJ loads as 663 halfedges (not a multiple of 6),
trips the sort.rs even-count assertion. See #1685 import work above.

### Hull vertex count (1)
`hull_tictac` (~66050 vs ~66038 verts).

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
