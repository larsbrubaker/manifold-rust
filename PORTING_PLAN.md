# Manifold C++ → Rust Porting Plan

This document tracks the incremental port of [Manifold](https://github.com/elalish/manifold) to Rust.
Every phase must pass all tests with **exact numerical match** to the C++ before the next phase begins.

## Guiding Principles

1. **Exact match** — Same floating-point results as C++. Instrument both implementations to compare.
2. **Phase by phase** — Complete each phase fully before starting the next.
3. **No stubs** — Every function implemented, no `todo!()`, no `unimplemented!()`.
4. **Tests first** — Port the relevant C++ tests to Rust before or alongside implementation.
5. **Dependency ordered** — Never implement a function before its dependencies exist.

## Key Architectural Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Linear algebra | Roll our own `linalg.rs` | Must control exact f64 semantics; no `glam` |
| Float precision | f64 throughout; f32 only at MeshGL boundary | Matches C++ double/float usage |
| `Vec<T>` / `SharedVec<T>` | `Vec<T>` / `Arc<Vec<T>>` | Rust ownership eliminates need for async dealloc |
| Parallelism | Sequential first; `rayon` behind `parallel` feature | Verify correctness before parallelizing |
| Radix sort | Custom impl in `src/sort.rs` | rayon has no radix sort; needed for Morton codes |
| Cross-section (Clipper2) | `clipper2-rust` crate | Our own Rust port; same API surface |
| Halfedge mesh | Index-based (`Vec<Halfedge>` + usize) | Rust ownership; same pattern as clipper2-rust |
| BVH/Collider | Direct port (radix tree + Morton codes) | Pure algorithmic; maps cleanly |
| WASM demo | Three.js + WebGL for 3D viewer | 8 demo pages, deployed to GitHub Pages |

## C++ Codebase Size Reference

| File | Lines |
|------|-------|
| `include/manifold/linalg.h` | 2590 |
| `src/vec.h` | 366 |
| `src/parallel.h` | 1171 |
| `src/impl.h` | 562 |
| `src/impl.cpp` | 919 |
| `src/iters.h` | 327 |
| `src/svd.h` | 304 |
| `src/collider.h` | 370 |
| `src/shared.h` | 211 |
| `src/utils.h` | 170 |
| `src/polygon.cpp` | 1006 |
| `src/sort.cpp` | 532 |
| `src/constructors.cpp` | 545 |
| `src/face_op.cpp` | 346 |
| `src/edge_op.cpp` | 972 |
| `src/properties.cpp` | 464 |
| `src/smoothing.cpp` | 996 |
| `src/boolean3.cpp` | 531 |
| `src/boolean_result.cpp` | 889 |
| `src/csg_tree.cpp` | 764 |
| `src/manifold.cpp` | 976 |
| `src/tree2d.cpp` / `tree2d.h` | 146 |
| `src/cross_section/cross_section.cpp` | 802 |
| `src/subdivision.cpp` | 811 |
| `src/sdf.cpp` | 538 |
| `src/minkowski.cpp` | 175 |
| `src/quickhull.cpp` + `quickhull.h` | 1131 |
| **Total src** | ~15,962 |
| **Test files** | ~6,490 |

## Phases

### Phase 0: Project Setup ✅
- [x] GitHub repo created: `larsbrubaker/manifold-rust`
- [x] C++ manifold added as git submodule at `cpp-reference/manifold`
- [x] Cargo workspace with `manifold-rust` lib + `demo/wasm` WASM crate
- [x] GitHub Actions workflow for demo deployment
- [x] CLAUDE.md with porting rules
- [x] PORTING_PLAN.md (this file)

---

### Phase 1: Linear Algebra & Vector Types
**C++ sources:** `include/manifold/linalg.h` (2590 lines), `src/vec.h` (366 lines), `include/manifold/vec_view.h` (149 lines)
**Rust module:** `src/linalg.rs`, `src/vec.rs`

Key types to port:
- `vec2`, `vec3`, `vec4` (float column vectors)
- `ivec2`, `ivec3`, `ivec4` (integer vectors)
- `bvec2`, `bvec3`, `bvec4` (bool vectors)
- `mat2x2`, `mat3x3`, `mat4x4`, `mat3x4` (matrices)
- All arithmetic operators, dot, cross, normalize, length, min, max, abs, clamp, mix
- `quat` (quaternion)
- `VecView<T>` / `VecDH<T>` — manifold's GPU-capable vector (port as `Vec<T>` wrapper)

Tests: port `polygon_test.cpp` math-related tests + add unit tests for each operator.

**Status:** ⬜ Not started

---

### Phase 2: Core Types & Common Definitions
**C++ sources:** `include/manifold/common.h` (594 lines), `include/manifold/polygon.h` (61 lines)
**Rust module:** `src/types.rs`

Key types:
- `Mesh`, `MeshGL`, `MeshGL64` — triangle mesh representations
- `Properties`, `BaryRef`, `MeshRelation`
- `Box`, `Smoothness`, `OpType`
- Error types, `Quality` settings

**Status:** ⬜ Not started

---

### Phase 3: Polygon Triangulation
**C++ sources:** `src/polygon.cpp` (1006 lines)
**Rust module:** `src/polygon.rs`

Key algorithms:
- Triangulation of simple polygons (ear-clipping / monotone subdivision)
- Polygon simplification
- `Triangulate()` — the main entry point

Tests: port `test/polygon_test.cpp` (134 lines)

**Numerical validation:** Compare triangulation output vertex-by-vertex with C++ output.

**Status:** ⬜ Not started

---

### Phase 4: Mesh Data Structure (Impl)
**C++ sources:** `src/impl.h` (562 lines), `src/impl.cpp` (919 lines), `src/iters.h` (327 lines), `src/shared.h` (211 lines), `src/utils.h` (170 lines)
**Rust module:** `src/impl_mesh.rs`, `src/iters.rs`

Key structures:
- `ManifoldImpl` — the core internal mesh representation
- Half-edge data structure
- Vertex, face, edge connectivity
- Bounding box computation

**Status:** ⬜ Not started

---

### Phase 5: Sort & Parallel Utilities
**C++ sources:** `src/sort.cpp` (532 lines), `src/parallel.h` (1171 lines)
**Rust module:** `src/sort.rs`

Note: Manifold uses Thrust (GPU) or TBB/OpenMP for parallelism. We port to:
- Sequential first, verified correct
- Then add `rayon` for parallel passes where beneficial
- All outputs must be deterministic and match C++ sequential results

**Status:** ⬜ Not started

---

### Phase 6: Constructors
**C++ sources:** `src/constructors.cpp` (545 lines)
**Rust module:** `src/constructors.rs`

Primitive mesh generators:
- `Sphere(radius, circular_segments)`
- `Cube(size, center)`
- `Cylinder(height, radius_low, radius_high, circular_segments, center)`
- `Tetrahedron()`
- `Extrude(cross_section, height, ...)`
- `Revolve(cross_section, ...)`
- `LevelSet()` / `SDF`-based

Tests: port `test/manifold_test.cpp` constructor tests.

**Status:** ⬜ Not started

---

### Phase 7: Face & Edge Operations
**C++ sources:** `src/face_op.cpp` (346 lines), `src/edge_op.cpp` (972 lines)
**Rust module:** `src/face_op.rs`, `src/edge_op.rs`

Operations:
- Face normal computation
- Edge collapse, split, flip
- Mesh cleanup and manifold repair

**Status:** ⬜ Not started

---

### Phase 8: Properties
**C++ sources:** `src/properties.cpp` (464 lines)
**Rust module:** `src/properties.rs`

Mesh property calculations:
- Volume, surface area
- Genus, is-manifold checks
- Bounding box
- Normal computation and propagation

Tests: port `test/properties_test.cpp` (302 lines)

**Status:** ⬜ Not started

---

### Phase 9: SVD & Smoothing
**C++ sources:** `src/svd.h` (304 lines), `src/smoothing.cpp` (996 lines)
**Rust module:** `src/svd.rs`, `src/smoothing.rs`

Algorithms:
- Singular Value Decomposition (3×3)
- Tangent space computation
- Smooth normal interpolation
- Halfedge tangent vectors

Tests: port `test/smooth_test.cpp` (390 lines)

**Numerical validation:** SVD results must match exactly — critical for smooth subdivision.

**Status:** ⬜ Not started

---

### Phase 10: Collider & Spatial Indexing
**C++ sources:** `src/collider.h` (370 lines), `src/tree2d.cpp` + `tree2d.h` (146 lines)
**Rust module:** `src/collider.rs`, `src/tree2d.rs`

Algorithms:
- Bounding volume hierarchy (BVH) for triangle collision detection
- 2D R-tree for cross-section operations
- Ray-triangle intersection

**Status:** ⬜ Not started

---

### Phase 11: Boolean Operations (Core)
**C++ sources:** `src/boolean3.cpp` (531 lines), `src/boolean3.h`, `src/boolean_result.cpp` (889 lines)
**Rust module:** `src/boolean3.rs`, `src/boolean_result.rs`

This is the heart of the library:
- 3D mesh intersection computation
- Half-edge boolean result assembly
- Union, Intersection, Difference

Tests: port `test/boolean_test.cpp` (744 lines), `test/boolean_complex_test.cpp` (1558 lines)

**Numerical validation:** Critical — compare intermediate intersection points with C++ output.

**Status:** ⬜ Not started

---

### Phase 12: CSG Tree
**C++ sources:** `src/csg_tree.cpp` (764 lines), `src/csg_tree.h`
**Rust module:** `src/csg_tree.rs`

- Lazy CSG tree evaluation
- Batch boolean operations
- Tree simplification and caching

**Status:** ⬜ Not started

---

### Phase 13: Public Manifold API
**C++ sources:** `src/manifold.cpp` (976 lines), `include/manifold/manifold.h` (545 lines)
**Rust module:** `src/manifold.rs` (public API)

The user-facing API:
- `Manifold::Boolean(other, op)`
- `Manifold::Translate(v)`, `Rotate(q)`, `Scale(v)`, `Transform(m)`
- `Manifold::GetMesh()`, `NumVert()`, `NumTri()`, etc.
- `Manifold::Decompose()`, `Compose()`
- `Manifold::Status()` (error reporting)

Tests: port `test/manifold_test.cpp` (1156 lines)

**Status:** ⬜ Not started

---

### Phase 14: Cross Section (2D)
**C++ sources:** `src/cross_section/cross_section.cpp` (802 lines), `include/manifold/cross_section.h` (184 lines)
**Rust module:** `src/cross_section.rs`

2D polygon operations (uses Clipper2 internally in C++ — we use clipper2-rust):
- Boolean operations on 2D polygons
- Offset / Minkowski sum
- Area, bounds

Tests: port `test/cross_section_test.cpp` (265 lines)

**Dependency:** `clipper2-rust` crate (already exists)

**Status:** ⬜ Not started

---

### Phase 15: Subdivision
**C++ sources:** `src/subdivision.cpp` (811 lines)
**Rust module:** `src/subdivision.rs`

- Loop subdivision
- Refine by edge length
- Warp functions

**Status:** ⬜ Not started

---

### Phase 16: SDF Mesh Generation
**C++ sources:** `src/sdf.cpp` (538 lines)
**Rust module:** `src/sdf.rs`

- Marching cubes / level-set meshing
- `LevelSet(sdf, bounds, edge_length)`

Tests: port `test/sdf_test.cpp` (214 lines)

**Status:** ⬜ Not started

---

### Phase 17: Minkowski & Quickhull
**C++ sources:** `src/minkowski.cpp` (175 lines), `src/quickhull.cpp` (842 lines), `src/quickhull.h` (289 lines)
**Rust module:** `src/minkowski.rs`, `src/quickhull.rs`

- Minkowski sum/difference of convex meshes
- 3D convex hull

Tests: port `test/hull_test.cpp` (264 lines)

**Status:** ⬜ Not started

---

### Phase 18: WASM Demo
**Target:** Interactive demo site at `https://larsbrubaker.github.io/manifold-rust/`

Demo pages (mirroring clipper2-rust demo style):
1. Basic boolean operations (union, intersection, difference of primitives)
2. CSG tree demo (complex multi-step operations)
3. Smooth subdivision showcase
4. SDF mesh generation
5. Cross-section extrude/revolve
6. Convex hull
7. Performance benchmark in browser
8. Interactive 3D viewer (Three.js or Babylon.js)

**Status:** ⬜ Not started (begins after Phase 13 at minimum)

---

## Numerical Validation Strategy

For each phase, we validate exact match by:

1. **Build C++ with instrumentation** — add `printf` / logging at key points in the C++ source
2. **Run C++ tests** — capture intermediate values
3. **Run Rust tests** — compare output
4. **Diff** — any mismatch is a bug in the Rust port, not in the test

### Building the C++ reference

```bash
cd cpp-reference/manifold
cmake -B build -DMANIFOLD_TEST=ON
cmake --build build
./build/test/manifold_test
```

### Comparing outputs

Key values to compare per test:
- Number of vertices, triangles
- Bounding box (min/max xyz)
- Volume and surface area
- Individual vertex positions (sorted)
- Individual triangle indices (canonicalized)

---

## Dependency Graph

```
Phase 1 (linalg)
    └─► Phase 2 (types)
            ├─► Phase 3 (polygon)
            └─► Phase 4 (impl)
                    ├─► Phase 5 (sort)
                    ├─► Phase 6 (constructors)
                    ├─► Phase 7 (face_op, edge_op)
                    ├─► Phase 8 (properties)
                    ├─► Phase 9 (smoothing)  ← needs svd
                    └─► Phase 10 (collider)
                            └─► Phase 11 (boolean3, boolean_result)
                                    └─► Phase 12 (csg_tree)
                                            └─► Phase 13 (manifold API) ✓ usable
                                                    ├─► Phase 14 (cross_section)
                                                    ├─► Phase 15 (subdivision)
                                                    ├─► Phase 16 (sdf)
                                                    ├─► Phase 17 (minkowski, quickhull)
                                                    └─► Phase 18 (WASM demo)
```

---

## Progress Summary

| Phase | Module(s) | C++ Lines | Status |
|-------|-----------|-----------|--------|
| 0 | Setup | — | ✅ Done |
| 1 | linalg, vec | ~3,100 | ⬜ |
| 2 | types, common | ~655 | ⬜ |
| 3 | polygon | ~1,100 | ⬜ |
| 4 | impl | ~2,000 | ⬜ |
| 5 | sort, parallel | ~1,700 | ⬜ |
| 6 | constructors | ~545 | ⬜ |
| 7 | face_op, edge_op | ~1,318 | ⬜ |
| 8 | properties | ~464 | ⬜ |
| 9 | svd, smoothing | ~1,300 | ⬜ |
| 10 | collider, tree2d | ~516 | ⬜ |
| 11 | boolean3, boolean_result | ~1,420 | ⬜ |
| 12 | csg_tree | ~764 | ⬜ |
| 13 | manifold API | ~1,521 | ⬜ |
| 14 | cross_section | ~986 | ⬜ |
| 15 | subdivision | ~811 | ⬜ |
| 16 | sdf | ~538 | ⬜ |
| 17 | minkowski, quickhull | ~1,306 | ⬜ |
| 18 | WASM demo | — | ⬜ |
