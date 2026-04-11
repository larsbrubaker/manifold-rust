# Manifold C++ → Rust Porting Plan

This document tracks the incremental port of [Manifold](https://github.com/elalish/manifold) to Rust.
Every phase must pass all tests with **exact numerical match** to the C++ before the next phase begins.

## Current Status: 331 tests passing, 0 failing, 15 ignored

~55% of C++ tests ported. ~90 remaining.

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
| Parallelism | Sequential first; `rayon` behind `parallel` feature | Verify correctness before parallelizing |
| Cross-section | `clipper2-rust` crate | Our own Rust port; same API surface |
| BVH/Collider | Radix tree + Morton codes (sorted internally) | Matches C++ algorithm exactly |

---

## Completed Phases (0–11, 14, QuickHull, WASM)

| Phase | Module(s) | Status |
|-------|-----------|--------|
| 0 | Project setup | ✅ |
| 1 | linalg, vec | ✅ |
| 2 | types, common | ✅ |
| 3 | polygon | ✅ |
| 4 | impl_mesh | ✅ |
| 5 | sort, parallel | ✅ |
| 6 | constructors | ✅ |
| 7 | face_op, edge_op | ✅ |
| 8 | properties | ✅ |
| 9 | svd, smoothing | ✅ |
| 10 | collider (BVH) | ✅ |
| 11 | boolean3, boolean_result | ✅ |
| 14 | cross_section | ✅ |
| — | quickhull | ✅ |
| — | disjoint_sets | ✅ |
| 18 | WASM demo (10 pages) | ✅ |

---

## Unported C++ Tests — The Roadmap

### Priority 1: Properties tests (11 unported)
These test existing functionality (min_gap, triangle distance) that should mostly work now.

- [ ] TriangleDistanceClosestPointsOnVertices
- [ ] TriangleDistanceClosestPointOnEdge
- [ ] TriangleDistanceClosestPointOnEdge2
- [ ] TriangleDistanceClosestPointOnFace
- [ ] TriangleDistanceOverlapping
- [ ] MinGapCubeSphereOverlapping
- [ ] MinGapAfterTransformationsOutOfBounds
- [ ] MingapStretchyBracelet

### Priority 2: Boolean tests (5 unported)
Fill gaps in boolean coverage.

- [ ] AlmostCoplanar
- [ ] Coplanar
- [ ] CornerUnion
- [ ] Perturb1
- [ ] TreeTransforms

### Priority 3: BooleanComplex tests (4 unported)
Complex boolean scenarios with OBJ meshes.

- [ ] GenericTwinBooleanTest7863
- [ ] InterpolatedNormals
- [ ] Ring

### Priority 4: Smooth tests (11 unported) — requires Phase 15 completion
Smooth subdivision with Bezier tangents.

- [ ] Csaszar, Manual, Normals, FacetedNormals, NormalTransform
- [ ] RefineQuads, SineSurface, Tetrahedron
- [ ] ToLength, Torus, TruncatedCone

### Priority 5: Manifold API tests (32 unported)
Covers MeshGL round-trip, mesh relations, decompose, invalid inputs.

- [ ] GetMeshGL, MeshGLRoundTrip, ObjRoundTrip
- [ ] MeshRelation, MeshRelationRefine, MeshRelationRefinePrecision, MeshRelationTransform
- [ ] FaceIDRoundTrip, MeshID, Merge, MergeDegenerates, MergeEmpty, MergeRefine
- [ ] Decompose, DecomposeProps
- [ ] Invalid, InvalidInput1–7
- [ ] MeshDeterminism, Warp2, WarpBatch
- [ ] Project, Slice, SliceEmptyObject
- [ ] OpenscadCrash

### Priority 6: CrossSection tests (6 unported)

- [ ] BevelOffset, FillRule, HullError, MirrorCheckAxis, NegativeOffset, Rect

### Priority 7: SDF tests (9 unported) — requires Phase 16 completion

- [ ] Blobs, Bounds, Bounds2, Bounds3, CubeVoid, Resize, SineSurface, SphereShell, Void

### Priority 8: Sample tests (12 unported) — integration tests

- [ ] Bracelet, CondensedMatter16/64, Frame, FrameReduced
- [ ] GyroidModule, Knot13/42, Scallop, Sponge1/4, TetPuzzle

---

## Incomplete Phases

### Phase 12: CSG Tree — ⚠️ Minimal
**Status:** Basic recursive `evaluate()` only.
**Missing:** Lazy transform propagation, caching, BatchBoolean with priority queue, BatchUnion with bbox partitioning, Compose for efficient disjoint merge.

### Phase 13: Public Manifold API — ⚠️ Core Working
**Status:** All constructors, transforms, booleans, smoothing, hull, split working.
**Missing:** `simplify()`, `set_properties()`, `Decompose`, per-mesh epsilon tracking, faceID in MeshGL.

### Phase 15: Subdivision — ⚠️ Midpoint Only
**Status:** Simple 4x midpoint subdivision, sphere from subdivided octahedron.
**Missing:** Partition class, Bezier vertex placement, curvature-based refine_to_tolerance, edge length-based refine_to_length, property interpolation, quad support.

### Phase 16: SDF — ❌ Not Started
**Missing:** Marching Tetrahedra on BCC grid, ITP root-finding, surface detection, sparse grid storage.

### Phase 17: Minkowski — ❌ Not Started
**Missing:** Convex-Convex (pairwise vertex sums → hull), Non-Convex (per-triangle sums → hull → BatchBoolean).

---

## Ignored Tests (15)

| Test | Reason |
|------|--------|
| 4× Minkowski | Phase 17 not started |
| hull_tictac | QuickHull edge case |
| hull_degenerate_2d | QuickHull 2D degenerate |
| properties_tolerance | Smooth subdivision needed |
| tolerance_sphere | Smooth subdivision needed |
| min_gap_transformed | Slow in debug (passes in release) |
| boolean_precision | Per-mesh epsilon not implemented |
| create_properties_slow | Slow sphere boolean in debug |
| complex_close | Slow 256-seg sphere ×10 boolean |
| perturb3 | BatchBoolean precision |
| craycloud | OBJ loader needs vertex merge |
| opposite_face | MeshGL degenerate face import |

---

## Future: Performance Parity

After all tests pass, validate Rust performance matches C++:
- Benchmark identical operations (sphere booleans, complex OBJ booleans)
- Profile hotspots (collider queries, boolean result assembly)
- Add rayon parallelism behind feature flag
- Target: within 2× of C++ for all operations
