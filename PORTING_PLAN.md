# manifold-rust — open work

The open work on [Manifold](https://github.com/elalish/manifold)'s Rust port and on the
robust boolean engine: for each item, what a caller would notice and what is needed. How
the code works is in its header and doc comments; history is in `git log`.

Every item is judged by CLAUDE.md's "Philosophy": **accuracy and utility over C++
parity** — the port's own test suite, not the C++ reference, is the bar, and a fix whose
only benefit is matching the C++ more closely is not pursued for that alone — and
**output changes are documented** at the code that differs and in the commit, for
manifold-sharp, which verifies against this tree bit-for-bit.

The C++ reference is the `cpp-reference/manifold` submodule (v3.5.4).

---

## Exact engine open work

- **`refine` with tangents may group coplanar faces differently from C++
  `Impl::Refine`** (`smoothing.cpp:1120-1128`). The C++ clears tangents, then branches —
  with tangents, `SetNormalsAndCoplanar` + `CalculateBBox`; without,
  `CalculateVertNormals` — and only *then* calls `SortGeometry`, with no `SetEpsilon`
  anywhere. Our `manifold_smooth.rs` tail runs `calculate_bbox` + `set_epsilon`
  unconditionally and computes normals *after* `sort_geometry`. The tangent-free leg
  gives the same values either way (the sort is only a Morton relabelling); on the
  `had_tangents` leg, running `set_normals_and_coplanar` after the sort can group faces
  into different `coplanar_id` sets, since the flood fill walks faces in storage order,
  and no test asserts those groupings. Needed: decide which order is right and pin the
  groupings with a test, **coordinated with manifold-sharp**, whose `FinishRefine`
  transcribes our order and whose oracle lane compares meshes with no slack.

## Known mismatches with the C++

Open differences in output from C++ Manifold. Deliberate ones are noted at the code that
differs; these are not deliberate, and any fix that moves output changes manifold-sharp
in step.

- **Chained booleans evaluate eagerly.** `Manifold` holds its impl directly, so
  `a + b + c` runs pairwise booleans and `batch_boolean` is a left fold, where C++
  builds a lazy `CsgOpNode` tree and flattens it into one `BatchBoolean`/`BatchUnion`.
  A caller sees different triangle counts from the same expression: on
  `TEST(Smooth, Fillet)`'s inputs C++ refines to 3770 triangles and this port to 3822
  (C++ with the first union forced also gives 3822). Needed: make `Manifold` hold a CSG
  node and evaluate lazily; `src/csg_tree.rs` already ports the tree, but only the
  Minkowski sums and convex dilation use it. A progress-reporting batch boolean needs
  the same change.
- **`smooth_by_normals` + `refine_to_tolerance` drift by a few ULPs** from C++ once the
  union order is matched: volume `0x40be40dfe1b35e91` here against `0x40be40dfe1b7737f`.
  Both meet the C++ test's own `EXPECT_NEAR` bounds. Needed: an instrumented trace
  against C++ to find the first differing operation.
- **`swap_degenerates` cuts solid beside a zero-thickness fin.** The exact union of
  `src/testdata/minkowski-641145-union-{a,b}.txt` loses solid of an operand; C++ v3.5.2
  does the same. `set_normals_and_coplanar`'s orientation-blind flood fill hands a sound
  triangle a reversed seed normal, and `recursive_edge_swap` then cuts a wedge. Two
  narrow fixes broke `test_cpp_simplify` and `test_cpp_nonconvex_convex_minkowski_sum`.
  Pinned, ignored, by
  `minkowski::union_regression_tests::exact_union_of_thingi641145_partial_unions_contains_both_operands`,
  shared with manifold-sharp. Needed: an orientation-aware seed for the coplanar flood
  fill that keeps both C++ tests passing.
- **`boolean_with_token` composes the disjoint operands of an Add** (the non-overlapping
  fast path in `src/boolean3.rs`). C++ `Boolean3` has no such shortcut, so when leaf
  boxes overlap but the meshes don't, the order of `vert_properties`, `tri_verts` and
  `face_id` differs from C++ (example:
  `manifold::tests::compose::test_batch_union_of_rotated_instanced_leaves_is_independent_of_scheduling`;
  the C++ hash is `0x5a5fb492e1484228`). Stays until manifold-sharp changes together
  with this port.

## Robust engine open work

`src/robust/mod.rs`'s header says how the engine is verified on Thingi10K
(`examples/thingi_sweep.rs`, `examples/volume_referee.rs`).

- **One Thingi10K mesh still times out.** A robust union of #252784 (889 k triangles,
  heavily self-intersecting) with its rotated copy runs past a 570 s budget in the
  arrangements and CDT, sequential or `parallel`. Levers: a cheaper exact tier for
  genuine degeneracies (a Shewchuk-style adaptive expansion; 55% of exact orient3d calls
  overflow the i128 tier through shared-scale inflation), less duplication in the
  per-triangle prims, chunking the piece materialization, walk-based CDT point location
  and a sweep-line crossings pass for triangles carrying hundreds of thousands of
  segments.
- **Teardown is slow on giant inputs.** Dropping millions of rational prims through the
  system allocator takes 8–44 s, paid after successful runs as well as after a cancel.
  Needed: arena allocation for the arrangement's rationals.
- **Cancel overshoots on giant inputs.** Under `parallel` the rayon inner loops don't
  check the token (#252784 ran 867 s wall at a 570 s budget); sequentially a
  later-phase stretch overshoots by ~100 s at that budget, and the CDT's exact-rational
  sort cannot be interrupted. Needed: token checks in the parallel loops and the
  uncovered phases.
- **`triangulate_pseudo` recurses** (`src/robust/cdt_constraints.rs`), so an
  arrangement of ~15 k points in one triangle risks the stack. Needed: an explicit work
  stack.

---

## Performance open work

- **Peak memory is ~10% above C++** on large booleans (1.47 GB against 1.34 GB peak
  working set on the 2 M-triangle sphere difference). The excess is transient:
  `face2tri` temporaries and `intersect12`'s per-edge `Vec` of `Vec`s, where C++ uses
  counted two-pass output. Measure with
  `MANIFOLD_TIMING=1 cargo run --release --example mem_profile [round]`.
- **More `parallel` sites.** Parallel stable sorts in `sort_geometry`, parallel
  `create_properties` barycentrics and a parallel collider build would speed up the
  `parallel` feature further; each must keep results bit-identical to the sequential
  build.

## Out of scope (do not re-litigate)

- ExecutionContext's C++ *shape* (#1663/#1669/#1699/#1704): `src/cancel.rs` and
  `src/progress.rs` cover cancellation and progress natively; the C++ ctx-aware
  factories and per-phase timing prints are not ported.
- Halfedge refactoring (#1709): structural, no behavior change.
- CrossSection backend selector (#1710), concurrent-const-access safety (#1636), the
  affinity_partitioner race (#1664): threading concerns that do not arise in this port.
