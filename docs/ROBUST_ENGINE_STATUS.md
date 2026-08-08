# Robust boolean engine — status

Snapshot for picking this work up on another machine. Current as of `1bf570c`
on `main`.

## What changed

The robust engine's classification was rewritten onto the **mesh-arrangement
formulation** (Zhou, Grinspun, Zorin, Jacobson 2016 — what libigl's
`mesh_boolean` uses). It replaces both the Barki 2015 paper's local Prop 2/3
ring walk and this engine's own per-component winding queries, neither of
which could handle the full Thingi10K corpus.

The pipeline is now:

```
intersection_graph  →  cells::build_cells  →  cells::windings  →  cells::extract  →  assemble
   (arrangement)        (cell complex)        (winding/cell)      (predicate)
```

Three properties carry the robustness, all in `src/robust/cells.rs`:

1. **Winding propagates combinatorially**, cell to cell. Every cell's winding
   derives from one traversal, so adjacent regions *cannot* disagree. The old
   engine's `NotClosed` failures were exactly two independent queries
   contradicting each other across a shared intersection segment.
2. **Orientation is derived from cell labels**, never inherited from the input
   face. A wall is emitted iff its two cells disagree about containment, wound
   so its normal points out of the solid. Inverted regions of self-intersecting
   scans therefore cost nothing.
3. **Coincident stacks add their winding steps.** A doubled sheet steps by two,
   a fold cancels to zero. There is no regularization pass; thin material
   cancels arithmetically.

Each operand's solid is `{w ≥ 1}` — a negative winding is inverted geometry and
never material. Every operation is one predicate on the winding vector, so
`Subtract` needs no operand flip.

`classify.rs` and `propagate.rs` were deleted; their jobs no longer exist.

## Current state

All 643 library tests pass. `cargo test --release --lib`.

Full-corpus sweep, **3781 of 9934 meshes** (run was still going):

| metric | value |
|---|---|
| robust errors | 32 (0.85%) — **30 are 60s timeouts, only 2 genuine** |
| volume agrees with exact engine | 2661 |
| volume disagrees | 568 (17.6% of compared) |
| not comparable (one engine errored) | 552 |
| import failures | 518 (13.7%) |

Against the previous engine on a 1076-mesh overlap, robust errors nearly
halved (36 → 19) with 728 agree→agree.

Performance is **~1.44×** the previous robust engine. The cost is diffuse;
two targeted optimizations after the first (skipping exact radial work on
ordinary manifold edges) bought ~1% each, so profiling is needed before
optimizing further rather than more guessing.

## Open items, highest value first

1. **568 volume disagreements (17.6%).** The largest unknown. These are *not*
   new — the pre-swap engine had a comparable rate on its overlap. Unresolved
   whether the robust or the exact engine is wrong; the exact engine is the
   nominal reference but these are self-intersecting inputs where its
   assumptions are shaky. Start by bucketing: pick a handful, check whether the
   robust result is closed and self-consistent (`cells::inconsistent_walls`),
   and compare against a third opinion.
2. **30 timeouts at 60s.** Sweep with a longer `--timeout-secs` to see whether
   they complete at all, then profile the worst.
3. **Performance.** 1.44× wants a profiler, not another guess.
4. **`is_strictly_non_delaunay` in `cdt.rs`** is dead — the one remaining
   compiler warning.

## Reproducing the data

Thingi10K lives at `C:\Development\rust-apps\Thingi10K` (sibling repo, meshes
under `meshes/`, zipped cache under `mesh-cache/`).

```bash
cargo run --release --example thingi_sweep -- --db out.db --timeout-secs 60
```

Useful flags: `--limit N`, `--skip N`, `--ids A,B,C`, `--notes TEXT`. Results
go to SQLite (`runs` + `results` tables) keyed by run, so two runs can be
diffed directly. `*.db` is gitignored.

The comparison that matters is the verdict transition matrix between two runs —
`agree`/`disagree`/`na` per mesh — since a raw disagree count moves when
meshes stop erroring and enter the comparison for the first time. Counting
totals instead of transitions hid a real regression once already.

**Gate:** volume + closure. Area and triangle counts are advisory — derived
orientation and one-representative-per-stack legitimately change surface
tessellation without changing the solid.

## Regression tests worth knowing about

- `robust::thingi_tests` — the two real-world `NotClosed` repros (#59082 ∪
  #1313535, #74660 ∪ #1147177), plus arrangement-consistency checks on both
  pairs via `cells::inconsistent_walls`.
- `robust::cells::tests` — cell counts, winding regions, doubled-shell
  multiplicity, end-to-end boolean volumes, and an inverted-operand case.
- `robust::nonmanifold_tests::soup_survives_non_axis_aligned_rotation` — soup
  impls carry no collider; rotating one used to crash.

## Two bugs the unit tests could not have found

Both were caught only by the corpus sweep, which is the argument for running it.

- **Collider crash** (pre-existing): `soupify` skips `sort_geometry`, where the
  collider is built, so soup impls carry a zero-leaf tree; a non-axis-aligned
  transform refit it and indexed out of bounds. The guarding `debug_assert` is
  compiled out in release.
- **`outer_cell`** (introduced during this work): it *deduced* the unbounded
  cell from the lexicographically extreme vertex. That is wrong when the chosen
  face has material on its outward side, which real scans have. Anchoring an
  interior cell at zero shifts every winding by a constant, and `{w ≥ 1}` then
  excludes nearly the whole model — two unions collapsed from 51372 and 5856
  triangles to 4 and 40. The output still looked closed and well-oriented.
  Every component is now anchored by an exact point query instead.

The second is the cautionary one: a 175-mesh subset reported zero regressions
for it. Structural changes need the full corpus before the numbers mean
anything.
