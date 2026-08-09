# Robust boolean engine — status

Snapshot for picking this work up on another machine. Current as of `9fbc40a`
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

All 644 library tests pass. `cargo test --release --lib`.

Full-corpus sweep, **3781 of 9934 meshes** (an earlier partial run on the
other machine; a fresh full baseline is accumulating as run 2 of
`thingi_sweep.db` on this one):

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

## The volume disagreements are (mostly) the exact engine's fault

The former top unknown is resolved in direction. `examples/volume_referee.rs`
arbitrates a disagreement independently of both engines: it Monte-Carlo
samples `vol({w_A ≥ 1} ∪ {w_B ≥ 1})` of the *inputs* with the exact winding
query — the solid both engines claim to compute — and reports an estimate
with a standard error. On the first 37 genuine volume disagreements of the
fresh sweep:

- **Robust right 9, exact right 0.** Every arbitrated large disagreement went
  robust's way, including one where exact reported 0.297 for a model whose
  `{w ≥ 1}` material is literally empty (referee: 0 hits in 100k samples).
- The exact engine's failure signature is 2–3× **overcounting** on
  self-overlapping scans: its volume integrates winding, so doubly-wound
  regions count twice. The referee's `overlapA` column (divergence volume ÷
  once-counted material of operand A) makes this visible per mesh.
- The remaining ~26 were sub-1% differences inside sampling noise —
  tessellation-level drift, not classification errors.

Consequence: exact-vs-robust `match_ok` on self-intersecting inputs measures
conformance to a broken reference. Treat volume mismatches as suspect-exact
until the referee says otherwise.

The sweep's `match_ok=0` bucket also conflates gates: of run 2's first 76
mismatches, only 37 differed in volume; 18 differed only in area and 21 only
in tri/vert counts (both advisory under derived orientation).

## One real robust-side bug found and fixed

Thingi10K #36088 (440 tris) produced 4 `inconsistent_walls`: a coincident
doubled sheet's two copies triangulated an exactly-cocircular square with
*opposite diagonals* (the CDT tie fell to construction history), so the stack
split into vi-distinct walls declaring `[1,0]` steps across a `[2,0]`
crossing. Windings survived only because BFS reached those cells by other
paths. Fixed in `9350880`: cocircular flips now resolve toward a canonical
coordinate-rank diagonal, making the triangulation a function of coordinates
— coincident copies agree by construction. Regression test:
`thingi_36088_arrangement_is_consistent`.

## Batch-mining findings (2026-08-08 session)

The 8-hour full sweep was replaced by exponentially growing random batches
(1, 2, 4, … 512 meshes; `--ids $(shuf …)` against `thingi_sweep.db`), mining
each for failures before the next. What that surfaced and fixed, in order:

- **Cocircular CDT tie fix validated** — it also cured the one real
  NotClosed of the sequential window (#36374, now a regression test).
- **Cancel was unresponsive** (#42211 ran 565 s past a 60 s deadline) —
  token now checked per triangle and inside the arrangement sweeps
  (`build_graph_with_token`); the same mesh stops within a second.
- **Timeout family root-caused**: heavily self-intersecting scans put
  hundreds of segments in single triangles, and the per-triangle quadratic
  sweeps paid an exact predicate per pair. Conservative-box prefilters plus
  a hashed segment registry cut #42211's robust pass from ~580 s to ~145 s
  (crossings 205→18.5 s, candidate points 200→20 s, on-seg 83→6.9 s,
  registries 72→27 s). Outputs unchanged. Remaining profile: arrangements
  ~87 s (~39 s of it interning/dedup wrapper cost), candidate points 20 s.
- **Referee tally through batch 256: robust right 31, exact right 0.**
  Recurring exact-engine failure families: inverted meshes (referee finds
  zero `{w≥1}` material, exact reports the divergence integral), and
  N-times-wound shells (exact = N× the physical volume; the #138xxx family
  is exactly 3×). Robust's occasional dramatic *simplifications* are also
  correct: #88147 emits 430 tris vs exact's 2222 with identical volume and
  area (coincident stacks cancelled arithmetically).
- Batches 1 through 256: **zero robust-side hard failures at HEAD** (no
  panic, no NotClosed, no wrong-volume verdicts). Only Cancelled on
  genuinely huge inputs (≥500 k tris).

## The one robust-wrong verdict, fixed (`61b2dea`)

#301921 ∪ rotated-self was the referee's only robust-wrong call (0.4033 vs
true 0.4381). Root cause chain: the extracted boundary legitimately touches
itself along edges (2k half-edges per id-edge); the generic import paired
those fans arbitrarily; `cleanup_topology`'s `dedupe_edge`/`update_vert`
then walked the crosslinked pairing out of the vertex star and physically
moved geometry. `robust/pairing.rs` now pairs multi-edge fans geometrically
(same solid wedge) and splits pinched vertices, making cleanup
volume-neutral; the residual 0.6% (swap_degenerates acting on sign-flipped
flood-filled normals) is resolved by skipping the swap for robust-assembled
outputs — the first entry in `docs/CPP_DIVERGENCES.md` under the new
documented-divergence policy. Robust volume is now extraction-exact on the
repro and both engines agree with the referee.

## Auto engine is clean-by-default (`7d86b25`)

Auto routed to Exact for any manifold-connectivity operand, which shipped
wrong volumes for manifold-but-self-overlapping inputs (gallery frame
#92068 ∪ #39926: 15% overcount, operand overlap ratio 3.0). Auto now picks
Exact only when both operands are manifold AND `has_self_intersections()`
is false (cached exact self-scan, ~0.2 ms with collider reuse; duplicates
count as contact). The demo badges self-intersecting operands under Auto.

## Full-corpus verdict (2026-08-09)

Every mesh in the local Thingi10K corpus has now been swept at
current-era code (9,905 distinct meshes; the pre-fix runs-1/2 window was
re-swept after the CDT/pairing fixes landed). Corpus-wide, latest result
per mesh:

| status | count |
|---|---|
| NoError | 8,823 |
| Cancelled (120 s, large self-intersecting inputs) | 43 |
| **NotClosed** | **0** (was 174 at pre-fix code) |
| import failures (unreadable/unimportable, no engine attempt) | 1,039 |

Every volume disagreement above the noise floor was arbitrated by the
Monte-Carlo referee: **robust right in every case but one** (the fixed
#301921), across ~150 arbitrated meshes — the exact engine's winding
integral over- or under-counts on manifold-but-self-overlapping input,
which is also why `Auto` now guards on `has_self_intersections()`. One
cosmetic arbiter artifact: #79187 reports "neither" because robust's
1e-5 residual sliver exceeds the referee's zero-hit absolute band —
robust is right there to 1e-5 absolute.

## Open items

1. **43 timeouts** on very large self-intersecting meshes. Roadmap from
   instrumentation: 55% of exact orient3d calls fall through to BigInt
   because one near-zero full-mantissa coordinate inflates the shared
   dyadic scale — the levers are fewer filter escalations or a
   Shewchuk-style adaptive-expansion tier, plus (architectural) a
   small-value-inline bigint backend to cut allocation churn. The i64
   predicate tier (native-neutral, big wasm win) is already in.
2. **Progress API** is in (library → FFI 0.3.1 → wasm → demo busy card).
   A progress-capable *batch* boolean needs `csg_tree.rs`.
3. `robust/intersection_graph.rs` is over the 800-line limit (~1,110)
   and needs splitting before further additions.

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

**Arbitrating a disagreement:**

```bash
cargo run --release --example volume_referee -- --ids A,B,C --samples 100000
```

Prints the sampled ground-truth volume ±σ, both engines' answers, the
arrangement's inconsistent-wall count, and operand A's self-overlap ratio.
`--from-run N` pulls every `match_ok=0` mesh of a sweep run instead.

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
