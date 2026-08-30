# Deliberate divergences from the C++ reference

Per CLAUDE.md, byte-matching `cpp-reference/manifold` is the default; the
entries below are the deliberate exceptions, each with the evidence that
justified it. Trace-diff debugging against the C++ must expect these.

Two kinds of entry live here, and they are not the same claim:

- **Justified divergence** — an accuracy fix, a real bug fix in the C++, or a
  measured improvement. CLAUDE.md's three qualifying reasons. Entry 1.
- **Inherited divergence, documented and scheduled for coordinated
  harmonization** — an output shape this port already shipped, which we would
  resolve toward the C++ on the merits but cannot change unilaterally, because a
  downstream consumer verifies against this tree bit-for-bit. These are
  *disclosures*, not justifications: the entry states what differs, why it
  cannot be fixed here alone, and what a coordinated fix would take. Entry 2.

The second kind is deliberately uncomfortable to write, which is the point — it
is a debt with a name attached, not a decision that ends the discussion. Nothing
belongs in either category for convenience.

## 1. Robust-engine outputs skip `swap_degenerates` (2026-08-08)

**What differs:** meshes assembled by the robust boolean engine
(`src/robust/assemble.rs`) do not run `edge_op::swap_degenerates` during
their import/simplification. The exact engine's pipeline is unchanged and
still byte-matches C++ (`Impl::SimplifyTopology`).

**Where:** `Manifold::from_mesh_gl64_robust_assembled`
(`src/manifold_meshgl.rs`, `pub(crate)`, used only by `robust::assemble`)
runs `cleanup_topology` + `collapse_short_edges` + `calculate_vert_normals`
in place of `edge_op::remove_degenerates`, and `assemble` composes the same
pieces plus `collapse_colinear_edges` in place of `simplify_topology`.
`edge_op`'s own functions are untouched; only the composition differs. Every
other import path, including user meshes arriving through
`from_mesh_gl_robust` / `from_mesh_gl64_robust`, is byte-identical to before.

**Why:** `face_op::set_normals_and_coplanar` (faithful port of C++
`Impl::SetNormalsAndCoplanar`, impl.cpp:214) flood-fills a seed triangle's
normal onto every coplanar neighbor with no orientation check. A boolean
result legitimately contains coplanar *antiparallel* adjacencies (material
below the plane on one side of a shared edge, above it on the other), so a
handful of large triangles receive sign-flipped normals. `swap_degenerates`
then misclassifies them as degenerate (`ccw <= 0`) and swaps large
non-planar quads — physically moving the surface. On Thingi10K #301921 ∪
rotated-self this moved the robust union volume by −2.5e-3 (and the exact
engine's own result by +5.6e-5; both engines' outputs are exposed, exact
merely got tessellation-lucky). Skipping the swap keeps robust volumes
extraction-exact. The cost is slightly larger tri/vert-count drift versus
the exact engine, which the sweep gate already treats as advisory.

**Evidence:** `robust::thingi_tests::thingi_301921_union_rotated_self_matches_exact`,
the staged volume traces in the session log (extraction 0.438010 preserved
through cleanup, lost only in swap), and the Monte-Carlo referee
(`examples/volume_referee.rs`).

## 2. The centered cylinder centers in place, not through `Translate().AsOriginal()` (2026-08-30)

**What differs:** `constructors::cylinder`'s `center` branch shifts `vert_pos.z`
in place and repairs the derived caches. The C++ re-centers by composing a
transform, which additionally assigns a fresh original ID, re-marks coplanar
faces, and recomputes every coordinate. Unlike entry 1 this is not an accuracy
fix — it is an **inherited** output shape being deliberately preserved, and the
one entry here that we would resolve toward the C++ if we could do it alone.

**Where:** `src/constructors.rs`, the `if center` block at the end of `cylinder`.
Reached by `Manifold::cylinder_centered(…, center: true)` and — via the recursive
`cylinder(height, radius_high, 0.0, n, true)` its cone branch starts from — every
cone, centered or not. `Manifold::cylinder` does *not* reach it:
`manifold_shape.rs:48` forwards to `cylinder_centered` with `center: false`.
C++ v3.5.2 is `src/constructors.cpp:155-157`:

```cpp
Manifold cylinder = Manifold::Extrude({circle}, height, 0, 0.0, vec2(scale));
if (center)
  cylinder = cylinder.Translate(vec3(0.0, 0.0, -height / 2.0)).AsOriginal();
return cylinder;
```

`AsOriginal` (`src/manifold.cpp:449-457`) copies the impl, then runs
`InitializeOriginal()` and `SetNormalsAndCoplanar()`.

**What is observably different.** Measured on `cylinder(4.0, 1.0, 1.0, 8, true)`
against the same mesh built uncentered and put through `transform`:

| | in place (ours) | through a transform (C++'s shape) |
|---|---|---|
| `mesh_relation.original_id` | `1` — extrude's | `-1` from the transform, then a fresh ID from `InitializeOriginal` |
| `epsilon` | `4e-12` | `3.999999999999999e-12` |
| `tolerance` | `4e-12` | `4e-12` |
| positions | 4 coordinates are `-0.0` | the same 4 are `+0.0`; **zero** value differences |

The signed zeros are on x and y, never z, and the mechanism is exactly that: the
in-place edit touches only z, so x and y keep the `-0.0` that `cosd`/`sind`
produced at the quarter angles, where a transform recomputes them as
`1.0*x + 0.0*y + 0.0*z + 0.0` and normalizes `-0.0` to `+0.0`. The epsilon gap is
one ULP, from `Impl::Transform` scaling epsilon by a spectral norm that comes back
`0.9999999999999998` rather than exactly `1.0` for a translation. Beyond the
table, `AsOriginal`'s `SetNormalsAndCoplanar()` re-marks coplanar faces, which we
do not re-run.

**Why we keep it.** This divergence **predates** the cache repair that sits at
the same lines: the in-place shim shipped without `AsOriginal`'s semantics, so
`cylinder_centered`'s originalID, epsilon and signed zeros have differed from the
C++ in every released 0.14.x. The repair deliberately fixed only the stale
caches — the face BVH was describing the pre-shift positions, and every boolean
against a centered cylinder tripped `pair_up`'s non-manifold assert — while
retaining the existing output semantics exactly. `sort_geometry` was the right
size of claim for that: it rebuilds the derived caches from the positions already
there, where re-centering through a transform would also have moved the positions.
`set_epsilon` is deliberately not re-run for the same reason, and that is what
keeps epsilon at `4e-12`.

Switching to the C++ form is therefore not a free correction but a breaking
behavioral change, and it has a downstream consumer that would break:
[manifold-sharp](https://github.com/larsbrubaker/manifold-sharp) is a pure C#
port of this crate whose test contract is bit-exactness against this tree, and
whose own ledger (`docs/RUST_DIVERGENCES.md`, entries 4 and 5) cites these
current outputs. Its oracle lane compares exported meshes row for row with no
slack, so a signed zero is a failure there.

**The two halves of `cylinder` disagree today, and we say so plainly.** The cone
branch a few lines above *does* finish with `cone.initialize_original()` and
`set_normals_and_coplanar(&mut cone)` — it models `AsOriginal` faithfully,
because it was written from the C++'s Mirror/Translate/AsOriginal chain. So a
cone reports a fresh original ID and a centered cylinder reports `extrude`'s.
That inconsistency is real, it is not defensible on its merits, and resolving it
belongs to a harmonization pass coordinated with manifold-sharp so both trees
move together. Until then, do not "fix" one half in isolation.

**Evidence:** the table above, measured with a scratch probe comparing the two
constructions coordinate by coordinate on raw f64 bits (4 signed-zero
differences, 0 value differences); the C++ sources cited above;
`constructors::tests::centered_cylinder_collider_matches_its_vertex_positions`
and `centered_cylinder_is_usable_in_a_boolean`, which pin the cache repair that
this entry's semantics were preserved *through*; and manifold-sharp's 34/34
oracle lane run against this tree's cdylib.
