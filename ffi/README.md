# manifold-ffi

C ABI shared library (`manifold_rs`) exposing the manifold-rust geometry kernel
to non-Rust callers — the immediate consumer is MatterCAD via P/Invoke, which
previously called the C++ `manifoldc` library.

The surface is intentionally small: import a triangle soup, run an n-ary
boolean, export the result. Op codes match `ManifoldOpType` from `manifoldc`
(0 = union, 1 = difference, 2 = intersection) so call sites port over
unchanged.

`manifold_rs.h` is the normative ABI description — ownership rules, sentinels
and array lifetimes are documented there. This file is the overview.

## Build

```sh
cargo build --release -p manifold-ffi
```

Output: `target/release/manifold_rs.dll` (plus `manifold_rs.dll.lib` on
Windows), `target/release/libmanifold_rs.so` on Linux,
`target/release/libmanifold_rs.dylib` on macOS.

The `parallel` feature of the core crate is enabled. It only parallelises
determinism-preserving sites, so results stay bit-identical to a sequential
build.

Tests:

```sh
cargo test --release -p manifold-ffi
```

## Functions

| Function | Purpose |
| --- | --- |
| `manifold_rs_version` | Static version string: FFI version plus core crate version. |
| `manifold_rs_from_mesh` | Build a manifold from interleaved vertex properties and triangle indices. |
| `manifold_rs_as_original` | Copy re-tagged as an original mesh (fresh mesh ID, no boolean history). |
| `manifold_rs_original_id` | Original mesh ID, or -1. |
| `manifold_rs_status` | Error code (see table below), or -1 for a NULL handle. |
| `manifold_rs_destroy` | Free a manifold handle. |
| `manifold_rs_batch_boolean` | N-ary union / difference / intersection through the CSG tree. |
| `manifold_rs_get_meshgl` | Export a manifold as a mesh handle. |
| `manifold_rs_meshgl_num_prop` | Properties per vertex (>= 3). |
| `manifold_rs_meshgl_vert_properties` | Borrowed float array, `num_prop` per vertex. |
| `manifold_rs_meshgl_tri_verts` | Borrowed index array, 3 per triangle. |
| `manifold_rs_meshgl_run_index` | Borrowed run start offsets plus end sentinel. |
| `manifold_rs_meshgl_run_original_id` | Borrowed source mesh ID per run. |
| `manifold_rs_meshgl_face_id` | Borrowed source face ID per triangle. |
| `manifold_rs_meshgl_destroy` | Free a mesh handle. |
| `manifold_rs_last_error` | Thread-local message for the most recent failure. |
| `manifold_rs_batch_boolean_ct` | As above, plus an optional cancellation token. |
| `manifold_rs_cancel_token_new` | Create a cancellation token. |
| `manifold_rs_cancel_token_cancel` | Request cancellation (callable from any thread). |
| `manifold_rs_cancel_token_is_cancelled` | Whether cancellation was requested. |
| `manifold_rs_cancel_token_destroy` | Free a token handle. |
| `manifold_rs_from_mesh64` | Build a manifold from `double` / `uint64_t` arrays. |
| `manifold_rs_get_meshgl64` | Export a manifold as a double-precision mesh handle. |
| `manifold_rs_meshgl64_num_prop` | Properties per vertex (>= 3). |
| `manifold_rs_meshgl64_vert_properties` | Borrowed `double` array, `num_prop` per vertex. |
| `manifold_rs_meshgl64_tri_verts` | Borrowed `uint64_t` index array, 3 per triangle. |
| `manifold_rs_meshgl64_run_index` | Borrowed `uint64_t` run start offsets plus end sentinel. |
| `manifold_rs_meshgl64_run_original_id` | Borrowed `uint32_t` source mesh ID per run. |
| `manifold_rs_meshgl64_face_id` | Borrowed `uint64_t` source face ID per triangle. |
| `manifold_rs_meshgl64_destroy` | Free a 64-bit mesh handle. |
| `manifold_rs_from_mesh_robust` | Non-manifold-tolerant import: closed orientable soup is retained for the robust engine. |
| `manifold_rs_from_mesh64_robust` | Double-precision counterpart of `manifold_rs_from_mesh_robust`. |
| `manifold_rs_set_boolean_engine` | Process-global engine default: 0 Exact, 1 Robust, 2 Auto. |
| `manifold_rs_get_boolean_engine` | The current engine default. |

## Robust (non-manifold) booleans

The default boolean engine is the exact C++-matching pipeline and requires
strictly manifold operands. `manifold_rs_from_mesh_robust` additionally
accepts *closed, orientable* but non-manifold geometry (shared edges or
vertices, disconnected shells, internal voids); such a handle reports status
0 and works with booleans once the engine default is `MANIFOLD_RS_ENGINE_AUTO`
(exact for manifold pairs, robust when a soup operand is involved) or
`MANIFOLD_RS_ENGINE_ROBUST` (always). Geometry that is not even closed
reports status 15, and pairing-dependent operations (`as_original`, …) on a
soup handle return empty results with status 2. The robust engine computes
with exact rational arithmetic; on manifold inputs it agrees with the exact
engine to near-f64 precision, but triangulation of results may differ.

## Cancellation

`manifold_rs_batch_boolean_ct` takes an optional `CancelTokenRs*`. Pass NULL to
get the uncancellable behaviour of `manifold_rs_batch_boolean` at zero cost —
the kernel does not read any flag on that path.

```c
CancelTokenRs* t = manifold_rs_cancel_token_new();
/* worker thread */ ManifoldRs* r = manifold_rs_batch_boolean_ct(parts, n, op, t);
/* UI thread    */ manifold_rs_cancel_token_cancel(t);
/* after the worker returns */
if (manifold_rs_status(r) == 14) { /* cancelled */ }
manifold_rs_cancel_token_destroy(t);
```

Three things to know:

- **Cancelling from another thread is the supported use**, not a hazard. The
  flag is atomic, and any number of threads may hold and cancel the same token.
  This is the one exception to the "don't use a handle concurrently" rule below.
- **A cancelled call returns a valid handle with status 14, not NULL.** NULL
  stays reserved for argument errors and caught panics, so a caller can tell
  "the user cancelled" from "something broke". The handle still needs
  destroying.
- **Destroying the token early is safe, but useless.** An operation clones the
  shared flag when it starts, so `manifold_rs_cancel_token_destroy` during a
  running call is memory-safe — no use-after-free, no leak. It just leaves
  nothing to cancel through, and the operation runs to completion. From C# that
  still means keeping the handle rooted for the lifetime of the operation
  (including across awaits) if you want a later cancel to land. The one hard
  rule is not to race destroy against `..._cancel` / `..._is_cancelled`, which
  read the handle itself.

Cancellation is cooperative: the kernel checks the flag at phase boundaries in
the CSG tree, the boolean stages and the triangulation, so a cancelled call
unwinds shortly after the request rather than instantly.

## Double-precision mesh path

`manifold_rs_from_mesh64` / `manifold_rs_get_meshgl64` and the
`manifold_rs_meshgl64_*` accessors mirror the single-precision family with the
element types of the core crate's `MeshGL64`: `double` vertex properties and
`uint64_t` indices, but `run_original_id` stays `uint32_t` (a mesh ID is not an
index). The kernel computes in double precision and this path is lossless end
to end: coordinates and tolerance go in and come out at full f64 precision,
with none of the f32 narrowing (or the FLT_EPSILON tolerance floor) the
single-precision export applies. Indices are the one exception — the kernel
indexes vertices with 32 bits, so `manifold_rs_from_mesh64` rejects indices
above `UINT32_MAX` rather than wrapping them. The `ManifoldRs*` it produces is
an ordinary handle and mixes freely with ones from `manifold_rs_from_mesh`.

## Memory ownership

- `ManifoldRs*`, `MeshGlRs*`, `MeshGl64Rs*` and `CancelTokenRs*` are opaque,
  caller-owned handles. Every non-NULL handle must be released with its
  matching destroy (`manifold_rs_destroy`, `manifold_rs_meshgl_destroy`,
  `manifold_rs_meshgl64_destroy`, `manifold_rs_cancel_token_destroy`). Every
  destroy ignores NULL.
- `manifold_rs_from_mesh` / `manifold_rs_from_mesh64` **copy** their input
  arrays; the caller may free them as soon as the call returns.
- `manifold_rs_batch_boolean` does **not** consume its operands — the caller
  still owns and must destroy each of them.
- The `meshgl` array accessors return **borrowed** pointers into the mesh
  handle. They are invalidated by `manifold_rs_meshgl_destroy`, so copy the
  data out first. They are never freed by the caller.
- `manifold_rs_version` returns a static string; never free it.
- For an empty array the accessors return a non-NULL pointer that must not be
  dereferenced — always check the length written to `out_len`.

## Failure handling

No function unwinds or aborts. The geometry kernel can panic on degenerate
input, so every exported function runs its body under `catch_unwind` and
converts a panic into its failure sentinel (`NULL`, `-1` or `0`), recording the
panic message for `manifold_rs_last_error`.

Validation failures are reported two different ways on purpose:

- Arguments that cannot describe a mesh at all (`num_prop < 3`, indivisible
  lengths, NULL arrays, unknown op code) return NULL and set the last-error
  message.
- A well-formed request that produces an invalid solid returns a **handle with
  a non-zero status**, matching the C++ behaviour. Check `manifold_rs_status`
  after importing. Exporting such a manifold is safe and yields an empty mesh.

Two consequences worth knowing before wiring this up:

- A boolean operand carrying a non-zero status is absorbed as empty geometry
  and the result may still report status 0, so a failed import surfaces as a
  part silently missing from the output rather than as an error. Check the
  status of every operand before combining them.
- `manifold_rs_as_original` on an empty or error manifold returns a plain copy
  and leaves `original_id` at -1, so don't assume a fresh ID was assigned.

### Note for .NET callers

Read `manifold_rs_last_error` on the **same thread** that made the failing
call, inside the same worker delegate. The slot is thread-local: if the failing
call is awaited and the continuation resumes on a different thread-pool thread,
the message is gone and you read an empty string. Keep the call and its error
read together in whatever you hand to `Task.Run`.

## Error codes (`manifold_rs_status`)

| Code | Meaning |
| --- | --- |
| -1 | NULL handle (not a real status) |
| 0 | No error |
| 1 | Non-finite vertex |
| 2 | Not manifold |
| 3 | Vertex out of bounds |
| 4 | Properties wrong length |
| 5 | Missing position properties |
| 6 | Merge vectors different lengths |
| 7 | Merge index out of bounds |
| 8 | Transform wrong length |
| 9 | Run index wrong length |
| 10 | Face ID wrong length |
| 11 | Invalid construction |
| 12 | Result too large |
| 13 | Invalid tangents |
| 14 | Cancelled (interrupted through a `CancelTokenRs`) |

## Threading

Distinct handles may be used concurrently from different threads. A single
handle must not be used concurrently with a call that frees it. The last-error
slot is thread-local, so read it on the thread that saw the failure.

`CancelTokenRs` is the deliberate exception: cancelling (or querying) a token
while another thread is inside an operation holding it is exactly what the type
is for. Only `manifold_rs_cancel_token_destroy` must be ordered after every
call that used the token.
