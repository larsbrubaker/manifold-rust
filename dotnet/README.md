# ManifoldRust (.NET binding)

Managed bindings for `manifold_rs`, the C ABI in `../ffi`. The surface mirrors
the ABI one-for-one: import a triangle soup, run an n-ary boolean, export the
result.

`../ffi/manifold_rs.h` is the normative description of the ABI. This binding
only adds .NET lifetime and error handling on top of it; when the two disagree,
the header is right.

## Build and test

The native library is produced by cargo, not by MSBuild, so build it first:

```sh
cargo build --release -p manifold-ffi      # from the repository root
```

Then, from this `dotnet` directory:

```sh
dotnet build ManifoldRust/ManifoldRust.csproj
dotnet test ManifoldRust.Tests/ManifoldRust.Tests.csproj
```

Run the .NET commands from this directory: `global.json` here opts `dotnet test`
into Microsoft.Testing.Platform, which TUnit requires on the .NET 10 SDK.

The test project copies `../target/release/manifold_rs.dll` (or the `.so` /
`.dylib`) next to its binaries and fails the build with a pointer at the cargo
command if that file is missing.

## Finding the native library

The assembly installs a `DllImportResolver` at module load, which probes in this
order:

1. `ManifoldNative.LibraryPath` — a **full path** set by the host in code.
2. The `MANIFOLD_RS_NATIVE` environment variable, also a full path to the
   `.dll` / `.so` / `.dylib`. Handy for running against a Rust build without
   copying anything.
3. Default runtime resolution — this is what a future NuGet package's
   `runtimes/<rid>/native/` folder will hit.
4. `AppContext.BaseDirectory`, using the platform's file name. A fallback for
   hosts whose native probing paths do not include it, such as a plugin loaded
   from its own folder.

The explicit overrides deliberately come first. Default resolution already
searches the application base directory, so an override checked after it could
never win in the common case of a library sitting next to the app — which is
exactly when you want to redirect it.

```csharp
// Before the first call into the library, and only useful for a host that
// cannot use the environment: the resolver slot is already taken by this
// assembly, so SetDllImportResolver is not available to callers.
ManifoldNative.LibraryPath = Path.Combine(myNativeFolder, "manifold_rs.dll");
```

Once the library is loaded the runtime caches the handle and never consults the
resolver again, so setting `LibraryPath` afterwards has no effect. A path that
fails to load is not an error; probing falls through to the next candidate.

The package carries natives for `win-x64`, `linux-x64`, `osx-arm64` and
`osx-x64` under `runtimes/<rid>/native/`, which step 3 above picks up with no
configuration. Any other runtime identifier needs one of the overrides.

### browser-wasm

A browser has no dynamic loading, so none of the four steps above apply — the
resolver stands down under `OperatingSystem.IsBrowser()`. The wasm build is a
static archive, `runtimes/browser-wasm/native/manifold_rs.a`, that has to be on
mono-wasm's emcc link line instead, which is what `build/ManifoldRust.props`
does: it adds the archive as a `NativeFileReference` and widens
`WasmOptConfigurationFlags` to cover the features rustc emits and the SDK's own
post-link `wasm-opt` run does not enable. That happens for any project the SDK
is building for `browser-wasm`, and does nothing at all for anyone else.

It is packed to both `build/` and `buildTransitive/`, which matters more than it
looks: NuGet imports `build/` only into projects that reference the package
*directly*, and the consumer that needs the link is usually the application at
the end of a chain — app → some geometry library → `ManifoldRust`. With
`build/` alone that application imports nothing, links nothing, and fails at its
first CSG call with no build-time hint that anything was missing.

Three consequences worth knowing:

- Linking it turns on `WasmBuildNative`, so the wasm module is relinked with
  emcc and the build takes correspondingly longer. Measured on the 10.0.400 SDK
  that happens on a plain `dotnet build` as well as on `dotnet publish` —
  publish only defers the build-time pass to its own nested one.
  `-p:ManifoldRustLinkBrowserWasm=false` opts out; the `DllImport`s are then
  unbacked and the first CSG call throws, which is the right failure for a build
  that asked not to link it.
- If the package was packed without the browser-wasm leg — a local `dotnet pack`
  on a machine that never built the archive does exactly that, with only a
  warning — the build fails with a message naming that opt-out, rather than
  somewhere inside the native build.
- The file name `manifold_rs.a` is load bearing. mono-wasm derives its P/Invoke
  module names from the file names it links, so `DllImport("manifold_rs")`
  resolves only against an archive named exactly that — cargo's own
  `libmanifold_rs.a` would not do.

## Version handshake

Because the native library can be redirected — that is the whole point of
`LibraryPath` and `MANIFOLD_RS_NATIVE` — nothing stops a stale one from being
loaded by a newer binding. A missing export raises
`EntryPointNotFoundException`, which is at least legible; the dangerous case is
an export that still exists but whose contract changed.

So the first call into the library compares `Manifold.NativeVersion` against the
ABI version this binding was written for and throws `ManifoldException` on a
mismatch. The check runs once per process and costs one call to
`manifold_rs_version` after that.

```
manifold-ffi 0.3.2 (manifold-rust 0.14.0)
             ^^^^ compared, to major.minor
```

Only major and minor are compared: a patch release of the FFI crate cannot
change the ABI, so pinning to it would reject a perfectly good native.
`Manifold.NativeVersion` itself never throws, so it is still usable to report
what actually loaded.

## API sketch

```csharp
using ManifoldRust;

// AsOriginal returns a copy, so the import it was called on is a separate
// manifold that has to be disposed too - hence the inner using.
static Manifold Import(ReadOnlySpan<float> verts, ReadOnlySpan<uint> tris)
{
    using Manifold imported = Manifold.FromMesh(verts, tris);
    return imported.AsOriginal();
}

using Manifold body = Import(vertProperties, triVerts);
using Manifold hole = Import(holeVerts, holeTris);

if (body.Status != ManifoldStatus.NoError || hole.Status != ManifoldStatus.NoError)
{
    // See the caveat below - do not just hand these to BatchBoolean.
}

using Manifold result = Manifold.BatchBoolean(
    new[] { body, hole },
    ManifoldOpType.Subtract);

MeshGL mesh = result.GetMeshGL();
// mesh.VertProperties, mesh.TriVerts, mesh.RunIndex, mesh.RunOriginalId, mesh.FaceId
```

| Member | Notes |
| --- | --- |
| `Manifold.FromMesh(ReadOnlySpan<float>, ReadOnlySpan<uint>, uint numProp = 3)` | Copies both arrays natively. Throws `ManifoldException` only when the arguments cannot describe a mesh; a mesh that merely fails validation comes back as a manifold with a non-zero `Status`. |
| `Manifold.AsOriginal()` | Returns a **new** manifold - dispose the one it was called on as well. Fresh mesh ID, boolean history dropped. On an empty or error manifold it is a plain copy and `OriginalId` stays -1. |
| `Manifold.OriginalId` | -1 until `AsOriginal` has been called. |
| `Manifold.Status` | `ManifoldStatus.NoError` or the reason the solid is empty. |
| `Manifold.FromMesh64(ReadOnlySpan<double>, ReadOnlySpan<ulong>, uint numProp = 3)` | Double-precision import. Same rules as `FromMesh`, plus: an index above `uint.MaxValue` is rejected rather than wrapped. See below. |
| `Manifold.FromMesh64(ReadOnlySpan<double>, ReadOnlySpan<uint>, uint numProp = 3)` | The same, widening 32-bit indices into a scratch buffer for you. |
| `Manifold.BatchBoolean(IReadOnlyList<Manifold>, ManifoldOpType)` | `Subtract` is the first operand minus the union of the rest. Operands are not consumed. |
| `Manifold.BatchBoolean(IReadOnlyList<Manifold>, ManifoldOpType, CancellationToken)` | Throws `OperationCanceledException` if the token is signalled before the boolean finishes. See below. |
| `Manifold.BatchBoolean(IReadOnlyList<Manifold>, ManifoldOpType, CancelToken)` | Lower level: reports cancellation as a result with `Status == Cancelled` instead of throwing. |
| `Manifold.Boolean(Manifold, Manifold, ManifoldOpType, BooleanEngine, ...)` | Binary boolean with the engine pinned, plus optional `WindingRule`, progress sink and cancellation. `Union` / `Subtract` / `Intersect` are the named spellings. Use `BatchBoolean` for three or more operands: only it runs the CSG tree. |
| `Manifold.RepairOrientation()` | Returns a **new** manifold with inside-out shells rewound so every body reads as solid under the `Positive` winding rule. A mesh that needs no repair comes back as a plain copy, so it is safe to call unconditionally. |
| `Manifold.RebuildSolid(WindingRule, ...)` | Returns a **new** manifold rebuilt from the winding numbers: self-intersections, doubled sheets and interior walls are removed, not just rewound. Costs what a boolean costs, so there are `CancelToken` and `CancellationToken` overloads. The result is re-triangulated — counts change and properties are re-interpolated. |
| `Manifold.HasSelfIntersections` | Whether two of the mesh's own triangles genuinely cross, overlap or coincide — shared edges and vertices do not count. Cached natively, so repeated reads are free. |
| `Manifold.GetMeshGL()` | Eagerly copies every array into managed memory and frees the native mesh before returning. |
| `Manifold.GetMeshGL64()` | The same, with `double[]` coordinates. |
| `Manifold.NativeVersion` | Version string of the loaded shared library. Reading it forces the load, so it is a cheap early check that the native library is findable. |
| `ManifoldNative.LibraryPath` | Host override for where the shared library lives. See above. |
| `MeshGL` | Plain data: `NumProp`, `VertProperties`, `TriVerts`, `RunIndex`, `RunOriginalId`, `FaceId`, plus `TriangleCount` / `VertexCount`. No unmanaged lifetime, nothing to dispose. |
| `MeshGL64` | The same fields, with `double[] VertProperties`. The index fields stay `uint[]`. |
| `CancelToken` | A cancellation flag that can be polled and reused across operations. `IDisposable`. |

`Manifold` is `IDisposable` and wraps a `SafeHandle`, so a forgotten `Dispose`
is collected rather than leaked, and `Dispose` may be called more than once. The
same is true of `CancelToken`.

## The f64 API

`FromMesh64` and `GetMeshGL64` exist for callers whose geometry is already
`double` — most mesh formats, and anything that came out of a CAD kernel — so
they do not have to narrow it by hand at the boundary.

**Coordinates are double precision end to end.** The kernel computes in
double precision, and this path feeds it and reads it back without narrowing
anywhere: a coordinate handed to `FromMesh64` that a boolean leaves untouched
comes back out of `GetMeshGL64` bit-identical, and booleans themselves run at
full f64 precision. (The `float`-based `FromMesh`/`GetMeshGL` remain the lossy
pair: their coordinates narrow at the boundary and the exported tolerance is
floored at `float` epsilon times the bounding-box scale, which the f64 export
does not do.)

Indices are the case where the narrowing is not merely lossy but *wrong*: the
kernel indexes vertices with 32 bits, and `as u32` wraps. An index of 2^32 would
silently become 0 and produce plausible-looking, incorrect geometry reporting
`NoError`. So `FromMesh64` rejects any index above `uint.MaxValue` outright, and
`MeshGL64` exposes the index arrays as `uint[]` narrowed with the same check
rather than as `ulong[]` that would only push the problem to the caller. In
practice this caps a mesh at about four billion vertices.

```csharp
// Coordinates stay double; indices are the common int-based kind, so the
// widening overload does the conversion.
using Manifold body = Manifold.FromMesh64(doubleVerts, intBasedTris);
using Manifold hole = Manifold.FromMesh64(holeVerts, holeTris);

using Manifold result = Manifold.BatchBoolean(
    new[] { body, hole },
    ManifoldOpType.Subtract);

MeshGL64 mesh = result.GetMeshGL64();
// mesh.VertProperties is double[]; mesh.TriVerts is uint[].
```

Manifolds from either import path are interchangeable and can be mixed in the
same boolean.

## Cancellation

A boolean over large geometry can take seconds, and the ABI supports
interrupting one from another thread. The `CancellationToken` overload is the
one most callers want:

```csharp
using CancellationTokenSource source = new();
cancelButton.Click += (_, _) => source.Cancel();

try
{
    using Manifold result = await Task.Run(
        () => Manifold.BatchBoolean(parts, ManifoldOpType.Add, source.Token),
        source.Token);

    Display(result.GetMeshGL());
}
catch (OperationCanceledException)
{
    // The partial result was disposed before this was thrown.
}
```

Cancellation is **cooperative** and checked at phase boundaries inside the
kernel, so the call returns shortly after `Cancel()` rather than instantly. It
is also **sticky**: an already-cancelled token makes the call return
immediately, and there is no un-cancel.

**Completion wins.** Cancelling is a request, not a guarantee. If the kernel
finishes before it notices the flag — the normal outcome when the cancel lands
near the end of the work — you get the finished result and no exception, even
though the token is by then signalled. Code that must treat a late cancel as a
cancellation has to check the token itself after the call returns.

For a caller that wants to poll the flag, or reuse one flag across several
operations, `CancelToken` maps straight onto the ABI:

```csharp
using CancelToken token = new();
using Manifold result = Manifold.BatchBoolean(parts, ManifoldOpType.Add, token);

if (result.Status == ManifoldStatus.Cancelled)
{
    // This overload reports cancellation rather than throwing. The result is a
    // valid, empty manifold and you still own it.
}
```

`Cancel()` and `IsCancelled` may be called from any thread at any time,
including while an operation holding the token is running — that is what the
type is for. Disposing a token while such an operation is still running is
memory safe (the operation keeps its own share of the flag), but it silently
loses the ability to cancel that operation. The rule that actually matters:
**do not race `Dispose` against `Cancel`/`IsCancelled` on the same token.** Keep
the token alive until every call that might cancel through it has returned.

## Progress reporting and winding rules

The binary overloads take the engine explicitly and accept an
`IProgress<(string Phase, double? Fraction)>`. A `null` `Fraction` is the ABI's
*indeterminate*: the phase is running but has no meaningful ratio, so show a
spinner rather than a bar.

```csharp
IProgress<(string Phase, double? Fraction)> report = new Progress<(string, double?)>(
    p => status.Text = p.Item2 is double f ? $"{p.Item1} {f:P0}" : p.Item1);

using Manifold result = Manifold.Subtract(
    body, hole, BooleanEngine.Auto, WindingRule.Positive, report, cancellationToken);
```

**The callback may run on a native worker thread**, not the one that made the
call — the pipeline is parallel. It is never re-entered concurrently, but a sink
that touches UI state has to marshal itself; `System.Progress<T>` already posts
to the synchronization context that created it, which is why the example uses
it. Reporting never changes the result: the same inputs produce the same
triangles with or without a sink attached, and a `null` sink is exactly the
un-instrumented pipeline. An exception thrown by the sink cannot unwind through
the native frame, so the first one is captured and rethrown from the call.

`WindingRule` decides which winding numbers count as solid, and only the robust
engine reads it (`Exact` ignores it; `Auto` resolves to `Robust` whenever the
rule is `Nonzero`):

- `Positive` — `{w >= 1}`, the default everywhere and the rule clean,
  consistently wound data wants. An inside-out shell encloses `w = -1` and so
  contributes nothing.
- `Nonzero` — `{w != 0}`, which keeps inside-out shells as material. For scans
  and third-party geometry wound inconsistently that cannot be repaired first.

`RepairOrientation()` is the alternative to `Nonzero`: it rewinds the shells once
so every later operation can keep using `Positive`, instead of changing what
"solid" means for the whole call. `RebuildSolid(rule)` is the heavy option for
when rewinding is not enough: it re-derives the surface from the winding numbers
under `rule`, so self-intersections, doubled sheets and interior walls go away
too — at the cost of re-triangulating, and of only accepting input that is
already closed and orientable (an open mesh imports as `NotClosed` and never gets
that far). `HasSelfIntersections` answers the other
question worth asking before a boolean — geometry that overlaps itself breaks the
exact engine's winding integral, which is why `BooleanEngine.Auto` routes it to
`Robust`.

## Caveat: check `Status` on every input before `BatchBoolean`

An operand carrying a non-zero `Status` is absorbed as **empty geometry** and
the result can still report `NoError`. A part that failed to import therefore
shows up as geometry silently missing from the output, not as an exception.
Check `Status` on every manifold you build before combining them, and decide
there what a failed import should mean.

Exporting an error-status manifold is safe: `GetMeshGL` returns empty arrays.
That is the path where the C++ binding crashes.

## Caveat: read errors on the thread that saw them

The native last-error slot is thread local and is what
`ManifoldException.NativeError` carries. The binding reads it immediately after
observing a failure, which is only correct while the call and the read stay on
one thread.

For a caller, this means keeping a whole operation inside a single
`Task.Run` delegate rather than awaiting in the middle of one: if a failing call
is awaited and the continuation resumes on a different thread-pool thread, any
error message read after that point is gone.

## Releasing

The package is published by `.github/workflows/release-nuget.yml`, which runs on
a pushed tag of the form `nuget-v<version>`:

1. Bump `<Version>` in `ManifoldRust/ManifoldRust.csproj`, and `version` in
   `../ffi/Cargo.toml` if the ABI changed — update `NativeVersionCheck.ExpectedPrefix`
   to match when that is a major.minor change, since the handshake above compares
   that prefix. The package version may also move on its own when the managed
   surface grows over an unchanged ABI (0.3.1 → 0.4.0 was exactly that), so the
   two numbers are not required to be equal — only the handshake prefix is
   binding.
2. Commit, then tag: `git tag nuget-v0.2.0 && git push origin nuget-v0.2.0`.
3. CI, in order: checks that the tag version matches the csproj version; runs
   the full Rust and .NET test suites **at the tagged commit**; builds
   `manifold_rs` for every runtime identifier; packs; pushes to nuget.org with
   the `NUGET_API_KEY` repository secret.

The version check is its own first job, so a mismatched tag costs a few seconds
rather than five cargo builds — or a yanked package. The tests are re-run at the
tag rather than trusted from the branch, so a release from code that does not
pass is not possible.

Adding a platform means editing two files: a matrix leg in the workflow *and* a
`ManifoldRustNative*` property, `None` item and `EnsureNativeIsPacked` check in
the csproj. Forgetting the second half would produce a package quietly missing
that platform, so the pack job fails on any staged runtime identifier the csproj
does not pack. A release also packs with `ManifoldRustRequireAllNatives=true`,
which turns the "missing platform" warnings into errors: every platform ships,
or none does. Finally, the pack job unzips the package it just built and fails
if any of those natives — or `build/ManifoldRust.props` — is not in it.

`browser-wasm` is built by its own job rather than a matrix leg: it produces a
static archive from a different cargo invocation, needs an Emscripten on `PATH`,
and gets `.github/scripts/rewrite-wasm-target-features.py` run over it
afterwards. That script drops the two feature names rustc's LLVM records that
the Emscripten in the .NET wasm-tools workload cannot parse; without it a
consumer's build dies in `wasm-opt` with `Unknown option
'--enable-bulk-memory-opt'` *after* a link that succeeded. The script's header
explains the whole failure, and it can be run by hand over any local wasm
archive — with `--check` to assert an archive is fit to ship without modifying
it. That check is an allowlist: it fails on any feature name this pipeline has
not shipped before, not just the two known-bad ones, because
`rust-toolchain@stable` is unpinned and the next rustc is free to invent a third.

Two follow-ups deliberately left for a later release, recorded here so they are
not rediscovered from scratch:

- **`runtimes/browser-wasm/nativeassets/` instead of `native/`.** SkiaSharp uses
  `nativeassets/` for exactly this kind of link-only wasm archive, so NuGet does
  not also treat it as a copy-to-output runtime asset. Ours is copied to the
  output directory for a RID-targeted build, which costs a copy but does not
  reach the publish payload — not worth the churn mid-release, worth doing when
  the layout is next touched.
- **Dropping the emsdk from this job.** A staticlib is archived by rustc's own
  LLVM and nothing links, so `llvm-tools` plus an `llvm-objcopy` may be all this
  job needs, which would cut the emsdk clone and install out of the release path
  entirely. Test it *before* the next release rather than during one: the failure
  mode if rustc does want a linker driver is a red release run.

To pack locally, refresh the staged native first — the csproj packs from
`ManifoldRust/native-staging/<rid>/`, not from `target/release`, because a cargo
build can relink that file mid-pack:

```sh
cargo build --release -p manifold-ffi
cp target/release/manifold_rs.dll dotnet/ManifoldRust/native-staging/win-x64/
dotnet pack dotnet/ManifoldRust/ManifoldRust.csproj -c Release
```

A local pack warns about the platforms it has no native for and produces a
win-x64-only package; only a missing win-x64 native is a hard error.
