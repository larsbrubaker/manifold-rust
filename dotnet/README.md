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
manifold-ffi 0.2.0 (manifold-rust 0.9.3)
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

**They are not more precision end to end.** The kernel currently computes in
single precision internally, so coordinates handed to `FromMesh64` round-trip
through `float` and come back out of `GetMeshGL64` having lost about 1e-7 of
relative precision. This is input fidelity and future-proofing, not a
double-precision pipeline. At any physical tolerance the difference is
irrelevant; if you are relying on more than seven significant digits, this API
does not give it to you today.

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
   `../ffi/Cargo.toml` if the ABI changed — the two are compared at run time by
   the version handshake above, so they move together. Update
   `NativeVersionCheck.ExpectedPrefix` to match.
2. Commit, then tag: `git tag nuget-v0.2.0 && git push origin nuget-v0.2.0`.
3. CI, in order: checks that the tag version matches the csproj version; runs
   the full Rust and .NET test suites **at the tagged commit**; builds
   `manifold_rs` for every runtime identifier; packs; pushes to nuget.org with
   the `NUGET_API_KEY` repository secret.

The version check is its own first job, so a mismatched tag costs a few seconds
rather than four cargo builds — or a yanked package. The tests are re-run at the
tag rather than trusted from the branch, so a release from code that does not
pass is not possible.

Adding a platform means editing two files: a matrix leg in the workflow *and* a
`ManifoldRustNative*` property, `None` item and `EnsureNativeIsPacked` check in
the csproj. Forgetting the second half would produce a package quietly missing
that platform, so the pack job fails on any staged runtime identifier the csproj
does not pack. A release also packs with `ManifoldRustRequireAllNatives=true`,
which turns the "missing platform" warnings into errors: every platform ships,
or none does.

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
