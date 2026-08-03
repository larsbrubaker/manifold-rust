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

The package does not yet carry natives; that arrives with the CI job that builds
every runtime identifier.

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
| `Manifold.BatchBoolean(IReadOnlyList<Manifold>, ManifoldOpType)` | `Subtract` is the first operand minus the union of the rest. Operands are not consumed. |
| `Manifold.GetMeshGL()` | Eagerly copies every array into managed memory and frees the native mesh before returning. |
| `Manifold.NativeVersion` | Version string of the loaded shared library. Reading it forces the load, so it is a cheap early check that the native library is findable. |
| `ManifoldNative.LibraryPath` | Host override for where the shared library lives. See above. |
| `MeshGL` | Plain data: `NumProp`, `VertProperties`, `TriVerts`, `RunIndex`, `RunOriginalId`, `FaceId`, plus `TriangleCount` / `VertexCount`. No unmanaged lifetime, nothing to dispose. |

`Manifold` is `IDisposable` and wraps a `SafeHandle`, so a forgotten `Dispose`
is collected rather than leaked, and `Dispose` may be called more than once.

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
