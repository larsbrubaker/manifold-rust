// Copyright 2026 Lars Brubaker
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// One-to-one declarations of the C ABI in ffi/manifold_rs.h, plus the resolver
// that finds the shared library. Nothing here interprets the ABI; the ownership
// and sentinel rules the header documents are enforced one layer up, in
// Manifold.cs and MeshGL.cs.
//
// Handles are passed as raw IntPtr rather than SafeHandle because the
// LibraryImport source generator does not marshal SafeHandle - callers are
// responsible for keeping the owning SafeHandle alive across the call
// (DangerousAddRef / DangerousRelease).

using System;
using System.IO;
using System.Reflection;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Text;

namespace ManifoldRust
{
	internal static unsafe partial class NativeMethods
	{
		/// <summary>
		/// Base name of the shared library: manifold_rs.dll, libmanifold_rs.so or
		/// libmanifold_rs.dylib depending on the platform.
		/// </summary>
		internal const string LibraryName = "manifold_rs";

		/// <summary>
		/// Environment variable holding a full path to the native library. See
		/// <see cref="ManifoldNative.LibraryPath"/> for why the overrides are probed
		/// before the runtime's own resolution rather than after it.
		/// </summary>
		internal const string NativePathVariable = "MANIFOLD_RS_NATIVE";

		// CA2255 warns that module initializers are usually the application's job.
		// Here it is the point: the resolver has to be installed before any P/Invoke
		// in this assembly runs, whatever entry point the consumer reaches first.
#pragma warning disable CA2255
		[ModuleInitializer]
		internal static void Initialize()
		{
			NativeLibrary.SetDllImportResolver(typeof(NativeMethods).Assembly, Resolve);
		}
#pragma warning restore CA2255

		/// <summary>
		/// Probes for the native library.
		/// </summary>
		/// <remarks>
		/// The explicit overrides come first, because they are only ever set by
		/// someone who has decided which build to load - and the runtime's default
		/// resolution already searches the application base directory, so an override
		/// checked after it could never win whenever a library happened to sit next
		/// to the app. Order is:
		/// <list type="number">
		/// <item><description><see cref="ManifoldNative.LibraryPath"/>.</description></item>
		/// <item><description>The <c>MANIFOLD_RS_NATIVE</c> environment variable.</description></item>
		/// <item><description>Default runtime resolution, which is what a package's runtimes/ folder relies on.</description></item>
		/// <item><description>The application base directory, for hosts whose native probing paths do not include it (a plugin loaded from its own folder, for example).</description></item>
		/// </list>
		/// </remarks>
		private static IntPtr Resolve(string libraryName, Assembly assembly, DllImportSearchPath? searchPath)
		{
			if (libraryName != LibraryName)
			{
				return IntPtr.Zero;
			}

			IntPtr handle;

			string? configuredPath = ManifoldNative.LibraryPath;
			if (!string.IsNullOrEmpty(configuredPath) && NativeLibrary.TryLoad(configuredPath, out handle))
			{
				return handle;
			}

			string? environmentPath = Environment.GetEnvironmentVariable(NativePathVariable);
			if (!string.IsNullOrEmpty(environmentPath) && NativeLibrary.TryLoad(environmentPath, out handle))
			{
				return handle;
			}

			// NativeLibrary.TryLoad does not re-enter this resolver, so this is the
			// runtime's own probing logic and not an infinite loop.
			if (NativeLibrary.TryLoad(libraryName, assembly, searchPath, out handle))
			{
				return handle;
			}

			foreach (string fileName in PlatformFileNames())
			{
				string candidate = Path.Combine(AppContext.BaseDirectory, fileName);
				if (File.Exists(candidate) && NativeLibrary.TryLoad(candidate, out handle))
				{
					return handle;
				}
			}

			// Zero lets the runtime raise its own DllNotFoundException, which names
			// the library and the platform's search path.
			return IntPtr.Zero;
		}

		private static string[] PlatformFileNames()
		{
			if (RuntimeInformation.IsOSPlatform(OSPlatform.Windows))
			{
				return new[] { LibraryName + ".dll" };
			}

			if (RuntimeInformation.IsOSPlatform(OSPlatform.OSX))
			{
				return new[] { "lib" + LibraryName + ".dylib" };
			}

			return new[] { "lib" + LibraryName + ".so" };
		}

		/// <summary>
		/// Reads the calling thread's last native failure message.
		/// </summary>
		/// <remarks>
		/// The slot is thread local and is not cleared by later successful calls, so
		/// this must be called on the same thread that saw the failure, before that
		/// thread makes any other call that could fail. Sized with a length-only
		/// query first so a long panic message is never truncated.
		/// </remarks>
		internal static string GetLastError()
		{
			nuint length = manifold_rs_last_error(null, 0);
			if (length == 0)
			{
				return string.Empty;
			}

			int size = checked((int)length);
			byte[] buffer = new byte[size];

			nuint written;
			fixed (byte* target = buffer)
			{
				written = manifold_rs_last_error(target, (nuint)size);
			}

			// written is the full message length, which can exceed the buffer if the
			// message somehow grew between the two calls.
			int copied = Math.Min(size, checked((int)written));
			return Encoding.UTF8.GetString(buffer, 0, copied);
		}

		[LibraryImport(LibraryName)]
		[UnmanagedCallConv(CallConvs = new[] { typeof(CallConvCdecl) })]
		internal static partial IntPtr manifold_rs_version();

		[LibraryImport(LibraryName)]
		[UnmanagedCallConv(CallConvs = new[] { typeof(CallConvCdecl) })]
		internal static partial IntPtr manifold_rs_from_mesh(
			float* vertProperties,
			nuint vertPropertiesLen,
			uint* triVerts,
			nuint triVertsLen,
			uint numProp);

		[LibraryImport(LibraryName)]
		[UnmanagedCallConv(CallConvs = new[] { typeof(CallConvCdecl) })]
		internal static partial IntPtr manifold_rs_as_original(IntPtr m);

		[LibraryImport(LibraryName)]
		[UnmanagedCallConv(CallConvs = new[] { typeof(CallConvCdecl) })]
		internal static partial int manifold_rs_original_id(IntPtr m);

		[LibraryImport(LibraryName)]
		[UnmanagedCallConv(CallConvs = new[] { typeof(CallConvCdecl) })]
		internal static partial int manifold_rs_status(IntPtr m);

		[LibraryImport(LibraryName)]
		[UnmanagedCallConv(CallConvs = new[] { typeof(CallConvCdecl) })]
		internal static partial void manifold_rs_destroy(IntPtr m);

		[LibraryImport(LibraryName)]
		[UnmanagedCallConv(CallConvs = new[] { typeof(CallConvCdecl) })]
		internal static partial IntPtr manifold_rs_batch_boolean(IntPtr* manifolds, nuint count, int op);

		// token may be IntPtr.Zero, which the ABI reads as "uncancellable" and makes
		// this call identical to manifold_rs_batch_boolean.
		[LibraryImport(LibraryName)]
		[UnmanagedCallConv(CallConvs = new[] { typeof(CallConvCdecl) })]
		internal static partial IntPtr manifold_rs_batch_boolean_ct(IntPtr* manifolds, nuint count, int op, IntPtr token);

		[LibraryImport(LibraryName)]
		[UnmanagedCallConv(CallConvs = new[] { typeof(CallConvCdecl) })]
		internal static partial IntPtr manifold_rs_cancel_token_new();

		[LibraryImport(LibraryName)]
		[UnmanagedCallConv(CallConvs = new[] { typeof(CallConvCdecl) })]
		internal static partial void manifold_rs_cancel_token_cancel(IntPtr t);

		[LibraryImport(LibraryName)]
		[UnmanagedCallConv(CallConvs = new[] { typeof(CallConvCdecl) })]
		internal static partial int manifold_rs_cancel_token_is_cancelled(IntPtr t);

		[LibraryImport(LibraryName)]
		[UnmanagedCallConv(CallConvs = new[] { typeof(CallConvCdecl) })]
		internal static partial void manifold_rs_cancel_token_destroy(IntPtr t);

		[LibraryImport(LibraryName)]
		[UnmanagedCallConv(CallConvs = new[] { typeof(CallConvCdecl) })]
		internal static partial IntPtr manifold_rs_get_meshgl(IntPtr m);

		[LibraryImport(LibraryName)]
		[UnmanagedCallConv(CallConvs = new[] { typeof(CallConvCdecl) })]
		internal static partial uint manifold_rs_meshgl_num_prop(IntPtr g);

		[LibraryImport(LibraryName)]
		[UnmanagedCallConv(CallConvs = new[] { typeof(CallConvCdecl) })]
		internal static partial float* manifold_rs_meshgl_vert_properties(IntPtr g, out nuint outLen);

		[LibraryImport(LibraryName)]
		[UnmanagedCallConv(CallConvs = new[] { typeof(CallConvCdecl) })]
		internal static partial uint* manifold_rs_meshgl_tri_verts(IntPtr g, out nuint outLen);

		[LibraryImport(LibraryName)]
		[UnmanagedCallConv(CallConvs = new[] { typeof(CallConvCdecl) })]
		internal static partial uint* manifold_rs_meshgl_run_index(IntPtr g, out nuint outLen);

		[LibraryImport(LibraryName)]
		[UnmanagedCallConv(CallConvs = new[] { typeof(CallConvCdecl) })]
		internal static partial uint* manifold_rs_meshgl_run_original_id(IntPtr g, out nuint outLen);

		[LibraryImport(LibraryName)]
		[UnmanagedCallConv(CallConvs = new[] { typeof(CallConvCdecl) })]
		internal static partial uint* manifold_rs_meshgl_face_id(IntPtr g, out nuint outLen);

		[LibraryImport(LibraryName)]
		[UnmanagedCallConv(CallConvs = new[] { typeof(CallConvCdecl) })]
		internal static partial void manifold_rs_meshgl_destroy(IntPtr g);

		// The 64-bit mesh family spells every element type out rather than reusing
		// the declarations above: the widths differ per field (run_original_id stays
		// uint32_t because a mesh ID is not an index), and getting one wrong would
		// misread the array instead of failing to compile.
		[LibraryImport(LibraryName)]
		[UnmanagedCallConv(CallConvs = new[] { typeof(CallConvCdecl) })]
		internal static partial IntPtr manifold_rs_from_mesh64(
			double* vertProperties,
			nuint vertPropertiesLen,
			ulong* triVerts,
			nuint triVertsLen,
			uint numProp);

		[LibraryImport(LibraryName)]
		[UnmanagedCallConv(CallConvs = new[] { typeof(CallConvCdecl) })]
		internal static partial IntPtr manifold_rs_get_meshgl64(IntPtr m);

		[LibraryImport(LibraryName)]
		[UnmanagedCallConv(CallConvs = new[] { typeof(CallConvCdecl) })]
		internal static partial uint manifold_rs_meshgl64_num_prop(IntPtr g);

		[LibraryImport(LibraryName)]
		[UnmanagedCallConv(CallConvs = new[] { typeof(CallConvCdecl) })]
		internal static partial double* manifold_rs_meshgl64_vert_properties(IntPtr g, out nuint outLen);

		[LibraryImport(LibraryName)]
		[UnmanagedCallConv(CallConvs = new[] { typeof(CallConvCdecl) })]
		internal static partial ulong* manifold_rs_meshgl64_tri_verts(IntPtr g, out nuint outLen);

		[LibraryImport(LibraryName)]
		[UnmanagedCallConv(CallConvs = new[] { typeof(CallConvCdecl) })]
		internal static partial ulong* manifold_rs_meshgl64_run_index(IntPtr g, out nuint outLen);

		[LibraryImport(LibraryName)]
		[UnmanagedCallConv(CallConvs = new[] { typeof(CallConvCdecl) })]
		internal static partial uint* manifold_rs_meshgl64_run_original_id(IntPtr g, out nuint outLen);

		[LibraryImport(LibraryName)]
		[UnmanagedCallConv(CallConvs = new[] { typeof(CallConvCdecl) })]
		internal static partial ulong* manifold_rs_meshgl64_face_id(IntPtr g, out nuint outLen);

		[LibraryImport(LibraryName)]
		[UnmanagedCallConv(CallConvs = new[] { typeof(CallConvCdecl) })]
		internal static partial void manifold_rs_meshgl64_destroy(IntPtr g);

		[LibraryImport(LibraryName)]
		[UnmanagedCallConv(CallConvs = new[] { typeof(CallConvCdecl) })]
		internal static partial nuint manifold_rs_last_error(byte* buf, nuint bufLen);
	}
}
