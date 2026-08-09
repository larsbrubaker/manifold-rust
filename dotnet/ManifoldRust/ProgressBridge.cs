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

// The managed side of the progress callback in ffi/manifold_rs.h: the static
// unmanaged entry point the kernel calls, the phase-id to name table, and the
// per-call state that carries the caller's IProgress across the boundary.
//
// Used only by Manifold.Boolean (ManifoldBoolean.cs). Nothing here is public:
// callers see an IProgress<(string, double?)> and never a phase id.

using System;
using System.Globalization;
using System.Runtime.CompilerServices;
using System.Runtime.ExceptionServices;
using System.Runtime.InteropServices;
using System.Threading;

namespace ManifoldRust
{
	/// <summary>
	/// Phase id to display name, as the native table defines them. The names are
	/// static native strings, so they are read once and cached: a callback that
	/// fires a hundred times a phase must not marshal the same string a hundred
	/// times.
	/// </summary>
	internal static class ProgressPhases
	{
		// Lazy rather than a static initializer: the table can only be read once
		// the native library is loaded, and Lazy also settles the race between two
		// booleans starting at the same time without a lock in the callback path.
		private static readonly Lazy<string[]> Table = new Lazy<string[]>(ReadTable);

		/// <summary>
		/// Reads the table if it has not been read yet. Called before the boolean
		/// starts so the first callback - which may arrive on a native worker
		/// thread - never has to make a native call of its own.
		/// </summary>
		internal static void Preload()
		{
			_ = Table.Value;
		}

		/// <summary>
		/// The display name for a phase id, or a <c>phase N</c> placeholder for an
		/// id this build of the native library does not name. Phase ids are append
		/// only, so a newer native can report an id the table does not cover; the
		/// header's rule is to fall back to the raw id rather than fail.
		/// </summary>
		internal static string Name(uint phaseId)
		{
			string[] table = Table.Value;
			return phaseId < (uint)table.Length
				? table[phaseId]
				: string.Create(CultureInfo.InvariantCulture, $"phase {phaseId}");
		}

		private static string[] ReadTable()
		{
			uint count = NativeMethods.manifold_rs_progress_phase_count();
			string[] table = new string[count];
			for (uint id = 0; id < count; id++)
			{
				IntPtr text = NativeMethods.manifold_rs_progress_phase_name(id);
				table[id] = text == IntPtr.Zero
					? string.Create(CultureInfo.InvariantCulture, $"phase {id}")
					: (Marshal.PtrToStringUTF8(text) ?? string.Empty);
			}

			return table;
		}
	}

	/// <summary>
	/// One boolean call's progress state: the caller's sink plus the first
	/// exception it threw. Handed to the native side as a <see cref="GCHandle"/>
	/// in the opaque <c>user</c> word.
	/// </summary>
	/// <remarks>
	/// The native side promises the callback is never re-entered concurrently but
	/// explicitly allows it to run on an internal worker thread, so this type is
	/// written for "one thread at a time, not necessarily the same one":
	/// <see cref="Failure"/> is published with a volatile write and read after the
	/// native call has returned, which is a full ordering point.
	/// </remarks>
	internal sealed unsafe class ProgressBridge
	{
		private readonly IProgress<(string Phase, double? Fraction)> sink;
		private ExceptionDispatchInfo? failure;

		internal ProgressBridge(IProgress<(string Phase, double? Fraction)> sink)
		{
			this.sink = sink;
		}

		/// <summary>
		/// The function pointer to hand to the native side. A static
		/// <c>[UnmanagedCallersOnly]</c> method rather than a marshalled delegate:
		/// there is no stub to keep rooted, so the callback cannot outlive its
		/// trampoline no matter which thread it fires on, and it is AOT safe.
		/// </summary>
		internal static delegate* unmanaged[Cdecl]<uint, double, void*, void> Callback => &Invoke;

		/// <summary>
		/// The first exception the caller's <see cref="IProgress{T}"/> threw, or
		/// null. Rethrown by the caller once the native call has returned - an
		/// exception must never unwind through
		/// <see cref="UnmanagedCallersOnlyAttribute"/>, which tears the process
		/// down.
		/// </summary>
		internal ExceptionDispatchInfo? Failure => Volatile.Read(ref this.failure);

		[UnmanagedCallersOnly(CallConvs = new[] { typeof(CallConvCdecl) })]
		private static void Invoke(uint phaseId, double fraction, void* user)
		{
			// Everything in here is inside an unmanaged frame, so nothing may throw.
			try
			{
				if (user is null)
				{
					return;
				}

				if (GCHandle.FromIntPtr((IntPtr)user).Target is ProgressBridge bridge)
				{
					bridge.Report(phaseId, fraction);
				}
			}
			catch
			{
				// Only reachable if the handle itself is bad, which would already be
				// a corrupted call. Swallowing beats crashing the process.
			}
		}

		private void Report(uint phaseId, double fraction)
		{
			try
			{
				// A negative fraction is the ABI's "indeterminate": the phase is
				// running but has no meaningful ratio, which is a spinner rather
				// than a bar. It maps to null, not to 0.
				this.sink.Report((ProgressPhases.Name(phaseId), fraction < 0.0 ? null : fraction));
			}
			catch (Exception exception)
			{
				if (this.failure is null)
				{
					Volatile.Write(ref this.failure, ExceptionDispatchInfo.Capture(exception));
				}
			}
		}
	}
}
