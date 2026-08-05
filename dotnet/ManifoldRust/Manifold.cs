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

// The public entry point: import a triangle soup, run an n-ary boolean, export
// the result. This is the layer that applies the rules ffi/manifold_rs.h states
// - a NULL return means failure and the reason is in the thread's last-error
// slot, a non-zero status means the solid is empty rather than broken - and
// turns them into ordinary .NET exceptions and properties.

using System;
using System.Buffers;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using System.Threading;

namespace ManifoldRust
{
	/// <summary>
	/// A solid built from a triangle mesh, ready to be combined with others through
	/// <see cref="BatchBoolean(IReadOnlyList{Manifold}, ManifoldOpType)"/>.
	/// </summary>
	/// <remarks>
	/// An import that fails validation still produces an instance; it carries a
	/// non-zero <see cref="Status"/> and behaves as empty geometry everywhere else.
	/// A boolean silently absorbs such an operand, so a part that failed to import
	/// goes missing from the output instead of raising an error - check
	/// <see cref="Status"/> on every operand before combining them.
	/// </remarks>
	public sealed class Manifold : IDisposable
	{
		private readonly ManifoldHandle handle;

		private Manifold(ManifoldHandle handle)
		{
			this.handle = handle;
		}

		/// <summary>
		/// Version string of the loaded native library, for example
		/// "manifold-ffi 0.2.0 (manifold-rust 0.9.3)". Reading it also forces the
		/// native library to load, which is a cheap way to fail early if it is missing.
		/// </summary>
		/// <remarks>
		/// This is the raw string and never throws, so it can report what actually
		/// loaded even when the version handshake every other entry point runs has
		/// rejected it.
		/// </remarks>
		public static string NativeVersion
		{
			get
			{
				IntPtr text = NativeMethods.manifold_rs_version();
				return text == IntPtr.Zero ? string.Empty : (Marshal.PtrToStringUTF8(text) ?? string.Empty);
			}
		}

		/// <summary>
		/// Builds a manifold from interleaved vertex properties and triangle indices.
		/// Both arrays are copied by the native side, so the spans need not stay alive.
		/// </summary>
		/// <param name="vertProperties">
		/// <paramref name="numProp"/> floats per vertex; the first three are x, y, z.
		/// </param>
		/// <param name="triVerts">
		/// Three vertex indices per triangle, counter-clockwise seen from outside.
		/// </param>
		/// <param name="numProp">Properties per vertex, at least 3.</param>
		/// <returns>
		/// A manifold, which may carry a non-zero <see cref="Status"/> if the mesh
		/// failed validation.
		/// </returns>
		/// <exception cref="ManifoldException">
		/// The arguments cannot describe a mesh at all (num_prop below 3, a length
		/// that is not evenly divisible), or the import panicked.
		/// </exception>
		public static unsafe Manifold FromMesh(ReadOnlySpan<float> vertProperties, ReadOnlySpan<uint> triVerts, uint numProp = 3)
		{
			NativeVersionCheck.Verify();

			IntPtr result;
			fixed (float* vertPointer = vertProperties)
			fixed (uint* triPointer = triVerts)
			{
				result = NativeMethods.manifold_rs_from_mesh(
					vertPointer,
					(nuint)vertProperties.Length,
					triPointer,
					(nuint)triVerts.Length,
					numProp);
			}

			if (result == IntPtr.Zero)
			{
				// Read the message before anything else runs on this thread: the slot
				// is thread local and belongs to the call that just failed.
				throw new ManifoldException("manifold_rs_from_mesh failed", null, NativeMethods.GetLastError());
			}

			return new Manifold(new ManifoldHandle(result));
		}

		/// <summary>
		/// <see cref="FromMesh"/> via the robust (non-manifold-tolerant) import.
		/// Manifold input behaves exactly like <see cref="FromMesh"/>. Closed,
		/// orientable but non-manifold input is retained (status
		/// <see cref="ManifoldStatus.NoError"/>) for use with
		/// <see cref="BooleanEngine.Robust"/> / <see cref="BooleanEngine.Auto"/>
		/// booleans, transforms, hulls, and mesh export; other operations on
		/// such a manifold return empty results with
		/// <see cref="ManifoldStatus.NotManifold"/>. Geometry that is not even
		/// closed reports <see cref="ManifoldStatus.NotClosed"/>.
		/// </summary>
		/// <inheritdoc cref="FromMesh"/>
		public static unsafe Manifold FromMeshRobust(ReadOnlySpan<float> vertProperties, ReadOnlySpan<uint> triVerts, uint numProp = 3)
		{
			NativeVersionCheck.Verify();

			IntPtr result;
			fixed (float* vertPointer = vertProperties)
			fixed (uint* triPointer = triVerts)
			{
				result = NativeMethods.manifold_rs_from_mesh_robust(
					vertPointer,
					(nuint)vertProperties.Length,
					triPointer,
					(nuint)triVerts.Length,
					numProp);
			}

			if (result == IntPtr.Zero)
			{
				throw new ManifoldException("manifold_rs_from_mesh_robust failed", null, NativeMethods.GetLastError());
			}

			return new Manifold(new ManifoldHandle(result));
		}

		/// <summary>
		/// The process-global default boolean engine. Every boolean —
		/// <see cref="BatchBoolean(System.Collections.Generic.IReadOnlyList{Manifold}, ManifoldOpType)"/>
		/// and friends — runs on this engine. The default,
		/// <see cref="BooleanEngine.Exact"/>, preserves the exact
		/// C++-matching pipeline byte for byte; set
		/// <see cref="BooleanEngine.Auto"/> to transparently handle
		/// non-manifold operands imported via <see cref="FromMeshRobust"/>.
		/// </summary>
		public static BooleanEngine DefaultBooleanEngine
		{
			get
			{
				NativeVersionCheck.Verify();
				return (BooleanEngine)NativeMethods.manifold_rs_get_boolean_engine();
			}
			set
			{
				NativeVersionCheck.Verify();
				if (NativeMethods.manifold_rs_set_boolean_engine((int)value) != 0)
				{
					throw new ManifoldException("manifold_rs_set_boolean_engine failed", null, NativeMethods.GetLastError());
				}
			}
		}

		/// <summary>
		/// Builds a manifold from double-precision vertex properties and 64-bit
		/// triangle indices. Both arrays are copied by the native side, so the spans
		/// need not stay alive.
		/// </summary>
		/// <remarks>
		/// The kernel computes in double precision and this path is lossless for
		/// coordinates: nothing narrows on import or export, so geometry survives at
		/// full <c>double</c> precision end to end. Indexing is the one limit — the
		/// kernel indexes vertices with 32 bits, so an index above
		/// <see cref="uint.MaxValue"/> is rejected outright rather than wrapped into
		/// wrong geometry that reports <see cref="ManifoldStatus.NoError"/>.
		/// </remarks>
		/// <param name="vertProperties">
		/// <paramref name="numProp"/> doubles per vertex; the first three are x, y, z.
		/// </param>
		/// <param name="triVerts">
		/// Three vertex indices per triangle, counter-clockwise seen from outside.
		/// </param>
		/// <param name="numProp">Properties per vertex, at least 3.</param>
		/// <returns>
		/// A manifold, which may carry a non-zero <see cref="Status"/> if the mesh
		/// failed validation.
		/// </returns>
		/// <exception cref="ManifoldException">
		/// The arguments cannot describe a mesh at all (num_prop below 3, a length
		/// that is not evenly divisible, an index above <see cref="uint.MaxValue"/>),
		/// or the import panicked.
		/// </exception>
		public static unsafe Manifold FromMesh64(ReadOnlySpan<double> vertProperties, ReadOnlySpan<ulong> triVerts, uint numProp = 3)
		{
			NativeVersionCheck.Verify();

			IntPtr result;
			fixed (double* vertPointer = vertProperties)
			fixed (ulong* triPointer = triVerts)
			{
				result = NativeMethods.manifold_rs_from_mesh64(
					vertPointer,
					(nuint)vertProperties.Length,
					triPointer,
					(nuint)triVerts.Length,
					numProp);
			}

			if (result == IntPtr.Zero)
			{
				// Read the message before anything else runs on this thread: the slot
				// is thread local and belongs to the call that just failed.
				throw new ManifoldException("manifold_rs_from_mesh64 failed", null, NativeMethods.GetLastError());
			}

			return new Manifold(new ManifoldHandle(result));
		}

		/// <summary>
		/// Double-precision counterpart of <see cref="FromMeshRobust"/> — the
		/// robust (non-manifold-tolerant) import. See that method for the
		/// semantics and
		/// <see cref="FromMesh64(ReadOnlySpan{double}, ReadOnlySpan{ulong}, uint)"/>
		/// for the precision notes.
		/// </summary>
		/// <inheritdoc cref="FromMesh64(ReadOnlySpan{double}, ReadOnlySpan{ulong}, uint)"/>
		public static unsafe Manifold FromMesh64Robust(ReadOnlySpan<double> vertProperties, ReadOnlySpan<ulong> triVerts, uint numProp = 3)
		{
			NativeVersionCheck.Verify();

			IntPtr result;
			fixed (double* vertPointer = vertProperties)
			fixed (ulong* triPointer = triVerts)
			{
				result = NativeMethods.manifold_rs_from_mesh64_robust(
					vertPointer,
					(nuint)vertProperties.Length,
					triPointer,
					(nuint)triVerts.Length,
					numProp);
			}

			if (result == IntPtr.Zero)
			{
				throw new ManifoldException("manifold_rs_from_mesh64_robust failed", null, NativeMethods.GetLastError());
			}

			return new Manifold(new ManifoldHandle(result));
		}

		/// <summary>
		/// <see cref="FromMesh64(ReadOnlySpan{double}, ReadOnlySpan{ulong}, uint)"/>
		/// for callers whose coordinates are already <c>double</c> but whose indices
		/// are 32 bit - the common case, since that is what every mesh format and
		/// graphics API uses. The indices are widened into a scratch buffer; the
		/// result is identical to passing them as <c>ulong</c>.
		/// </summary>
		/// <inheritdoc cref="FromMesh64(ReadOnlySpan{double}, ReadOnlySpan{ulong}, uint)"/>
		public static Manifold FromMesh64(ReadOnlySpan<double> vertProperties, ReadOnlySpan<uint> triVerts, uint numProp = 3)
		{
			// A cube is 36 indices, so small meshes - which is what most test and
			// primitive geometry is - never touch the pool. 512 ulongs is 4 KB of
			// stack, well inside what a frame can spare.
			const int StackallocLimit = 512;

			ulong[]? rented = null;
			Span<ulong> widened = triVerts.Length <= StackallocLimit
				? stackalloc ulong[StackallocLimit]
				: (rented = ArrayPool<ulong>.Shared.Rent(triVerts.Length));

			try
			{
				// The rented array can be longer than asked for, and the stackalloc
				// always is, so slice before filling.
				widened = widened.Slice(0, triVerts.Length);
				for (int i = 0; i < triVerts.Length; i++)
				{
					widened[i] = triVerts[i];
				}

				return FromMesh64(vertProperties, widened, numProp);
			}
			finally
			{
				if (rented is not null)
				{
					ArrayPool<ulong>.Shared.Return(rented);
				}
			}
		}

		/// <summary>
		/// A copy re-tagged as an original mesh: it gets a fresh mesh ID and its
		/// boolean history is dropped, so results derived from it report that ID in
		/// <see cref="MeshGL.RunOriginalId"/>.
		/// </summary>
		/// <remarks>
		/// An empty manifold - which includes every manifold with a non-zero
		/// <see cref="Status"/> - comes back as a plain copy with
		/// <see cref="OriginalId"/> still -1, so do not assume an ID was assigned.
		/// </remarks>
		/// <exception cref="ObjectDisposedException">This manifold has been disposed.</exception>
		public Manifold AsOriginal()
		{
			this.ThrowIfDisposed();

			bool taken = false;
			try
			{
				this.handle.DangerousAddRef(ref taken);
				IntPtr result = NativeMethods.manifold_rs_as_original(this.handle.DangerousGetHandle());
				if (result == IntPtr.Zero)
				{
					throw new ManifoldException("manifold_rs_as_original failed", null, NativeMethods.GetLastError());
				}

				return new Manifold(new ManifoldHandle(result));
			}
			finally
			{
				if (taken)
				{
					this.handle.DangerousRelease();
				}
			}
		}

		/// <summary>
		/// Original mesh ID, or -1 when this is not an original mesh. Imports are not
		/// originals until <see cref="AsOriginal"/> is called.
		/// </summary>
		public int OriginalId
		{
			get
			{
				this.ThrowIfDisposed();

				bool taken = false;
				try
				{
					this.handle.DangerousAddRef(ref taken);
					return NativeMethods.manifold_rs_original_id(this.handle.DangerousGetHandle());
				}
				finally
				{
					if (taken)
					{
						this.handle.DangerousRelease();
					}
				}
			}
		}

		/// <summary>
		/// Validation state of this solid. Anything but <see cref="ManifoldStatus.NoError"/>
		/// means it is empty geometry.
		/// </summary>
		public ManifoldStatus Status
		{
			get
			{
				this.ThrowIfDisposed();

				bool taken = false;
				try
				{
					this.handle.DangerousAddRef(ref taken);
					return ReadStatus(this.handle.DangerousGetHandle());
				}
				finally
				{
					if (taken)
					{
						this.handle.DangerousRelease();
					}
				}
			}
		}

		/// <summary>
		/// Runs an n-ary boolean over the operands. The operands are not consumed and
		/// stay usable afterwards.
		/// </summary>
		/// <remarks>
		/// <see cref="ManifoldOpType.Subtract"/> is the first operand minus the union
		/// of all the others. Prefer one call with n operands over n-1 pairwise calls:
		/// the native side runs the whole set through the CSG tree, which is what makes
		/// large unions tractable.
		/// <para>
		/// An operand with a non-zero <see cref="Status"/> is absorbed as empty
		/// geometry and the result can still report
		/// <see cref="ManifoldStatus.NoError"/>, so a failed import shows up as a part
		/// silently missing from the output. Check every operand first.
		/// </para>
		/// </remarks>
		/// <exception cref="ArgumentException">The list is empty or contains a null entry.</exception>
		/// <exception cref="ObjectDisposedException">One of the operands has been disposed.</exception>
		/// <exception cref="ManifoldException">The operation failed or panicked natively.</exception>
		public static Manifold BatchBoolean(IReadOnlyList<Manifold> manifolds, ManifoldOpType op)
		{
			// This ends up calling manifold_rs_batch_boolean_ct with a NULL token,
			// which is not an approximation of the uncancellable entry point: that
			// entry point is itself implemented by delegating to _ct with NULL.
			return BatchBooleanCore(manifolds, op, null);
		}

		/// <summary>
		/// Runs an n-ary boolean that can be interrupted through
		/// <paramref name="cancellationToken"/>.
		/// </summary>
		/// <remarks>
		/// <para>
		/// Cancellation is cooperative and checked at phase boundaries inside the
		/// kernel, so the call returns shortly after cancellation is requested rather
		/// than instantly. Everything
		/// <see cref="BatchBoolean(IReadOnlyList{Manifold}, ManifoldOpType)"/> says
		/// about operands applies here too.
		/// </para>
		/// <para>
		/// <b>Completion wins.</b> Cancelling is a request, not a guarantee: if the
		/// kernel finishes before it observes the flag - which is the normal outcome
		/// when the cancel arrives near the end of the work - the completed result is
		/// returned and nothing is thrown, even though
		/// <paramref name="cancellationToken"/> is by then signalled. A caller that
		/// must treat a late cancel as a cancellation has to check the token itself
		/// after this returns.
		/// </para>
		/// </remarks>
		/// <exception cref="OperationCanceledException">
		/// <paramref name="cancellationToken"/> was signalled and the kernel observed
		/// it before finishing. The partial result is discarded before this is thrown.
		/// </exception>
		/// <exception cref="ArgumentException">The list is empty or contains a null entry.</exception>
		/// <exception cref="ObjectDisposedException">One of the operands has been disposed.</exception>
		/// <exception cref="ManifoldException">The operation failed or panicked natively.</exception>
		public static Manifold BatchBoolean(IReadOnlyList<Manifold> manifolds, ManifoldOpType op, CancellationToken cancellationToken)
		{
			// CancellationToken.None and any other token that can never be signalled
			// have nothing to cancel, so they should not pay for a native token and a
			// registration. This is also what makes the overload safe to use
			// unconditionally in a wrapper that may or may not have a real token.
			if (!cancellationToken.CanBeCanceled)
			{
				return BatchBooleanCore(manifolds, op, null);
			}

			// Argument checking happens before the token is allocated: a caller that
			// passed an empty list should get the ArgumentException without a native
			// allocation and free in between.
			ValidateOperands(manifolds);

			// Declaration order matters: `using` disposes in reverse, so the
			// registration is torn down - waiting for any callback already running -
			// before the token it cancels through is destroyed. That is the one
			// ordering the ABI leaves to the caller.
			using CancelToken token = new CancelToken();
			using CancellationTokenRegistration registration = cancellationToken.Register(
				static state => ((CancelToken)state!).Cancel(),
				token);

			Manifold result = BatchBoolean(manifolds, op, token);

			if (result.Status == ManifoldStatus.Cancelled)
			{
				result.Dispose();
				throw new OperationCanceledException(cancellationToken);
			}

			return result;
		}

		/// <summary>
		/// Runs an n-ary boolean against a <see cref="CancelToken"/> directly, for
		/// callers that want to poll the token or reuse one across several operations.
		/// </summary>
		/// <remarks>
		/// Unlike the <see cref="CancellationToken"/> overload this does not throw on
		/// cancellation: a cancelled operation comes back as a valid, empty manifold
		/// whose <see cref="Status"/> is <see cref="ManifoldStatus.Cancelled"/>, which
		/// the caller still owns and must dispose. A token that is already cancelled
		/// makes the call return immediately.
		/// </remarks>
		/// <exception cref="ArgumentException">The list is empty or contains a null entry.</exception>
		/// <exception cref="ObjectDisposedException">An operand or the token has been disposed.</exception>
		/// <exception cref="ManifoldException">The operation failed or panicked natively.</exception>
		public static Manifold BatchBoolean(IReadOnlyList<Manifold> manifolds, ManifoldOpType op, CancelToken cancelToken)
		{
			if (cancelToken is null)
			{
				throw new ArgumentNullException(nameof(cancelToken));
			}

			return BatchBooleanCore(manifolds, op, cancelToken);
		}

		/// <summary>
		/// The single implementation behind every <c>BatchBoolean</c> overload. A null
		/// <paramref name="cancelToken"/> becomes a NULL token natively, which the ABI
		/// reads as "uncancellable" and costs nothing.
		/// </summary>
		private static unsafe Manifold BatchBooleanCore(IReadOnlyList<Manifold> manifolds, ManifoldOpType op, CancelToken? cancelToken)
		{
			NativeVersionCheck.Verify();
			ValidateOperands(manifolds);

			int count = manifolds.Count;
			IntPtr[] pointers = new IntPtr[count];
			ManifoldHandle[] referenced = new ManifoldHandle[count];
			int taken = 0;
			bool tokenTaken = false;

			try
			{
				// The token is reference counted for the call exactly like the operands
				// are. Destroying it mid-call would be memory safe - the kernel clones
				// the flag on entry - but it would silently lose the ability to cancel,
				// which is worse than an exception.
				IntPtr tokenPointer = IntPtr.Zero;
				if (cancelToken is not null)
				{
					cancelToken.AddRef(ref tokenTaken);
					tokenPointer = cancelToken.Handle.DangerousGetHandle();
				}

				for (int i = 0; i < count; i++)
				{
					Manifold operand = manifolds[i]
						?? throw new ArgumentException($"Operand {i} is null.", nameof(manifolds));

					operand.ThrowIfDisposed();

					// Every operand stays reference counted for the whole native call,
					// so a concurrent Dispose or a finalizer cannot free one out from
					// under the kernel.
					bool added = false;
					operand.handle.DangerousAddRef(ref added);
					if (!added)
					{
						throw new ObjectDisposedException(nameof(Manifold));
					}

					referenced[taken++] = operand.handle;
					pointers[i] = operand.handle.DangerousGetHandle();
				}

				IntPtr result;
				fixed (IntPtr* operands = pointers)
				{
					result = NativeMethods.manifold_rs_batch_boolean_ct(operands, (nuint)count, (int)op, tokenPointer);
				}

				// Cancellation comes back as a valid handle with a Cancelled status, not
				// as NULL, so this really is only argument errors and caught panics.
				if (result == IntPtr.Zero)
				{
					throw new ManifoldException($"manifold_rs_batch_boolean_ct ({op}, {count} operands) failed", null, NativeMethods.GetLastError());
				}

				return new Manifold(new ManifoldHandle(result));
			}
			finally
			{
				for (int i = 0; i < taken; i++)
				{
					referenced[i].DangerousRelease();
				}

				if (tokenTaken)
				{
					cancelToken!.Handle.DangerousRelease();
				}
			}
		}

		/// <summary>
		/// Exports the solid as a managed mesh. A manifold with a non-zero
		/// <see cref="Status"/> exports as empty arrays rather than failing - the C++
		/// binding crashes on that path, this one does not.
		/// </summary>
		/// <exception cref="ObjectDisposedException">This manifold has been disposed.</exception>
		/// <exception cref="ManifoldException">The export panicked natively.</exception>
		public MeshGL GetMeshGL()
		{
			this.ThrowIfDisposed();

			bool taken = false;
			try
			{
				this.handle.DangerousAddRef(ref taken);
				IntPtr raw = this.handle.DangerousGetHandle();

				IntPtr mesh = NativeMethods.manifold_rs_get_meshgl(raw);
				if (mesh == IntPtr.Zero)
				{
					// Read the last error first: it belongs to the export call, and any
					// later native call on this thread could overwrite it.
					string nativeError = NativeMethods.GetLastError();
					throw new ManifoldException("manifold_rs_get_meshgl failed", ReadStatus(raw), nativeError);
				}

				try
				{
					return MeshGL.CopyFrom(mesh);
				}
				finally
				{
					// The arrays MeshGL copied borrow from this handle, so it can only be
					// destroyed once the copies are made.
					NativeMethods.manifold_rs_meshgl_destroy(mesh);
				}
			}
			finally
			{
				if (taken)
				{
					this.handle.DangerousRelease();
				}
			}
		}

		/// <summary>
		/// Exports the solid as a managed mesh with double-precision coordinates. A
		/// manifold with a non-zero <see cref="Status"/> exports as empty arrays rather
		/// than failing.
		/// </summary>
		/// <remarks>
		/// Coordinates come back at full <c>double</c> precision straight from the
		/// f64 kernel — see <see cref="MeshGL64"/>, which also explains why the
		/// index arrays still come back as <c>uint[]</c>.
		/// </remarks>
		/// <exception cref="ObjectDisposedException">This manifold has been disposed.</exception>
		/// <exception cref="ManifoldException">
		/// The export panicked natively, or an exported index did not fit in 32 bits.
		/// </exception>
		public MeshGL64 GetMeshGL64()
		{
			this.ThrowIfDisposed();

			bool taken = false;
			try
			{
				this.handle.DangerousAddRef(ref taken);
				IntPtr raw = this.handle.DangerousGetHandle();

				IntPtr mesh = NativeMethods.manifold_rs_get_meshgl64(raw);
				if (mesh == IntPtr.Zero)
				{
					// Read the last error first: it belongs to the export call, and any
					// later native call on this thread could overwrite it.
					string nativeError = NativeMethods.GetLastError();
					throw new ManifoldException("manifold_rs_get_meshgl64 failed", ReadStatus(raw), nativeError);
				}

				try
				{
					return MeshGL64.CopyFrom(mesh);
				}
				finally
				{
					// The arrays MeshGL64 copied borrow from this handle, so it can only
					// be destroyed once the copies are made.
					NativeMethods.manifold_rs_meshgl64_destroy(mesh);
				}
			}
			finally
			{
				if (taken)
				{
					this.handle.DangerousRelease();
				}
			}
		}

		/// <summary>
		/// Releases the native solid. Safe to call more than once; the underlying
		/// SafeHandle also frees it from the finalizer if this is never called.
		/// </summary>
		public void Dispose()
		{
			this.handle.Dispose();
		}

		/// <summary>
		/// The argument checks every <c>BatchBoolean</c> overload shares. Split out so
		/// the cancellable path can run them before it allocates a native token, and
		/// so a bad argument never costs an allocate/free pair.
		/// </summary>
		/// <remarks>
		/// Only the list itself is checked here. Whether an individual operand is
		/// disposed cannot be settled without reference counting it, so that stays
		/// where the reference counting happens.
		/// </remarks>
		private static void ValidateOperands(IReadOnlyList<Manifold> manifolds)
		{
			if (manifolds is null)
			{
				throw new ArgumentNullException(nameof(manifolds));
			}

			if (manifolds.Count == 0)
			{
				throw new ArgumentException("A boolean needs at least one operand.", nameof(manifolds));
			}
		}

		private static ManifoldStatus ReadStatus(IntPtr raw)
		{
			int status = NativeMethods.manifold_rs_status(raw);
			if (status < 0)
			{
				// -1 means a NULL handle, which cannot happen for a live wrapper.
				throw new ManifoldException($"manifold_rs_status returned {status} for a live handle");
			}

			return (ManifoldStatus)status;
		}

		private void ThrowIfDisposed()
		{
			ObjectDisposedException.ThrowIf(this.handle.IsClosed || this.handle.IsInvalid, this);
		}
	}
}
