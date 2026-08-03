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
using System.Collections.Generic;
using System.Runtime.InteropServices;

namespace ManifoldRust
{
	/// <summary>
	/// A solid built from a triangle mesh, ready to be combined with others through
	/// <see cref="BatchBoolean"/>.
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
		/// "manifold-ffi 0.1.0 (manifold-rust 0.9.3)". Reading it also forces the
		/// native library to load, which is a cheap way to fail early if it is missing.
		/// </summary>
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
		public static unsafe Manifold BatchBoolean(IReadOnlyList<Manifold> manifolds, ManifoldOpType op)
		{
			if (manifolds is null)
			{
				throw new ArgumentNullException(nameof(manifolds));
			}

			int count = manifolds.Count;
			if (count == 0)
			{
				throw new ArgumentException("A boolean needs at least one operand.", nameof(manifolds));
			}

			IntPtr[] pointers = new IntPtr[count];
			ManifoldHandle[] referenced = new ManifoldHandle[count];
			int taken = 0;

			try
			{
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
					result = NativeMethods.manifold_rs_batch_boolean(operands, (nuint)count, (int)op);
				}

				if (result == IntPtr.Zero)
				{
					throw new ManifoldException($"manifold_rs_batch_boolean ({op}, {count} operands) failed", null, NativeMethods.GetLastError());
				}

				return new Manifold(new ManifoldHandle(result));
			}
			finally
			{
				for (int i = 0; i < taken; i++)
				{
					referenced[i].DangerousRelease();
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
		/// Releases the native solid. Safe to call more than once; the underlying
		/// SafeHandle also frees it from the finalizer if this is never called.
		/// </summary>
		public void Dispose()
		{
			this.handle.Dispose();
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
