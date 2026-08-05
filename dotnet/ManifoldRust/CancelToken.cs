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

// The cooperative cancellation flag. This is the one handle in the ABI that is
// meant to be touched from two threads at once - one inside the boolean, one
// calling Cancel - so it is also the one place where the general "a handle
// belongs to one thread at a time" rule does not apply.
//
// Most callers should use the CancellationToken overload of
// Manifold.BatchBoolean instead; this type is for callers that want to poll a
// token or reuse one across several operations.

using System;

namespace ManifoldRust
{
	/// <summary>
	/// A cancellation flag that can interrupt a long-running operation from another
	/// thread. Hand it to
	/// <see cref="Manifold.BatchBoolean(System.Collections.Generic.IReadOnlyList{Manifold}, ManifoldOpType, CancelToken)"/>
	/// and call <see cref="Cancel"/> from wherever the user's Cancel button lives.
	/// </summary>
	/// <remarks>
	/// <para>
	/// Cancellation is <b>sticky</b>: a cancelled token stays cancelled and every
	/// later operation given it returns immediately as
	/// <see cref="ManifoldStatus.Cancelled"/>. Create a new token for work that
	/// should be allowed to run.
	/// </para>
	/// <para>
	/// <see cref="Cancel"/> and <see cref="IsCancelled"/> may be called from any
	/// thread at any time, including while an operation holding this token is
	/// running - that is what the type is for.
	/// </para>
	/// <para>
	/// Disposing a token while an operation that was given it is still running is
	/// memory safe (the operation keeps its own share of the underlying flag), but
	/// it means the operation can no longer be cancelled. The hard rule is the
	/// other one: <b>do not race <see cref="Dispose"/> against <see cref="Cancel"/>
	/// or <see cref="IsCancelled"/> on the same token</b>. Keep the token alive
	/// until every call that might cancel through it has returned.
	/// </para>
	/// </remarks>
	public sealed class CancelToken : IDisposable
	{
		private readonly CancelTokenHandle handle;

		/// <summary>
		/// Creates a fresh, uncancelled token.
		/// </summary>
		/// <exception cref="ManifoldException">The native allocation failed.</exception>
		public CancelToken()
		{
			NativeVersionCheck.Verify();

			IntPtr raw = NativeMethods.manifold_rs_cancel_token_new();
			if (raw == IntPtr.Zero)
			{
				throw new ManifoldException("manifold_rs_cancel_token_new failed", null, NativeMethods.GetLastError());
			}

			this.handle = new CancelTokenHandle(raw);
		}

		/// <summary>
		/// Requests cancellation. Callable from any thread; calling it more than once
		/// is a no-op.
		/// </summary>
		/// <exception cref="ObjectDisposedException">This token has been disposed.</exception>
		public void Cancel()
		{
			this.ThrowIfDisposed();

			bool taken = false;
			try
			{
				this.handle.DangerousAddRef(ref taken);
				NativeMethods.manifold_rs_cancel_token_cancel(this.handle.DangerousGetHandle());
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
		/// Whether cancellation has been requested.
		/// </summary>
		/// <exception cref="ObjectDisposedException">This token has been disposed.</exception>
		public bool IsCancelled
		{
			get
			{
				this.ThrowIfDisposed();

				bool taken = false;
				try
				{
					this.handle.DangerousAddRef(ref taken);
					return NativeMethods.manifold_rs_cancel_token_is_cancelled(this.handle.DangerousGetHandle()) != 0;
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
		/// Releases the native token. Safe to call more than once; the underlying
		/// SafeHandle also frees it from the finalizer if this is never called.
		/// </summary>
		public void Dispose()
		{
			this.handle.Dispose();
		}

		/// <summary>
		/// The raw handle for a native call, or <see cref="IntPtr.Zero"/> for no
		/// token. The caller must have reference counted the handle first - see
		/// <see cref="AddRef"/>.
		/// </summary>
		internal CancelTokenHandle Handle => this.handle;

		/// <summary>
		/// Reference counts the handle for the duration of a native call, so a
		/// concurrent Dispose or a finalizer cannot free it out from under the kernel.
		/// </summary>
		internal void AddRef(ref bool taken)
		{
			this.ThrowIfDisposed();

			this.handle.DangerousAddRef(ref taken);
			if (!taken)
			{
				throw new ObjectDisposedException(nameof(CancelToken));
			}
		}

		private void ThrowIfDisposed()
		{
			ObjectDisposedException.ThrowIf(this.handle.IsClosed || this.handle.IsInvalid, this);
		}
	}
}
