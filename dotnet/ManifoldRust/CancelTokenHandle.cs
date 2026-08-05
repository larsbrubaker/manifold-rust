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

// Lifetime of a CancelTokenRs* handle. Same shape as ManifoldHandle, kept in its
// own type so there is exactly one place that can free a token.

using System;
using System.Runtime.InteropServices;

namespace ManifoldRust
{
	/// <summary>
	/// Owns a native <c>CancelTokenRs*</c>.
	/// </summary>
	/// <remarks>
	/// Destroying a token while an operation that was given it is still running is
	/// memory safe - the operation cloned the shared flag when it started - so the
	/// SafeHandle reference counting here is about the other direction: it keeps
	/// <c>cancel</c> and <c>is_cancelled</c>, which read the handle itself, from
	/// racing a <see cref="System.IDisposable.Dispose"/> or a finalizer.
	/// </remarks>
	internal sealed class CancelTokenHandle : SafeHandle
	{
		internal CancelTokenHandle(IntPtr existingHandle)
			: base(IntPtr.Zero, ownsHandle: true)
		{
			this.SetHandle(existingHandle);
		}

		public override bool IsInvalid => this.handle == IntPtr.Zero;

		protected override bool ReleaseHandle()
		{
			NativeMethods.manifold_rs_cancel_token_destroy(this.handle);
			return true;
		}
	}
}
