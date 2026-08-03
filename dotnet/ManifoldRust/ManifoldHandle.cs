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

// Lifetime of a ManifoldRs* handle. Kept separate from Manifold so the only
// place that can free native memory is one short, obvious type.

using System;
using System.Runtime.InteropServices;

namespace ManifoldRust
{
	/// <summary>
	/// Owns a native <c>ManifoldRs*</c>. Every handle the FFI hands back non-NULL
	/// must reach <c>manifold_rs_destroy</c> exactly once, which is what the
	/// SafeHandle reference counting guarantees even if a finalizer runs while a
	/// native call is in flight.
	/// </summary>
	internal sealed class ManifoldHandle : SafeHandle
	{
		internal ManifoldHandle(IntPtr existingHandle)
			: base(IntPtr.Zero, ownsHandle: true)
		{
			this.SetHandle(existingHandle);
		}

		public override bool IsInvalid => this.handle == IntPtr.Zero;

		protected override bool ReleaseHandle()
		{
			NativeMethods.manifold_rs_destroy(this.handle);
			return true;
		}
	}
}
