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

// The two mesh-health entry points of ffi/manifold_rs.h: shell orientation
// repair and the self-intersection scan. Both are about deciding whether
// geometry is fit for the exact boolean engine before handing it over, so they
// sit together, apart from the import/export surface in Manifold.cs.

using System;

namespace ManifoldRust
{
	public sealed partial class Manifold
	{
		/// <summary>
		/// A copy with inside-out shells rewound so every body reads as solid
		/// material under the robust engine's <see cref="WindingRule.Positive"/>
		/// semantics.
		/// </summary>
		/// <remarks>
		/// <para>
		/// Outermost shells end up wound outward (winding +1), legitimate cavity
		/// shells stay - or become - inward wound, and coincident or doubled
		/// sheets are left alone. It works on both strict and
		/// <see cref="FromMeshRobust"/> imports, before and independent of any
		/// boolean, so it is the fix-up to run on scanned or third-party geometry
		/// that comes in wound inconsistently.
		/// </para>
		/// <para>
		/// A mesh that needs no repair comes back as a plain copy, so this is safe
		/// to call unconditionally. Repairing is the alternative to
		/// <see cref="WindingRule.Nonzero"/>: it corrects the data once instead of
		/// changing what every later operation means by "solid".
		/// </para>
		/// </remarks>
		/// <returns>A new manifold; the original is unchanged and still owned by the caller.</returns>
		/// <exception cref="ObjectDisposedException">This manifold has been disposed.</exception>
		/// <exception cref="ManifoldException">The repair panicked natively.</exception>
		public Manifold RepairOrientation()
		{
			this.ThrowIfDisposed();

			bool taken = false;
			try
			{
				this.handle.DangerousAddRef(ref taken);
				IntPtr result = NativeMethods.manifold_rs_repair_orientation(this.handle.DangerousGetHandle());
				if (result == IntPtr.Zero)
				{
					throw new ManifoldException("manifold_rs_repair_orientation failed", null, NativeMethods.GetLastError());
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
		/// Whether two of this solid's own triangles genuinely intersect - they
		/// cross, they overlap, or they are coincident surface (a doubled or
		/// multiply wound sheet).
		/// </summary>
		/// <remarks>
		/// <para>
		/// Merely sharing edges and vertices, as every closed mesh does, is not a
		/// self-intersection. A mesh with non-finite positions answers
		/// <c>true</c>: that is the safe verdict for geometry no exact predicate
		/// can evaluate.
		/// </para>
		/// <para>
		/// A mesh can be topologically manifold and still answer <c>true</c>. Such
		/// input breaks the assumptions of <see cref="BooleanEngine.Exact"/> - its
		/// winding integral over-counts doubly wound material - which is why
		/// <see cref="BooleanEngine.Auto"/> routes it to
		/// <see cref="BooleanEngine.Robust"/> instead.
		/// </para>
		/// <para>
		/// The verdict is cached natively on the mesh, so this is a property
		/// rather than a method: repeated reads are free after the first.
		/// </para>
		/// </remarks>
		/// <exception cref="ObjectDisposedException">This manifold has been disposed.</exception>
		/// <exception cref="ManifoldException">The scan panicked natively.</exception>
		public bool HasSelfIntersections
		{
			get
			{
				this.ThrowIfDisposed();

				bool taken = false;
				try
				{
					this.handle.DangerousAddRef(ref taken);
					int verdict = NativeMethods.manifold_rs_has_self_intersections(this.handle.DangerousGetHandle());
					if (verdict < 0)
					{
						// -1 is a NULL handle or a caught panic; a live wrapper rules
						// out the first, so this is the second.
						throw new ManifoldException("manifold_rs_has_self_intersections failed", null, NativeMethods.GetLastError());
					}

					return verdict != 0;
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
	}
}
