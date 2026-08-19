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

// The mesh-health entry points of ffi/manifold_rs.h: shell orientation repair,
// the full winding-number rebuild, and the self-intersection scan. All three
// are about deciding whether geometry is fit for the exact boolean engine -
// and making it fit - before handing it over, so they sit together, apart from
// the import/export surface in Manifold.cs.

using System;
using System.Threading;

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
		/// A fresh, properly paired solid enclosing the same material this mesh
		/// bounds under <paramref name="rule"/>, rebuilt from the winding numbers
		/// rather than repaired in place.
		/// </summary>
		/// <remarks>
		/// <para>
		/// This is the heavy sibling of <see cref="RepairOrientation"/>, and they
		/// fix different things. <see cref="RepairOrientation"/> changes nothing
		/// but triangle winding: it is fast, exact, and hands back the same
		/// triangles, vertices and properties. <c>RebuildSolid</c> runs the full
		/// robust pipeline on this one mesh - exact self-intersection, arrangement,
		/// winding classification, reassembly - so it also removes what winding
		/// alone cannot: self-intersections, doubled or coincident sheets, and
		/// interior walls with material on both sides. The price is that the result
		/// is re-triangulated. Triangle and vertex counts change and per-vertex
		/// properties are re-interpolated, so nothing that indexes into the
		/// original mesh survives.
		/// </para>
		/// <para>
		/// Reach for <see cref="RepairOrientation"/> first when the geometry is
		/// sound and only the winding is wrong; reach for this when the geometry
		/// itself is wrong.
		/// </para>
		/// <para>
		/// The input must already be closed and orientable, which is what
		/// <see cref="FromMeshRobust"/> admits. An open mesh never gets that far -
		/// it imports as an empty manifold with
		/// <see cref="ManifoldStatus.NotClosed"/> - and rebuilding it is a no-op.
		/// </para>
		/// </remarks>
		/// <param name="rule">
		/// Which winding numbers count as solid.
		/// <see cref="WindingRule.Positive"/> (<c>{w &gt;= 1}</c>) discards an
		/// inside-out body as enclosing nothing; <see cref="WindingRule.Nonzero"/>
		/// (<c>{w != 0}</c>) reads it as material and rewinds it outward, which is
		/// what scans and CAD exports with arbitrarily wound shells want.
		/// </param>
		/// <returns>A new manifold; the original is unchanged and still owned by the caller.</returns>
		/// <exception cref="ObjectDisposedException">This manifold has been disposed.</exception>
		/// <exception cref="ManifoldException">
		/// The rebuild failed or panicked natively, which includes an unrecognised
		/// <paramref name="rule"/> value.
		/// </exception>
		public Manifold RebuildSolid(WindingRule rule)
		{
			return this.RebuildSolidCore(rule, null);
		}

		/// <summary>
		/// <see cref="RebuildSolid(WindingRule)"/> against a
		/// <see cref="CancelToken"/> directly, for callers that want to poll the
		/// token or reuse one across several operations.
		/// </summary>
		/// <param name="rule">Which winding numbers count as solid, exactly as in <see cref="RebuildSolid(WindingRule)"/>.</param>
		/// <param name="cancelToken">
		/// The cancellation flag, or null for an uncancellable call. Unlike the
		/// <see cref="CancellationToken"/> overload this does not throw: a
		/// cancelled rebuild comes back as a valid, empty manifold whose
		/// <see cref="Status"/> is <see cref="ManifoldStatus.Cancelled"/>.
		/// </param>
		/// <inheritdoc cref="RebuildSolid(WindingRule)"/>
		public Manifold RebuildSolid(WindingRule rule, CancelToken? cancelToken)
		{
			return this.RebuildSolidCore(rule, cancelToken);
		}

		/// <summary>
		/// <see cref="RebuildSolid(WindingRule)"/> under a
		/// <see cref="CancellationToken"/>. A rebuild costs what a boolean against
		/// a partner costs, so an interactive caller wants one of the cancellable
		/// forms.
		/// </summary>
		/// <exception cref="OperationCanceledException">
		/// <paramref name="cancellationToken"/> was signalled and the kernel
		/// observed it before finishing. The partial result is disposed first.
		/// </exception>
		/// <inheritdoc cref="RebuildSolid(WindingRule)"/>
		public Manifold RebuildSolid(WindingRule rule, CancellationToken cancellationToken)
		{
			// A token that can never be signalled has nothing to cancel, so it must
			// not pay for a native token and a registration.
			if (!cancellationToken.CanBeCanceled)
			{
				return this.RebuildSolidCore(rule, null);
			}

			// Declaration order matters: `using` disposes in reverse, so the
			// registration is torn down - waiting for any callback already running -
			// before the token it cancels through is destroyed.
			using CancelToken token = new CancelToken();
			using CancellationTokenRegistration registration = cancellationToken.Register(
				static state => ((CancelToken)state!).Cancel(),
				token);

			Manifold result = this.RebuildSolidCore(rule, token);

			if (result.Status == ManifoldStatus.Cancelled)
			{
				result.Dispose();
				throw new OperationCanceledException(cancellationToken);
			}

			return result;
		}

		/// <summary>
		/// The one place the rebuild crosses the ABI. Both handles are ref-counted
		/// for exactly the native call, and a null token is passed as
		/// <see cref="IntPtr.Zero"/>, which the ABI documents as the uncancellable
		/// path.
		/// </summary>
		private Manifold RebuildSolidCore(WindingRule rule, CancelToken? cancelToken)
		{
			this.ThrowIfDisposed();

			bool taken = false;
			bool tokenTaken = false;
			try
			{
				IntPtr tokenPointer = IntPtr.Zero;
				if (cancelToken is not null)
				{
					cancelToken.AddRef(ref tokenTaken);
					tokenPointer = cancelToken.Handle.DangerousGetHandle();
				}

				this.handle.DangerousAddRef(ref taken);
				if (!taken)
				{
					throw new ObjectDisposedException(nameof(Manifold));
				}

				IntPtr result = NativeMethods.manifold_rs_rebuild_solid_ct(
					this.handle.DangerousGetHandle(),
					(int)rule,
					tokenPointer);

				// Cancellation comes back as a valid handle with a Cancelled status,
				// not as NULL, so this really is only argument errors and panics.
				if (result == IntPtr.Zero)
				{
					throw new ManifoldException($"manifold_rs_rebuild_solid ({rule}) failed", null, NativeMethods.GetLastError());
				}

				return new Manifold(new ManifoldHandle(result));
			}
			finally
			{
				if (taken)
				{
					this.handle.DangerousRelease();
				}

				if (tokenTaken)
				{
					cancelToken!.Handle.DangerousRelease();
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
