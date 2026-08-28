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

// The binary boolean of ffi/manifold_rs.h - manifold_rs_boolean_progress_rule -
// and the Union / Subtract / Intersect spellings of it. This is the entry point
// that takes the engine, the winding rule and a progress sink explicitly, rather
// than reading the process-global engine the way BatchBoolean in Manifold.cs
// does. Use BatchBoolean for three or more operands: it runs the CSG tree, which
// is what makes large unions tractable, and this one cannot.

using System;
using System.Runtime.InteropServices;
using System.Threading;

namespace ManifoldRust
{
	public sealed partial class Manifold
	{
		/// <summary>
		/// Combines two solids on an explicitly chosen engine, under the default
		/// <see cref="WindingRule.Positive"/> rule.
		/// </summary>
		/// <param name="a">The first operand; for <see cref="ManifoldOpType.Subtract"/> the one being cut.</param>
		/// <param name="b">The second operand.</param>
		/// <param name="op">Which boolean to run.</param>
		/// <param name="engine">
		/// Which implementation runs it. This is taken explicitly rather than from
		/// <see cref="DefaultBooleanEngine"/>, because a caller reaching for this
		/// overload usually wants the engine pinned.
		/// </param>
		/// <inheritdoc cref="Boolean(Manifold, Manifold, ManifoldOpType, BooleanEngine, WindingRule, IProgress{ValueTuple{string, double?}}, CancelToken)"/>
		public static Manifold Boolean(Manifold a, Manifold b, ManifoldOpType op, BooleanEngine engine)
		{
			return BooleanCore(a, b, op, engine, WindingRule.Positive, null, null);
		}

		/// <summary>
		/// Combines two solids on an explicitly chosen engine and winding rule.
		/// </summary>
		/// <inheritdoc cref="Boolean(Manifold, Manifold, ManifoldOpType, BooleanEngine, WindingRule, IProgress{ValueTuple{string, double?}}, CancelToken)"/>
		public static Manifold Boolean(Manifold a, Manifold b, ManifoldOpType op, BooleanEngine engine, WindingRule windingRule)
		{
			return BooleanCore(a, b, op, engine, windingRule, null, null);
		}

		/// <summary>
		/// Combines two solids, reporting pipeline progress to
		/// <paramref name="progress"/> as it runs.
		/// </summary>
		/// <inheritdoc cref="Boolean(Manifold, Manifold, ManifoldOpType, BooleanEngine, WindingRule, IProgress{ValueTuple{string, double?}}, CancelToken)"/>
		public static Manifold Boolean(
			Manifold a,
			Manifold b,
			ManifoldOpType op,
			BooleanEngine engine,
			WindingRule windingRule,
			IProgress<(string Phase, double? Fraction)>? progress)
		{
			return BooleanCore(a, b, op, engine, windingRule, progress, null);
		}

		/// <summary>
		/// Combines two solids with progress reporting and cooperative
		/// cancellation, throwing <see cref="OperationCanceledException"/> if the
		/// kernel observes the cancellation before it finishes.
		/// </summary>
		/// <remarks>
		/// <para>
		/// <b>Completion wins.</b> Cancelling is a request: if the kernel finishes
		/// before it observes the flag - the normal outcome when the cancel arrives
		/// near the end - the completed result is returned and nothing is thrown,
		/// even though <paramref name="cancellationToken"/> is by then signalled.
		/// This mirrors
		/// <see cref="BatchBoolean(System.Collections.Generic.IReadOnlyList{Manifold}, ManifoldOpType, CancellationToken)"/>.
		/// </para>
		/// </remarks>
		/// <exception cref="OperationCanceledException">
		/// <paramref name="cancellationToken"/> was signalled and the kernel
		/// observed it before finishing. The partial result is disposed first.
		/// </exception>
		/// <inheritdoc cref="Boolean(Manifold, Manifold, ManifoldOpType, BooleanEngine, WindingRule, IProgress{ValueTuple{string, double?}}, CancelToken)"/>
		public static Manifold Boolean(
			Manifold a,
			Manifold b,
			ManifoldOpType op,
			BooleanEngine engine,
			WindingRule windingRule,
			IProgress<(string Phase, double? Fraction)>? progress,
			CancellationToken cancellationToken)
		{
			// A token that can never be signalled has nothing to cancel, so it must
			// not pay for a native token and a registration.
			if (!cancellationToken.CanBeCanceled)
			{
				return BooleanCore(a, b, op, engine, windingRule, progress, null);
			}

			// Declaration order matters: `using` disposes in reverse, so the
			// registration is torn down - waiting for any callback already running -
			// before the token it cancels through is destroyed.
			using CancelToken token = new CancelToken();
			using CancellationTokenRegistration registration = cancellationToken.Register(
				static state => ((CancelToken)state!).Cancel(),
				token);

			Manifold result = BooleanCore(a, b, op, engine, windingRule, progress, token);

			if (result.Status == ManifoldStatus.Cancelled)
			{
				result.Dispose();
				throw new OperationCanceledException(cancellationToken);
			}

			return result;
		}

		/// <summary>
		/// Combines two solids with progress reporting against a
		/// <see cref="CancelToken"/> directly, for callers that want to poll the
		/// token or reuse one across several operations.
		/// </summary>
		/// <param name="a">The first operand; for <see cref="ManifoldOpType.Subtract"/> the one being cut.</param>
		/// <param name="b">The second operand.</param>
		/// <param name="op">Which boolean to run.</param>
		/// <param name="engine">
		/// Which implementation runs it, taken explicitly rather than from
		/// <see cref="DefaultBooleanEngine"/>.
		/// </param>
		/// <param name="windingRule">
		/// Which winding numbers count as solid.
		/// <see cref="WindingRule.Positive"/> is the default and the rule clean,
		/// consistently wound data wants; <see cref="WindingRule.Nonzero"/> keeps
		/// inside-out shells as material. The rule is a robust-engine semantic:
		/// <see cref="BooleanEngine.Exact"/> ignores it, and
		/// <see cref="BooleanEngine.Auto"/> resolves to
		/// <see cref="BooleanEngine.Robust"/> whenever it is
		/// <see cref="WindingRule.Nonzero"/>.
		/// </param>
		/// <param name="progress">
		/// Receives the phase name and its completion fraction as the pipeline
		/// runs, or null for no reporting. A null <c>Fraction</c> means the phase
		/// is indeterminate - show a spinner rather than a bar.
		/// </param>
		/// <param name="cancelToken">
		/// The cancellation flag, or null for an uncancellable call. Unlike the
		/// <see cref="CancellationToken"/> overload this does not throw: a
		/// cancelled operation comes back as a valid, empty manifold whose
		/// <see cref="Status"/> is <see cref="ManifoldStatus.Cancelled"/>.
		/// </param>
		/// <remarks>
		/// <para>
		/// <b>The progress callback may run on a native worker thread</b>, not the
		/// thread that made this call - the pipeline is parallel. It is never
		/// re-entered concurrently, but a sink that touches UI state must marshal
		/// itself (<see cref="System.Progress{T}"/> already posts to the
		/// synchronization context that created it). Reporting never changes the
		/// result: the same inputs produce the same triangles with or without a
		/// sink attached.
		/// </para>
		/// <para>
		/// An exception thrown by <paramref name="progress"/> cannot unwind through
		/// the native frame, so the first one is captured and rethrown here once
		/// the boolean has returned, with the result disposed.
		/// </para>
		/// </remarks>
		/// <returns>A new manifold; both operands are unchanged and still owned by the caller.</returns>
		/// <exception cref="ArgumentNullException"><paramref name="a"/> or <paramref name="b"/> is null.</exception>
		/// <exception cref="ObjectDisposedException">An operand or the token has been disposed.</exception>
		/// <exception cref="ManifoldException">
		/// The operation failed or panicked natively, which includes an
		/// unrecognised <paramref name="op"/>, <paramref name="engine"/> or
		/// <paramref name="windingRule"/> value.
		/// </exception>
		public static Manifold Boolean(
			Manifold a,
			Manifold b,
			ManifoldOpType op,
			BooleanEngine engine,
			WindingRule windingRule,
			IProgress<(string Phase, double? Fraction)>? progress,
			CancelToken? cancelToken)
		{
			return BooleanCore(a, b, op, engine, windingRule, progress, cancelToken);
		}

		/// <summary>
		/// The union of two solids on an explicitly chosen engine, under the
		/// default <see cref="WindingRule.Positive"/> rule.
		/// </summary>
		/// <inheritdoc cref="Boolean(Manifold, Manifold, ManifoldOpType, BooleanEngine, WindingRule, IProgress{ValueTuple{string, double?}}, CancelToken)"/>
		public static Manifold Union(Manifold a, Manifold b, BooleanEngine engine)
		{
			return BooleanCore(a, b, ManifoldOpType.Add, engine, WindingRule.Positive, null, null);
		}

		/// <summary>
		/// The union of two solids, with an explicit winding rule, progress
		/// reporting and cancellation.
		/// </summary>
		/// <inheritdoc cref="Boolean(Manifold, Manifold, ManifoldOpType, BooleanEngine, WindingRule, IProgress{ValueTuple{string, double?}}, CancellationToken)"/>
		public static Manifold Union(
			Manifold a,
			Manifold b,
			BooleanEngine engine,
			WindingRule windingRule,
			IProgress<(string Phase, double? Fraction)>? progress,
			CancellationToken cancellationToken = default)
		{
			return Boolean(a, b, ManifoldOpType.Add, engine, windingRule, progress, cancellationToken);
		}

		/// <summary>
		/// <paramref name="a"/> minus <paramref name="b"/> on an explicitly chosen
		/// engine, under the default <see cref="WindingRule.Positive"/> rule.
		/// </summary>
		/// <inheritdoc cref="Boolean(Manifold, Manifold, ManifoldOpType, BooleanEngine, WindingRule, IProgress{ValueTuple{string, double?}}, CancelToken)"/>
		public static Manifold Subtract(Manifold a, Manifold b, BooleanEngine engine)
		{
			return BooleanCore(a, b, ManifoldOpType.Subtract, engine, WindingRule.Positive, null, null);
		}

		/// <summary>
		/// <paramref name="a"/> minus <paramref name="b"/>, with an explicit
		/// winding rule, progress reporting and cancellation.
		/// </summary>
		/// <inheritdoc cref="Boolean(Manifold, Manifold, ManifoldOpType, BooleanEngine, WindingRule, IProgress{ValueTuple{string, double?}}, CancellationToken)"/>
		public static Manifold Subtract(
			Manifold a,
			Manifold b,
			BooleanEngine engine,
			WindingRule windingRule,
			IProgress<(string Phase, double? Fraction)>? progress,
			CancellationToken cancellationToken = default)
		{
			return Boolean(a, b, ManifoldOpType.Subtract, engine, windingRule, progress, cancellationToken);
		}

		/// <summary>
		/// The intersection of two solids on an explicitly chosen engine, under the
		/// default <see cref="WindingRule.Positive"/> rule.
		/// </summary>
		/// <inheritdoc cref="Boolean(Manifold, Manifold, ManifoldOpType, BooleanEngine, WindingRule, IProgress{ValueTuple{string, double?}}, CancelToken)"/>
		public static Manifold Intersect(Manifold a, Manifold b, BooleanEngine engine)
		{
			return BooleanCore(a, b, ManifoldOpType.Intersect, engine, WindingRule.Positive, null, null);
		}

		/// <summary>
		/// The intersection of two solids, with an explicit winding rule, progress
		/// reporting and cancellation.
		/// </summary>
		/// <inheritdoc cref="Boolean(Manifold, Manifold, ManifoldOpType, BooleanEngine, WindingRule, IProgress{ValueTuple{string, double?}}, CancellationToken)"/>
		public static Manifold Intersect(
			Manifold a,
			Manifold b,
			BooleanEngine engine,
			WindingRule windingRule,
			IProgress<(string Phase, double? Fraction)>? progress,
			CancellationToken cancellationToken = default)
		{
			return Boolean(a, b, ManifoldOpType.Intersect, engine, windingRule, progress, cancellationToken);
		}

		/// <summary>
		/// The single implementation behind every binary boolean overload. A null
		/// <paramref name="cancelToken"/> becomes a NULL token natively ("
		/// uncancellable", which costs nothing) and a null
		/// <paramref name="progress"/> becomes a NULL callback, which the ABI
		/// documents as exactly the un-instrumented pipeline.
		/// </summary>
		private static unsafe Manifold BooleanCore(
			Manifold a,
			Manifold b,
			ManifoldOpType op,
			BooleanEngine engine,
			WindingRule windingRule,
			IProgress<(string Phase, double? Fraction)>? progress,
			CancelToken? cancelToken)
		{
			NativeVersionCheck.Verify();

			if (a is null)
			{
				throw new ArgumentNullException(nameof(a));
			}

			if (b is null)
			{
				throw new ArgumentNullException(nameof(b));
			}

			a.ThrowIfDisposed();
			b.ThrowIfDisposed();

			bool takenA = false;
			bool takenB = false;
			bool tokenTaken = false;

			// The bridge is rooted by the GCHandle for exactly the native call, so
			// the callback can fire from a worker thread without a GC hole, and the
			// handle is freed on every path out.
			ProgressBridge? bridge = null;
			GCHandle bridgeHandle = default;

			try
			{
				IntPtr tokenPointer = IntPtr.Zero;
				if (cancelToken is not null)
				{
					cancelToken.AddRef(ref tokenTaken);
					tokenPointer = cancelToken.Handle.DangerousGetHandle();
				}

				a.handle.DangerousAddRef(ref takenA);
				if (!takenA)
				{
					throw new ObjectDisposedException(nameof(Manifold));
				}

				b.handle.DangerousAddRef(ref takenB);
				if (!takenB)
				{
					throw new ObjectDisposedException(nameof(Manifold));
				}

				// Held as IntPtr rather than as the function-pointer type it really is,
				// because that is how the import declares it - see
				// NativeMethods.manifold_rs_boolean_progress_rule for the mono-wasm
				// reason. IntPtr.Zero is the same NULL the ABI reads as "no callback".
				IntPtr callback = IntPtr.Zero;
				void* user = null;
				if (progress is not null)
				{
					// Read the phase-name table before the kernel starts, so the first
					// callback - which may already be on a worker thread - does no
					// native work of its own.
					ProgressPhases.Preload();
					bridge = new ProgressBridge(progress);
					bridgeHandle = GCHandle.Alloc(bridge);
					callback = (IntPtr)ProgressBridge.Callback;
					user = (void*)GCHandle.ToIntPtr(bridgeHandle);
				}

				IntPtr result = NativeMethods.manifold_rs_boolean_progress_rule(
					a.handle.DangerousGetHandle(),
					b.handle.DangerousGetHandle(),
					(int)op,
					(int)engine,
					(int)windingRule,
					tokenPointer,
					callback,
					user);

				// Cancellation comes back as a valid handle with a Cancelled status,
				// not as NULL, so this really is only argument errors and panics.
				if (result == IntPtr.Zero)
				{
					throw new ManifoldException(
						$"manifold_rs_boolean_progress_rule ({op}, {engine}, {windingRule}) failed",
						null,
						NativeMethods.GetLastError());
				}

				Manifold combined = new Manifold(new ManifoldHandle(result));

				// A sink that threw is the caller's bug, but it has to surface as an
				// exception on their thread rather than as a silently ignored one - so
				// the result is dropped and the original exception rethrown with its
				// stack intact.
				if (bridge?.Failure is System.Runtime.ExceptionServices.ExceptionDispatchInfo failure)
				{
					combined.Dispose();
					failure.Throw();
				}

				return combined;
			}
			finally
			{
				if (bridgeHandle.IsAllocated)
				{
					bridgeHandle.Free();
				}

				if (takenA)
				{
					a!.handle.DangerousRelease();
				}

				if (takenB)
				{
					b!.handle.DangerousRelease();
				}

				if (tokenTaken)
				{
					cancelToken!.Handle.DangerousRelease();
				}
			}
		}
	}
}
