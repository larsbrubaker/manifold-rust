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

// Tests for cancellation, at both levels: the CancelToken that maps one-to-one
// onto the ABI, and the CancellationToken overload that turns a cancelled result
// into the exception .NET callers expect. The interesting case - cancelling from
// a different thread while the kernel is running - is the whole reason the
// feature exists, so it is tested against real work rather than a pre-set flag.

using System;
using System.Diagnostics;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using TUnit.Assertions;
using TUnit.Assertions.Extensions;
using TUnit.Core;
using TUnit.Core.Exceptions;

namespace ManifoldRust.Tests
{
	public class CancellationTests
	{
		[Test]
		public async Task CancelTokenRoundTripsItsFlagAndToleratesRepeatedDispose()
		{
			CancelToken token = new CancelToken();

			await Assert.That(token.IsCancelled).IsFalse();
			token.Cancel();
			await Assert.That(token.IsCancelled).IsTrue();

			// Sticky: cancelling twice is a no-op, not an error, and there is no
			// un-cancel.
			token.Cancel();
			await Assert.That(token.IsCancelled).IsTrue();

			token.Dispose();
			token.Dispose();

			await Assert.That(() => token.Cancel()).Throws<ObjectDisposedException>();
			await Assert.That(() => { _ = token.IsCancelled; }).Throws<ObjectDisposedException>();
		}

		[Test]
		public async Task ACancelledCancelTokenYieldsAnEmptyResultRatherThanAnException()
		{
			using Manifold a = TestMeshes.OriginalCube(0, 0, 0, 1);
			using Manifold b = TestMeshes.OriginalCube(0.5f, 0, 0, 1);

			using CancelToken token = new CancelToken();
			token.Cancel();

			// The low-level overload reports cancellation the way the ABI does: a
			// valid, empty handle carrying the Cancelled status. Only the
			// CancellationToken overload turns that into an exception.
			using Manifold result = Manifold.BatchBoolean(new[] { a, b }, ManifoldOpType.Add, token);

			await Assert.That(result.Status).IsEqualTo(ManifoldStatus.Cancelled);
			await Assert.That(result.GetMeshGL().TriangleCount).IsEqualTo(0);
		}

		[Test]
		public async Task AnAlreadyCancelledCancellationTokenThrowsAndLeavesTheLibraryUsable()
		{
			using Manifold a = TestMeshes.OriginalCube(0, 0, 0, 1);
			using Manifold b = TestMeshes.OriginalCube(0.5f, 0, 0, 1);

			using CancellationTokenSource source = new CancellationTokenSource();
			source.Cancel();

			OperationCanceledException? caught = null;
			try
			{
				using Manifold unexpected = Manifold.BatchBoolean(new[] { a, b }, ManifoldOpType.Add, source.Token);
			}
			catch (OperationCanceledException exception)
			{
				caught = exception;
			}

			await Assert.That(caught).IsNotNull();
			await Assert.That(caught!.CancellationToken).IsEqualTo(source.Token);

			// The cancelled result handle is disposed on the way out. Nothing asserts
			// that directly - there is no managed hook for it - but a leak or a double
			// free would take the next operation down with it, so this is the check
			// that has teeth.
			using Manifold union = Manifold.BatchBoolean(new[] { a, b }, ManifoldOpType.Add);
			await Assert.That(union.Status).IsEqualTo(ManifoldStatus.NoError);
			await Assert.That(union.GetMeshGL().TriangleCount).IsGreaterThan(0);
		}

		[Test]
		public async Task ADisposedOperandStillReleasesTheTokenOnTheWayOut()
		{
			// The partial-failure path: the operand loop throws part way through, so
			// the finally block has to release the token reference it took before the
			// loop started. A leaked reference would not fail the throwing call - it
			// would wedge the token's SafeHandle - so the token is reused afterwards.
			using Manifold live = TestMeshes.OriginalCube(0, 0, 0, 1);
			Manifold dead = TestMeshes.CubeManifold(0.5f, 0, 0, 1);
			dead.Dispose();

			using CancelToken token = new CancelToken();

			await Assert.That(() => { _ = Manifold.BatchBoolean(new[] { live, dead }, ManifoldOpType.Add, token); })
				.Throws<ObjectDisposedException>();

			// The token survived the failed call: still usable, and still uncancelled,
			// because it was only ever borrowed.
			await Assert.That(token.IsCancelled).IsFalse();
			using Manifold union = Manifold.BatchBoolean(new[] { live, live }, ManifoldOpType.Add, token);
			await Assert.That(union.Status).IsEqualTo(ManifoldStatus.NoError);
		}

		[Test]
		public async Task CancellationTokenNoneMatchesTheUncancellableOverload()
		{
			using Manifold a = TestMeshes.OriginalCube(0, 0, 0, 1);
			using Manifold b = TestMeshes.OriginalCube(0.5f, 0, 0, 1);
			Manifold[] operands = { a, b };

			using Manifold plain = Manifold.BatchBoolean(operands, ManifoldOpType.Add);
			using Manifold withToken = Manifold.BatchBoolean(operands, ManifoldOpType.Add, CancellationToken.None);

			await Assert.That(withToken.Status).IsEqualTo(plain.Status);

			MeshGL fromPlain = plain.GetMeshGL();
			MeshGL fromToken = withToken.GetMeshGL();

			await Assert.That(fromToken.VertProperties.SequenceEqual(fromPlain.VertProperties)).IsTrue();
			await Assert.That(fromToken.TriVerts.SequenceEqual(fromPlain.TriVerts)).IsTrue();
			await Assert.That(fromToken.RunIndex.SequenceEqual(fromPlain.RunIndex)).IsTrue();
		}

		// NotInParallel because this is the one test in the assembly that measures
		// time. TUnit runs classes concurrently, so on a two-core runner a neighbour's
		// sphere union would inflate the cancelled measurement and break the ratio -
		// a scheduling artefact reported as a cancellation bug.
		[Test]
		[NotInParallel]
		public async Task CancellingFromAnotherThreadCutsASlowUnionShort()
		{
			using Manifold a = TestMeshes.OriginalSphere(0, 0, 0, 1);
			using Manifold b = TestMeshes.OriginalSphere(0.5f, 0, 0, 1);
			await Assert.That(a.Status).IsEqualTo(ManifoldStatus.NoError);
			await Assert.That(b.Status).IsEqualTo(ManifoldStatus.NoError);

			Manifold[] operands = { a, b };

			// Baseline measured in this test rather than assumed, so the assertion
			// below is a ratio and not an absolute millisecond count that a loaded
			// build agent would fail.
			Stopwatch uncancelledWatch = Stopwatch.StartNew();
			using (Manifold baseline = Manifold.BatchBoolean(operands, ManifoldOpType.Add))
			{
				uncancelledWatch.Stop();
				await Assert.That(baseline.Status).IsEqualTo(ManifoldStatus.NoError);
			}

			TimeSpan uncancelled = uncancelledWatch.Elapsed;
			if (uncancelled <= TimeSpan.FromMilliseconds(20))
			{
				// A precondition of the test, not a claim about the library: with no
				// measurable window there is nothing to cancel in the middle of. A
				// machine fast enough to hit this deserves a skip, not a red build -
				// the correctness of cancellation is covered by the pre-cancelled
				// tests above, which have no timing in them at all.
				throw new SkipTestException(
					$"the uncancelled union took only {uncancelled.TotalMilliseconds:F1} ms, too fast to cancel mid-flight");
			}

			using CancellationTokenSource source = new CancellationTokenSource();
			using ManualResetEventSlim started = new ManualResetEventSlim(false);
			Stopwatch cancelledWatch = new Stopwatch();

			Task<OperationCanceledException?> work = Task.Run(() =>
			{
				started.Set();
				cancelledWatch.Start();
				try
				{
					using Manifold unexpected = Manifold.BatchBoolean(operands, ManifoldOpType.Add, source.Token);
					return null;
				}
				catch (OperationCanceledException exception)
				{
					return exception;
				}
				finally
				{
					cancelledWatch.Stop();
				}
			});

			// The handshake pins the start so the cancel lands on work in flight
			// rather than racing the task scheduler, and the wait after it is a small
			// fraction of the measured baseline so the ratio stays stable on a
			// machine where both numbers inflate. There is no event to wait on here:
			// "the kernel is part way through" is not something the ABI signals.
			started.Wait();
			await Task.Delay(TimeSpan.FromMilliseconds(Math.Max(1.0, uncancelled.TotalMilliseconds / 16.0)));
			source.Cancel();

			OperationCanceledException? caught = await work;

			await Assert.That(caught).IsNotNull();
			await Assert.That(cancelledWatch.Elapsed * 2).IsLessThan(uncancelled);
		}
	}
}
