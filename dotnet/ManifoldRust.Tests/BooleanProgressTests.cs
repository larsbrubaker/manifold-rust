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

// Tests for the binary boolean overloads: progress reporting and the winding
// rule. They mirror ffi/src/tests_progress.rs and
// boolean_progress_rule_keeps_inside_out_geometry in ffi/src/tests_robust.rs, so
// the managed layer is held to the same observable behaviour as the ABI.

using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using TUnit.Assertions;
using TUnit.Assertions.Extensions;
using TUnit.Core;

namespace ManifoldRust.Tests
{
	public class BooleanProgressTests
	{
		/// <summary>
		/// Collects every report. Locked because the ABI promises only that the
		/// callback is not re-entered concurrently, not that it runs on the thread
		/// that started the boolean.
		/// </summary>
		private sealed class Collector : IProgress<(string Phase, double? Fraction)>
		{
			private readonly List<(string Phase, double? Fraction)> events = new List<(string, double?)>();

			public void Report((string Phase, double? Fraction) value)
			{
				lock (this.events)
				{
					this.events.Add(value);
				}
			}

			public List<(string Phase, double? Fraction)> Snapshot()
			{
				lock (this.events)
				{
					return new List<(string Phase, double? Fraction)>(this.events);
				}
			}
		}

		[Test]
		public async Task AProgressSinkSeesOrderedPhasesAndTheSameResultAsThePlainCall()
		{
			using Manifold a = TestMeshes.CubeManifold(0, 0, 0, 1);
			using Manifold b = TestMeshes.CubeManifold(0.5f, 0.25f, 0.25f, 1);

			Collector collector = new Collector();
			using Manifold watched = Manifold.Union(a, b, BooleanEngine.Robust, WindingRule.Positive, collector);
			using Manifold plain = Manifold.Union(a, b, BooleanEngine.Robust);

			await Assert.That(watched.Status).IsEqualTo(ManifoldStatus.NoError);

			List<(string Phase, double? Fraction)> events = collector.Snapshot();
			await Assert.That(events.Count).IsGreaterThan(0);

			// Every report names a phase and carries either a fraction in [0, 1] or
			// the indeterminate null.
			foreach ((string phase, double? fraction) in events)
			{
				await Assert.That(phase).IsNotNullOrEmpty();
				if (fraction is double value)
				{
					await Assert.That(value).IsGreaterThanOrEqualTo(0.0);
					await Assert.That(value).IsLessThanOrEqualTo(1.0);
				}
			}

			// Ordered: the kernel walks the phase table forwards, so a phase that has
			// been left is never reported again. Checking for contiguous runs is the
			// managed equivalent of the ABI test's monotonic phase ids.
			List<string> order = new List<string>();
			foreach ((string phase, double? _) in events)
			{
				if (order.Count == 0 || order[order.Count - 1] != phase)
				{
					await Assert.That(order.Contains(phase)).IsFalse();
					order.Add(phase);
				}
			}

			// Reporting must not change the geometry.
			MeshGL fromWatched = watched.GetMeshGL();
			MeshGL fromPlain = plain.GetMeshGL();
			await Assert.That(fromWatched.VertProperties.SequenceEqual(fromPlain.VertProperties)).IsTrue();
			await Assert.That(fromWatched.TriVerts.SequenceEqual(fromPlain.TriVerts)).IsTrue();
			await Assert.That(fromWatched.RunIndex.SequenceEqual(fromPlain.RunIndex)).IsTrue();
		}

		[Test]
		public async Task TheExactEngineReportsOneIndeterminatePhase()
		{
			using Manifold a = TestMeshes.CubeManifold(0, 0, 0, 1);
			using Manifold b = TestMeshes.CubeManifold(0.5f, 0.25f, 0.25f, 1);

			Collector collector = new Collector();
			using Manifold result = Manifold.Boolean(a, b, ManifoldOpType.Add, BooleanEngine.Exact, WindingRule.Positive, collector);

			await Assert.That(result.Status).IsEqualTo(ManifoldStatus.NoError);

			List<(string Phase, double? Fraction)> events = collector.Snapshot();
			await Assert.That(events.Count).IsEqualTo(1);
			await Assert.That(events[0].Phase).IsEqualTo("exact boolean");

			// Indeterminate: the exact pipeline has no work total to divide by, so it
			// reports null rather than a fraction of zero.
			await Assert.That(events[0].Fraction).IsNull();
		}

		[Test]
		public async Task AnExceptionFromTheSinkSurfacesOnTheCallingThread()
		{
			// The callback runs inside a native frame, where an escaping exception
			// would take the process down; it has to come out of the call instead.
			using Manifold a = TestMeshes.CubeManifold(0, 0, 0, 1);
			using Manifold b = TestMeshes.CubeManifold(0.5f, 0.25f, 0.25f, 1);

			InvalidOperationException thrown = new InvalidOperationException("sink failed");
			IProgress<(string, double?)> sink = new DelegateProgress(_ => throw thrown);

			await Assert.That(() => Manifold.Union(a, b, BooleanEngine.Robust, WindingRule.Positive, sink))
				.Throws<InvalidOperationException>();

			// The library is still usable afterwards: nothing was left reference
			// counted or half destroyed.
			using Manifold union = Manifold.Union(a, b, BooleanEngine.Robust);
			await Assert.That(union.Status).IsEqualTo(ManifoldStatus.NoError);
		}

		[Test]
		public async Task TheNonzeroRuleKeepsAnInsideOutOperand()
		{
			// An inside-out 2-cube at the origin unioned with a correctly wound
			// 2-cube offset by one along each axis. Under the positive rule only the
			// wound cube is material (volume 8); under the nonzero rule the inverted
			// one counts too (8 + 8 - 1 of overlap = 15), with its emitted
			// orientation corrected so the signed volume stays positive.
			(float[] verts, uint[] tris) = TestMeshes.Cube(0, 0, 0, 2);
			using Manifold inverted = Manifold.FromMesh(verts, TestMeshes.Inverted(tris));
			using Manifold wound = TestMeshes.CubeManifold(1, 1, 1, 2);
			await Assert.That(inverted.Status).IsEqualTo(ManifoldStatus.NoError);
			await Assert.That(wound.Status).IsEqualTo(ManifoldStatus.NoError);

			using Manifold positive = Manifold.Union(inverted, wound, BooleanEngine.Robust, WindingRule.Positive, null);
			using Manifold nonzero = Manifold.Union(inverted, wound, BooleanEngine.Robust, WindingRule.Nonzero, null);

			await Assert.That(TestMeshes.SignedVolume(positive.GetMeshGL())).IsEqualTo(8.0).Within(1e-6);
			await Assert.That(TestMeshes.SignedVolume(nonzero.GetMeshGL())).IsEqualTo(15.0).Within(1e-6);
		}

		[Test]
		public async Task AnUnknownWindingRuleIsAnArgumentErrorNotACrash()
		{
			using Manifold a = TestMeshes.CubeManifold(0, 0, 0, 1);
			using Manifold b = TestMeshes.CubeManifold(0.5f, 0, 0, 1);

			await Assert.That(() => Manifold.Boolean(a, b, ManifoldOpType.Add, BooleanEngine.Robust, (WindingRule)7))
				.Throws<ManifoldException>();
		}

		[Test]
		public async Task ACancelledBinaryBooleanReportsCancelledAndThrowsForACancellationToken()
		{
			using Manifold a = TestMeshes.CubeManifold(0, 0, 0, 1);
			using Manifold b = TestMeshes.CubeManifold(0.5f, 0, 0, 1);

			using CancelToken token = new CancelToken();
			token.Cancel();

			// The CancelToken overload reports cancellation the way the ABI does.
			using (Manifold cancelled = Manifold.Boolean(a, b, ManifoldOpType.Add, BooleanEngine.Robust, WindingRule.Positive, null, token))
			{
				await Assert.That(cancelled.Status).IsEqualTo(ManifoldStatus.Cancelled);
			}

			using CancellationTokenSource source = new CancellationTokenSource();
			source.Cancel();

			await Assert.That(() => Manifold.Union(a, b, BooleanEngine.Robust, WindingRule.Positive, null, source.Token))
				.Throws<OperationCanceledException>();
		}

		[Test]
		public async Task NullOperandsAreRejected()
		{
			using Manifold a = TestMeshes.CubeManifold(0, 0, 0, 1);

			await Assert.That(() => Manifold.Union(null!, a, BooleanEngine.Exact)).Throws<ArgumentNullException>();
			await Assert.That(() => Manifold.Union(a, null!, BooleanEngine.Exact)).Throws<ArgumentNullException>();
		}

		/// <summary>
		/// An <see cref="IProgress{T}"/> over a lambda.
		/// <see cref="System.Progress{T}"/> is unusable here: it posts to the
		/// synchronization context that created it, which would swallow the
		/// exception this test is about instead of letting it reach the binding.
		/// </summary>
		private sealed class DelegateProgress : IProgress<(string, double?)>
		{
			private readonly Action<(string, double?)> action;

			internal DelegateProgress(Action<(string, double?)> action)
			{
				this.action = action;
			}

			public void Report((string, double?) value)
			{
				this.action(value);
			}
		}
	}
}
