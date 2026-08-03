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

// Tests for the managed binding. Every one of them goes through the public API
// and the real manifold_rs shared library - there is no fake, because the whole
// point of this assembly is the marshalling between the two.

using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using TUnit.Assertions;
using TUnit.Assertions.Extensions;
using TUnit.Core;

namespace ManifoldRust.Tests
{
	public class ManifoldTests
	{
		[Test]
		public async Task ImportedCubeBecomesAnOriginalWithNoError()
		{
			using Manifold cube = TestMeshes.CubeManifold(0, 0, 0, 1);
			using Manifold original = cube.AsOriginal();

			await Assert.That(original.Status).IsEqualTo(ManifoldStatus.NoError);
			await Assert.That(original.OriginalId).IsGreaterThanOrEqualTo(0);

			MeshGL mesh = original.GetMeshGL();
			await Assert.That(mesh.NumProp).IsEqualTo(3u);
			await Assert.That(mesh.TriVerts.Length).IsEqualTo(36);
			await Assert.That(mesh.FaceId.Length).IsEqualTo(12);
			await Assert.That(mesh.VertProperties.Length % 3).IsEqualTo(0);
		}

		[Test]
		public async Task UnionOfTwoCubesKeepsBothOriginalIds()
		{
			using Manifold a = TestMeshes.OriginalCube(0, 0, 0, 1);
			using Manifold b = TestMeshes.OriginalCube(0.5f, 0, 0, 1);

			uint idA = (uint)a.OriginalId;
			uint idB = (uint)b.OriginalId;

			using Manifold union = Manifold.BatchBoolean(new[] { a, b }, ManifoldOpType.Add);
			await Assert.That(union.Status).IsEqualTo(ManifoldStatus.NoError);

			MeshGL mesh = union.GetMeshGL();
			await Assert.That(mesh.TriVerts.Length).IsGreaterThan(0);
			await Assert.That(mesh.TriVerts.Length % 3).IsEqualTo(0);

			// One run per source mesh, plus the trailing end sentinel in run_index.
			await Assert.That(mesh.RunOriginalId.Length).IsEqualTo(2);
			await Assert.That(mesh.RunOriginalId).Contains(idA);
			await Assert.That(mesh.RunOriginalId).Contains(idB);
			await Assert.That(mesh.RunIndex.Length).IsEqualTo(3);
			await Assert.That(mesh.RunIndex[^1]).IsEqualTo((uint)mesh.TriVerts.Length);
		}

		[Test]
		public async Task SubtractAppliesEveryCutterToTheFirstOperand()
		{
			using Manifold body = TestMeshes.OriginalCube(0, 0, 0, 1);
			using Manifold cutterA = TestMeshes.OriginalCube(-0.1f, -0.1f, -0.1f, 0.3f);
			using Manifold cutterB = TestMeshes.OriginalCube(0.8f, 0.8f, 0.8f, 0.3f);

			using Manifold carved = Manifold.BatchBoolean(
				new[] { body, cutterA, cutterB },
				ManifoldOpType.Subtract);

			await Assert.That(carved.Status).IsEqualTo(ManifoldStatus.NoError);

			MeshGL mesh = carved.GetMeshGL();
			await Assert.That(mesh.TriVerts.Length).IsGreaterThan(0);

			// Both corners really were removed: a plain cube is 12 triangles, and
			// cutting a corner off can only add more.
			await Assert.That(mesh.TriangleCount).IsGreaterThan(12);
		}

		[Test]
		public async Task IntersectOfOverlappingCubesIsNonEmpty()
		{
			using Manifold a = TestMeshes.CubeManifold(0, 0, 0, 1);
			using Manifold b = TestMeshes.CubeManifold(0.5f, 0, 0, 1);

			using Manifold overlap = Manifold.BatchBoolean(new[] { a, b }, ManifoldOpType.Intersect);

			await Assert.That(overlap.Status).IsEqualTo(ManifoldStatus.NoError);

			MeshGL mesh = overlap.GetMeshGL();
			await Assert.That(mesh.TriVerts.Length).IsGreaterThan(0);
			await Assert.That(mesh.TriVerts.Length % 3).IsEqualTo(0);
		}

		[Test]
		public async Task NonManifoldInputReportsStatusAndExportsEmpty()
		{
			(float[] verts, uint[] tris) = TestMeshes.NonManifoldTetrahedron();

			using Manifold bad = Manifold.FromMesh(verts, tris);

			// A validation failure still yields a usable handle, the way the C++ API
			// returns an error-status manifold.
			await Assert.That(bad.Status).IsEqualTo(ManifoldStatus.NotManifold);

			// The C++ binding crashes when exporting an error-status manifold. This
			// one must simply hand back empty arrays.
			MeshGL mesh = bad.GetMeshGL();
			await Assert.That(mesh.TriVerts.Length).IsEqualTo(0);
			await Assert.That(mesh.VertProperties.Length).IsEqualTo(0);
		}

		[Test]
		public async Task ABadOperandIsAbsorbedAsEmptyGeometryWithoutAnyError()
		{
			// This is the documented footgun, pinned here because it is the one a
			// caller cannot see coming: a part that failed to import does not make the
			// boolean fail, it just quietly contributes nothing. Anything that feeds
			// user geometry into BatchBoolean has to check Status itself.
			(float[] badVerts, uint[] badTris) = TestMeshes.NonManifoldTetrahedron();

			using Manifold good = TestMeshes.OriginalCube(0, 0, 0, 1);
			using Manifold bad = Manifold.FromMesh(badVerts, badTris);
			await Assert.That(bad.Status).IsEqualTo(ManifoldStatus.NotManifold);

			using Manifold union = Manifold.BatchBoolean(new[] { good, bad }, ManifoldOpType.Add);

			// No error is reported anywhere - the failed operand has simply vanished.
			await Assert.That(union.Status).IsEqualTo(ManifoldStatus.NoError);

			MeshGL fromUnion = union.GetMeshGL();
			MeshGL fromGoodAlone = good.GetMeshGL();
			await Assert.That(fromUnion.TriVerts.SequenceEqual(fromGoodAlone.TriVerts)).IsTrue();
			await Assert.That(fromUnion.VertProperties.SequenceEqual(fromGoodAlone.VertProperties)).IsTrue();
		}

		[Test]
		public async Task ArgumentsThatCannotDescribeAMeshThrowWithTheNativeMessage()
		{
			(float[] verts, uint[] tris) = TestMeshes.Cube(0, 0, 0, 1);

			// num_prop below the three position slots. The native side returns NULL
			// and records why, which is the only path that exercises the last-error
			// read: its two-call sizing protocol and the UTF-8 decode of the result.
			ManifoldException tooFewProperties = await ThrowsManifoldException(
				() => Manifold.FromMesh(verts, tris, numProp: 2));

			await Assert.That(tooFewProperties.NativeError).IsNotNull();
			await Assert.That(tooFewProperties.NativeError!).Contains("num_prop 2 < 3");
			await Assert.That(tooFewProperties.Message).Contains("manifold_rs_from_mesh");
		}

		[Test]
		public async Task VertexArrayThatIsNotAWholeNumberOfVerticesThrowsWithTheNativeMessage()
		{
			(float[] verts, uint[] tris) = TestMeshes.Cube(0, 0, 0, 1);

			// One float short of a whole vertex, so the length check rejects it.
			float[] truncated = verts[..^1];

			ManifoldException ragged = await ThrowsManifoldException(
				() => Manifold.FromMesh(truncated, tris));

			await Assert.That(ragged.NativeError).IsNotNull();
			await Assert.That(ragged.NativeError!).Contains($"vert_properties_len {truncated.Length} is not a multiple of num_prop 3");
		}

		[Test]
		public async Task EmptyInputProducesAnEmptyManifoldRatherThanAFailure()
		{
			// Empty spans give the native side a null pointer with length 0, which it
			// routes around from_raw_parts rather than dereferencing.
			using Manifold empty = Manifold.FromMesh(ReadOnlySpan<float>.Empty, ReadOnlySpan<uint>.Empty);

			await Assert.That(empty.Status).IsEqualTo(ManifoldStatus.NoError);

			MeshGL mesh = empty.GetMeshGL();
			await Assert.That(mesh.TriVerts.Length).IsEqualTo(0);
			await Assert.That(mesh.VertProperties.Length).IsEqualTo(0);
		}

		[Test]
		public async Task BatchBooleanRejectsAnEmptyOperandList()
		{
			await Assert.That(() => { _ = Manifold.BatchBoolean(Array.Empty<Manifold>(), ManifoldOpType.Add); })
				.Throws<ArgumentException>();
		}

		[Test]
		public async Task BatchBooleanRejectsADisposedOperand()
		{
			using Manifold live = TestMeshes.CubeManifold(0, 0, 0, 1);
			Manifold dead = TestMeshes.CubeManifold(0.5f, 0, 0, 1);
			dead.Dispose();

			List<Manifold> operands = new() { live, dead };

			await Assert.That(() => { _ = Manifold.BatchBoolean(operands, ManifoldOpType.Add); })
				.Throws<ObjectDisposedException>();
		}

		[Test]
		public async Task DisposeIsIdempotentAndExportAfterDisposeThrows()
		{
			Manifold cube = TestMeshes.CubeManifold(0, 0, 0, 1);
			cube.Dispose();
			cube.Dispose();

			await Assert.That(() => { _ = cube.GetMeshGL(); }).Throws<ObjectDisposedException>();
			await Assert.That(() => { _ = cube.Status; }).Throws<ObjectDisposedException>();
		}

		[Test]
		public async Task RepeatedUnionsProduceIdenticalGeometry()
		{
			MeshGL first = UnionOfTwoCubes();
			MeshGL second = UnionOfTwoCubes();

			// The kernel is deterministic - including the parallel feature the FFI
			// crate enables - so two identical inputs must give bit-identical output.
			// SequenceEqual rather than a set comparison: order is part of the claim.
			await Assert.That(second.VertProperties.SequenceEqual(first.VertProperties)).IsTrue();
			await Assert.That(second.TriVerts.SequenceEqual(first.TriVerts)).IsTrue();
		}

		/// <summary>
		/// Runs <paramref name="action"/> and hands back the <see cref="ManifoldException"/>
		/// it threw, so a test can assert on the native message it carries. Failing
		/// the assertion here rather than letting the exception escape keeps the
		/// reason for the failure readable.
		/// </summary>
		private static async Task<ManifoldException> ThrowsManifoldException(Func<Manifold> action)
		{
			ManifoldException? caught = null;

			try
			{
				using Manifold unexpected = action();
			}
			catch (ManifoldException exception)
			{
				caught = exception;
			}

			await Assert.That(caught).IsNotNull();
			return caught!;
		}

		private static MeshGL UnionOfTwoCubes()
		{
			using Manifold a = TestMeshes.OriginalCube(0, 0, 0, 1);
			using Manifold b = TestMeshes.OriginalCube(0.5f, 0, 0, 1);
			using Manifold union = Manifold.BatchBoolean(new[] { a, b }, ManifoldOpType.Add);
			return union.GetMeshGL();
		}
	}
}
