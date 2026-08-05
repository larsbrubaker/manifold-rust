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

// Tests for the double-precision entry points. The claim being pinned is not
// "f64 is more accurate" - the kernel narrows internally, so it is not - but
// that the wide path is a faithful alternative spelling of the narrow one: the
// same geometry in, the same geometry out, with provenance intact.

using System;
using System.Linq;
using System.Threading.Tasks;
using TUnit.Assertions;
using TUnit.Assertions.Extensions;
using TUnit.Core;

namespace ManifoldRust.Tests
{
	public class Mesh64Tests
	{
		[Test]
		public async Task Mesh64CubeRoundTripsThroughAUnionWithBothOriginalIds()
		{
			using Manifold a = TestMeshes.OriginalCube64(0, 0, 0, 1);
			using Manifold b = TestMeshes.OriginalCube64(0.5, 0, 0, 1);

			await Assert.That(a.Status).IsEqualTo(ManifoldStatus.NoError);
			await Assert.That(b.Status).IsEqualTo(ManifoldStatus.NoError);

			uint idA = (uint)a.OriginalId;
			uint idB = (uint)b.OriginalId;

			using Manifold union = Manifold.BatchBoolean(new[] { a, b }, ManifoldOpType.Add);
			await Assert.That(union.Status).IsEqualTo(ManifoldStatus.NoError);

			MeshGL64 mesh = union.GetMeshGL64();
			await Assert.That(mesh.NumProp).IsEqualTo(3u);
			await Assert.That(mesh.TriVerts.Length % 3).IsEqualTo(0);
			await Assert.That(mesh.TriangleCount).IsGreaterThan(0);

			// One run per source mesh, plus the trailing end sentinel in run_index.
			await Assert.That(mesh.RunOriginalId.Length).IsEqualTo(2);
			await Assert.That(mesh.RunOriginalId).Contains(idA);
			await Assert.That(mesh.RunOriginalId).Contains(idB);
			await Assert.That(mesh.RunIndex.Length).IsEqualTo(3);
			await Assert.That(mesh.RunIndex[^1]).IsEqualTo((uint)mesh.TriVerts.Length);
			await Assert.That(mesh.FaceId.Length).IsEqualTo(mesh.TriangleCount);
		}

		[Test]
		public async Task TheWidePathAgreesWithTheNarrowOneOnTheSameGeometry()
		{
			using Manifold wideA = TestMeshes.OriginalCube64(0, 0, 0, 1);
			using Manifold wideB = TestMeshes.OriginalCube64(0.5, 0, 0, 1);
			using Manifold wideUnion = Manifold.BatchBoolean(new[] { wideA, wideB }, ManifoldOpType.Add);

			using Manifold narrowA = TestMeshes.OriginalCube(0, 0, 0, 1);
			using Manifold narrowB = TestMeshes.OriginalCube(0.5f, 0, 0, 1);
			using Manifold narrowUnion = Manifold.BatchBoolean(new[] { narrowA, narrowB }, ManifoldOpType.Add);

			MeshGL64 wide = wideUnion.GetMeshGL64();
			MeshGL narrow = narrowUnion.GetMeshGL();

			// Same topology, because the kernel saw the same numbers: the doubles are
			// exactly representable as floats and the indices are the same values.
			await Assert.That(wide.TriangleCount).IsEqualTo(narrow.TriangleCount);
			await Assert.That(wide.TriVerts.SequenceEqual(narrow.TriVerts)).IsTrue();
			await Assert.That(wide.VertexCount).IsEqualTo(narrow.VertexCount);

			// And the coordinates match once the float side is widened - which is the
			// honest statement of what the f64 path buys today: the same numbers in a
			// wider container, since the kernel computed them in single precision.
			await Assert.That(wide.VertProperties.SequenceEqual(narrow.VertProperties.Select(v => (double)v))).IsTrue();
		}

		[Test]
		public async Task WideningOverloadIsIdenticalToPassingUlongIndices()
		{
			(double[] verts, ulong[] wideTris) = TestMeshes.Cube64(0, 0, 0, 1);
			(_, uint[] narrowTris) = TestMeshes.Cube(0, 0, 0, 1);

			using Manifold fromWide = Manifold.FromMesh64(verts, wideTris);
			using Manifold fromNarrow = Manifold.FromMesh64(verts, narrowTris);

			await Assert.That(fromNarrow.Status).IsEqualTo(fromWide.Status);

			MeshGL64 wide = fromWide.GetMeshGL64();
			MeshGL64 narrow = fromNarrow.GetMeshGL64();

			await Assert.That(narrow.NumProp).IsEqualTo(wide.NumProp);
			await Assert.That(narrow.VertProperties.SequenceEqual(wide.VertProperties)).IsTrue();
			await Assert.That(narrow.TriVerts.SequenceEqual(wide.TriVerts)).IsTrue();
			await Assert.That(narrow.RunIndex.SequenceEqual(wide.RunIndex)).IsTrue();
			await Assert.That(narrow.FaceId.SequenceEqual(wide.FaceId)).IsTrue();
		}

		[Test]
		public async Task AnIndexTooLargeForTheKernelIsRejectedRatherThanWrapped()
		{
			// The one case the extra width really cannot carry: the kernel indexes
			// with 32 bits, and narrowing wraps. Silently building geometry from
			// index 0 instead of 2^32 would report NoError, so the FFI refuses.
			(double[] verts, _) = TestMeshes.Cube64(0, 0, 0, 1);
			ulong[] tris = { 0, 1, (ulong)uint.MaxValue + 1 };

			ManifoldException? caught = null;
			try
			{
				using Manifold unexpected = Manifold.FromMesh64(verts, tris);
			}
			catch (ManifoldException exception)
			{
				caught = exception;
			}

			await Assert.That(caught).IsNotNull();
			await Assert.That(caught!.NativeError).IsNotNull();
			await Assert.That(caught.NativeError!).Contains("exceeds u32::MAX");
		}

		[Test]
		public async Task NarrowingAnExportedArrayRefusesToWrap()
		{
			// The export side has the same hazard, but reaching it would need a mesh
			// with more than four billion vertices, so the guard is tested where it
			// lives instead of through the public API.
			ulong[] fits = { 0, 1, uint.MaxValue };
			await Assert.That(NativeArrays.NarrowChecked(fits, "tri_verts").SequenceEqual(new uint[] { 0, 1, uint.MaxValue })).IsTrue();

			ulong[] overflows = { 0, 1, (ulong)uint.MaxValue + 1 };
			ManifoldException? caught = null;
			try
			{
				_ = NativeArrays.NarrowChecked(overflows, "tri_verts");
			}
			catch (ManifoldException exception)
			{
				caught = exception;
			}

			await Assert.That(caught).IsNotNull();
			await Assert.That(caught!.Message).Contains("tri_verts[2]");

			await Assert.That(NativeArrays.NarrowChecked(Array.Empty<ulong>(), "tri_verts").Length).IsEqualTo(0);
		}
	}
}
