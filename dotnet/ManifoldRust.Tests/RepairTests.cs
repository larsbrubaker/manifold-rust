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

// Managed-side tests for Manifold.RepairOrientation and
// Manifold.HasSelfIntersections, mirroring repair_orientation_rewinds_inverted_cube
// and has_self_intersections_flags_doubled_surface in ffi/src/tests_robust.rs so
// the binding is checked against the same geometry as the ABI.

using System;
using System.Threading.Tasks;
using TUnit.Assertions;
using TUnit.Assertions.Extensions;
using TUnit.Core;

namespace ManifoldRust.Tests
{
	public class RepairTests
	{
		[Test]
		public async Task RepairOrientationRewindsAnInvertedCube()
		{
			// A cube with every triangle reversed: a perfectly valid manifold that is
			// wound inside out, so it encloses no {w >= 1} material at all.
			(float[] verts, uint[] tris) = TestMeshes.Cube(0, 0, 0, 2);
			using Manifold broken = Manifold.FromMesh(verts, TestMeshes.Inverted(tris));
			await Assert.That(broken.Status).IsEqualTo(ManifoldStatus.NoError);
			await Assert.That(TestMeshes.SignedVolume(broken.GetMeshGL())).IsLessThan(0.0);

			using Manifold repaired = broken.RepairOrientation();
			await Assert.That(repaired.Status).IsEqualTo(ManifoldStatus.NoError);
			await Assert.That(TestMeshes.SignedVolume(repaired.GetMeshGL())).IsEqualTo(8.0).Within(1e-6);
		}

		[Test]
		public async Task RepairOrientationLeavesAWellWoundCubeAlone()
		{
			// The no-op path is the one callers will hit most: repairing has to be
			// safe to run unconditionally, geometry unchanged.
			using Manifold cube = TestMeshes.CubeManifold(0, 0, 0, 2);
			using Manifold repaired = cube.RepairOrientation();

			await Assert.That(repaired.Status).IsEqualTo(ManifoldStatus.NoError);
			await Assert.That(TestMeshes.SignedVolume(repaired.GetMeshGL())).IsEqualTo(8.0).Within(1e-6);
		}

		[Test]
		public async Task RepairOrientationOnADisposedManifoldThrows()
		{
			Manifold cube = TestMeshes.CubeManifold(0, 0, 0, 2);
			cube.Dispose();

			await Assert.That(() => cube.RepairOrientation()).Throws<ObjectDisposedException>();
			await Assert.That(() => { _ = cube.HasSelfIntersections; }).Throws<ObjectDisposedException>();
		}

		[Test]
		public async Task ACleanCubeHasNoSelfIntersections()
		{
			// Shared edges and vertices are not intersections; a closed mesh is made
			// of them.
			using Manifold cube = TestMeshes.CubeManifold(0, 0, 0, 2);
			await Assert.That(cube.HasSelfIntersections).IsFalse();

			// The verdict is cached natively, so reading twice must agree.
			await Assert.That(cube.HasSelfIntersections).IsFalse();
		}

		[Test]
		public async Task ADoubledSheetHasSelfIntersections()
		{
			// The same cube with every triangle duplicated in reverse winding: a
			// coincident (doubled) sheet. Pairing fails, so it comes in as soup.
			(float[] verts, uint[] tris) = TestMeshes.Cube(0, 0, 0, 2);
			(float[] soupVerts, uint[] soupTris) = TestMeshes.AsSoup(verts, tris);

			float[] doubledVerts = new float[soupVerts.Length * 2];
			soupVerts.CopyTo(doubledVerts, 0);
			for (int t = 0; t < soupTris.Length; t += 3)
			{
				// Reverse winding: the second copy of each triangle faces the other way.
				uint[] order = { soupTris[t], soupTris[t + 2], soupTris[t + 1] };
				for (int k = 0; k < 3; k++)
				{
					int source = (int)order[k] * 3;
					int target = soupVerts.Length + ((t + k) * 3);
					doubledVerts[target + 0] = soupVerts[source + 0];
					doubledVerts[target + 1] = soupVerts[source + 1];
					doubledVerts[target + 2] = soupVerts[source + 2];
				}
			}

			uint[] doubledTris = new uint[soupTris.Length * 2];
			for (uint i = 0; i < doubledTris.Length; i++)
			{
				doubledTris[i] = i;
			}

			using Manifold doubled = Manifold.FromMeshRobust(doubledVerts, doubledTris);
			await Assert.That(doubled.Status).IsEqualTo(ManifoldStatus.NoError);
			await Assert.That(doubled.HasSelfIntersections).IsTrue();
		}
	}
}
