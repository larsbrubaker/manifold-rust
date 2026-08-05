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

// Managed-side tests for the robust (non-manifold) import and the boolean
// engine selector, mirroring ffi/src/tests_robust.rs.

using System;
using System.Collections.Generic;
using System.Threading.Tasks;
using TUnit.Assertions;
using TUnit.Assertions.Extensions;
using TUnit.Core;

namespace ManifoldRust.Tests
{
	public class RobustEngineTests
	{
		/// <summary>
		/// A mesh re-expressed as fully disconnected triangle soup (three
		/// duplicated verts per triangle) — geometrically identical, guaranteed
		/// to fail strict halfedge pairing.
		/// </summary>
		private static (float[] VertProperties, uint[] TriVerts) AsSoup(float[] verts, uint[] tris)
		{
			float[] soupVerts = new float[tris.Length * 3];
			for (int i = 0; i < tris.Length; i++)
			{
				int o = (int)tris[i] * 3;
				soupVerts[(i * 3) + 0] = verts[o + 0];
				soupVerts[(i * 3) + 1] = verts[o + 1];
				soupVerts[(i * 3) + 2] = verts[o + 2];
			}

			uint[] soupTris = new uint[tris.Length];
			for (uint i = 0; i < tris.Length; i++)
			{
				soupTris[i] = i;
			}

			return (soupVerts, soupTris);
		}

		[Test]
		public async Task RobustImportAcceptsSoupThatStrictRejects()
		{
			(float[] verts, uint[] tris) = TestMeshes.Cube(0, 0, 0, 2);
			(float[] soupVerts, uint[] soupTris) = AsSoup(verts, tris);

			using Manifold strict = Manifold.FromMesh(soupVerts, soupTris);
			await Assert.That(strict.Status).IsEqualTo(ManifoldStatus.NotManifold);

			using Manifold robust = Manifold.FromMeshRobust(soupVerts, soupTris);
			await Assert.That(robust.Status).IsEqualTo(ManifoldStatus.NoError);
		}

		[Test]
		public async Task RobustImportRejectsOpenMeshWithNotClosed()
		{
			(float[] verts, uint[] tris) = TestMeshes.Cube(0, 0, 0, 2);
			(float[] soupVerts, uint[] soupTris) = AsSoup(verts, tris);
			uint[] open = new uint[soupTris.Length - 3];
			Array.Copy(soupTris, open, open.Length);

			using Manifold m = Manifold.FromMeshRobust(soupVerts, open);
			await Assert.That(m.Status).IsEqualTo(ManifoldStatus.NotClosed);
		}

		[Test]
		public async Task EngineSelectorRoundTripsAndAutoHandlesSoup()
		{
			// The selector is process-global — keep every interaction inside
			// one test body and restore the default on all paths.
			await Assert.That(Manifold.DefaultBooleanEngine).IsEqualTo(BooleanEngine.Exact);
			try
			{
				Manifold.DefaultBooleanEngine = BooleanEngine.Auto;
				await Assert.That(Manifold.DefaultBooleanEngine).IsEqualTo(BooleanEngine.Auto);

				// Edge-kissing cubes as one soup + a manifold cutter: under
				// Auto the difference succeeds where Exact would refuse.
				(float[] v1, uint[] t1) = TestMeshes.Cube(0, 0, 0, 2);
				(float[] v2, uint[] t2) = TestMeshes.Cube(2, 2, 0, 2);
				(float[] s1, uint[] _) = AsSoup(v1, t1);
				(float[] s2, uint[] _) = AsSoup(v2, t2);
				float[] merged = new float[s1.Length + s2.Length];
				s1.CopyTo(merged, 0);
				s2.CopyTo(merged, s1.Length);
				uint[] mergedTris = new uint[merged.Length / 3];
				for (uint i = 0; i < mergedTris.Length; i++)
				{
					mergedTris[i] = i;
				}

				using Manifold soup = Manifold.FromMeshRobust(merged, mergedTris);
				await Assert.That(soup.Status).IsEqualTo(ManifoldStatus.NoError);

				using Manifold cutter = TestMeshes.CubeManifold(0.5f, 0.5f, 0.5f, 1.0f);
				using Manifold diff = Manifold.BatchBoolean(
					new List<Manifold> { soup, cutter },
					ManifoldOpType.Subtract);
				await Assert.That(diff.Status).IsEqualTo(ManifoldStatus.NoError);
				MeshGL mesh = diff.GetMeshGL();
				await Assert.That(mesh.TriVerts.Length).IsGreaterThan(0);
			}
			finally
			{
				Manifold.DefaultBooleanEngine = BooleanEngine.Exact;
			}
		}
	}
}
