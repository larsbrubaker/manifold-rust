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

// Hand-built triangle soups for the binding tests, transliterated from
// ffi/src/tests.rs so the managed tests feed the FFI exactly what the Rust
// tests do. Written out by hand on purpose: the tests must not depend on a
// mesh generator that lives on the other side of the boundary they are testing.

namespace ManifoldRust.Tests
{
	internal static class TestMeshes
	{
		/// <summary>
		/// Axis-aligned cube as interleaved positions plus counter-clockwise-outward
		/// triangles.
		/// </summary>
		internal static (float[] VertProperties, uint[] TriVerts) Cube(float x, float y, float z, float size)
		{
			float s = size;

			float[] verts =
			{
				x, y, z,
				x + s, y, z,
				x + s, y + s, z,
				x, y + s, z,
				x, y, z + s,
				x + s, y, z + s,
				x + s, y + s, z + s,
				x, y + s, z + s,
			};

			uint[] tris =
			{
				0, 2, 1, 0, 3, 2, // -Z
				4, 5, 6, 4, 6, 7, // +Z
				0, 1, 5, 0, 5, 4, // -Y
				1, 2, 6, 1, 6, 5, // +X
				2, 3, 7, 2, 7, 6, // +Y
				3, 0, 4, 3, 4, 7, // -X
			};

			return (verts, tris);
		}

		/// <summary>
		/// Imports a cube through the public API. The result is not re-tagged as an
		/// original; call <see cref="Manifold.AsOriginal"/> when the test cares about
		/// run IDs.
		/// </summary>
		internal static Manifold CubeManifold(float x, float y, float z, float size)
		{
			(float[] verts, uint[] tris) = Cube(x, y, z, size);
			return Manifold.FromMesh(verts, tris);
		}

		/// <summary>
		/// A cube re-tagged as an original mesh, so booleans made from it report its
		/// ID in their runs. The intermediate import is disposed here rather than
		/// left to the finalizer.
		/// </summary>
		internal static Manifold OriginalCube(float x, float y, float z, float size)
		{
			using Manifold imported = CubeManifold(x, y, z, size);
			return imported.AsOriginal();
		}

		/// <summary>
		/// A tetrahedron with one face wound the wrong way: four verts and four tris,
		/// so it clears the size checks, but two half-edges then run in the same
		/// direction and the topology check rejects it as not manifold.
		/// </summary>
		internal static (float[] VertProperties, uint[] TriVerts) NonManifoldTetrahedron()
		{
			float[] verts =
			{
				0.0f, 0.0f, 0.0f,
				1.0f, 0.0f, 0.0f,
				0.0f, 1.0f, 0.0f,
				0.0f, 0.0f, 1.0f,
			};

			uint[] tris = { 0, 1, 2, 0, 1, 3, 0, 3, 2, 1, 2, 3 };

			return (verts, tris);
		}
	}
}
