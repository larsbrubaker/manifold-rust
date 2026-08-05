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

using System;

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
		/// The same cube as <see cref="Cube"/> in double precision, with 64-bit
		/// indices. The values are identical, so a test can compare the two import
		/// paths on geometry that is the same by construction rather than by
		/// coincidence.
		/// </summary>
		internal static (double[] VertProperties, ulong[] TriVerts) Cube64(double x, double y, double z, double size)
		{
			(float[] verts, uint[] tris) = Cube((float)x, (float)y, (float)z, (float)size);

			double[] wideVerts = new double[verts.Length];
			for (int i = 0; i < verts.Length; i++)
			{
				wideVerts[i] = verts[i];
			}

			ulong[] wideTris = new ulong[tris.Length];
			for (int i = 0; i < tris.Length; i++)
			{
				wideTris[i] = tris[i];
			}

			return (wideVerts, wideTris);
		}

		/// <summary>
		/// Imports a cube through the double-precision entry point, re-tagged as an
		/// original mesh so booleans made from it report its ID in their runs.
		/// </summary>
		internal static Manifold OriginalCube64(double x, double y, double z, double size)
		{
			(double[] verts, ulong[] tris) = Cube64(x, y, z, size);
			using Manifold imported = Manifold.FromMesh64(verts, tris);
			return imported.AsOriginal();
		}

		/// <summary>
		/// A closed UV sphere: one vertex per pole and <paramref name="slices"/> per
		/// intermediate ring, with the seam column sharing the ring's first vertex, so
		/// the triangle soup is manifold by index rather than by coincident positions.
		/// </summary>
		/// <remarks>
		/// Written here rather than taken from the kernel's own sphere primitive for
		/// the same reason as every other mesh in this file: the tests must not depend
		/// on the side of the boundary they are testing. It exists so the cancellation
		/// test has geometry whose union takes long enough to be interrupted.
		/// </remarks>
		internal static (float[] VertProperties, uint[] TriVerts) UvSphere(float cx, float cy, float cz, float radius, int slices, int stacks)
		{
			// Below these there is no closed surface to build: fewer than three slices
			// gives degenerate fans, and fewer than two stacks leaves no ring between
			// the poles. Caught here rather than as an IndexOutOfRangeException from
			// the middle of a fill loop, where the cause would not be obvious.
			if (slices <= 2)
			{
				throw new ArgumentOutOfRangeException(nameof(slices), slices, "A sphere needs at least 3 slices.");
			}

			if (stacks <= 1)
			{
				throw new ArgumentOutOfRangeException(nameof(stacks), stacks, "A sphere needs at least 2 stacks.");
			}

			int ringCount = stacks - 1;
			float[] verts = new float[(2 + (ringCount * slices)) * 3];

			int next = 0;
			void Emit(float x, float y, float z)
			{
				verts[next++] = cx + x;
				verts[next++] = cy + y;
				verts[next++] = cz + z;
			}

			const uint NorthPole = 0;
			Emit(0, 0, radius);

			for (int ring = 0; ring < ringCount; ring++)
			{
				double theta = Math.PI * (ring + 1) / stacks;
				double sinTheta = Math.Sin(theta);
				double cosTheta = Math.Cos(theta);

				for (int slice = 0; slice < slices; slice++)
				{
					double phi = 2.0 * Math.PI * slice / slices;
					Emit(
						(float)(radius * sinTheta * Math.Cos(phi)),
						(float)(radius * sinTheta * Math.Sin(phi)),
						(float)(radius * cosTheta));
				}
			}

			uint southPole = (uint)(1 + (ringCount * slices));
			Emit(0, 0, -radius);

			uint[] tris = new uint[((slices * 2) + ((ringCount - 1) * slices * 2)) * 3];
			int triNext = 0;
			void EmitTriangle(uint a, uint b, uint c)
			{
				tris[triNext++] = a;
				tris[triNext++] = b;
				tris[triNext++] = c;
			}

			// Ring vertices start at index 1; ring 0 is the one nearest the north pole.
			uint Ring(int ring, int slice) => (uint)(1 + (ring * slices) + (slice % slices));

			for (int slice = 0; slice < slices; slice++)
			{
				EmitTriangle(NorthPole, Ring(0, slice), Ring(0, slice + 1));
			}

			for (int ring = 0; ring < ringCount - 1; ring++)
			{
				for (int slice = 0; slice < slices; slice++)
				{
					uint upper = Ring(ring, slice);
					uint upperNext = Ring(ring, slice + 1);
					uint lower = Ring(ring + 1, slice);
					uint lowerNext = Ring(ring + 1, slice + 1);

					EmitTriangle(upper, lower, lowerNext);
					EmitTriangle(upper, lowerNext, upperNext);
				}
			}

			for (int slice = 0; slice < slices; slice++)
			{
				// Reversed relative to the north fan, so this cap also faces outward.
				EmitTriangle(southPole, Ring(ringCount - 1, slice + 1), Ring(ringCount - 1, slice));
			}

			return (verts, tris);
		}

		/// <summary>
		/// A sphere imported and re-tagged as an original mesh. Big enough by default
		/// that a union of two of them is slow enough to cancel mid-flight.
		/// </summary>
		internal static Manifold OriginalSphere(float cx, float cy, float cz, float radius, int slices = 400, int stacks = 200)
		{
			(float[] verts, uint[] tris) = UvSphere(cx, cy, cz, radius, slices, stacks);
			using Manifold imported = Manifold.FromMesh(verts, tris);
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
