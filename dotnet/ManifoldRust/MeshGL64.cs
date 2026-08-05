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

// The double-precision counterpart of MeshGL, with the same no-native-lifetime
// design: the native MeshGl64Rs is created, drained and destroyed inside
// Manifold.GetMeshGL64.
//
// Only the coordinates stay wide. The index-shaped fields the ABI hands back as
// uint64_t come out of here as uint[], because the kernel that produced them
// indexes with 32 bits and every consumer of a triangle list wants uint anyway.
// The narrowing is checked rather than cast - see NativeArrays.NarrowChecked.

using System;

namespace ManifoldRust
{
	/// <summary>
	/// A triangle mesh exported from a <see cref="Manifold"/> with
	/// double-precision coordinates. All arrays are managed copies; the mesh has no
	/// native resources and needs no disposal.
	/// </summary>
	/// <remarks>
	/// The kernel computes in double precision, and this path is lossless:
	/// coordinates that went in through
	/// <see cref="Manifold.FromMesh64(ReadOnlySpan{double}, ReadOnlySpan{ulong}, uint)"/>
	/// and were left untouched by an operation come back bit-identical, with no
	/// <c>float</c> round-trip anywhere.
	/// </remarks>
	public sealed class MeshGL64
	{
		private MeshGL64(
			uint numProp,
			double[] vertProperties,
			uint[] triVerts,
			uint[] runIndex,
			uint[] runOriginalId,
			uint[] faceId)
		{
			this.NumProp = numProp;
			this.VertProperties = vertProperties;
			this.TriVerts = triVerts;
			this.RunIndex = runIndex;
			this.RunOriginalId = runOriginalId;
			this.FaceId = faceId;
		}

		/// <summary>
		/// Properties per vertex in <see cref="VertProperties"/>. Always at least 3;
		/// slots 0-2 are the position.
		/// </summary>
		public uint NumProp { get; }

		/// <summary>Interleaved vertex properties, <see cref="NumProp"/> doubles per vertex.</summary>
		public double[] VertProperties { get; }

		/// <summary>Three vertex indices per triangle, counter-clockwise seen from outside.</summary>
		public uint[] TriVerts { get; }

		/// <summary>
		/// Start offsets into <see cref="TriVerts"/>, one per run plus a trailing end
		/// sentinel, so this is one longer than <see cref="RunOriginalId"/>.
		/// </summary>
		public uint[] RunIndex { get; }

		/// <summary>Source mesh ID for each run, matching <see cref="Manifold.OriginalId"/> of the input it came from.</summary>
		public uint[] RunOriginalId { get; }

		/// <summary>The source face each triangle came from, one entry per triangle.</summary>
		public uint[] FaceId { get; }

		/// <summary>Number of triangles in the mesh.</summary>
		public int TriangleCount => this.TriVerts.Length / 3;

		/// <summary>Number of vertices in the mesh.</summary>
		public int VertexCount => this.NumProp == 0 ? 0 : this.VertProperties.Length / (int)this.NumProp;

		/// <summary>
		/// Copies every array out of a live native mesh handle. The caller still owns
		/// the handle and must destroy it; nothing in the returned object points at it.
		/// </summary>
		internal static unsafe MeshGL64 CopyFrom(IntPtr mesh)
		{
			uint numProp = NativeMethods.manifold_rs_meshgl64_num_prop(mesh);

			double* vertProperties = NativeMethods.manifold_rs_meshgl64_vert_properties(mesh, out nuint vertPropertiesLen);
			ulong* triVerts = NativeMethods.manifold_rs_meshgl64_tri_verts(mesh, out nuint triVertsLen);
			ulong* runIndex = NativeMethods.manifold_rs_meshgl64_run_index(mesh, out nuint runIndexLen);
			uint* runOriginalId = NativeMethods.manifold_rs_meshgl64_run_original_id(mesh, out nuint runOriginalIdLen);
			ulong* faceId = NativeMethods.manifold_rs_meshgl64_face_id(mesh, out nuint faceIdLen);

			return new MeshGL64(
				numProp,
				NativeArrays.Copy(vertProperties, vertPropertiesLen),
				NativeArrays.NarrowChecked(triVerts, triVertsLen, "tri_verts"),
				NativeArrays.NarrowChecked(runIndex, runIndexLen, "run_index"),
				// run_original_id is uint32_t on both paths: a mesh ID is not an index
				// into the mesh, so there is nothing to narrow.
				NativeArrays.Copy(runOriginalId, runOriginalIdLen),
				NativeArrays.NarrowChecked(faceId, faceIdLen, "face_id"));
		}
	}
}
