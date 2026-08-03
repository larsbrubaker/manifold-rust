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

// The managed copy of an exported mesh. Deliberately not a handle wrapper: the
// native MeshGlRs is created, drained and destroyed inside Manifold.GetMeshGL,
// so nothing a caller holds has an unmanaged lifetime and there is no way to
// read through a pointer the FFI has already invalidated.

using System;

namespace ManifoldRust
{
	/// <summary>
	/// A triangle mesh exported from a <see cref="Manifold"/>. All arrays are
	/// managed copies; the mesh has no native resources and needs no disposal.
	/// </summary>
	public sealed class MeshGL
	{
		private MeshGL(
			uint numProp,
			float[] vertProperties,
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

		/// <summary>Interleaved vertex properties, <see cref="NumProp"/> floats per vertex.</summary>
		public float[] VertProperties { get; }

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
		internal static unsafe MeshGL CopyFrom(IntPtr mesh)
		{
			uint numProp = NativeMethods.manifold_rs_meshgl_num_prop(mesh);

			float* vertProperties = NativeMethods.manifold_rs_meshgl_vert_properties(mesh, out nuint vertPropertiesLen);
			uint* triVerts = NativeMethods.manifold_rs_meshgl_tri_verts(mesh, out nuint triVertsLen);
			uint* runIndex = NativeMethods.manifold_rs_meshgl_run_index(mesh, out nuint runIndexLen);
			uint* runOriginalId = NativeMethods.manifold_rs_meshgl_run_original_id(mesh, out nuint runOriginalIdLen);
			uint* faceId = NativeMethods.manifold_rs_meshgl_face_id(mesh, out nuint faceIdLen);

			return new MeshGL(
				numProp,
				Copy(vertProperties, vertPropertiesLen),
				Copy(triVerts, triVertsLen),
				Copy(runIndex, runIndexLen),
				Copy(runOriginalId, runOriginalIdLen),
				Copy(faceId, faceIdLen));
		}

		/// <summary>
		/// Copies <paramref name="length"/> elements out of a borrowed native array.
		/// An empty array comes back as a non-NULL pointer that must not be
		/// dereferenced, so the length is what decides, not the pointer.
		/// </summary>
		private static unsafe T[] Copy<T>(T* source, nuint length)
			where T : unmanaged
		{
			if (source == null || length == 0)
			{
				return Array.Empty<T>();
			}

			T[] result = new T[checked((int)length)];
			new ReadOnlySpan<T>(source, result.Length).CopyTo(result);
			return result;
		}
	}
}
