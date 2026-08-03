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

// The error codes a manifold handle can carry. Values and order come from the
// table in ffi/manifold_rs.h, which mirrors the C++ Manifold::Error enum.

namespace ManifoldRust
{
	/// <summary>
	/// Validation state of a <see cref="Manifold"/>. Anything other than
	/// <see cref="NoError"/> means the solid is empty: it exports as an empty mesh
	/// and is silently absorbed as empty geometry by
	/// <see cref="Manifold.BatchBoolean"/>, so check it after every import.
	/// </summary>
	public enum ManifoldStatus
	{
		/// <summary>The manifold is valid.</summary>
		NoError = 0,

		/// <summary>A vertex position was NaN or infinite.</summary>
		NonFiniteVertex = 1,

		/// <summary>The triangle soup does not describe a closed, orientable surface.</summary>
		NotManifold = 2,

		/// <summary>A triangle referenced a vertex index past the end of the vertex array.</summary>
		VertexOutOfBounds = 3,

		/// <summary>The property array length is not a whole number of vertices.</summary>
		PropertiesWrongLength = 4,

		/// <summary>Fewer than three properties per vertex, so there is no position.</summary>
		MissingPositionProperties = 5,

		/// <summary>The merge-from and merge-to vectors had different lengths.</summary>
		MergeVectorsDifferentLengths = 6,

		/// <summary>A merge entry referenced a vertex index that does not exist.</summary>
		MergeIndexOutOfBounds = 7,

		/// <summary>A transform array was not the expected length.</summary>
		TransformWrongLength = 8,

		/// <summary>The run index array was not the expected length.</summary>
		RunIndexWrongLength = 9,

		/// <summary>The face ID array was not one entry per triangle.</summary>
		FaceIdWrongLength = 10,

		/// <summary>The requested construction could not produce a valid solid.</summary>
		InvalidConstruction = 11,

		/// <summary>The result exceeded the kernel's size limits.</summary>
		ResultTooLarge = 12,

		/// <summary>Supplied halfedge tangents were invalid.</summary>
		InvalidTangents = 13,
	}
}
