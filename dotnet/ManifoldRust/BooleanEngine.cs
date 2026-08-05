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

// Boolean engine codes. The numbering matches the MANIFOLD_RS_ENGINE_*
// defines in ffi/manifold_rs.h.

namespace ManifoldRust
{
	/// <summary>
	/// Which 3D boolean implementation the library runs; see
	/// <see cref="Manifold.DefaultBooleanEngine"/>.
	/// </summary>
	public enum BooleanEngine
	{
		/// <summary>
		/// The ported exact pipeline (default). Byte-identical results to the
		/// C++ reference; requires strictly manifold operands.
		/// </summary>
		Exact = 0,

		/// <summary>
		/// The robust engine: exact rational arithmetic that accepts closed,
		/// orientable but non-manifold geometry (imported via
		/// <see cref="Manifold.FromMeshRobust"/>). Slower; triangulation of
		/// results may differ from <see cref="Exact"/>.
		/// </summary>
		Robust = 1,

		/// <summary>
		/// <see cref="Exact"/> unless an operand carries non-manifold soup
		/// geometry, then <see cref="Robust"/>.
		/// </summary>
		Auto = 2,
	}
}
