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

// Boolean op codes. The numbering matches C++ ManifoldOpType so call sites
// ported from manifoldc keep working unchanged.

namespace ManifoldRust
{
	/// <summary>
	/// The boolean operation
	/// <see cref="Manifold.BatchBoolean(System.Collections.Generic.IReadOnlyList{Manifold}, ManifoldOpType)"/>
	/// applies.
	/// </summary>
	public enum ManifoldOpType
	{
		/// <summary>Union of every operand.</summary>
		Add = 0,

		/// <summary>The first operand minus the union of all the others.</summary>
		Subtract = 1,

		/// <summary>Intersection of every operand.</summary>
		Intersect = 2,
	}
}
