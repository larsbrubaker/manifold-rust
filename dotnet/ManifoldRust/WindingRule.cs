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

// Winding rules for the robust engine. The numbering matches the
// MANIFOLD_RS_WINDING_* defines in ffi/manifold_rs.h.

namespace ManifoldRust
{
	/// <summary>
	/// Which winding numbers the robust boolean engine counts as solid material.
	/// </summary>
	/// <remarks>
	/// <para>
	/// The rule is a robust-engine semantic. <see cref="BooleanEngine.Exact"/>
	/// ignores it, and <see cref="BooleanEngine.Auto"/> resolves to
	/// <see cref="BooleanEngine.Robust"/> whenever the rule is
	/// <see cref="Nonzero"/>.
	/// </para>
	/// <para>
	/// The choice only matters for geometry that is wound inconsistently — for a
	/// mesh whose shells are all wound outward the two rules agree exactly.
	/// </para>
	/// </remarks>
	public enum WindingRule
	{
		/// <summary>
		/// <c>{w &gt;= 1}</c>: only material the surface encloses with a positive
		/// winding number is solid. This is the default everywhere in the library
		/// and the rule that matches the exact engine, so it is what clean,
		/// consistently wound data should use. An inside-out shell encloses
		/// <c>w = -1</c> and therefore contributes nothing.
		/// </summary>
		Positive = 0,

		/// <summary>
		/// <c>{w != 0}</c>: any non-zero winding number is solid, so an inside-out
		/// shell is kept as material instead of vanishing. For scans and other
		/// imported geometry whose shells are wound inconsistently and cannot be
		/// repaired first — see <see cref="Manifold.RepairOrientation"/>, which
		/// fixes the winding itself and lets <see cref="Positive"/> keep working.
		/// </summary>
		Nonzero = 1,
	}
}
