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

// Host configuration for locating the shared library. This exists because the
// assembly installs its own DllImportResolver at module load: a host cannot
// call NativeLibrary.SetDllImportResolver for this assembly afterwards - the
// slot is already taken and a second registration throws - so an application
// with its own native layout needs a way in that is not an environment variable.

namespace ManifoldRust
{
	/// <summary>
	/// Settings for how the native manifold_rs library is located.
	/// </summary>
	public static class ManifoldNative
	{
		/// <summary>
		/// Full path to the native library (<c>manifold_rs.dll</c>,
		/// <c>libmanifold_rs.so</c> or <c>libmanifold_rs.dylib</c>), or null to let
		/// the normal probing decide.
		/// </summary>
		/// <remarks>
		/// This is checked before the <c>MANIFOLD_RS_NATIVE</c> environment variable
		/// and before the runtime's own resolution, so a host that lays out its
		/// natives somewhere unusual always wins.
		/// <para>
		/// It must be set before the first call into the library. Once the operating
		/// system has loaded the library the runtime caches the handle and never
		/// consults the resolver again, so a later change has no effect.
		/// </para>
		/// <para>
		/// A path that fails to load is not an error here; probing simply falls
		/// through to the next candidate.
		/// </para>
		/// </remarks>
		/// <seealso cref="Manifold.NativeVersion"/>
		public static string? LibraryPath { get; set; }
	}
}
