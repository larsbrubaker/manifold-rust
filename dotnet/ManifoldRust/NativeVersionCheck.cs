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

// The managed/native version handshake. The binding and the cdylib are shipped
// as one package but a consumer can point MANIFOLD_RS_NATIVE or
// ManifoldNative.LibraryPath at any library on disk, so nothing stops an old
// native from being loaded by a new binding. That combination is not a missing
// export - the runtime raises EntryPointNotFoundException for those - it is the
// silent case: an export that exists but whose contract changed. Checking the
// version string once turns it into one clear exception.

using System;

namespace ManifoldRust
{
	internal static class NativeVersionCheck
	{
		/// <summary>
		/// Prefix the native version string must start with. Only major.minor is
		/// compared: a patch release of the FFI crate cannot change the ABI, so
		/// pinning to it would reject a perfectly good native.
		/// </summary>
		/// <remarks>
		/// The trailing dot is load bearing. Without it "manifold-ffi 0.2" also
		/// prefixes "manifold-ffi 0.20.0", which is a different ABI that would then
		/// be accepted silently - the exact failure this check exists to prevent.
		/// </remarks>
		internal const string ExpectedPrefix = "manifold-ffi 0.3.";

		// Lazy caches the exception from a failing factory and re-throws it on every
		// later access, so a mismatched native costs exactly one native call for the
		// life of the process and still fails every entry point.
		private static readonly Lazy<bool> Checked = new Lazy<bool>(Check);

		/// <summary>
		/// Ensures the loaded native library speaks the ABI this binding was written
		/// against. Cheap after the first call; the check itself runs once per process.
		/// </summary>
		/// <exception cref="ManifoldException">
		/// The native library reports a version this binding does not support.
		/// </exception>
		internal static void Verify()
		{
			_ = Checked.Value;
		}

		/// <summary>
		/// Whether <paramref name="version"/> - a string in the shape
		/// <c>manifold-ffi 0.2.0 (manifold-rust 0.9.3)</c> - is an ABI this binding
		/// can talk to. Separate from <see cref="Check"/> so the comparison can be
		/// tested without a native library.
		/// </summary>
		internal static bool IsCompatible(string version)
		{
			return version.StartsWith(ExpectedPrefix, StringComparison.Ordinal);
		}

		/// <summary>
		/// The message thrown for an incompatible native. Separate from
		/// <see cref="Check"/> for the same reason as <see cref="IsCompatible"/>: this
		/// text is the entire value of the handshake to whoever hits it, so it is
		/// worth being able to test it without arranging for a stale library.
		/// </summary>
		internal static string MismatchMessage(string version)
		{
			return $"The loaded manifold_rs library reports \"{version}\", but this binding requires \"{ExpectedPrefix}x\". " +
				"Update the native library, or check ManifoldNative.LibraryPath and the MANIFOLD_RS_NATIVE environment variable for a stale one.";
		}

		private static bool Check()
		{
			string version = Manifold.NativeVersion;

			if (!IsCompatible(version))
			{
				throw new ManifoldException(MismatchMessage(version));
			}

			return true;
		}
	}
}
