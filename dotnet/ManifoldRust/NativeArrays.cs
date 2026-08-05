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

// Draining borrowed native arrays into managed ones. Shared by MeshGL and
// MeshGL64, which follow the same accessor contract from ffi/manifold_rs.h: the
// pointer borrows from the mesh handle and the length, not the pointer, is what
// says whether there is anything to read.

using System;

namespace ManifoldRust
{
	internal static class NativeArrays
	{
		/// <summary>
		/// Copies <paramref name="length"/> elements out of a borrowed native array.
		/// An empty array comes back as a non-NULL pointer that must not be
		/// dereferenced, so the length is what decides, not the pointer.
		/// </summary>
		internal static unsafe T[] Copy<T>(T* source, nuint length)
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

		/// <summary>
		/// Copies a borrowed native <c>uint64_t</c> array into a managed
		/// <c>uint[]</c>, the width the rest of this binding - and every graphics API
		/// a caller will hand the indices to - uses.
		/// </summary>
		/// <remarks>
		/// The kernel narrows indices to 32 bits internally, so nothing it exports
		/// can exceed <see cref="uint.MaxValue"/> today. The check is here anyway
		/// because the alternative to checking is <em>wrapping</em>: an index of 2^32
		/// would silently become 0 and produce geometry that looks fine and is wrong.
		/// The FFI rejects the same case on import for the same reason.
		/// </remarks>
		internal static unsafe uint[] NarrowChecked(ulong* source, nuint length, string field)
		{
			if (source == null || length == 0)
			{
				return Array.Empty<uint>();
			}

			return NarrowChecked(new ReadOnlySpan<ulong>(source, checked((int)length)), field);
		}

		/// <summary>
		/// Narrows every element to <c>uint</c>, refusing to wrap. See the pointer
		/// overload for why this is checked rather than a cast.
		/// </summary>
		/// <exception cref="ManifoldException">
		/// An element does not fit in 32 bits; the message names the field and the
		/// offending element.
		/// </exception>
		internal static uint[] NarrowChecked(ReadOnlySpan<ulong> source, string field)
		{
			if (source.Length == 0)
			{
				return Array.Empty<uint>();
			}

			uint[] result = new uint[source.Length];
			for (int i = 0; i < source.Length; i++)
			{
				ulong value = source[i];
				if (value > uint.MaxValue)
				{
					throw new ManifoldException(
						$"{field}[{i}] = {value} does not fit in 32 bits, so it cannot be narrowed without wrapping");
				}

				result[i] = (uint)value;
			}

			return result;
		}
	}
}
