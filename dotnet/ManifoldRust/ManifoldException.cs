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

// The exception every native failure sentinel is turned into. The native side
// never throws and never aborts; it returns NULL / -1 / 0 and records a message
// in a thread-local slot, so this type is where that message surfaces.

using System;

namespace ManifoldRust
{
	/// <summary>
	/// Thrown when a native manifold_rs call reports failure. <see cref="NativeError"/>
	/// carries the message the native side recorded for the failing call, and
	/// <see cref="Status"/> the validation state of the manifold involved when one
	/// was available.
	/// </summary>
	public sealed class ManifoldException : Exception
	{
		/// <summary>
		/// Creates an exception describing a native failure.
		/// </summary>
		/// <param name="message">What the caller was trying to do.</param>
		/// <param name="status">Status of the manifold involved, when known.</param>
		/// <param name="nativeError">Message read from manifold_rs_last_error, when any.</param>
		public ManifoldException(string message, ManifoldStatus? status = null, string? nativeError = null)
			: base(BuildMessage(message, status, nativeError))
		{
			this.Status = status;
			this.NativeError = nativeError;
		}

		/// <summary>Status of the manifold involved in the failure, if one was known.</summary>
		public ManifoldStatus? Status { get; }

		/// <summary>The native last-error message, or null if the slot was empty.</summary>
		public string? NativeError { get; }

		private static string BuildMessage(string message, ManifoldStatus? status, string? nativeError)
		{
			string text = message;

			if (status.HasValue)
			{
				text += $" (manifold status: {status.Value})";
			}

			if (!string.IsNullOrEmpty(nativeError))
			{
				text += $": {nativeError}";
			}

			return text;
		}
	}
}
