#!/usr/bin/env python3
# Copyright 2026 Lars Brubaker
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Drop wasm target features a consumer's Binaryen cannot parse from a static archive.

Why this exists
---------------
`cargo rustc -p manifold-ffi --target wasm32-unknown-emscripten --crate-type staticlib`
produces an archive of wasm objects, and every object carries a `target_features`
custom section naming the features its code assumes. Whoever links the archive - for
manifold_rs that is mono-wasm, through the .NET wasm-tools workload - unions those
sections into the linked module, and emcc turns the union into `--enable-<feature>`
flags for `wasm-opt`.

rustc's LLVM is newer than the Emscripten the workload ships (LLVM 21 against
Emscripten 3.1.56 = LLVM 19 / Binaryen 116), and LLVM 20 introduced two feature *names*
that Binaryen has never heard of, both of them implied refinements of features already
in the list:

    bulk-memory-opt          (implied by bulk-memory)
    call-indirect-overlong   (implied by reference-types)

wasm-opt then dies with "Unknown option '--enable-bulk-memory-opt'" *after* a link that
succeeded, which is what makes the failure so confusing. RUSTFLAGS cannot fix it: the
precompiled rust-std rlibs pulled into the archive carry the same section and are not
rebuilt. Deleting the section outright cannot fix it either - it also names the features
the module genuinely needs, so emcc stops passing `--enable-bulk-memory` and wasm-opt
fails validating the `memory.copy` instructions instead.

So the section is *rewritten*: every feature Binaryen understands is kept and only the
unknown names are dropped. This belongs to the artifact rather than to whoever links it,
which is why it runs in this repo's release pipeline (.github/workflows/release-nuget.yml)
and not in the consumer.

`llvm-objcopy` applies one --remove-section/--add-section pair to every member of an
archive in a single pass, so one call fixes the whole thing. The members do not all carry
the same feature set, so the payload written to all of them is their union minus the
dropped names - see the comment on that code for why the union is the safe flattening.

Revisit the whole thing once the wasm-tools workload moves past Emscripten 3.1.56.

Section payload format (WebAssembly tool conventions, "target_features"):

    LEB128 feature count, then per feature:
      1 byte prefix: '+' (used) or '-' (disallowed)
      LEB128 name length
      ASCII name
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile

# The two LLVM 20 names Binaryen 116 rejects. Overridable with --drop so a newer rustc
# that adds a third one can be handled without editing this file.
DEFAULT_DROPPED = ["bulk-memory-opt", "call-indirect-overlong"]

# What --check accepts, and it is an allowlist rather than a denylist on purpose. The
# release workflow builds with rust-toolchain@stable, which is unpinned by design, so the
# next rustc could start emitting a *third* name Binaryen has never heard of - and a check
# that only looked for the two known-bad names would wave it through and break every
# consumer's wasm-opt run instead.
#
# The list is exactly the set the wasm build of manifold-ffi has been observed to produce
# and link successfully: the seven "+" features that survive the rewrite (which are also
# the five WasmOptConfigurationFlags in build/ManifoldRust.props plus the bulk-memory and
# exception-handling the .NET SDK already enables for its own sources), and the
# "-shared-mem" the precompiled rust-std objects carry because this is a build without
# pthreads. It is not "everything Binaryen 116 understands" - Binaryen knows plenty of
# features rustc has never emitted here, and admitting them untested would defeat the
# point.
#
# So a new name failing this check is not necessarily a bug: it means a toolchain changed
# and someone has to decide whether the consumer's Binaryen knows the name (add it here)
# or does not (add it to DEFAULT_DROPPED). Both are one-line changes; neither should
# happen silently in a release.
DEFAULT_ALLOWED = [
    "bulk-memory",
    "exception-handling",
    "multivalue",
    "mutable-globals",
    "nontrapping-fptoint",
    "reference-types",
    "shared-mem",
    "sign-ext",
]

# A wasm custom section is: id 0, LEB size, LEB name length, name. "target_features" is
# 15 bytes, so a 0x0F immediately before the name pins a section start without having to
# walk each member's section list - and scanning the archive as one byte string finds the
# sections of every member at once.
MARKER = b"\x0ftarget_features"


def read_leb128(data, at):
    """Return (value, position after the value) for the unsigned LEB128 at `at`."""
    result = 0
    shift = 0
    while True:
        byte = data[at]
        at += 1
        result |= (byte & 0x7F) << shift
        if not byte & 0x80:
            return result, at
        shift += 7


def write_leb128(value):
    """Encode `value` as unsigned LEB128."""
    out = bytearray()
    while True:
        byte = value & 0x7F
        value >>= 7
        if value:
            out.append(byte | 0x80)
        else:
            out.append(byte)
            return bytes(out)


def parse_sections(data):
    """Every target_features payload in `data`, as a list of feature-string lists.

    A feature string is its prefix byte followed by its name, e.g. "+bulk-memory", so a
    '+' and a '-' of the same name never compare equal.
    """
    sections = []
    at = data.find(MARKER)
    while at >= 0:
        position = at + len(MARKER)
        count, position = read_leb128(data, position)
        features = []
        for _ in range(count):
            prefix = data[position:position + 1]
            length, position = read_leb128(data, position + 1)
            features.append((prefix + data[position:position + length]).decode("ascii"))
            position += length
        sections.append(features)
        at = data.find(MARKER, position)
    return sections


def encode_section(features):
    """The payload bytes for `features`, in the order given."""
    payload = bytearray(write_leb128(len(features)))
    for feature in features:
        name = feature[1:].encode("ascii")
        payload += feature[0].encode("ascii") + write_leb128(len(name)) + name
    return bytes(payload)


def describe(sections):
    """One line per distinct feature set, with how many objects carry it."""
    distinct = {}
    for features in sections:
        distinct.setdefault(tuple(features), 0)
        distinct[tuple(features)] += 1
    return "\n".join(
        "  {0} object(s): {1}".format(count, " ".join(features))
        for features, count in distinct.items()
    )


def find_objcopy(explicit):
    """The llvm-objcopy to use: the one asked for, else one on PATH, else the workload's.

    A plain `objcopy` is deliberately not accepted - GNU objcopy does not understand wasm
    objects, so falling back to it would fail in a much less obvious way.
    """
    if explicit:
        return explicit
    found = shutil.which("llvm-objcopy")
    if found:
        return found
    for variable in ("EMSDK", "DOTNET_EMSCRIPTEN_LLVM_ROOT"):
        root = os.environ.get(variable)
        if not root:
            continue
        for candidate in (
            os.path.join(root, "llvm-objcopy"),
            os.path.join(root, "bin", "llvm-objcopy"),
            os.path.join(root, "upstream", "bin", "llvm-objcopy"),
        ):
            if os.path.exists(candidate):
                return candidate
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("archive", help="the wasm static archive to rewrite, in place")
    parser.add_argument(
        "--drop",
        action="append",
        metavar="FEATURE",
        help="a feature name to remove; repeatable. Defaults to " + ", ".join(DEFAULT_DROPPED),
    )
    parser.add_argument(
        "--allow",
        action="append",
        metavar="FEATURE",
        help="a feature name --check accepts; repeatable. Defaults to " + ", ".join(DEFAULT_ALLOWED),
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="do not rewrite; exit non-zero if the archive names any feature outside --allow",
    )
    parser.add_argument(
        "--objcopy",
        help="path to llvm-objcopy (default: PATH, then EMSDK / DOTNET_EMSCRIPTEN_LLVM_ROOT)",
    )
    arguments = parser.parse_args()

    dropped = arguments.drop if arguments.drop else DEFAULT_DROPPED
    allowed = arguments.allow if arguments.allow else DEFAULT_ALLOWED

    with open(arguments.archive, "rb") as handle:
        data = handle.read()

    sections = parse_sections(data)
    if not sections:
        # Not a warning: something that carries no target_features is not a wasm archive
        # built by rustc, and rewriting the wrong file silently is the failure to avoid.
        sys.stderr.write(
            "error: no target_features section in '{0}' - is it a wasm archive?\n".format(
                arguments.archive
            )
        )
        return 1

    print("{0}: {1} object(s) with a target_features section".format(arguments.archive, len(sections)))
    print(describe(sections))

    present = sorted({
        feature[1:]
        for features in sections
        for feature in features
        if feature[1:] in dropped
    })

    if arguments.check:
        unexpected = sorted({
            feature[1:]
            for features in sections
            for feature in features
            if feature[1:] not in allowed
        })
        if unexpected:
            # Two different failures wearing the same face, and the fix differs, so say
            # which one this is rather than making the reader diff two lists by eye.
            known_bad = [name for name in unexpected if name in dropped]
            surprises = [name for name in unexpected if name not in dropped]
            if known_bad:
                sys.stderr.write(
                    "error: '{0}' still names {1}, which the consuming Binaryen cannot parse. "
                    "The rewrite did not run, or ran against a different file.\n".format(
                        arguments.archive, ", ".join(known_bad)
                    )
                )
            if surprises:
                sys.stderr.write(
                    "error: '{0}' names {1}, which this pipeline has never shipped. A toolchain "
                    "changed. If the consuming Binaryen understands the name, add it to "
                    "DEFAULT_ALLOWED; if it does not, add it to DEFAULT_DROPPED so the rewrite "
                    "strips it.\n".format(arguments.archive, ", ".join(surprises))
                )
            return 1
        print("check: every feature named is one of {0}".format(", ".join(allowed)))
        return 0

    if not present:
        # Idempotent, and correct for a future toolchain that never emits these names.
        print("nothing to drop; leaving '{0}' untouched".format(arguments.archive))
        return 0

    # llvm-objcopy writes one payload to every member, and the members do NOT all agree:
    # a cargo-built archive has the rust-std objects carrying "-shared-mem" that the
    # crate's own objects do not. The union is the safe flattening, because it is exactly
    # what wasm-ld would compute anyway - "+" means the output requires the feature and
    # "-" means no object may use it, and both of those already hold for the whole link
    # the moment one member says so. Taking any single member's list instead would either
    # invent or discard a constraint.
    kept = []
    for features in sections:
        for feature in features:
            if feature[1:] not in dropped and feature not in kept:
                kept.append(feature)

    # No sane payload exists if one object requires a feature another forbids, and
    # picking either way would turn a link-time error into a silent miscompile.
    conflicting = sorted({
        feature[1:] for feature in kept if ("+" + feature[1:] in kept and "-" + feature[1:] in kept)
    })
    if conflicting:
        sys.stderr.write(
            "error: '{0}' has objects that require {1} and objects that forbid it. There is no "
            "single payload llvm-objcopy could write to every member.\n".format(
                arguments.archive, ", ".join(conflicting)
            )
        )
        return 1

    # Required features first, then forbidden ones, each alphabetical - the order LLVM
    # itself writes, so a rewritten section reads like an unrewritten one.
    kept.sort(key=lambda feature: (feature[0] != "+", feature[1:]))
    print("dropping {0}".format(", ".join(present)))
    print("keeping  {0}".format(" ".join(kept)))

    objcopy = find_objcopy(arguments.objcopy)
    if not objcopy:
        sys.stderr.write(
            "error: no llvm-objcopy found. Pass --objcopy, or set EMSDK / "
            "DOTNET_EMSCRIPTEN_LLVM_ROOT.\n"
        )
        return 1

    with tempfile.TemporaryDirectory() as directory:
        payload = os.path.join(directory, "target_features.bin")
        with open(payload, "wb") as handle:
            handle.write(encode_section(kept))

        subprocess.run(
            [
                objcopy,
                "--remove-section=target_features",
                "--add-section=target_features=" + payload,
                arguments.archive,
            ],
            check=True,
        )

    with open(arguments.archive, "rb") as handle:
        rewritten = parse_sections(handle.read())

    if len(rewritten) != len(sections):
        sys.stderr.write(
            "error: '{0}' had {1} target_features sections before the rewrite and {2} after. "
            "llvm-objcopy did not touch every member.\n".format(
                arguments.archive, len(sections), len(rewritten)
            )
        )
        return 1

    still_there = sorted({
        feature[1:] for features in rewritten for feature in features if feature[1:] in dropped
    })
    if still_there:
        sys.stderr.write(
            "error: '{0}' still names {1} after the rewrite.\n".format(
                arguments.archive, ", ".join(still_there)
            )
        )
        return 1

    print("rewrote {0} object(s) in '{1}'".format(len(rewritten), arguments.archive))
    return 0


if __name__ == "__main__":
    sys.exit(main())
