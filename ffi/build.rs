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

// Publishes the core crate's version to the FFI build as MANIFOLD_RUST_VERSION.
// Cargo only exposes CARGO_PKG_VERSION for the package being compiled, and
// manifold_rs_version() must report the geometry kernel version — that is the
// number a caller needs when triaging a result mismatch.

use std::path::Path;

fn main() {
    let manifest = Path::new(env!("CARGO_MANIFEST_DIR")).join("..").join("Cargo.toml");
    println!("cargo:rerun-if-changed={}", manifest.display());

    let text = std::fs::read_to_string(&manifest).expect("read core crate Cargo.toml");
    // The root manifest opens with [workspace], so start after [package] — and
    // stop at the next section header, or a `version` belonging to a later
    // table (a dependency, say) would be reported as the kernel version.
    // Both failures panic: a wrong version string here is a lie the consumer
    // would carry into bug reports, so a broken build is preferable.
    let package = text
        .split("[package]")
        .nth(1)
        .unwrap_or_else(|| panic!("no [package] section in {}", manifest.display()));
    let version = package
        .lines()
        .take_while(|line| !line.trim_start().starts_with('['))
        .find_map(|line| line.strip_prefix("version = "))
        .map(|v| v.trim().trim_matches('"').to_string())
        .unwrap_or_else(|| {
            panic!(
                "no `version = \"...\"` line in the [package] section of {}",
                manifest.display()
            )
        });

    println!("cargo:rustc-env=MANIFOLD_RUST_VERSION={version}");
}
