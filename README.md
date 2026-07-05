# manifold-rust

[![Demo](README_HERO.png)](https://larsbrubaker.github.io/manifold-rust/)

## Support the Project

<a href="https://buymeacoffee.com/larsbrubaker"><img src="https://cdn.buymeacoffee.com/buttons/v2/default-yellow.png" alt="Buy Me A Coffee" height="50" width="210"></a>

manifold-rust is open-source and free to use, maintained in spare time as a labor of love. Friends James Smith and Dan Ruskin help out from time to time too.

If you find it useful, here are a few ways to help keep development going:

- **Donations:** [Buy Me a Coffee](https://buymeacoffee.com/larsbrubaker) — every coffee helps.
- **Star the repo:** Costs nothing and helps others find the project.
- **Report issues:** [Open an issue](https://github.com/larsbrubaker/manifold-rust/issues) for bugs or feature ideas.
- **Contribute:** PRs welcome — open an issue first to discuss larger changes.

Pure Rust port of [Manifold](https://github.com/elalish/manifold) — a geometry library for 3D boolean operations on triangle meshes.

> Part of the [rust-apps](https://github.com/larsbrubaker/rust-apps) suite — a collection of Rust graphics and geometry libraries by Lars Brubaker.

> **Status: Port complete.** All 18 phases of the C++ engine (v3.5.0) are implemented and
> every C++ test is ported or covered — 520 tests passing, 0 failing (the handful of
> `#[ignore]`d tests are debug-build-speed only and pass in release). Heavy boolean/CSG
> workloads run at parity with the sequential C++ build, and the optional `parallel`
> feature roughly doubles them. See [PORTING_PLAN.md](PORTING_PLAN.md) for the full record.

## What is Manifold?

Manifold is a high-performance C++ library for 3D solid modeling. It supports:

- Boolean operations (union, intersection, difference) on triangle meshes
- Mesh constructors (sphere, cube, cylinder, extrude, revolve)
- Cross-section (2D polygon) operations
- Smooth subdivision and SDF-based mesh generation
- Convex hull
- Minkowski sum/difference

This Rust port targets **exact numerical match** with the C++ implementation — same algorithms, same floating-point results, same triangle topology. Exactness is validated by instrumented, boolean-by-boolean trace comparison against a locally built C++ reference (see `validate-reference.ps1`), down to the tie-breaking order of symbolic-perturbation predicates.

## Why

[MatterHackers](https://www.matterhackers.com) uses 3D mesh boolean operations extensively in production for 3D printing workflows. A pure Rust implementation avoids FFI overhead and integrates cleanly with Rust tooling including WASM compilation.

## Usage

```rust
use manifold_rust::manifold::Manifold;
use manifold_rust::linalg::Vec3;

// Constructors
let cube = Manifold::cube(Vec3::new(1.0, 1.0, 1.0), true);
let sphere = Manifold::sphere(0.6, 32);

// Guaranteed-manifold booleans (also available as + - operators)
let difference = cube.difference(&sphere);
assert_eq!(difference.status(), manifold_rust::types::Error::NoError);

// Measure
println!("volume = {}", difference.volume());
println!("area   = {}", difference.surface_area());
println!("genus  = {}", difference.genus());

// Mesh I/O via MeshGL
let mesh = difference.get_mesh_gl(0);
let round_tripped = Manifold::from_mesh_gl(&mesh);
```

Enable parallel execution (results stay bit-identical to the sequential build —
only determinism-preserving sites are parallelized):

```toml
[dependencies]
manifold-rust = { version = "0.9", features = ["parallel"] }
```

## Demo

An interactive WASM demo is live at <https://larsbrubaker.github.io/manifold-rust/> —
booleans, extrude/revolve with twist, convex hull, subdivision, a Menger sponge, and more,
all running the Rust engine compiled to WebAssembly.

## Building

```bash
cargo build
cargo test
```

Restore the upstream C++ reference before doing exact-match validation:

```bash
git submodule update --init --recursive
```

Build and compare against the C++ reference with:

```powershell
./validate-reference.ps1
```

You can also run a narrower validation slice by phase, for example:

```powershell
./validate-reference.ps1 -Phase phase8
```

This configures and builds `cpp-reference/manifold`, runs the matching C++
reference tests, and then runs the corresponding Rust tests for the selected
phase.

For the WASM demo:
```bash
cd demo
bun run build:wasm
bun run dev
```

## Architecture

The port follows the C++ module structure:

| Rust module | C++ source | Description |
|-------------|-----------|-------------|
| `vec` / `linalg` | `linalg.h`, `vec.h` | Vector math, linear algebra |
| `polygon` | `polygon.cpp` | 2D polygon triangulation |
| `impl` | `impl.cpp`, `impl.h` | Core mesh data structure |
| `constructors` | `constructors.cpp` | Primitive mesh constructors |
| `boolean3` | `boolean3.cpp` | 3D boolean operations |
| `boolean_result` | `boolean_result.cpp` | Boolean output assembly |
| `csg_tree` | `csg_tree.cpp` | CSG tree evaluation |
| `edge_op` | `edge_op.cpp` | Edge manipulation |
| `face_op` | `face_op.cpp` | Face manipulation |
| `smoothing` | `smoothing.cpp` | Smooth normals and subdivision |
| `subdivision` | `subdivision.cpp` | Mesh subdivision |
| `properties` | `properties.cpp` | Mesh properties |
| `sdf` | `sdf.cpp` | SDF-based mesh generation |
| `quickhull` | `quickhull.cpp` | Convex hull |
| `minkowski` | `minkowski.cpp` | Minkowski operations |
| `cross_section` | `cross_section/` | 2D cross section |

## License

Apache-2.0 — matching the original Manifold library.

## Credits

- **[Emmett Lalish](https://github.com/elalish)** — Author of the original Manifold library.
- **[Lars Brubaker](https://github.com/larsbrubaker)** — Port author.
- **[MatterHackers](https://www.matterhackers.com)** — Sponsor.
