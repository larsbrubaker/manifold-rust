// About page: project story and methodology

export function init(container: HTMLElement): void {
  container.innerHTML = `
    <div class="home-page">
      <h2>About This Project</h2>

      <div class="about-content">
        <h3>What is Manifold?</h3>
        <p>
          <a href="https://github.com/elalish/manifold" target="_blank">Manifold</a> is a high-performance
          C++ library for 3D boolean mesh operations (union, intersection, difference) created by Emmett Lalish.
          It is used in 3D printing, CAD, and geometry processing pipelines where robust, watertight
          mesh operations are critical.
        </p>

        <h3>Why Port to Rust?</h3>
        <p>
          <a href="https://www.matterhackers.com" target="_blank">MatterHackers</a> needs native Rust
          for 3D printing workflows — no FFI overhead, full WASM compatibility for browser-based
          slicers, and Rust's safety guarantees for production geometry code.
        </p>

        <h3>Porting Methodology</h3>
        <ul>
          <li><strong>Exact numerical match</strong> — Same IEEE 754 f64 semantics, same algorithms,
          identical results on identical inputs. Not "close enough" — bit-for-bit identical.</li>
          <li><strong>Phase by phase</strong> — 18 phases from linear algebra up through boolean operations
          and WASM demos. Each phase fully tested before the next begins.</li>
          <li><strong>Test-first</strong> — C++ tests are ported to Rust alongside the implementation.
          Every function must reproduce the C++ output exactly.</li>
          <li><strong>No stubs</strong> — Every function is complete and production-ready.
          No <code>todo!()</code>, no <code>unimplemented!()</code>.</li>
        </ul>

        <h3>Architecture</h3>
        <table class="about-table">
          <tr><th>Decision</th><th>Choice</th><th>Rationale</th></tr>
          <tr><td>Linear algebra</td><td>Custom <code>linalg.rs</code></td><td>Must control exact f64 semantics; no glam</td></tr>
          <tr><td>Float precision</td><td>f64 throughout</td><td>f32 only at MeshGL boundary</td></tr>
          <tr><td>Cross-section</td><td><code>clipper2-rust</code></td><td>Our own Rust port of Clipper2</td></tr>
          <tr><td>Parallelism</td><td>Sequential first</td><td>Rayon behind feature flag later</td></tr>
          <tr><td>3D viewer</td><td>Three.js + WebGL</td><td>Industry-standard for browser 3D</td></tr>
        </table>

        <h3>Current Status</h3>
        <p>
          Phases 0–10 and 14 are complete (163 tests passing, ~12,000 lines of Rust).
          Phase 11 (Boolean Operations) is the critical next step — porting the full
          edge-face intersection algorithm with symbolic perturbation.
        </p>

        <h3>Related Projects</h3>
        <ul>
          <li><a href="https://larsbrubaker.github.io/agg-rust/" target="_blank">agg-rust</a> — Pure Rust port of AGG 2.6 (2D vector graphics, 903 tests, 64 demos)</li>
          <li><a href="https://larsbrubaker.github.io/clipper2-rust/" target="_blank">clipper2-rust</a> — Pure Rust port of Clipper2 (polygon clipping, 444 tests, 8 demos)</li>
          <li><a href="https://larsbrubaker.github.io/tess2-rust/" target="_blank">tess2-rust</a> — Pure Rust port of libtess2 (polygon tessellation, 11 tests, 6 demos)</li>
        </ul>
      </div>
    </div>
  `;
}
