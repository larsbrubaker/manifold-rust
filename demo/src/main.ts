// Main entry point - SPA router and WASM initialization

import { initWasm } from './wasm.ts';

// Demo page modules (lazy loaded)
type DemoInit = (container: HTMLElement) => (() => void) | void;
const demoModules: Record<string, () => Promise<{ init: DemoInit }>> = {
  'extrude-revolve': () => import('./demos/extrude-revolve.ts' as any),
  'convex-hull': () => import('./demos/convex-hull.ts' as any),
  'boolean-gallery': () => import('./demos/boolean-gallery.ts' as any),
  'menger-sponge': () => import('./demos/menger-sponge.ts' as any),
  'extrude-twist': () => import('./demos/extrude-twist.ts' as any),
  'revolve-partial': () => import('./demos/revolve-partial.ts' as any),
  'smooth-shapes': () => import('./demos/smooth-shapes.ts' as any),
  'properties': () => import('./demos/properties.ts' as any),
  'test-gallery': () => import('./demos/test-gallery.ts' as any),
  'about': () => import('./demos/about.ts' as any),
};

let currentCleanup: (() => void) | null = null;

// Mobile sidebar toggle
const menuToggle = document.getElementById('menu-toggle')!;
const sidebar = document.getElementById('sidebar')!;
const sidebarOverlay = document.getElementById('sidebar-overlay')!;

function openSidebar() {
  sidebar.classList.add('open');
  menuToggle.classList.add('open');
  sidebarOverlay.classList.add('visible');
}

function closeSidebar() {
  sidebar.classList.remove('open');
  menuToggle.classList.remove('open');
  sidebarOverlay.classList.remove('visible');
}

menuToggle.addEventListener('click', () => {
  if (sidebar.classList.contains('open')) closeSidebar();
  else openSidebar();
});

sidebarOverlay.addEventListener('click', closeSidebar);

document.querySelectorAll('.nav-link').forEach(link => {
  link.addEventListener('click', closeSidebar);
});

function getRoute(): string {
  const hash = window.location.hash.slice(2) || '';
  return hash || 'home';
}

function updateNav(route: string) {
  document.querySelectorAll('.nav-link').forEach(el => {
    const r = (el as HTMLElement).dataset.route;
    el.classList.toggle('active', r === route);
  });
}

function renderHome(container: HTMLElement) {
  container.innerHTML = `
    <div class="home-page">
      <div class="github-badge">
        <a href="https://github.com/larsbrubaker/manifold-rust" target="_blank" class="github-badge-link">
          <svg height="20" viewBox="0 0 16 16" width="20" fill="currentColor"><path d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0016 8c0-4.42-3.58-8-8-8z"/></svg>
          <span>larsbrubaker/manifold-rust</span>
        </a>
      </div>
      <div class="hero">
        <h1>Manifold <span>for Rust</span></h1>
        <p>
          A pure Rust port of the Manifold 3D geometry library. Explore interactive demos
          showcasing mesh primitives, boolean operations, extrusion, convex hulls, and more
          &mdash; all running in your browser via WebAssembly.
        </p>
      </div>
      <div class="feature-grid">
        <a href="#/extrude-revolve" class="feature-card">
          <div class="card-thumb" style="background-image:url(public/thumbs/extrude-revolve.jpg)"></div>
          <h3>Extrude &amp; Revolve</h3>
          <p>Turn 2D cross-sections into 3D solids by extrusion or revolution with adjustable parameters.</p>
        </a>
        <a href="#/convex-hull" class="feature-card">
          <div class="card-thumb" style="background-image:url(public/thumbs/convex-hull.jpg)"></div>
          <h3>Convex Hull</h3>
          <p>Generate convex hulls from random 3D point clouds using the QuickHull algorithm.</p>
        </a>
        <a href="#/boolean-gallery" class="feature-card">
          <div class="card-thumb" style="background-image:url(public/thumbs/boolean-gallery.jpg)"></div>
          <h3>Boolean Gallery</h3>
          <p>Mix and match cubes, spheres, and cylinders with union, intersection, and difference in 3D.</p>
        </a>
        <a href="#/menger-sponge" class="feature-card">
          <div class="card-thumb" style="background-image:url(public/thumbs/menger-sponge.jpg)"></div>
          <h3>Menger Sponge</h3>
          <p>Recursive boolean subtraction creates this classic fractal at adjustable depth levels.</p>
        </a>
        <a href="#/extrude-twist" class="feature-card">
          <div class="card-thumb" style="background-image:url(public/thumbs/extrude-twist.jpg)"></div>
          <h3>Extrude + Twist</h3>
          <p>Extrude circular profiles with twist rotation and taper scaling for spiral shapes.</p>
        </a>
        <a href="#/revolve-partial" class="feature-card">
          <div class="card-thumb" style="background-image:url(public/thumbs/revolve-partial.jpg)"></div>
          <h3>Partial Revolve</h3>
          <p>Revolve 2D profiles by partial angles to create torus arcs, rings, and curved solids.</p>
        </a>
        <a href="#/smooth-shapes" class="feature-card">
          <div class="card-thumb" style="background-image:url(public/thumbs/smooth-shapes.jpg)"></div>
          <h3>Subdivision</h3>
          <p>Subdivide base shapes with adjustable refinement levels for smoother, denser meshes.</p>
        </a>
        <a href="#/properties" class="feature-card">
          <div class="card-thumb" style="background-image:url(public/thumbs/properties.jpg)"></div>
          <h3>Properties</h3>
          <p>Measure volume, surface area, vertex and triangle counts for any generated mesh.</p>
        </a>
        <a href="#/test-gallery" class="feature-card">
          <div class="card-thumb" style="background-image:url(public/thumbs/test-gallery.jpg)"></div>
          <h3>Test Gallery</h3>
          <p>Browse WASM visualizations of key tests: mirror, split, warp, spiral, boolean precision, and more.</p>
        </a>
        <a href="#/about" class="feature-card" style="border-color: var(--accent); background: var(--accent-light);">
          <span class="card-icon">&#128214;</span>
          <h3>About This Project</h3>
          <p>How we ported 16,000 lines of C++ to Rust with exact numerical matching.</p>
        </a>
      </div>

      <div class="about-section">
        <h2>About This Project</h2>
        <p>
          This is a pure Rust port of Emmett Lalish's
          <a href="https://github.com/elalish/manifold" target="_blank">Manifold</a> C++ library,
          implementing 3D boolean mesh operations with exact numerical matching.
          The demos above run entirely in your browser via WebAssembly compiled from the Rust source.
        </p>
        <p style="margin-top: 12px">
          Ported by <strong>Lars Brubaker</strong>, sponsored by
          <a href="https://www.matterhackers.com" target="_blank">MatterHackers</a>.
        </p>
        <div class="stats-row">
          <div class="stat">
            <div class="stat-value">12K+</div>
            <div class="stat-label">Lines of Rust</div>
          </div>
          <div class="stat">
            <div class="stat-value">228</div>
            <div class="stat-label">Tests Passing</div>
          </div>
          <div class="stat">
            <div class="stat-value">f64</div>
            <div class="stat-label">Precision</div>
          </div>
          <div class="stat">
            <div class="stat-value">WASM</div>
            <div class="stat-label">Browser Ready</div>
          </div>
        </div>
      </div>
    </div>
  `;
}

async function navigate(route: string) {
  const container = document.getElementById('main-content')!;

  if (currentCleanup) {
    currentCleanup();
    currentCleanup = null;
  }

  updateNav(route);

  if (route === 'home') {
    renderHome(container);
    return;
  }

  const loader = demoModules[route];
  if (!loader) {
    container.innerHTML = `<div class="home-page"><h2>Page not found</h2><p>Unknown route: ${route}</p></div>`;
    return;
  }

  container.innerHTML = `<div class="home-page" style="display:flex;align-items:center;justify-content:center;height:80vh;"><p style="color:var(--text-muted)">Loading demo...</p></div>`;

  try {
    await initWasm();
    const mod = await loader();
    container.innerHTML = '';
    const cleanup = mod.init(container);
    if (cleanup) currentCleanup = cleanup;
  } catch (e) {
    console.error('Failed to load demo:', e);
    container.innerHTML = `<div class="home-page"><h2>Error loading demo</h2><pre style="color:#c44">${e}</pre></div>`;
  }
}

window.addEventListener('hashchange', () => navigate(getRoute()));
navigate(getRoute());

// ===================== Thumbnail Capture (localhost only) =====================

const DEMO_ROUTES = [
  'extrude-revolve', 'convex-hull', 'boolean-gallery', 'menger-sponge',
  'extrude-twist', 'revolve-partial', 'smooth-shapes', 'properties', 'test-gallery',
];

function isLocalhost(): boolean {
  return location.hostname === 'localhost' || location.hostname === '127.0.0.1';
}

async function captureAllThumbnails() {
  const btn = document.getElementById('capture-thumbs-btn')!;
  btn.classList.add('capturing');
  btn.setAttribute('disabled', 'true');

  const originalRoute = getRoute();
  const results: string[] = [];

  for (const route of DEMO_ROUTES) {
    btn.textContent = `Capturing ${route}...`;
    window.location.hash = `#/${route}`;
    await navigate(route);

    // Wait for WASM + Three.js to render
    await new Promise(r => setTimeout(r, 2000));

    const canvas = document.querySelector('#viewer-container canvas') as HTMLCanvasElement | null;
    if (!canvas) {
      results.push(`${route}: no canvas found`);
      continue;
    }

    // Capture at 320x240
    const tmp = document.createElement('canvas');
    tmp.width = 320;
    tmp.height = 240;
    const ctx = tmp.getContext('2d')!;
    ctx.drawImage(canvas, 0, 0, 320, 240);
    const dataUrl = tmp.toDataURL('image/jpeg', 0.85);

    try {
      const resp = await fetch('/api/save-thumb', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ name: route, data: dataUrl }),
      });
      const json = await resp.json();
      results.push(`${route}: ${json.size} bytes`);
    } catch (e) {
      results.push(`${route}: upload failed - ${e}`);
    }
  }

  // Navigate back to original route
  window.location.hash = originalRoute === 'home' ? '#/' : `#/${originalRoute}`;
  await navigate(originalRoute);

  btn.classList.remove('capturing');
  btn.removeAttribute('disabled');
  btn.textContent = 'Update Thumbnails';

  console.log('Thumbnail capture complete:', results);
}

// Show button only on localhost
if (isLocalhost()) {
  const btn = document.getElementById('capture-thumbs-btn')!;
  btn.style.display = 'block';
  btn.addEventListener('click', captureAllThumbnails);
}
