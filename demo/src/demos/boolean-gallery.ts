// Boolean Gallery: booleans on either built-in primitives or random
// Thingi10K mesh pairs, with a choice of boolean engine (Exact / Robust /
// Auto). A Source dropdown switches between the two mesh sources; every
// other control — operation, engine, offsets, animate, wireframe — applies
// to whichever source is active. Every operation's inputs are captured for
// the Copy Debug Info button so failures can be reproduced from a pasted
// report.

import { ThreeViewer } from '../three-viewer.ts';
import { createSlider, createDropdown, createCheckbox, createButton, createReadout, updateReadout } from '../controls.ts';
import {
  booleanGalleryMeshEngine, importTriangleSoup, importedBoolean,
  type MeshData, type ImportedMesh, type BooleanEngine,
} from '../wasm.ts';
import { loadSetting, saveSetting } from '../settings.ts';
import { pickRandomModel, pairOperandKinds, fetchMesh, meshZipUrl, type ThingiModel, type PairKind } from '../thingi.ts';

const DEMO = 'boolean-gallery';

const SOURCES = [
  { value: 'builtin', text: 'Built In' },
  { value: 'thingi', text: 'Thingi10K' },
];

const SHAPES = [
  { value: '0', text: 'Cube' },
  { value: '1', text: 'Sphere' },
  { value: '2', text: 'Cylinder' },
  { value: '3', text: 'Spiky Dodecahedron' },
];

const OPS = [
  { value: '0', text: 'Union' },
  { value: '1', text: 'Intersection' },
  { value: '2', text: 'Difference' },
];

const ENGINES = [
  { value: '0', text: 'Exact (C++-matching)' },
  { value: '1', text: 'Robust (non-manifold)' },
  { value: '2', text: 'Auto' },
];

const PAIR_KINDS = [
  { value: 'mm', text: 'Manifold + Manifold' },
  { value: 'mn', text: 'Manifold + Non-manifold' },
  { value: 'nn', text: 'Non-manifold + Non-manifold' },
];

const OP_COLORS = [0x4488cc, 0x44aa44, 0xcc4444];
const OP_NAMES = ['union', 'intersection', 'difference'];
const ENGINE_NAMES = ['exact', 'robust', 'auto'];
const SHAPE_NAMES = ['cube', 'sphere', 'cylinder', 'spiky-dodecahedron'];

interface ThingiOperand {
  model: ThingiModel;
  handle: ImportedMesh;
  normalization: { center: [number, number, number]; scale: number };
}

export function init(container: HTMLElement): () => void {
  container.innerHTML = `
    <div class="demo-page">
      <div class="demo-header">
        <h2>Boolean Gallery</h2>
        <p>Combine meshes with union, intersection, and difference. Pick the source —
        built-in primitives or random real-world Thingi10K pairs (including non-manifold
        ones) — and the engine: the exact C++-matching pipeline or the robust engine that
        accepts non-manifold input.</p>
      </div>
      <div class="demo-layout">
        <div class="demo-canvas-area" id="viewer-container"></div>
        <div class="demo-controls" id="controls"></div>
      </div>
    </div>
  `;

  const viewerEl = document.getElementById('viewer-container')!;
  const controlsEl = document.getElementById('controls')!;
  const viewer = new ThreeViewer(viewerEl);

  let source = loadSetting(DEMO, 'source', 'builtin') as 'builtin' | 'thingi';
  let shapeA = loadSetting(DEMO, 'shapeA', 3);
  let shapeB = loadSetting(DEMO, 'shapeB', 3);
  let op = loadSetting(DEMO, 'op', 0);
  let engine = loadSetting(DEMO, 'engine', 0) as BooleanEngine;
  let offsetX = loadSetting(DEMO, 'offsetX', 0.3);
  let offsetY = loadSetting(DEMO, 'offsetY', 0.0);
  let offsetZ = loadSetting(DEMO, 'offsetZ', 0.0);
  let wireframe = loadSetting(DEMO, 'wireframe', false);
  let animate = loadSetting(DEMO, 'animate', true);
  let animating = false;
  let animId = 0;
  let rotX = 0, rotY = 0, rotZ = 0;
  const ROT_SPEED_X = 0.7;  // degrees per frame
  const ROT_SPEED_Y = 1.5;
  const ROT_SPEED_Z = 0.3;

  // Thingi10K state
  let pairKind = loadSetting(DEMO, 'pairKind', 'mn') as PairKind;
  let thingiA: ThingiOperand | null = null;
  let thingiB: ThingiOperand | null = null;
  let skippedImports: { id: number; status: string }[] = [];
  let loadingPair = false;
  let disposed = false;

  let lastDebugInfo: any = null;
  let frameCount = 0;

  // Thingi mode is active only once a pair is actually loaded; while a pair
  // is in flight the built-in shapes stay on screen.
  const inThingi = () => source === 'thingi' && !!(thingiA && thingiB);

  const readout = createReadout();
  const errorBox = document.createElement('div');
  errorBox.className = 'demo-note';
  errorBox.style.display = 'none';
  const thingiInfo = document.createElement('div');
  thingiInfo.className = 'demo-note';
  thingiInfo.style.display = 'none';

  function describeOperand(o: ThingiOperand) {
    return {
      id: o.model.id,
      name: o.model.name,
      url: meshZipUrl(o.model),
      faces: o.model.faces,
      edge_manifold: o.model.edge_manifold,
      vertex_manifold: o.model.vertex_manifold,
      imported: { num_tri: o.handle.num_tri, is_soup: o.handle.is_soup, status: o.handle.status },
      normalization: o.normalization,
    };
  }

  // Capture the full operation before running it, so even a crash/hang
  // leaves a reproducible record (also persisted to localStorage).
  function captureDebugInfo(): any {
    const info: any = {
      demo: DEMO,
      timestamp: new Date().toISOString(),
      frame: frameCount++,
      source,
      engine: ENGINE_NAMES[engine],
      op: OP_NAMES[op] || op,
      offset: [offsetX, offsetY, offsetZ],
      rotation_degrees: [rotX, rotY, rotZ],
      result: null,
      error: null,
    };
    if (inThingi()) {
      info.pair_kind = pairKind;
      info.model_a = describeOperand(thingiA!);
      info.model_b = describeOperand(thingiB!);
      if (skippedImports.length) info.skipped_imports = skippedImports;
    } else {
      info.shape_a = SHAPE_NAMES[shapeA] || shapeA;
      info.shape_b = SHAPE_NAMES[shapeB] || shapeB;
    }
    lastDebugInfo = info;
    try { localStorage.setItem('boolean-gallery-debug', JSON.stringify(info)); } catch { /* ignore */ }
    return info;
  }

  function showReadout(data: MeshData) {
    errorBox.style.display = 'none';
    updateReadout(readout, [
      { label: 'Vertices', value: String(data.num_vert) },
      { label: 'Triangles', value: String(data.num_tri) },
      { label: 'Volume', value: data.volume.toFixed(4) },
      { label: 'Surface Area', value: data.surface_area.toFixed(4) },
    ]);
  }

  function update(silent = false) {
    if (source === 'thingi' && !inThingi()) return; // pair still loading
    const info = captureDebugInfo();
    try {
      let data: MeshData;
      if (inThingi()) {
        data = importedBoolean(thingiA!.handle, thingiB!.handle, op, engine, offsetX, offsetY, offsetZ, rotX, rotY, rotZ);
      } else {
        data = booleanGalleryMeshEngine(shapeA, shapeB, op, offsetX, offsetY, offsetZ, rotX, rotY, rotZ, engine);
      }
      info.result = { num_tri: data.num_tri, num_vert: data.num_vert, volume: data.volume };
      try { localStorage.setItem('boolean-gallery-debug', JSON.stringify(info)); } catch { /* ignore */ }
      viewer.setMesh(data);
      viewer.setColor(OP_COLORS[op] || 0x4488cc);
      errorBox.style.display = 'none';
      showReadout(data);
    } catch (e) {
      info.error = String(e);
      try { localStorage.setItem('boolean-gallery-debug', JSON.stringify(info)); } catch { /* ignore */ }
      console.error('Boolean op failed:', info, e);
      if (!silent) {
        errorBox.style.display = 'block';
        errorBox.innerHTML = `<strong>Boolean operation failed:</strong> ${String(e)}<br>` +
          `Use <em>Copy Debug Info</em> to capture the inputs for a bug report.`;
        updateReadout(readout, []);
      }
    }
  }

  // Booleans run synchronously on the main thread; animating a pair whose
  // boolean takes seconds per frame would freeze the tab. Auto-pause when a
  // frame blows this budget — the user can re-enable Animate any time.
  const ANIMATE_BUDGET_MS = 1500;

  function animateStep() {
    if (!animating) return;
    rotX = (rotX + ROT_SPEED_X) % 360;
    rotY = (rotY + ROT_SPEED_Y) % 360;
    rotZ = (rotZ + ROT_SPEED_Z) % 360;
    const t0 = performance.now();
    update(true);
    const dt = performance.now() - t0;
    if (dt > ANIMATE_BUDGET_MS) {
      toggleAnimate(false);
      const box = animateBox.querySelector('input') as HTMLInputElement | null;
      if (box) box.checked = false;
      errorBox.style.display = 'block';
      errorBox.innerHTML = `<strong>Animation paused:</strong> this boolean takes ` +
        `${(dt / 1000).toFixed(1)} s per frame. Re-check Animate to continue anyway.`;
      return;
    }
    animId = requestAnimationFrame(animateStep);
  }

  function toggleAnimate(on: boolean) {
    saveSetting(DEMO, 'animate', on);
    animate = on;
    if (on) {
      if (!animating) {
        animating = true;
        animateStep();
      }
    } else {
      animating = false;
      cancelAnimationFrame(animId);
      // Stay at current rotation — don't reset
    }
  }

  function freeThingiPair() {
    thingiA?.handle.free();
    thingiB?.handle.free();
    thingiA = thingiB = null;
  }

  function setThingiInfo(html: string | null) {
    if (html) {
      thingiInfo.innerHTML = html;
      thingiInfo.style.display = 'block';
    } else {
      thingiInfo.style.display = 'none';
    }
  }

  function operandLabel(o: ThingiOperand): string {
    // Topological kind per the dataset flags (what the pools guarantee);
    // "soup" marks the rare case where even welding could not pair it.
    let kind = o.model.edge_manifold && o.model.vertex_manifold ? 'manifold' : 'non-manifold';
    if (o.handle.is_soup) kind += ', soup';
    return `#${o.model.id} &ldquo;${o.model.name}&rdquo; (${o.model.faces} tris, ${kind})`;
  }

  // Some dataset models flagged "closed" still fail the robust importer's
  // stricter closed check — skip those and re-roll, keeping a record of the
  // skips for the debug info.
  async function loadOperand(kind: 'm' | 'n', exclude?: ThingiModel): Promise<ThingiOperand> {
    const MAX_ATTEMPTS = 5;
    for (let attempt = 0; attempt < MAX_ATTEMPTS; attempt++) {
      const model = await pickRandomModel(kind, exclude);
      setThingiInfo(`Fetching #${model.id} &ldquo;${model.name}&rdquo;…`);
      const parsed = await fetchMesh(model);
      if (disposed) throw new Error('demo closed');
      const handle = importTriangleSoup(parsed.positions, parsed.indices);
      if (handle.ok) {
        return { model, handle, normalization: { center: parsed.center, scale: parsed.scale } };
      }
      skippedImports.push({ id: model.id, status: handle.status });
      handle.free();
    }
    throw new Error(`No importable ${kind === 'm' ? 'manifold' : 'non-manifold'} mesh found in ` +
      `${MAX_ATTEMPTS} attempts (skipped: ${skippedImports.map(s => `#${s.id} ${s.status}`).join(', ')})`);
  }

  async function loadRandomPair() {
    if (loadingPair) return;
    loadingPair = true;
    loadBtn.disabled = true;
    loadBtn.textContent = 'Loading pair…';
    errorBox.style.display = 'none';
    skippedImports = [];
    try {
      const [kindA, kindB] = pairOperandKinds(pairKind);
      const a = await loadOperand(kindA);
      let b: ThingiOperand;
      try {
        b = await loadOperand(kindB, kindA === kindB ? a.model : undefined);
      } catch (e) {
        a.handle.free();
        throw e;
      }
      freeThingiPair();
      thingiA = a;
      thingiB = b;
      // Random static orientation for B so every pair meets differently;
      // if Animate is on, rotation keeps advancing from here.
      rotX = Math.floor(Math.random() * 360);
      rotY = Math.floor(Math.random() * 360);
      rotZ = Math.floor(Math.random() * 360);
      setThingiInfo(`A: ${operandLabel(thingiA)}<br>B: ${operandLabel(thingiB)}<br>` +
        `Loaded rotation: [${rotX}, ${rotY}, ${rotZ}]&deg;`);
      update();
    } catch (e) {
      console.error('Thingi10K pair load failed:', e);
      setThingiInfo(null);
      errorBox.style.display = 'block';
      errorBox.innerHTML = `<strong>Failed to load Thingi10K pair:</strong> ${String(e)}`;
    } finally {
      loadingPair = false;
      loadBtn.disabled = false;
      loadBtn.textContent = 'Load Random Pair';
    }
  }

  async function copyDebugInfo() {
    const text = JSON.stringify(lastDebugInfo ?? { note: 'no operation recorded yet' }, null, 2);
    try {
      await navigator.clipboard.writeText(text);
      copyBtn.textContent = 'Copied!';
    } catch {
      // Clipboard API can be unavailable (non-secure context); fall back.
      const ta = document.createElement('textarea');
      ta.value = text;
      document.body.appendChild(ta);
      ta.select();
      document.execCommand('copy');
      ta.remove();
      copyBtn.textContent = 'Copied!';
    }
    setTimeout(() => { copyBtn.textContent = 'Copy Debug Info'; }, 1500);
  }

  // ---- Controls layout: Source → source-specific selectors → shared ----

  function setSource(s: 'builtin' | 'thingi') {
    source = s;
    saveSetting(DEMO, 'source', s);
    builtinSection.style.display = s === 'builtin' ? '' : 'none';
    thingiSection.style.display = s === 'thingi' ? '' : 'none';
    if (s === 'thingi') {
      if (inThingi()) {
        setThingiInfo(`A: ${operandLabel(thingiA!)}<br>B: ${operandLabel(thingiB!)}`);
        update();
      } else {
        loadRandomPair();
      }
    } else {
      update();
    }
  }

  controlsEl.appendChild(createDropdown('Source', SOURCES, source, v => setSource(v as 'builtin' | 'thingi')));

  const builtinSection = document.createElement('div');
  builtinSection.appendChild(createDropdown('Shape A', SHAPES, String(shapeA), v => { shapeA = parseInt(v); saveSetting(DEMO, 'shapeA', shapeA); update(); }));
  builtinSection.appendChild(createDropdown('Shape B', SHAPES, String(shapeB), v => { shapeB = parseInt(v); saveSetting(DEMO, 'shapeB', shapeB); update(); }));
  controlsEl.appendChild(builtinSection);

  const thingiSection = document.createElement('div');
  thingiSection.appendChild(createDropdown('Pair Type', PAIR_KINDS, pairKind, v => { pairKind = v as PairKind; saveSetting(DEMO, 'pairKind', pairKind); }));
  const loadBtn = createButton('Load Random Pair', () => { loadRandomPair(); });
  thingiSection.appendChild(loadBtn);
  thingiSection.appendChild(thingiInfo);
  controlsEl.appendChild(thingiSection);

  controlsEl.appendChild(createDropdown('Operation', OPS, String(op), v => { op = parseInt(v); saveSetting(DEMO, 'op', op); update(); }));
  controlsEl.appendChild(createDropdown('Engine', ENGINES, String(engine), v => { engine = parseInt(v) as BooleanEngine; saveSetting(DEMO, 'engine', engine); update(); }));
  controlsEl.appendChild(createSlider('Offset X ', -1.5, 1.5, offsetX, 0.1, v => { offsetX = v; saveSetting(DEMO, 'offsetX', v); update(); }));
  controlsEl.appendChild(createSlider('Offset Y ', -1.5, 1.5, offsetY, 0.1, v => { offsetY = v; saveSetting(DEMO, 'offsetY', v); update(); }));
  controlsEl.appendChild(createSlider('Offset Z ', -1.5, 1.5, offsetZ, 0.1, v => { offsetZ = v; saveSetting(DEMO, 'offsetZ', v); update(); }));
  const animateBox = createCheckbox('Animate', animate, toggleAnimate);
  controlsEl.appendChild(animateBox);
  controlsEl.appendChild(createCheckbox('Wireframe', wireframe, v => { saveSetting(DEMO, 'wireframe', v); viewer.setWireframe(v); }));
  const copyBtn = createButton('Copy Debug Info', () => { copyDebugInfo(); });
  controlsEl.appendChild(copyBtn);
  controlsEl.appendChild(errorBox);
  controlsEl.appendChild(readout);

  builtinSection.style.display = source === 'builtin' ? '' : 'none';
  thingiSection.style.display = source === 'thingi' ? '' : 'none';

  if (source === 'thingi') {
    loadRandomPair();
  } else {
    update();
  }
  if (wireframe) viewer.setWireframe(true);
  toggleAnimate(animate);

  return () => {
    disposed = true;
    animating = false;
    cancelAnimationFrame(animId);
    freeThingiPair();
    viewer.dispose();
  };
}
