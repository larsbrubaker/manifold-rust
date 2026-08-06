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

/// Rebuild a ThingiModel from a Copy-Debug-Info operand record: the URL
/// carries the repo and format, so a pasted report reproduces without the
/// metadata index (whose pools might have changed since the capture).
function modelFromDebugRecord(rec: any): ThingiModel {
  if (typeof rec?.id !== 'number' || typeof rec?.url !== 'string') {
    throw new Error('operand record is missing id/url');
  }
  const repoMatch = /Thingi10K-meshes-(\d+)@/.exec(rec.url);
  const formatMatch = /\.(\w+)\.zip$/.exec(rec.url);
  if (!repoMatch) throw new Error(`unrecognized mesh URL: ${rec.url}`);
  return {
    id: rec.id,
    thing_id: 0,
    name: rec.name ?? `#${rec.id}`,
    format: formatMatch ? formatMatch[1] : 'stl',
    repo: parseInt(repoMatch[1]),
    closed: true,
    edge_manifold: !!rec.edge_manifold,
    vertex_manifold: !!rec.vertex_manifold,
    faces: rec.faces ?? 0,
    vertices: 0,
  };
}

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
  // Wall time of the most recent boolean evaluation (not mesh upload), shown
  // in the info panel so slow pairs are visible at a glance.
  let lastFrameMs: number | null = null;

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

  function formatMs(ms: number): string {
    return ms >= 1000 ? `${(ms / 1000).toFixed(2)} s` : `${ms.toFixed(1)} ms`;
  }

  function showReadout(data: MeshData) {
    errorBox.style.display = 'none';
    const rows = [
      { label: 'Vertices', value: String(data.num_vert) },
      { label: 'Triangles', value: String(data.num_tri) },
      { label: 'Volume', value: data.volume.toFixed(4) },
      { label: 'Surface Area', value: data.surface_area.toFixed(4) },
    ];
    if (lastFrameMs !== null) {
      rows.push({ label: 'Last Frame', value: formatMs(lastFrameMs) });
    }
    updateReadout(readout, rows);
  }

  function update(silent = false) {
    if (source === 'thingi' && !inThingi()) return; // pair still loading
    const info = captureDebugInfo();
    try {
      let data: MeshData;
      const t0 = performance.now();
      if (inThingi()) {
        data = importedBoolean(thingiA!.handle, thingiB!.handle, op, engine, offsetX, offsetY, offsetZ, rotX, rotY, rotZ);
      } else {
        data = booleanGalleryMeshEngine(shapeA, shapeB, op, offsetX, offsetY, offsetZ, rotX, rotY, rotZ, engine);
      }
      lastFrameMs = performance.now() - t0;
      info.elapsed_ms = Math.round(lastFrameMs * 10) / 10;
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

  const sourceCtl = createDropdown('Source', SOURCES, source, v => setSource(v as 'builtin' | 'thingi'));
  controlsEl.appendChild(sourceCtl);

  const builtinSection = document.createElement('div');
  const shapeACtl = createDropdown('Shape A', SHAPES, String(shapeA), v => { shapeA = parseInt(v); saveSetting(DEMO, 'shapeA', shapeA); update(); });
  const shapeBCtl = createDropdown('Shape B', SHAPES, String(shapeB), v => { shapeB = parseInt(v); saveSetting(DEMO, 'shapeB', shapeB); update(); });
  builtinSection.appendChild(shapeACtl);
  builtinSection.appendChild(shapeBCtl);
  controlsEl.appendChild(builtinSection);

  const thingiSection = document.createElement('div');
  const pairKindCtl = createDropdown('Pair Type', PAIR_KINDS, pairKind, v => { pairKind = v as PairKind; saveSetting(DEMO, 'pairKind', pairKind); });
  thingiSection.appendChild(pairKindCtl);
  const loadBtn = createButton('Load Random Pair', () => { loadRandomPair(); });
  thingiSection.appendChild(loadBtn);
  thingiSection.appendChild(thingiInfo);
  controlsEl.appendChild(thingiSection);

  const opCtl = createDropdown('Operation', OPS, String(op), v => { op = parseInt(v); saveSetting(DEMO, 'op', op); update(); });
  const engineCtl = createDropdown('Engine', ENGINES, String(engine), v => { engine = parseInt(v) as BooleanEngine; saveSetting(DEMO, 'engine', engine); update(); });
  const offXCtl = createSlider('Offset X ', -1.5, 1.5, offsetX, 0.1, v => { offsetX = v; saveSetting(DEMO, 'offsetX', v); update(); });
  const offYCtl = createSlider('Offset Y ', -1.5, 1.5, offsetY, 0.1, v => { offsetY = v; saveSetting(DEMO, 'offsetY', v); update(); });
  const offZCtl = createSlider('Offset Z ', -1.5, 1.5, offsetZ, 0.1, v => { offsetZ = v; saveSetting(DEMO, 'offsetZ', v); update(); });
  controlsEl.appendChild(opCtl);
  controlsEl.appendChild(engineCtl);
  controlsEl.appendChild(offXCtl);
  controlsEl.appendChild(offYCtl);
  controlsEl.appendChild(offZCtl);
  const animateBox = createCheckbox('Animate', animate, toggleAnimate);
  controlsEl.appendChild(animateBox);
  controlsEl.appendChild(createCheckbox('Wireframe', wireframe, v => { saveSetting(DEMO, 'wireframe', v); viewer.setWireframe(v); }));
  const copyBtn = createButton('Copy Debug Info', () => { copyDebugInfo(); });
  controlsEl.appendChild(copyBtn);
  const useBtn = createButton('Use Debug Info', () => { useDebugInfo(); });
  controlsEl.appendChild(useBtn);
  controlsEl.appendChild(errorBox);
  controlsEl.appendChild(readout);

  // ---- Use Debug Info: restore the exact captured state from a paste ----

  function setDropdownValue(ctl: HTMLElement, value: string) {
    const sel = ctl.querySelector('select');
    if (sel) sel.value = value;
  }

  function setSliderValue(ctl: HTMLElement, v: number) {
    const input = ctl.querySelector('input');
    if (input) input.value = String(v);
    const span = ctl.querySelector('.slider-value');
    if (span) span.textContent = String(v);
  }

  function forceAnimateOff() {
    if (animate || animating) toggleAnimate(false);
    const box = animateBox.querySelector('input') as HTMLInputElement | null;
    if (box) box.checked = false;
  }

  function showUseError(msg: string) {
    errorBox.style.display = 'block';
    errorBox.innerHTML = `<strong>Use Debug Info:</strong> ${msg}`;
  }

  async function loadOperandFromRecord(rec: any): Promise<ThingiOperand> {
    const model = modelFromDebugRecord(rec);
    setThingiInfo(`Fetching #${model.id} &ldquo;${model.name}&rdquo;…`);
    const parsed = await fetchMesh(model);
    if (disposed) throw new Error('demo closed');
    const handle = importTriangleSoup(parsed.positions, parsed.indices);
    if (!handle.ok) {
      const status = handle.status;
      handle.free();
      throw new Error(`#${model.id} failed to import: ${status}`);
    }
    return { model, handle, normalization: { center: parsed.center, scale: parsed.scale } };
  }

  async function useDebugInfo() {
    errorBox.style.display = 'none';
    let text: string;
    try {
      text = await navigator.clipboard.readText();
    } catch (e) {
      showUseError(`could not read the clipboard (${String(e)}). The browser may need ` +
        `clipboard permission, or this page must run on a secure context (https or localhost).`);
      return;
    }
    let info: any;
    try {
      info = JSON.parse(text);
    } catch {
      showUseError('the clipboard does not contain debug info (not valid JSON). ' +
        'Press <em>Copy Debug Info</em> on a gallery frame first, then paste it back here.');
      return;
    }
    if (!info || info.demo !== DEMO || !Array.isArray(info.rotation_degrees) || !Array.isArray(info.offset)) {
      showUseError('the clipboard JSON is not Boolean Gallery debug info ' +
        '(expected the object produced by <em>Copy Debug Info</em>).');
      return;
    }
    if (loadingPair) {
      showUseError('a Thingi10K pair is still loading — try again when it finishes.');
      return;
    }

    // Land exactly on the captured frame: freeze animation first.
    forceAnimateOff();

    const opIdx = OP_NAMES.indexOf(info.op);
    if (opIdx >= 0) { op = opIdx; saveSetting(DEMO, 'op', op); setDropdownValue(opCtl, String(op)); }
    const engIdx = ENGINE_NAMES.indexOf(info.engine);
    if (engIdx >= 0) { engine = engIdx as BooleanEngine; saveSetting(DEMO, 'engine', engine); setDropdownValue(engineCtl, String(engine)); }
    [offsetX, offsetY, offsetZ] = info.offset.map(Number) as [number, number, number];
    saveSetting(DEMO, 'offsetX', offsetX);
    saveSetting(DEMO, 'offsetY', offsetY);
    saveSetting(DEMO, 'offsetZ', offsetZ);
    setSliderValue(offXCtl, offsetX);
    setSliderValue(offYCtl, offsetY);
    setSliderValue(offZCtl, offsetZ);
    [rotX, rotY, rotZ] = info.rotation_degrees.map(Number) as [number, number, number];

    if (info.model_a && info.model_b) {
      source = 'thingi';
      saveSetting(DEMO, 'source', source);
      setDropdownValue(sourceCtl, source);
      builtinSection.style.display = 'none';
      thingiSection.style.display = '';
      if (typeof info.pair_kind === 'string') {
        pairKind = info.pair_kind as PairKind;
        saveSetting(DEMO, 'pairKind', pairKind);
        setDropdownValue(pairKindCtl, pairKind);
      }
      loadingPair = true;
      loadBtn.disabled = true;
      useBtn.disabled = true;
      useBtn.textContent = 'Loading pair…';
      skippedImports = [];
      try {
        const a = await loadOperandFromRecord(info.model_a);
        let b: ThingiOperand;
        try {
          b = await loadOperandFromRecord(info.model_b);
        } catch (e) {
          a.handle.free();
          throw e;
        }
        freeThingiPair();
        thingiA = a;
        thingiB = b;
        setThingiInfo(`A: ${operandLabel(thingiA)}<br>B: ${operandLabel(thingiB)}<br>` +
          `Restored rotation: [${rotX}, ${rotY}, ${rotZ}]&deg;`);
        update();
      } catch (e) {
        setThingiInfo(inThingi() ? `A: ${operandLabel(thingiA!)}<br>B: ${operandLabel(thingiB!)}` : null);
        showUseError(`failed to load the captured pair: ${String(e)}`);
      } finally {
        loadingPair = false;
        loadBtn.disabled = false;
        useBtn.disabled = false;
        useBtn.textContent = 'Use Debug Info';
      }
    } else {
      source = 'builtin';
      saveSetting(DEMO, 'source', source);
      setDropdownValue(sourceCtl, source);
      builtinSection.style.display = '';
      thingiSection.style.display = 'none';
      const sa = SHAPE_NAMES.indexOf(info.shape_a);
      const sb = SHAPE_NAMES.indexOf(info.shape_b);
      if (sa >= 0) { shapeA = sa; saveSetting(DEMO, 'shapeA', shapeA); setDropdownValue(shapeACtl, String(shapeA)); }
      if (sb >= 0) { shapeB = sb; saveSetting(DEMO, 'shapeB', shapeB); setDropdownValue(shapeBCtl, String(shapeB)); }
      update();
    }
  }

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
