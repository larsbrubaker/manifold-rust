// Boolean Gallery: booleans on either built-in primitives or random
// Thingi10K mesh pairs, with a choice of boolean engine (Exact / Robust /
// Auto). A Source dropdown switches between the two mesh sources; every
// other control — operation, engine, offsets, animate, wireframe — applies
// to whichever source is active. Every operation's inputs are captured for
// the Copy Debug Info button so failures can be reproduced from a pasted
// report.

import { ThreeViewer } from '../three-viewer.ts';
import { createSlider, createDropdown, createCheckbox, createButton, createNumberInput, createReadout, updateReadout } from '../controls.ts';
import { type MeshData, type BooleanEngine } from '../wasm.ts';
import { BooleanRunner, type ImportInfo, type OperandArrays, type RunParams } from '../boolean-runner.ts';
import { loadSetting, saveSetting } from '../settings.ts';
import { pickRandomModel, pairOperandKinds, fetchMesh, meshZipUrl, findModelById, type ThingiModel, type PairKind } from '../thingi.ts';

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

/// Pair Type is really a mode selector: the three kinds above roll a random
/// pair, while the two below reveal a sub-control that names an exact pair
/// instead — a curated trouble case, or model ids typed in by hand.
type PairMode = PairKind | 'trouble' | 'pick';

const PAIR_MODES = PAIR_KINDS.concat([
  { value: 'trouble', text: 'Trouble Cases…' },
  { value: 'pick', text: 'Pick Models…' },
]);

const isRandomMode = (m: PairMode): m is PairKind =>
  m === 'mm' || m === 'mn' || m === 'nn';

/// Known Thingi10K trouble cases — the pairs the perf and bug hunts keep
/// returning to. Selecting one loads the exact configuration (models by id,
/// op, engine, offsets, rotation) through the same path as Use Debug Info.
const TROUBLE_CASES: {
  value: string; text: string; a: number; b: number;
  op: string; engine: string; offset: [number, number, number];
  rot: [number, number, number]; pairKind: PairKind;
}[] = [
  {
    value: '1663774-51334',
    text: '1663774 ∪ 51334 — heavy fins (slowest known)',
    a: 1663774, b: 51334, op: 'union', engine: 'robust',
    offset: [0.3, 0, 0], rot: [231.39999999999753, 124, 273.6000000000049], pairKind: 'nn',
  },
  {
    value: '91946-61459',
    text: '91946 ∪ 61459 — doubled-surface windows (perf)',
    a: 91946, b: 61459, op: 'union', engine: 'robust',
    offset: [0.3, 0, 0], rot: [236, 231, 42], pairKind: 'nn',
  },
  {
    value: '939888-93557',
    text: '939888 ∪ 93557 — self-overlapping operand',
    a: 939888, b: 93557, op: 'union', engine: 'robust',
    offset: [0.3, 0, 0], rot: [356, 140, 322], pairKind: 'nn',
  },
  {
    value: '1075458-91115',
    text: '1075458 − 91115 — CDT constraint recovery (fixed)',
    a: 1075458, b: 91115, op: 'difference', engine: 'robust',
    offset: [0.7, -0.2, 0.4], rot: [311, 55, 345], pairKind: 'mn',
  },
  {
    value: '92068-39926',
    text: '92068 ∪ 39926 — tripled facets (fixed NotClosed)',
    a: 92068, b: 39926, op: 'union', engine: 'robust',
    offset: [0.3, 0, 0], rot: [0, 0, 0], pairKind: 'nn',
  },
];

const OP_COLORS = [0x4488cc, 0x44aa44, 0xcc4444];
const OP_NAMES = ['union', 'intersection', 'difference'];
const ENGINE_NAMES = ['exact', 'robust', 'auto'];
const SHAPE_NAMES = ['cube', 'sphere', 'cylinder', 'spiky-dodecahedron'];

interface ThingiOperand {
  model: ThingiModel;
  /** Raw parsed soup, kept so a cancelled worker can be re-fed. */
  arrays: OperandArrays;
  /** Import result reported by the worker's wasm instance. */
  info: ImportInfo;
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
  let repairOrientation = loadSetting(DEMO, 'repairOrientation', false);
  let offsetX = loadSetting(DEMO, 'offsetX', 0.3);
  let offsetY = loadSetting(DEMO, 'offsetY', 0.0);
  let offsetZ = loadSetting(DEMO, 'offsetZ', 0.0);
  let wireframe = loadSetting(DEMO, 'wireframe', false);
  let xray = loadSetting(DEMO, 'xray', false);
  let animate = loadSetting(DEMO, 'animate', true);
  let animating = false;
  let animId = 0;
  let rotX = 0, rotY = 0, rotZ = 0;
  const ROT_SPEED_X = 0.7;  // degrees per frame
  const ROT_SPEED_Y = 1.5;
  const ROT_SPEED_Z = 0.3;

  // Thingi10K state
  let pairKind = loadSetting(DEMO, 'pairKind', 'mn') as PairKind;
  // Mode defaults to the saved random kind, so existing users land where they
  // left off. Pick Models ids: 0 means "empty".
  let pairMode = loadSetting(DEMO, 'pairMode', pairKind) as PairMode;
  let pickA = loadSetting(DEMO, 'pickA', 0);
  let pickB = loadSetting(DEMO, 'pickB', 0);
  let troubleCase = loadSetting(DEMO, 'troubleCase', '');
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
      imported: { num_tri: o.info.num_tri, is_soup: o.info.is_soup, status: o.info.status },
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
      repair_orientation: repairOrientation,
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

  // ---- Worker-driven evaluation: the boolean runs off the main thread ----

  const runner = new BooleanRunner();

  function update(silent = false) {
    if (source === 'thingi' && (!inThingi() || loadingPair)) return; // pair still loading
    const info = captureDebugInfo();
    const params: RunParams = inThingi()
      ? {
          source: 'thingi', op, engine, repair: repairOrientation,
          ox: offsetX, oy: offsetY, oz: offsetZ, rx: rotX, ry: rotY, rz: rotZ,
          tag: { info, silent },
        }
      : {
          source: 'builtin', shapeA, shapeB, op, engine, repair: repairOrientation,
          ox: offsetX, oy: offsetY, oz: offsetZ, rx: rotX, ry: rotY, rz: rotZ,
          tag: { info, silent },
        };
    runner.requestRun(params);
  }

  runner.onResult = ({ data, elapsedMs }, params) => {
    const { info } = params.tag;
    lastFrameMs = elapsedMs;
    info.elapsed_ms = Math.round(elapsedMs * 10) / 10;
    info.result = { num_tri: data.num_tri, num_vert: data.num_vert, volume: data.volume };
    try { localStorage.setItem('boolean-gallery-debug', JSON.stringify(info)); } catch { /* ignore */ }
    viewer.setMesh(data as MeshData);
    viewer.setColor(OP_COLORS[params.op] || 0x4488cc);
    errorBox.style.display = 'none';
    showReadout(data as MeshData);
    scheduleNextAnimateFrame(elapsedMs);
  };

  runner.onError = (error, params) => {
    const { info, silent } = params.tag;
    info.error = error;
    try { localStorage.setItem('boolean-gallery-debug', JSON.stringify(info)); } catch { /* ignore */ }
    console.error('Boolean op failed:', info, error);
    if (!silent) {
      errorBox.style.display = 'block';
      errorBox.innerHTML = `<strong>Boolean operation failed:</strong> ${error}<br>` +
        `Use <em>Copy Debug Info</em> to capture the inputs for a bug report.`;
      updateReadout(readout, []);
    }
    scheduleNextAnimateFrame(0);
  };

  // With the boolean off the main thread the tab never freezes, but pairs
  // that take seconds per frame still make a useless slideshow — auto-pause
  // and let the user re-enable Animate if they really want it.
  const ANIMATE_BUDGET_MS = 1500;

  function scheduleNextAnimateFrame(elapsedMs: number) {
    if (!animating) return;
    if (elapsedMs > ANIMATE_BUDGET_MS) {
      toggleAnimate(false);
      const box = animateBox.querySelector('input') as HTMLInputElement | null;
      if (box) box.checked = false;
      errorBox.style.display = 'block';
      errorBox.innerHTML = `<strong>Animation paused:</strong> this boolean takes ` +
        `${(elapsedMs / 1000).toFixed(1)} s per frame. Re-check Animate to continue anyway.`;
      return;
    }
    animId = requestAnimationFrame(animateStep);
  }

  function animateStep() {
    if (!animating) return;
    rotX = (rotX + ROT_SPEED_X) % 360;
    rotY = (rotY + ROT_SPEED_Y) % 360;
    rotZ = (rotZ + ROT_SPEED_Z) % 360;
    // The next frame is scheduled when this one's result arrives.
    update(true);
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

  // The worker owns the wasm-side meshes; the main thread only drops its
  // references. A stale worker slot is overwritten by the next import.
  function freeThingiPair() {
    thingiA = thingiB = null;
  }

  /// Re-feed the committed pair to the worker — after a failed load attempt
  /// or a cancel replaced worker slots with candidate meshes, this restores
  /// the state the on-screen labels describe.
  function restoreWorkerOperands() {
    if (thingiA) runner.importOperand('a', thingiA.arrays);
    if (thingiB) runner.importOperand('b', thingiB.arrays);
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
    if (o.info.is_soup) kind += ', soup';
    return `#${o.model.id} &ldquo;${o.model.name}&rdquo; (${o.model.faces} tris, ${kind})`;
  }

  // Some dataset models flagged "closed" still fail the robust importer's
  // stricter closed check — skip those and re-roll, keeping a record of the
  // skips for the debug info.
  async function loadOperand(kind: 'm' | 'n', slot: 'a' | 'b', exclude?: ThingiModel): Promise<ThingiOperand> {
    const MAX_ATTEMPTS = 5;
    for (let attempt = 0; attempt < MAX_ATTEMPTS; attempt++) {
      const model = await pickRandomModel(kind, exclude);
      setThingiInfo(`Fetching #${model.id} &ldquo;${model.name}&rdquo;…`);
      const parsed = await fetchMesh(model);
      if (disposed) throw new Error('demo closed');
      const arrays: OperandArrays = { positions: parsed.positions, indices: parsed.indices };
      const info = await runner.importOperand(slot, arrays);
      if (info.ok) {
        return { model, arrays, info, normalization: { center: parsed.center, scale: parsed.scale } };
      }
      skippedImports.push({ id: model.id, status: info.status });
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
    let loaded = false;
    try {
      const [kindA, kindB] = pairOperandKinds(pairKind);
      const a = await loadOperand(kindA, 'a');
      const b = await loadOperand(kindB, 'b', kindA === kindB ? a.model : undefined);
      freeThingiPair();
      thingiA = a;
      thingiB = b;
      loaded = true;
      // Random static orientation for B so every pair meets differently;
      // if Animate is on, rotation keeps advancing from here.
      rotX = Math.floor(Math.random() * 360);
      rotY = Math.floor(Math.random() * 360);
      rotZ = Math.floor(Math.random() * 360);
      setThingiInfo(`A: ${operandLabel(thingiA)}<br>B: ${operandLabel(thingiB)}<br>` +
        `Loaded rotation: [${rotX}, ${rotY}, ${rotZ}]&deg;`);
    } catch (e) {
      console.error('Thingi10K pair load failed:', e);
      restoreWorkerOperands();
      setThingiInfo(inThingi()
        ? `<strong>Load failed — showing previous pair:</strong><br>` +
          `A: ${operandLabel(thingiA!)}<br>B: ${operandLabel(thingiB!)}`
        : null);
      errorBox.style.display = 'block';
      errorBox.innerHTML = `<strong>Failed to load Thingi10K pair:</strong> ${String(e)}`;
    } finally {
      loadingPair = false;
      loadBtn.disabled = false;
      loadBtn.textContent = 'Load Random Pair';
    }
    // Re-render only on success (update() is a no-op while a pair loads).
    // On failure the previous result stays on screen and — crucially — the
    // error stays visible: a successful update() would hide the error box.
    if (loaded && inThingi()) update();
  }

  /// First pair for a fresh Thingi session, honoring the saved Pair Type:
  /// a saved Pick Models pair, the saved Trouble Case, or a random roll of
  /// the saved kind. Never silently rolls a random pair while the visible
  /// mode says otherwise.
  function loadInitialPair() {
    if (pairMode === 'pick' && pickA && pickB) {
      loadPickedModels();
    } else if (pairMode === 'trouble') {
      if (troubleCase && TROUBLE_CASES.some(c => c.value === troubleCase)) {
        loadTroubleCase(troubleCase);
      } else {
        setThingiInfo('Pick a known case from the <em>Trouble Cases</em> list.');
      }
    } else if (pairMode === 'pick') {
      setThingiInfo('Enter two model ids and press <em>Load These Models</em>.');
    } else {
      loadRandomPair();
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
        loadInitialPair();
      }
    } else {
      update();
    }
  }

  /// Switch the Thingi10K pair mode, showing only that mode's sub-control:
  /// Load Random Pair for the three random kinds, the curated list for
  /// Trouble Cases, the id inputs for Pick Models.
  ///
  /// `userInitiated` marks a change made through the dropdown: picking a
  /// random kind then immediately rolls a pair of that kind, so the display
  /// always matches what the Pair Type says (previously the old pair — of a
  /// different kind — stayed on screen until Load Random Pair was pressed).
  function setPairMode(mode: PairMode, userInitiated = false) {
    const changed = mode !== pairMode;
    pairMode = mode;
    saveSetting(DEMO, 'pairMode', mode);
    const random = isRandomMode(mode);
    if (random) {
      pairKind = mode;
      saveSetting(DEMO, 'pairKind', pairKind);
    }
    loadBtn.style.display = random ? '' : 'none';
    troubleCtl.style.display = mode === 'trouble' ? '' : 'none';
    pickSection.style.display = mode === 'pick' ? '' : 'none';
    if (userInitiated && changed && source === 'thingi') {
      if (random) {
        loadRandomPair();
      } else if (mode === 'trouble' && troubleCase && TROUBLE_CASES.some(c => c.value === troubleCase)) {
        loadTroubleCase(troubleCase);
      }
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
  const pairModeCtl = createDropdown('Pair Type', PAIR_MODES, pairMode, v => setPairMode(v as PairMode, true));
  thingiSection.appendChild(pairModeCtl);
  const loadBtn = createButton('Load Random Pair', () => { loadRandomPair(); });
  loadBtn.style.marginTop = '10px';
  thingiSection.appendChild(loadBtn);
  const troubleOptions = [{ value: '', text: 'Pick a known case…' }]
    .concat(TROUBLE_CASES.map(c => ({ value: c.value, text: c.text })));
  const troubleCtl = createDropdown('Trouble Cases', troubleOptions,
    TROUBLE_CASES.some(c => c.value === troubleCase) ? troubleCase : '', v => {
    if (v) {
      troubleCase = v;
      saveSetting(DEMO, 'troubleCase', v);
      loadTroubleCase(v);
    }
  });
  thingiSection.appendChild(troubleCtl);
  const pickSection = document.createElement('div');
  pickSection.className = 'control-substack';
  const pickACtl = createNumberInput('Model A id', pickA || null, v => {
    pickA = v ?? 0;
    saveSetting(DEMO, 'pickA', pickA);
  });
  const pickBCtl = createNumberInput('Model B id', pickB || null, v => {
    pickB = v ?? 0;
    saveSetting(DEMO, 'pickB', pickB);
  });
  const pickBtn = createButton('Load These Models', () => { loadPickedModels(); });
  pickSection.appendChild(pickACtl);
  pickSection.appendChild(pickBCtl);
  pickSection.appendChild(pickBtn);
  thingiSection.appendChild(pickSection);
  thingiSection.appendChild(thingiInfo);
  controlsEl.appendChild(thingiSection);

  // Engine first: it is the major choice (which pipeline runs at all), and
  // the operation is a detail within it.
  const engineCtl = createDropdown('Engine', ENGINES, String(engine), v => { engine = parseInt(v) as BooleanEngine; saveSetting(DEMO, 'engine', engine); update(); });
  // Repair sits with the engine choice: it changes what the operands *mean*
  // (inside-out bodies become solid material) before any boolean runs.
  const repairBox = createCheckbox('Repair orientation', repairOrientation, v => {
    repairOrientation = v;
    saveSetting(DEMO, 'repairOrientation', v);
    update();
  });
  const opCtl = createDropdown('Operation', OPS, String(op), v => { op = parseInt(v); saveSetting(DEMO, 'op', op); update(); });
  const offXCtl = createSlider('Offset X ', -1.5, 1.5, offsetX, 0.1, v => { offsetX = v; saveSetting(DEMO, 'offsetX', v); update(); });
  const offYCtl = createSlider('Offset Y ', -1.5, 1.5, offsetY, 0.1, v => { offsetY = v; saveSetting(DEMO, 'offsetY', v); update(); });
  const offZCtl = createSlider('Offset Z ', -1.5, 1.5, offsetZ, 0.1, v => { offsetZ = v; saveSetting(DEMO, 'offsetZ', v); update(); });
  controlsEl.appendChild(engineCtl);
  controlsEl.appendChild(repairBox);
  controlsEl.appendChild(opCtl);
  controlsEl.appendChild(offXCtl);
  controlsEl.appendChild(offYCtl);
  controlsEl.appendChild(offZCtl);
  const animateBox = createCheckbox('Animate', animate, toggleAnimate);
  controlsEl.appendChild(animateBox);
  controlsEl.appendChild(createCheckbox('Wireframe', wireframe, v => { saveSetting(DEMO, 'wireframe', v); viewer.setWireframe(v); }));
  controlsEl.appendChild(createCheckbox('X-Ray (depth peeled)', xray, v => { xray = v; saveSetting(DEMO, 'xray', v); viewer.setXRay(v); }));
  const copyBtn = createButton('Copy Debug Info', () => { copyDebugInfo(); });
  controlsEl.appendChild(copyBtn);
  const useBtn = createButton('Use Debug Info', () => { useDebugInfo(); });
  controlsEl.appendChild(useBtn);

  // Cancel appears only when a boolean has been computing for a moment
  // (avoids flicker on fast frames). Cancelling kills the worker outright —
  // the runner respawns it and re-imports the operands.
  const cancelBtn = createButton('Cancel Computation', () => {
    runner.cancel();
  });
  cancelBtn.style.display = 'none';
  cancelBtn.style.marginTop = '10px';
  controlsEl.appendChild(cancelBtn);
  let cancelShowTimer = 0;
  runner.onBusyChange = busy => {
    clearTimeout(cancelShowTimer);
    if (busy) {
      cancelShowTimer = window.setTimeout(() => { cancelBtn.style.display = ''; }, 400);
    } else {
      cancelBtn.style.display = 'none';
    }
  };
  runner.onCancelled = () => {
    errorBox.style.display = 'block';
    errorBox.innerHTML = '<strong>Computation cancelled.</strong> The previous result stays on screen.';
  };

  controlsEl.appendChild(errorBox);
  controlsEl.appendChild(readout);

  // ---- Use Debug Info: restore the exact captured state from a paste ----

  function setDropdownValue(ctl: HTMLElement, value: string) {
    const sel = ctl.querySelector('select');
    if (sel) sel.value = value;
  }

  function setCheckboxValue(ctl: HTMLElement, v: boolean) {
    const input = ctl.querySelector('input');
    if (input) input.checked = v;
  }

  function setSliderValue(ctl: HTMLElement, v: number) {
    const input = ctl.querySelector('input');
    if (input) input.value = String(v);
    const span = ctl.querySelector('.slider-value');
    if (span) span.textContent = String(v);
  }

  function setNumberValue(ctl: HTMLElement, v: number) {
    const input = ctl.querySelector('input');
    if (input) input.value = v ? String(v) : '';
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

  async function loadOperandFromRecord(rec: any, slot: 'a' | 'b'): Promise<ThingiOperand> {
    const model = modelFromDebugRecord(rec);
    setThingiInfo(`Fetching #${model.id} &ldquo;${model.name}&rdquo;…`);
    const parsed = await fetchMesh(model);
    if (disposed) throw new Error('demo closed');
    const arrays: OperandArrays = { positions: parsed.positions, indices: parsed.indices };
    const info = await runner.importOperand(slot, arrays);
    if (!info.ok) {
      throw new Error(`#${model.id} failed to import: ${info.status}`);
    }
    return { model, arrays, info, normalization: { center: parsed.center, scale: parsed.scale } };
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
    await applyDebugInfo(info);
  }

  /// Resolve two model ids in the metadata index and restore the gallery to
  /// that pair through the shared debug-info path. Backs both the Trouble
  /// Cases list and the Pick Models id inputs.
  async function loadPairByIds(
    aId: number, bId: number,
    cfg: {
      op: string; engine: string;
      offset: [number, number, number]; rot: [number, number, number];
      pairKind?: PairKind;
    },
    modeAfter: PairMode,
  ) {
    if (loadingPair) return;
    errorBox.style.display = 'none';
    // The index gives repo/name/format, which the shared restore path needs
    // to build a fetchable URL.
    let ma: ThingiModel | null = null;
    let mb: ThingiModel | null = null;
    try {
      [ma, mb] = await Promise.all([findModelById(aId), findModelById(bId)]);
    } catch (e) {
      showUseError(`could not load the Thingi10K index: ${String(e)}`);
      return;
    }
    if (!ma || !mb) {
      showUseError(`model #${!ma ? aId : bId} is not in the Thingi10K index ` +
        `(it must be an STL entry of the Thingi10K subset this demo mirrors).`);
      return;
    }
    const rec = (m: ThingiModel) => ({
      id: m.id, name: m.name, url: meshZipUrl(m), faces: m.faces,
      edge_manifold: m.edge_manifold, vertex_manifold: m.vertex_manifold,
    });
    await applyDebugInfo({
      demo: DEMO,
      op: cfg.op,
      engine: cfg.engine,
      offset: cfg.offset,
      rotation_degrees: cfg.rot,
      pair_kind: cfg.pairKind,
      model_a: rec(ma),
      model_b: rec(mb),
    }, modeAfter);
  }

  async function loadTroubleCase(value: string) {
    const c = TROUBLE_CASES.find(tc => tc.value === value);
    if (!c) return;
    await loadPairByIds(
      c.a, c.b,
      { op: c.op, engine: c.engine, offset: c.offset, rot: c.rot, pairKind: c.pairKind },
      'trouble',
    );
  }

  /// Load the two ids typed into Pick Models, keeping the current op, engine,
  /// offsets and rotation — so a pair can be swapped under a fixed setup.
  async function loadPickedModels() {
    if (!pickA || !pickB) {
      showUseError('enter both a Model A id and a Model B id first.');
      return;
    }
    await loadPairByIds(
      pickA, pickB,
      {
        op: OP_NAMES[op] || String(op),
        engine: ENGINE_NAMES[engine],
        offset: [offsetX, offsetY, offsetZ],
        rot: [rotX, rotY, rotZ],
      },
      'pick',
    );
  }

  /// Restore the gallery to the exact state a debug-info record describes.
  /// Shared by Use Debug Info (clipboard), the Trouble Cases dropdown and the
  /// Pick Models inputs. `modeAfter` is the Pair Type the restored state
  /// belongs to: a pasted report lands in Pick Models (its ids are now the
  /// typed-in pair), while the curated list stays on Trouble Cases.
  async function applyDebugInfo(info: any, modeAfter: PairMode = 'pick') {
    // Land exactly on the captured frame: freeze animation first.
    forceAnimateOff();

    const opIdx = OP_NAMES.indexOf(info.op);
    if (opIdx >= 0) { op = opIdx; saveSetting(DEMO, 'op', op); setDropdownValue(opCtl, String(op)); }
    const engIdx = ENGINE_NAMES.indexOf(info.engine);
    if (engIdx >= 0) { engine = engIdx as BooleanEngine; saveSetting(DEMO, 'engine', engine); setDropdownValue(engineCtl, String(engine)); }
    if (typeof info.repair_orientation === 'boolean') {
      repairOrientation = info.repair_orientation;
      saveSetting(DEMO, 'repairOrientation', repairOrientation);
      setCheckboxValue(repairBox, repairOrientation);
    }
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
      // Keep the random-pair kind in sync when the record names one, but the
      // visible mode is whatever brought us here.
      if (typeof info.pair_kind === 'string' && isRandomMode(info.pair_kind as PairMode)) {
        pairKind = info.pair_kind as PairKind;
        saveSetting(DEMO, 'pairKind', pairKind);
      }
      // Surface the restored pair in the Pick Models inputs so it can be
      // re-run, tweaked one id at a time, or read off at a glance.
      if (typeof info.model_a.id === 'number' && typeof info.model_b.id === 'number') {
        pickA = info.model_a.id;
        pickB = info.model_b.id;
        saveSetting(DEMO, 'pickA', pickA);
        saveSetting(DEMO, 'pickB', pickB);
        setNumberValue(pickACtl, pickA);
        setNumberValue(pickBCtl, pickB);
      }
      setPairMode(modeAfter);
      setDropdownValue(pairModeCtl, modeAfter);
      loadingPair = true;
      loadBtn.disabled = true;
      useBtn.disabled = true;
      useBtn.textContent = 'Loading pair…';
      skippedImports = [];
      let loaded = false;
      try {
        const a = await loadOperandFromRecord(info.model_a, 'a');
        const b = await loadOperandFromRecord(info.model_b, 'b');
        freeThingiPair();
        thingiA = a;
        thingiB = b;
        loaded = true;
        setThingiInfo(`A: ${operandLabel(thingiA)}<br>B: ${operandLabel(thingiB)}<br>` +
          `Restored rotation: [${rotX}, ${rotY}, ${rotZ}]&deg;`);
      } catch (e) {
        restoreWorkerOperands();
        setThingiInfo(inThingi()
          ? `<strong>Load failed — showing previous pair:</strong><br>` +
            `A: ${operandLabel(thingiA!)}<br>B: ${operandLabel(thingiB!)}`
          : null);
        showUseError(`failed to load the requested pair: ${String(e)}`);
      } finally {
        loadingPair = false;
        loadBtn.disabled = false;
        useBtn.disabled = false;
        useBtn.textContent = 'Use Debug Info';
      }
      // Re-render only on success (update() is a no-op while a pair loads).
      // On failure the previous result stays on screen and — crucially — the
      // error stays visible: a successful update() would hide the error box.
      if (loaded && inThingi()) update();
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
  setPairMode(pairMode);

  if (source === 'thingi') {
    loadInitialPair();
  } else {
    update();
  }
  if (wireframe) viewer.setWireframe(true);
  if (xray) viewer.setXRay(true);
  toggleAnimate(animate);

  return () => {
    disposed = true;
    animating = false;
    cancelAnimationFrame(animId);
    freeThingiPair();
    runner.dispose();
    viewer.dispose();
  };
}
