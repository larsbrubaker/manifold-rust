// Boolean Gallery: booleans on either built-in primitives or random
// Thingi10K mesh pairs, with a choice of boolean engine (Exact / Robust /
// Auto). The sidebar (boolean-panel.ts) is banded Input → Operation → View →
// Result; this file owns the state behind it, the worker-driven evaluation
// loop, and the debug-info capture/restore path, so any failure can be
// reproduced from a pasted report.

import { ThreeViewer } from '../three-viewer.ts';
import { BooleanReadout } from './boolean-readout.ts';
import { buildPanel, type GalleryPanel, type Source } from './boolean-panel.ts';
import { type MeshData, type BooleanEngine } from '../wasm.ts';
import { BooleanRunner, type ImportInfo, type OperandArrays, type RunParams } from '../boolean-runner.ts';
import { loadSetting, saveSetting } from '../settings.ts';
import { pickRandomModel, pairOperandKinds, fetchMesh, meshZipUrl, findModelById, type ThingiModel, type PairKind } from '../thingi.ts';
import {
  ENGINE_NAMES, OP_COLORS, OP_NAMES, SHAPE_NAMES, TROUBLE_CASES,
  isRandomMode, modelFromDebugRecord, type PairMode,
} from './boolean-gallery-data.ts';

const DEMO = 'boolean-gallery';

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
      <div class="demo-layout bool-layout">
        <div class="demo-canvas-area" id="viewer-container"></div>
        <div id="controls"></div>
      </div>
    </div>
  `;

  const viewerEl = document.getElementById('viewer-container')!;
  const controlsEl = document.getElementById('controls')!;
  const viewer = new ThreeViewer(viewerEl);

  let source = loadSetting(DEMO, 'source', 'builtin') as Source;
  let shapeA = loadSetting(DEMO, 'shapeA', 3);
  let shapeB = loadSetting(DEMO, 'shapeB', 3);
  let op = loadSetting(DEMO, 'op', 0);
  let engine = loadSetting(DEMO, 'engine', 0) as BooleanEngine;
  let repairOrientation = loadSetting(DEMO, 'repairOrientation', false);
  let nonzeroWinding = loadSetting(DEMO, 'nonzeroWinding', false);
  let offsetX = loadSetting(DEMO, 'offsetX', 0.3);
  let offsetY = loadSetting(DEMO, 'offsetY', 0.0);
  let offsetZ = loadSetting(DEMO, 'offsetZ', 0.0);
  let offsetsOpen = loadSetting(DEMO, 'offsetsOpen', false);
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

  // Thingi mode is active only once a pair is actually loaded; while a pair
  // is in flight the built-in shapes stay on screen.
  const inThingi = () => source === 'thingi' && !!(thingiA && thingiB);

  // ---- Sidebar ----

  const troubleOptions = [{ value: '', text: 'Pick a known case…' }]
    .concat(TROUBLE_CASES.map(c => ({ value: c.value, text: c.text })));

  const panel: GalleryPanel = buildPanel(controlsEl, {
    source, shapeA, shapeB, pairMode, troubleOptions,
    troubleCase: TROUBLE_CASES.some(c => c.value === troubleCase) ? troubleCase : '',
    pickA, pickB, op, engine, repair: repairOrientation, nonzero: nonzeroWinding,
    offset: [offsetX, offsetY, offsetZ], offsetsOpen, wireframe, xray, animate,
  }, {
    onSource: s => setSource(s),
    onShape: (slot, value) => {
      if (slot === 'a') { shapeA = value; saveSetting(DEMO, 'shapeA', shapeA); }
      else { shapeB = value; saveSetting(DEMO, 'shapeB', shapeB); }
      update();
    },
    onPairMode: mode => setPairMode(mode, true),
    onLoadRandom: () => { loadRandomPair(); },
    onSwap: () => { swapOperands(); },
    onTrouble: v => {
      if (!v) return;
      troubleCase = v;
      saveSetting(DEMO, 'troubleCase', v);
      loadTroubleCase(v);
    },
    onPickId: (slot, value) => {
      if (slot === 'a') { pickA = value; saveSetting(DEMO, 'pickA', pickA); }
      else { pickB = value; saveSetting(DEMO, 'pickB', pickB); }
    },
    onLoadPicked: () => { loadPickedModels(); },
    onOp: v => { op = v; saveSetting(DEMO, 'op', op); update(); },
    onEngine: v => {
      engine = v as BooleanEngine;
      saveSetting(DEMO, 'engine', engine);
      // Switching to Auto is the moment the self-intersection verdict becomes
      // meaningful, so re-feed the operands to have the worker measure it.
      if (wantSelfIntersecting() && inThingi()) void measureSelfIntersection();
      update();
    },
    // Repair sits with the engine choice: it changes what the operands *mean*
    // (inside-out bodies become solid material) before any boolean runs.
    onRepair: v => { repairOrientation = v; saveSetting(DEMO, 'repairOrientation', v); update(); },
    // Where repair rewrites the operands, this changes what the robust engine
    // counts as solid: {winding != 0} instead of {winding >= 1}. It is the only
    // thing that saves a single shell wound correctly in one region and
    // inside-out in another — exactly what repair (per shell) cannot fix.
    onNonzero: v => { nonzeroWinding = v; saveSetting(DEMO, 'nonzeroWinding', v); update(); },
    onOffset: (axis, v) => {
      if (axis === 'x') { offsetX = v; saveSetting(DEMO, 'offsetX', v); }
      else if (axis === 'y') { offsetY = v; saveSetting(DEMO, 'offsetY', v); }
      else { offsetZ = v; saveSetting(DEMO, 'offsetZ', v); }
      panel.setOffsets([offsetX, offsetY, offsetZ]);
      update();
    },
    onOffsetsToggle: open => { offsetsOpen = open; saveSetting(DEMO, 'offsetsOpen', open); },
    onWireframe: v => { wireframe = v; saveSetting(DEMO, 'wireframe', v); viewer.setWireframe(v); },
    onXRay: v => { xray = v; saveSetting(DEMO, 'xray', v); viewer.setXRay(v); },
    onAnimate: toggleAnimate,
    onCopy: () => { copyDebugInfo(); },
    onPaste: () => { useDebugInfo(); },
  });

  // Owns the stats grid *and* the "Last frame" boolean timing, which has to
  // keep updating (pending clock, failed, cancelled) between results.
  const stats = new BooleanReadout(panel.readoutEl);

  function describeOperand(o: ThingiOperand) {
    return {
      id: o.model.id,
      name: o.model.name,
      url: meshZipUrl(o.model),
      faces: o.model.faces,
      edge_manifold: o.model.edge_manifold,
      vertex_manifold: o.model.vertex_manifold,
      imported: {
        num_tri: o.info.num_tri,
        is_soup: o.info.is_soup,
        self_intersecting: o.info.self_intersecting,
        status: o.info.status,
      },
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
      nonzero_winding: nonzeroWinding,
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

  // ---- Worker-driven evaluation: the boolean runs off the main thread ----

  const runner = new BooleanRunner();
  // Busy overlay sits at the top-left of the 3D pane, over the scene; its
  // Cancel button is the one and only cancel affordance.
  runner.attachBusyIndicator(viewerEl);

  function update(silent = false) {
    if (source === 'thingi' && (!inThingi() || loadingPair)) return; // pair still loading
    const info = captureDebugInfo();
    const params: RunParams = inThingi()
      ? {
          source: 'thingi', op, engine, repair: repairOrientation, nonzero: nonzeroWinding,
          ox: offsetX, oy: offsetY, oz: offsetZ, rx: rotX, ry: rotY, rz: rotZ,
          tag: { info, silent },
        }
      : {
          source: 'builtin', shapeA, shapeB, op, engine, repair: repairOrientation, nonzero: nonzeroWinding,
          ox: offsetX, oy: offsetY, oz: offsetZ, rx: rotX, ry: rotY, rz: rotZ,
          tag: { info, silent },
        };
    runner.requestRun(params);
  }

  runner.onResult = ({ data, elapsedMs }, params) => {
    const { info } = params.tag;
    stats.setResult(elapsedMs);
    info.elapsed_ms = Math.round(elapsedMs * 10) / 10;
    info.result = { num_tri: data.num_tri, num_vert: data.num_vert, volume: data.volume };
    try { localStorage.setItem('boolean-gallery-debug', JSON.stringify(info)); } catch { /* ignore */ }
    viewer.setMesh(data as MeshData);
    viewer.setColor(OP_COLORS[params.op] || 0x4488cc);
    panel.setNotice(null);
    stats.setMesh(data as MeshData);
    scheduleNextAnimateFrame(elapsedMs);
  };

  runner.onError = (error, params) => {
    const { info, silent } = params.tag;
    stats.setFailed();
    info.error = error;
    try { localStorage.setItem('boolean-gallery-debug', JSON.stringify(info)); } catch { /* ignore */ }
    console.error('Boolean op failed:', info, error);
    if (!silent) {
      panel.setNotice(`<strong>Boolean operation failed:</strong> ${error}<br>` +
        `Use <em>Copy debug</em> to capture the inputs for a bug report.`);
      // The stats rows described a result that is gone; the timing row stays
      // and now reads "failed".
      stats.clear();
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
      panel.setAnimate(false);
      panel.setNotice(`<strong>Animation paused:</strong> this boolean takes ` +
        `${(elapsedMs / 1000).toFixed(1)} s per frame. Re-enable Animate to continue anyway.`);
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

  /// The self-intersection scan is only worth its (~1us per triangle) cost
  /// when the Auto engine will consult it — that is also when the operand
  /// card needs it to explain why Auto picked the robust engine.
  function wantSelfIntersecting(): boolean {
    return ENGINE_NAMES[engine] === 'auto';
  }

  /// Re-feed the committed pair to the worker — after a failed load attempt
  /// or a cancel replaced worker slots with candidate meshes, this restores
  /// the state the on-screen labels describe.
  function restoreWorkerOperands() {
    if (thingiA) runner.importOperand('a', thingiA.arrays, wantSelfIntersecting());
    if (thingiB) runner.importOperand('b', thingiB.arrays, wantSelfIntersecting());
  }

  /// Switching to Auto is the moment the self-intersection verdict starts
  /// mattering (it decides which engine runs), so re-feed the committed pair
  /// with the scan requested and adopt the fresh ImportInfo, which is what
  /// the operand card reads.
  async function measureSelfIntersection() {
    const pair: [ThingiOperand | null, 'a' | 'b'][] = [[thingiA, 'a'], [thingiB, 'b']];
    for (const [o, slot] of pair) {
      if (!o || o.info.self_intersecting !== null) continue;
      const info = await runner.importOperand(slot, o.arrays, true);
      if (info.ok && (slot === 'a' ? thingiA : thingiB) === o) o.info = info;
    }
    if (inThingi()) showOperands();
  }

  /// `#849728 · 702 tris · manifold` — the mono line under each mesh name.
  function operandMeta(o: ThingiOperand): { text: string; caution: boolean } {
    // Topological kind per the dataset flags (what the pools guarantee);
    // "soup" marks the rare case where even welding could not pair it.
    const manifold = o.model.edge_manifold && o.model.vertex_manifold;
    let kind = manifold ? 'manifold' : 'non-manifold';
    if (o.info.is_soup) kind += ', soup';
    // Self-intersecting operands force Auto onto the robust engine. Only
    // known when the scan was requested (Auto); null means "not measured".
    if (o.info.self_intersecting) kind += ', self-intersecting';
    return {
      text: `#${o.model.id} · ${o.model.faces} tris · ${kind}`,
      caution: !manifold || !!o.info.is_soup || !!o.info.self_intersecting,
    };
  }

  /// Repaint the A/B card from the committed pair (plus an optional status
  /// line under it). Rotation is shown as loaded/restored; while Animate runs
  /// it keeps advancing past the printed value, as before.
  function showOperands(message: string | null = null) {
    if (!thingiA || !thingiB) {
      panel.operands.setLines([], null);
      panel.operands.setMessage(message);
      return;
    }
    const a = operandMeta(thingiA);
    const b = operandMeta(thingiB);
    panel.operands.setLines([
      { slot: 'a', title: thingiA.model.name, thingiId: thingiA.model.id, meta: a.text, caution: a.caution },
      { slot: 'b', title: thingiB.model.name, thingiId: thingiB.model.id, meta: b.text, caution: b.caution },
    ], `rot [${Math.round(rotX)}, ${Math.round(rotY)}, ${Math.round(rotZ)}]°`);
    panel.operands.setMessage(message);
  }

  /// Trade which mesh is A and which is B, so Difference subtracts the other
  /// way round. The slots keep their identity — swatch colours, the offset and
  /// the rotation all belong to the slot, so B's transform now applies to the
  /// mesh that just moved into B. The worker's operand slots have to follow
  /// (its 'a'/'b' meshes are what the boolean reads); the messages are queued
  /// before the run below, and the worker processes them in order.
  function swapOperands() {
    if (loadingPair) return;
    if (source === 'thingi') {
      if (!thingiA || !thingiB) return;
      [thingiA, thingiB] = [thingiB, thingiA];
      void runner.importOperand('a', thingiA.arrays, wantSelfIntersecting());
      void runner.importOperand('b', thingiB.arrays, wantSelfIntersecting());
      showOperands();
    } else {
      [shapeA, shapeB] = [shapeB, shapeA];
      saveSetting(DEMO, 'shapeA', shapeA);
      saveSetting(DEMO, 'shapeB', shapeB);
      panel.setShape('a', shapeA);
      panel.setShape('b', shapeB);
    }
    // The debug record is rebuilt from this state, so it simply describes the
    // new assignment.
    update();
  }

  // Some dataset models flagged "closed" still fail the robust importer's
  // stricter closed check — skip those and re-roll, keeping a record of the
  // skips for the debug info.
  async function loadOperand(kind: 'm' | 'n', slot: 'a' | 'b', exclude?: ThingiModel): Promise<ThingiOperand> {
    const MAX_ATTEMPTS = 5;
    for (let attempt = 0; attempt < MAX_ATTEMPTS; attempt++) {
      const model = await pickRandomModel(kind, exclude);
      panel.operands.setMessage(`Fetching #${model.id} “${model.name}”…`);
      const parsed = await fetchMesh(model);
      if (disposed) throw new Error('demo closed');
      const arrays: OperandArrays = { positions: parsed.positions, indices: parsed.indices };
      const info = await runner.importOperand(slot, arrays, wantSelfIntersecting());
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
    panel.setLoading(true);
    panel.setNotice(null);
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
      showOperands();
    } catch (e) {
      console.error('Thingi10K pair load failed:', e);
      restoreWorkerOperands();
      showOperands(inThingi() ? 'Load failed — showing the previous pair.' : null);
      panel.setNotice(`<strong>Failed to load Thingi10K pair:</strong> ${String(e)}`);
    } finally {
      loadingPair = false;
      panel.setLoading(false);
    }
    // Re-render only on success (update() is a no-op while a pair loads).
    // On failure the previous result stays on screen and — crucially — the
    // error stays visible: a successful update() would hide the notice.
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
        showOperands('Pick a known case from the Trouble Cases list.');
      }
    } else if (pairMode === 'pick') {
      showOperands('Enter two model ids and press Load.');
    } else {
      loadRandomPair();
    }
  }

  async function copyDebugInfo() {
    const text = JSON.stringify(lastDebugInfo ?? { note: 'no operation recorded yet' }, null, 2);
    try {
      await navigator.clipboard.writeText(text);
    } catch {
      // Clipboard API can be unavailable (non-secure context); fall back.
      const ta = document.createElement('textarea');
      ta.value = text;
      document.body.appendChild(ta);
      ta.select();
      document.execCommand('copy');
      ta.remove();
    }
    panel.setCopyLabel('Copied!');
    setTimeout(() => { panel.setCopyLabel('Copy debug'); }, 1500);
  }

  // ---- Source / pair mode ----

  function setSource(s: Source) {
    source = s;
    saveSetting(DEMO, 'source', s);
    panel.setSource(s);
    if (s === 'thingi') {
      if (inThingi()) {
        showOperands();
        update();
      } else {
        loadInitialPair();
      }
    } else {
      update();
    }
  }

  /// Switch the Thingi10K pair mode, showing only that mode's sub-control:
  /// the ⟳ Load Random Pair button for the three random kinds, the curated
  /// list for Trouble Cases, the id fields for Pick Models.
  ///
  /// `userInitiated` marks a change made through the select: picking a random
  /// kind then immediately rolls a pair of that kind, so the display always
  /// matches what the Pair Type says.
  function setPairMode(mode: PairMode, userInitiated = false) {
    const changed = mode !== pairMode;
    pairMode = mode;
    saveSetting(DEMO, 'pairMode', mode);
    const random = isRandomMode(mode);
    if (random) {
      pairKind = mode;
      saveSetting(DEMO, 'pairKind', pairKind);
    }
    panel.setPairMode(mode);
    if (userInitiated && changed && source === 'thingi') {
      if (random) {
        loadRandomPair();
      } else if (mode === 'trouble' && troubleCase && TROUBLE_CASES.some(c => c.value === troubleCase)) {
        loadTroubleCase(troubleCase);
      }
    }
  }

  runner.onBusyChange = busy => {
    // Keep the "Last frame" row honest while the boolean runs: it shows a
    // live pending clock instead of the previous run's stale duration.
    stats.setBusy(busy);
  };
  runner.onCancelled = () => {
    stats.setCancelled();
    panel.setNotice('<strong>Computation cancelled.</strong> The previous result stays on screen.');
  };

  // ---- Paste debug info: restore the exact captured state ----

  function forceAnimateOff() {
    if (animate || animating) toggleAnimate(false);
    panel.setAnimate(false);
  }

  function showUseError(msg: string) {
    panel.setNotice(`<strong>Paste debug info:</strong> ${msg}`);
  }

  async function loadOperandFromRecord(rec: any, slot: 'a' | 'b'): Promise<ThingiOperand> {
    const model = modelFromDebugRecord(rec);
    panel.operands.setMessage(`Fetching #${model.id} “${model.name}”…`);
    const parsed = await fetchMesh(model);
    if (disposed) throw new Error('demo closed');
    const arrays: OperandArrays = { positions: parsed.positions, indices: parsed.indices };
    const info = await runner.importOperand(slot, arrays, wantSelfIntersecting());
    if (!info.ok) {
      throw new Error(`#${model.id} failed to import: ${info.status}`);
    }
    return { model, arrays, info, normalization: { center: parsed.center, scale: parsed.scale } };
  }

  async function useDebugInfo() {
    panel.setNotice(null);
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
        'Press <em>Copy debug</em> on a gallery frame first, then paste it back here.');
      return;
    }
    if (!info || info.demo !== DEMO || !Array.isArray(info.rotation_degrees) || !Array.isArray(info.offset)) {
      showUseError('the clipboard JSON is not Boolean Gallery debug info ' +
        '(expected the object produced by <em>Copy debug</em>).');
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
  /// Cases list and the Pick Models id fields.
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
    panel.setNotice(null);
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
  /// Shared by the Paste chip (clipboard), the Trouble Cases list and the
  /// Pick Models fields. `modeAfter` is the Pair Type the restored state
  /// belongs to: a pasted report lands in Pick Models (its ids are now the
  /// typed-in pair), while the curated list stays on Trouble Cases.
  async function applyDebugInfo(info: any, modeAfter: PairMode = 'pick') {
    // Land exactly on the captured frame: freeze animation first.
    forceAnimateOff();

    const opIdx = OP_NAMES.indexOf(info.op);
    if (opIdx >= 0) { op = opIdx; saveSetting(DEMO, 'op', op); panel.setOp(op); }
    const engIdx = ENGINE_NAMES.indexOf(info.engine);
    if (engIdx >= 0) { engine = engIdx as BooleanEngine; saveSetting(DEMO, 'engine', engine); panel.setEngine(engine); }
    if (typeof info.repair_orientation === 'boolean') {
      repairOrientation = info.repair_orientation;
      saveSetting(DEMO, 'repairOrientation', repairOrientation);
      panel.setRepair(repairOrientation);
    }
    if (typeof info.nonzero_winding === 'boolean') {
      nonzeroWinding = info.nonzero_winding;
      saveSetting(DEMO, 'nonzeroWinding', nonzeroWinding);
      panel.setNonzero(nonzeroWinding);
    }
    [offsetX, offsetY, offsetZ] = info.offset.map(Number) as [number, number, number];
    saveSetting(DEMO, 'offsetX', offsetX);
    saveSetting(DEMO, 'offsetY', offsetY);
    saveSetting(DEMO, 'offsetZ', offsetZ);
    panel.setOffsets([offsetX, offsetY, offsetZ]);
    [rotX, rotY, rotZ] = info.rotation_degrees.map(Number) as [number, number, number];

    if (info.model_a && info.model_b) {
      source = 'thingi';
      saveSetting(DEMO, 'source', source);
      panel.setSource(source);
      // Keep the random-pair kind in sync when the record names one, but the
      // visible mode is whatever brought us here.
      if (typeof info.pair_kind === 'string' && isRandomMode(info.pair_kind as PairMode)) {
        pairKind = info.pair_kind as PairKind;
        saveSetting(DEMO, 'pairKind', pairKind);
      }
      // Surface the restored pair in the Pick Models fields so it can be
      // re-run, tweaked one id at a time, or read off at a glance.
      if (typeof info.model_a.id === 'number' && typeof info.model_b.id === 'number') {
        pickA = info.model_a.id;
        pickB = info.model_b.id;
        saveSetting(DEMO, 'pickA', pickA);
        saveSetting(DEMO, 'pickB', pickB);
        panel.setPickId('a', pickA);
        panel.setPickId('b', pickB);
      }
      setPairMode(modeAfter);
      loadingPair = true;
      panel.setLoading(true);
      skippedImports = [];
      let loaded = false;
      try {
        const a = await loadOperandFromRecord(info.model_a, 'a');
        const b = await loadOperandFromRecord(info.model_b, 'b');
        freeThingiPair();
        thingiA = a;
        thingiB = b;
        loaded = true;
        showOperands();
      } catch (e) {
        restoreWorkerOperands();
        showOperands(inThingi() ? 'Load failed — showing the previous pair.' : null);
        showUseError(`failed to load the requested pair: ${String(e)}`);
      } finally {
        loadingPair = false;
        panel.setLoading(false);
      }
      // Re-render only on success (update() is a no-op while a pair loads).
      // On failure the previous result stays on screen and — crucially — the
      // error stays visible: a successful update() would hide the notice.
      if (loaded && inThingi()) update();
    } else {
      source = 'builtin';
      saveSetting(DEMO, 'source', source);
      panel.setSource(source);
      const sa = SHAPE_NAMES.indexOf(info.shape_a);
      const sb = SHAPE_NAMES.indexOf(info.shape_b);
      if (sa >= 0) { shapeA = sa; saveSetting(DEMO, 'shapeA', shapeA); panel.setShape('a', shapeA); }
      if (sb >= 0) { shapeB = sb; saveSetting(DEMO, 'shapeB', shapeB); panel.setShape('b', shapeB); }
      update();
    }
  }

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
    stats.dispose();
    runner.dispose();
    viewer.dispose();
  };
}
