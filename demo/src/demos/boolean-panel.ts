// The Boolean Gallery's sidebar ("layout 2a"): four bands in the order the
// user works — Input → Operation → View → Result.
//
// This module owns the panel's DOM only. Every control reports through the
// `PanelHandlers` callbacks and is re-rendered through the returned
// `GalleryPanel` setters, so all gallery state (and the persistence of it)
// stays in boolean-gallery.ts. The Result band hands its readout element to
// BooleanReadout, which paints the stats grid and the Last-frame row.
//
// Widgets come from gallery-widgets.ts; the offset sliders are the shared
// createSlider from controls.ts, unchanged, just tucked inside a disclosure.

import { createSlider } from '../controls.ts';
import {
  createSegmented, createOpGrid, createPill, createChip, createIconButton,
  createCheckRow, createSelect, createDisclosure, createOperandCard,
  type Toggleable, type CheckRow, type Disclosure, type OperandCard, type OperandLine,
} from './gallery-widgets.ts';
import { OP_CHOICES, PAIR_MODES, SHAPES, type PairMode } from './boolean-gallery-data.ts';

export type Source = 'builtin' | 'thingi';
export type Axis = 'x' | 'y' | 'z';

/** Engine segment values, indexed by the BooleanEngine enum value. */
type EngineKey = 'exact' | 'robust' | 'auto';
const ENGINE_KEYS: EngineKey[] = ['exact', 'robust', 'auto'];

export interface PanelState {
  source: Source;
  shapeA: number;
  shapeB: number;
  pairMode: PairMode;
  troubleCase: string;
  troubleOptions: { value: string; text: string }[];
  pickA: number;
  pickB: number;
  op: number;
  /** 0 exact, 1 robust, 2 auto (auto selects neither segment). */
  engine: number;
  repair: boolean;
  nonzero: boolean;
  offset: [number, number, number];
  offsetsOpen: boolean;
  wireframe: boolean;
  xray: boolean;
  animate: boolean;
}

export interface PanelHandlers {
  onSource(s: Source): void;
  onShape(slot: 'a' | 'b', value: number): void;
  onPairMode(mode: PairMode): void;
  onLoadRandom(): void;
  onTrouble(value: string): void;
  onPickId(slot: 'a' | 'b', value: number): void;
  onLoadPicked(): void;
  onOp(op: number): void;
  onEngine(engine: number): void;
  onRepair(v: boolean): void;
  onNonzero(v: boolean): void;
  onOffset(axis: Axis, value: number): void;
  onOffsetsToggle(open: boolean): void;
  onWireframe(v: boolean): void;
  onXRay(v: boolean): void;
  onAnimate(v: boolean): void;
  onCopy(): void;
  onPaste(): void;
}

/** Formatted offsets for the closed disclosure row: `0.1, 0, 0`, or `none`. */
export function offsetSummary(o: [number, number, number]): { text: string; muted: boolean } {
  if (o.every(v => v === 0)) return { text: 'none', muted: true };
  return { text: o.map(v => String(v)).join(', '), muted: false };
}

function band(root: HTMLElement, heading: string, extraClass = ''): HTMLElement {
  const el = document.createElement('div');
  el.className = extraClass ? `bg-band ${extraClass}` : 'bg-band';
  const h = document.createElement('div');
  h.className = 'bg-band-h';
  h.textContent = heading;
  el.appendChild(h);
  root.appendChild(el);
  return el;
}

function row(cls: string, ...children: HTMLElement[]): HTMLElement {
  const el = document.createElement('div');
  el.className = cls;
  for (const c of children) el.appendChild(c);
  return el;
}

/// Inline numeric field (a Thingi10K model id). Commits on Enter or blur so a
/// half-typed id never starts a fetch; 0 means "empty".
function numberField(placeholder: string, value: number, onCommit: (v: number) => void): HTMLInputElement {
  const input = document.createElement('input');
  input.type = 'number';
  input.className = 'bgw-num';
  input.placeholder = placeholder;
  input.title = placeholder;
  input.value = value ? String(value) : '';
  const commit = () => {
    const parsed = parseInt(input.value.trim(), 10);
    onCommit(Number.isFinite(parsed) ? parsed : 0);
  };
  input.addEventListener('change', commit);
  input.addEventListener('keydown', e => {
    if (e.key === 'Enter') { e.preventDefault(); commit(); }
  });
  return input;
}

export interface GalleryPanel {
  /** Target element for BooleanReadout (the Result band's stats area). */
  readoutEl: HTMLElement;
  operands: OperandCard;
  setSource(s: Source): void;
  setShape(slot: 'a' | 'b', value: number): void;
  setPairMode(mode: PairMode): void;
  setTrouble(value: string): void;
  setPickId(slot: 'a' | 'b', value: number): void;
  setOp(op: number): void;
  setEngine(engine: number): void;
  setRepair(v: boolean): void;
  setNonzero(v: boolean): void;
  setOffsets(o: [number, number, number]): void;
  setAnimate(on: boolean): void;
  /** Block the pair-loading controls while a fetch is in flight. */
  setLoading(loading: boolean): void;
  /** Transient label on the Copy chip ("Copied!"). */
  setCopyLabel(text: string): void;
  /** Error / notice banner at the foot of the panel; null hides it. */
  setNotice(html: string | null): void;
}

export function buildPanel(root: HTMLElement, state: PanelState, h: PanelHandlers): GalleryPanel {
  root.className = 'bg-panel';
  root.replaceChildren();

  // ---- Band 1: Input ----
  const inputBand = band(root, 'Input');
  const sourceSeg = createSegmented<Source>(
    [
      { value: 'builtin', text: 'Primitives', title: 'Built-in primitive shapes' },
      { value: 'thingi', text: 'Thingi10K', title: 'Real-world Thingi10K mesh pairs' },
    ],
    state.source, h.onSource);
  inputBand.appendChild(sourceSeg.el);

  const pairSelect = createSelect(PAIR_MODES, state.pairMode,
    v => h.onPairMode(v as PairMode), 'Pair type');
  const shapeSelects = {
    a: createSelect(SHAPES, String(state.shapeA), v => h.onShape('a', parseInt(v)), 'Shape A'),
    b: createSelect(SHAPES, String(state.shapeB), v => h.onShape('b', parseInt(v)), 'Shape B'),
  };
  const loadBtn = createIconButton('⟳', 'Load random pair', h.onLoadRandom);
  const pairRow = row('bgw-row', pairSelect.el, shapeSelects.a.el, shapeSelects.b.el, loadBtn);
  inputBand.appendChild(pairRow);

  const troubleSelect = createSelect(state.troubleOptions, state.troubleCase, h.onTrouble,
    'Known trouble cases');
  inputBand.appendChild(troubleSelect.el);

  const pickInputs = {
    a: numberField('Model A id', state.pickA, v => h.onPickId('a', v)),
    b: numberField('Model B id', state.pickB, v => h.onPickId('b', v)),
  };
  const pickBtn = createChip('Load', 'Load these two model ids', h.onLoadPicked);
  const pickRow = row('bgw-row', pickInputs.a, pickInputs.b, pickBtn);
  inputBand.appendChild(pickRow);

  const operands = createOperandCard();
  inputBand.appendChild(operands.el);

  // ---- Band 2: Operation ----
  const opBand = band(root, 'Operation');
  const opGrid = createOpGrid(OP_CHOICES, state.op, h.onOp);
  opBand.appendChild(opGrid.el);

  const engineLabel = document.createElement('span');
  engineLabel.className = 'bgw-inline-label';
  engineLabel.textContent = 'Engine';
  // 3-up rather than the spec's 2-up: Auto is real plumbing (it picks the
  // engine per input, and drives the self-intersection scan), so it stays a
  // first-class choice instead of an unselectable persisted state.
  const engineSeg = createSegmented<EngineKey>(
    [
      { value: 'exact', text: 'Exact', title: 'Exact (C++ matching)' },
      { value: 'robust', text: 'Robust', title: 'Robust (non-manifold)' },
      { value: 'auto', text: 'Auto', title: 'Auto (Exact when provably safe, else Robust)' },
    ],
    ENGINE_KEYS[state.engine] ?? 'exact',
    v => h.onEngine(ENGINE_KEYS.indexOf(v)),
    'compact');
  opBand.appendChild(row('bgw-inline', engineLabel, engineSeg.el));

  const repairRow = createCheckRow('Repair orientation',
    'Rewind inside-out shells of both operands before the boolean (robust engine only)',
    state.repair, h.onRepair);
  const nonzeroRow = createCheckRow('Keep inside-out geometry',
    'Nonzero winding: treat {winding != 0} as solid instead of {winding >= 1} (robust engine only)',
    state.nonzero, h.onNonzero);
  const flags = row('bgw-flags', repairRow.el, nonzeroRow.el);
  opBand.appendChild(flags);

  const offsets = createDisclosure('Offset B', state.offsetsOpen, h.onOffsetsToggle);
  const sliders: Record<Axis, HTMLElement> = {
    x: createSlider('Offset X ', -1.5, 1.5, state.offset[0], 0.1, v => h.onOffset('x', v)),
    y: createSlider('Offset Y ', -1.5, 1.5, state.offset[1], 0.1, v => h.onOffset('y', v)),
    z: createSlider('Offset Z ', -1.5, 1.5, state.offset[2], 0.1, v => h.onOffset('z', v)),
  };
  offsets.body.appendChild(sliders.x);
  offsets.body.appendChild(sliders.y);
  offsets.body.appendChild(sliders.z);
  opBand.appendChild(offsets.el);

  // ---- Band 3: View ----
  const viewBand = band(root, 'View', 'bg-band-view');
  const pills = {
    wireframe: createPill('Wireframe', 'Wireframe', state.wireframe, h.onWireframe),
    xray: createPill('X-Ray', 'X-Ray (depth peeled)', state.xray, h.onXRay),
    animate: createPill('Animate', 'Animate the rotation of mesh B', state.animate, h.onAnimate),
  };
  viewBand.appendChild(row('bgw-pills', pills.wireframe.el, pills.xray.el, pills.animate.el));

  // ---- Band 4: Result ----
  const resultBand = band(root, 'Result', 'bg-band-result');
  const heading = resultBand.firstElementChild as HTMLElement;
  const copyChip = createChip('Copy debug', 'Copy debug info for this frame to the clipboard', h.onCopy);
  const pasteChip = createChip('Paste', 'Restore the state from debug info in the clipboard', h.onPaste);
  const headRow = row('bgw-result-head', heading, row('bgw-chips', copyChip, pasteChip));
  resultBand.prepend(headRow);

  const readoutEl = document.createElement('div');
  readoutEl.className = 'bgr';
  resultBand.appendChild(readoutEl);

  const notice = document.createElement('div');
  notice.className = 'bgw-notice';
  notice.style.display = 'none';
  resultBand.appendChild(notice);

  // Primitives mode replaces the pair select with the two shape pickers and
  // hides every Thingi10K-only affordance, so only one of the two lives in
  // the pair row at a time.
  let currentSource: Source = state.source;
  let currentMode: PairMode = state.pairMode;

  function applyRows() {
    const thingi = currentSource === 'thingi';
    const random = currentMode === 'mm' || currentMode === 'mn' || currentMode === 'nn';
    pairSelect.el.style.display = thingi ? '' : 'none';
    shapeSelects.a.el.style.display = thingi ? 'none' : '';
    shapeSelects.b.el.style.display = thingi ? 'none' : '';
    loadBtn.style.display = thingi && random ? '' : 'none';
    troubleSelect.el.style.display = thingi && currentMode === 'trouble' ? '' : 'none';
    pickRow.style.display = thingi && currentMode === 'pick' ? '' : 'none';
    operands.setVisible(thingi);
  }

  function applyPairMode(mode: PairMode) {
    currentMode = mode;
    pairSelect.select.value = mode;
    applyRows();
  }

  function applySource(s: Source) {
    currentSource = s;
    applyRows();
  }

  function applyEngine(engine: number) {
    engineSeg.set(ENGINE_KEYS[engine] ?? 'exact');
    // Both flags are robust-engine features; Exact ignores them. Auto keeps
    // them live — nonzero winding there forces the robust engine, by design.
    const active = engine !== 0;
    repairRow.setEnabled(active);
    nonzeroRow.setEnabled(active);
  }

  function setSliderValue(ctl: HTMLElement, v: number) {
    const input = ctl.querySelector('input');
    if (input) input.value = String(v);
    const span = ctl.querySelector('.slider-value');
    if (span) span.textContent = String(v);
  }

  applySource(state.source);
  applyPairMode(state.pairMode);
  applyEngine(state.engine);
  const summary = offsetSummary(state.offset);
  offsets.setSummary(summary.text, summary.muted);

  return {
    readoutEl,
    operands,
    setSource(s) { sourceSeg.set(s); applySource(s); },
    setShape(slot, value) { shapeSelects[slot].select.value = String(value); },
    setPairMode(mode) { applyPairMode(mode); },
    setTrouble(value) { troubleSelect.select.value = value; },
    setPickId(slot, value) { pickInputs[slot].value = value ? String(value) : ''; },
    setOp(op) { opGrid.set(op); },
    setEngine(engine) { applyEngine(engine); },
    setRepair(v) { repairRow.set(v); },
    setNonzero(v) { nonzeroRow.set(v); },
    setOffsets(o) {
      setSliderValue(sliders.x, o[0]);
      setSliderValue(sliders.y, o[1]);
      setSliderValue(sliders.z, o[2]);
      const s = offsetSummary(o);
      offsets.setSummary(s.text, s.muted);
    },
    setAnimate(on) { pills.animate.set(on); },
    setLoading(loading) {
      loadBtn.disabled = loading;
      loadBtn.classList.toggle('is-busy', loading);
      pickBtn.disabled = loading;
      pasteChip.disabled = loading;
    },
    setCopyLabel(text) { copyChip.textContent = text; },
    setNotice(html) {
      if (html === null) {
        notice.style.display = 'none';
        notice.innerHTML = '';
      } else {
        notice.innerHTML = html;
        notice.style.display = '';
      }
    },
  };
}

export type { OperandLine, Toggleable, CheckRow, Disclosure, OperandCard };
