// Smoke tests for the overall-progress mapping (progress-model.ts) and for the
// bar it drives (busy-indicator.ts), run with `bun test` from demo/.
//
// The indicator half uses the same fake-DOM pattern as the other DOM smoke
// checks: a handful of stub element objects is enough, because the card only
// ever touches textContent, style.width and classList.

import { expect, test, beforeEach } from 'bun:test';
import { ProgressModel, phaseLabel } from './progress-model.ts';
import { BusyIndicator } from './busy-indicator.ts';

/** The robust engine's phase order, as the kernel emits it. */
const ROBUST_PHASES = [
  'narrow phase',
  'self intersections',
  'candidate points',
  'registries',
  'arrangements',
  'cells',
  'winding',
  'assemble',
];

test('overall fraction never decreases across a full robust run', () => {
  const model = new ProgressModel();
  let last = -1;
  for (const phase of ROBUST_PHASES) {
    for (const f of [0, 0.25, 0.5, 0.75, 1]) {
      const { fraction } = model.update(7, phase, f);
      expect(fraction).not.toBeNull();
      expect(fraction!).toBeGreaterThanOrEqual(last);
      last = fraction!;
    }
  }
  // The last phase completing means the whole operation is complete.
  expect(last).toBeCloseTo(1, 12);
});

test('a phase restarting at 0 does not rewind the bar', () => {
  const model = new ProgressModel();
  model.update(1, 'narrow phase', 1);
  const afterNarrow = model.update(1, 'self intersections', 0).fraction!;
  expect(afterNarrow).toBeGreaterThan(0);
  const mid = model.update(1, 'self intersections', 0.5).fraction!;
  // Entering the next phase at 0 must not undo the previous phase's share.
  const next = model.update(1, 'candidate points', 0).fraction!;
  expect(next).toBeGreaterThanOrEqual(mid);
});

test('out-of-order fractions inside one phase are clamped monotonic', () => {
  const model = new ProgressModel();
  const high = model.update(3, 'arrangements', 0.9).fraction!;
  const low = model.update(3, 'arrangements', 0.2).fraction!;
  expect(low).toBe(high);
});

test('phase ids map to the same result as phase names', () => {
  const byName = new ProgressModel();
  const byId = new ProgressModel();
  for (let i = 0; i < ROBUST_PHASES.length; i++) {
    const a = byName.update(1, ROBUST_PHASES[i]!, 0.5);
    const b = byId.update(1, i, 0.5);
    expect(b.fraction).toBe(a.fraction);
    expect(b.label).toBe(a.label);
  }
});

test('label carries the phase name and the overall percentage', () => {
  const model = new ProgressModel();
  // narrow phase = 5/100 of the work, half done -> 3%.
  expect(model.update(1, 'narrow phase', 0.5).label).toBe('narrow phase — 3%');
  // arrangements starts at 5+20+10+10 = 45 and is 30 wide; half of it -> 60%.
  expect(model.update(1, 'arrangements', 0.5).label).toBe('arrangements — 60%');
});

test('the exact engine stays indeterminate', () => {
  const model = new ProgressModel();
  const out = model.update(2, 'exact boolean', null);
  expect(out.fraction).toBeNull();
  expect(out.label).toBe('exact boolean');
});

test('an unknown phase keeps the bar where it was', () => {
  const model = new ProgressModel();
  const known = model.update(1, 'cells', 0.5).fraction!;
  const unknown = model.update(1, 99, 0.5);
  expect(unknown.fraction).toBe(known);
  // cells starts at 75/100 and is 15 wide, so half of it is 82.5% -> 83%.
  expect(unknown.label).toBe('phase 99 — 83%');
});

test('a new seq restarts the bar, the same seq continues it', () => {
  const model = new ProgressModel();
  const late = model.update(1, 'winding', 1).fraction!;
  expect(late).toBeGreaterThan(0.9);
  const fresh = model.update(2, 'narrow phase', 0).fraction!;
  expect(fresh).toBe(0);
  expect(phaseLabel(0)).toBe('narrow phase');
});

// --- fake DOM -------------------------------------------------------------

class FakeClassList {
  private set = new Set<string>();
  add(...c: string[]) { for (const x of c) this.set.add(x); }
  remove(...c: string[]) { for (const x of c) this.set.delete(x); }
  contains(c: string) { return this.set.has(c); }
}

class FakeElement {
  children: FakeElement[] = [];
  parentElement: FakeElement | null = null;
  textContent = '';
  type = '';
  style: Record<string, string> = {};
  classList = new FakeClassList();
  private _className = '';
  get className() { return this._className; }
  set className(v: string) {
    this._className = v;
    for (const c of v.split(/\s+/).filter(Boolean)) this.classList.add(c);
  }
  setAttribute() {}
  addEventListener() {}
  appendChild(child: FakeElement) {
    child.parentElement = this;
    this.children.push(child);
    return child;
  }
  remove() {
    this.parentElement = null;
  }
  /** Depth-first search by the single class name the indicator uses. */
  find(cls: string): FakeElement | null {
    if (this.classList.contains(cls)) return this;
    for (const c of this.children) {
      const hit = c.find(cls);
      if (hit) return hit;
    }
    return null;
  }
}

let body: FakeElement;

beforeEach(() => {
  body = new FakeElement();
  (globalThis as any).document = {
    createElement: () => new FakeElement(),
    querySelector: () => null,
    body,
  };
  (globalThis as any).window = {
    setTimeout: () => 1,
    clearTimeout: () => {},
    setInterval: () => 2,
    clearInterval: () => {},
  };
});

test('the indicator renders a monotonic fill and the phase labels', () => {
  const indicator = new BusyIndicator(0);
  indicator.show('Computing… Robust union');
  const model = new ProgressModel();

  const card = body.children[0]!;
  const track = card.find('busy-track')!;
  const bar = card.find('busy-bar')!;
  const label = card.find('busy-phase')!;

  expect(track.classList.contains('indeterminate')).toBe(true);

  let lastWidth = -1;
  const seenLabels: string[] = [];
  for (const phase of ROBUST_PHASES) {
    for (const f of [0, 0.5, 1]) {
      const out = model.update(11, phase, f);
      indicator.setPhase(out.label, out.fraction);
      seenLabels.push(label.textContent);
      expect(track.classList.contains('indeterminate')).toBe(false);
      const width = parseFloat(bar.style.width);
      expect(width).toBeGreaterThanOrEqual(lastWidth);
      lastWidth = width;
    }
  }
  expect(lastWidth).toBe(100);
  for (const phase of ROBUST_PHASES) {
    expect(seenLabels.some(l => l.startsWith(phase))).toBe(true);
  }
});

test('an exact-engine run leaves the indicator indeterminate', () => {
  const indicator = new BusyIndicator(0);
  indicator.show('Computing… Exact union');
  const model = new ProgressModel();
  const card = body.children[0]!;
  const track = card.find('busy-track')!;

  const out = model.update(4, 'exact boolean', null);
  indicator.setPhase(out.label, out.fraction);
  expect(track.classList.contains('indeterminate')).toBe(true);
  expect(card.find('busy-phase')!.textContent).toBe('exact boolean');
});
