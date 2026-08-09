// Regression tests for the Boolean Gallery info panel's "Last Frame" row
// (boolean-readout.ts), run with `bun test` from demo/.
//
// The bug: "Last Frame" is the boolean's wall time and was only ever written
// when a run *finished*, so during a multi-second computation — and after a
// failure or a cancel — the panel kept displaying the previous run's number
// as if it still described what was on screen.
//
// The panel only touches innerHTML, so a one-field stub stands in for the DOM.

import { expect, test } from 'bun:test';
import { BooleanReadout, frameStatRow, formatMs, frameIsSlow } from './boolean-readout.ts';
import type { MeshData } from '../wasm.ts';

function fakeEl(): HTMLElement {
  return { innerHTML: '' } as unknown as HTMLElement;
}

const MESH = {
  num_vert: 8, num_tri: 12, volume: 1, surface_area: 6,
} as unknown as MeshData;

/** Drives the injected clock so the pending row is deterministic. */
function harness() {
  let now = 0;
  const el = fakeEl();
  const readout = new BooleanReadout(el, () => now);
  return {
    el, readout,
    advance(ms: number) { now += ms; },
    text() { return el.innerHTML; },
  };
}

test('pending row ticks instead of freezing on the previous duration', () => {
  const h = harness();
  h.readout.setBusy(true);
  h.readout.setResult(120);
  h.readout.setMesh(MESH);
  h.readout.setBusy(false);
  expect(h.text()).toContain('120.0 ms');

  // A long run starts: the stale 120 ms must not be presented as the answer.
  h.readout.setBusy(true);
  h.advance(2500);
  h.readout.setBusy(true); // re-render at the tick's timing
  expect(h.text()).toContain('computing…');

  h.readout.dispose();
});

test('pending clock advances with wall time', () => {
  const stat = { status: 'computing' as const, lastMs: 120, startedMs: 1000 };
  expect(frameStatRow(stat, 1000)!.value).toBe('computing… 0.0 ms (last 120.0 ms)');
  expect(frameStatRow(stat, 3400)!.value).toBe('computing… 2.40 s (last 120.0 ms)');
});

test('nothing is claimed before the first run', () => {
  expect(frameStatRow({ status: 'idle', lastMs: null, startedMs: 0 }, 0)).toBeNull();
});

test('a failed run is not reported as a frame time', () => {
  const h = harness();
  h.readout.setBusy(true);
  h.readout.setResult(120);
  h.readout.setMesh(MESH);
  h.readout.setBusy(false);

  h.readout.setBusy(true);
  h.readout.setBusy(false); // runner idles before reporting the error
  h.readout.setFailed();
  expect(h.text()).toContain('failed');
  expect(h.text()).toContain('(last 120.0 ms)');
  h.readout.dispose();
});

test('a cancelled run says so and stops the clock', () => {
  const h = harness();
  h.readout.setBusy(true);
  h.readout.setBusy(false); // BooleanRunner.cancel() clears busy first
  h.readout.setCancelled();
  expect(h.text()).toContain('cancelled');
  h.readout.dispose();
});

test('coalesced runs restart the pending clock rather than showing the old value', () => {
  const h = harness();
  h.readout.setBusy(true);   // run 1 starts at t=0
  h.advance(500);
  h.readout.setResult(500);  // run 2 already queued: no busy transition
  h.advance(200);
  h.readout.setMesh(MESH);   // forces a re-render at t=700
  expect(h.text()).toContain('computing… 200.0 ms (last 500.0 ms)');
  h.readout.dispose();
});

test('mesh rows survive timing-row updates', () => {
  const h = harness();
  h.readout.setMesh(MESH);
  h.readout.setBusy(true);
  expect(h.text()).toContain('Triangles');
  expect(h.text()).toContain('computing…');
  h.readout.dispose();
});

test('formatMs switches to seconds at 1000 ms', () => {
  expect(formatMs(999.94)).toBe('999.9 ms');
  expect(formatMs(1000)).toBe('1.00 s');
});

// ---- Result band presentation (layout 2a) ----

test('counts are printed with thousands separators', () => {
  const h = harness();
  h.readout.setMesh({ num_vert: 8649, num_tri: 17294, volume: 1234.5, surface_area: 6 } as unknown as MeshData);
  expect(h.text()).toContain('8,649');
  expect(h.text()).toContain('17,294');
  expect(h.text()).toContain('1,234.5000');
  h.readout.dispose();
});

test('abbreviated labels keep the full word as a tooltip', () => {
  const h = harness();
  h.readout.setMesh(MESH);
  expect(h.text()).toContain('title="Triangles">Tris<');
  expect(h.text()).toContain('title="Surface area">Area<');
  h.readout.dispose();
});

test('a slow frame is flagged for the warning colour, a fast one is not', () => {
  const idle = (ms: number) => ({ status: 'idle' as const, lastMs: ms, startedMs: 0 });
  expect(frameIsSlow(idle(120), 0)).toBe(false);
  expect(frameIsSlow(idle(501), 0)).toBe(true);
  // A run that is *already* over the budget warns while it is still running.
  expect(frameIsSlow({ status: 'computing', lastMs: null, startedMs: 0 }, 400)).toBe(false);
  expect(frameIsSlow({ status: 'computing', lastMs: null, startedMs: 0 }, 900)).toBe(true);
  // Failure and cancellation are never "a normal frame time".
  expect(frameIsSlow({ status: 'failed', lastMs: 10, startedMs: 0 }, 0)).toBe(true);
  expect(frameIsSlow({ status: 'cancelled', lastMs: 10, startedMs: 0 }, 0)).toBe(true);
});

test('the warning class only reaches the DOM for slow frames', () => {
  const h = harness();
  h.readout.setBusy(true);
  h.readout.setResult(120);
  h.readout.setBusy(false);
  expect(h.text()).not.toContain('is-warn');

  h.readout.setBusy(true);
  h.readout.setResult(880);
  h.readout.setBusy(false);
  expect(h.text()).toContain('is-warn');
  h.readout.dispose();
});
