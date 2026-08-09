// Info panel for the Boolean Gallery: the Vertices / Triangles / Volume /
// Surface Area rows plus the "Last Frame" timing row.
//
// Split out of boolean-gallery.ts because the timing row has a lifecycle of
// its own. "Last Frame" is the *boolean's* wall time (BooleanRunner reports
// elapsedMs per run), not the render frame time, so it only has a new value
// when a run finishes. Rendering it from the result handler alone made it
// freeze on a stale number for the whole (possibly multi-second) run, and
// leave that number standing after a failed or cancelled run — the panel
// claimed a frame time that no longer described anything on screen.
//
// This class owns that state instead: while a run is in flight the row shows
// a ticking "computing…" clock, a failure or cancellation is labelled as
// such, and only a completed run prints a duration. The mesh rows are
// re-rendered from the cached MeshData so the timing row can update on its
// own schedule.

import type { MeshData } from '../wasm.ts';
import { updateReadout } from '../controls.ts';

export type FrameStatus = 'idle' | 'computing' | 'failed' | 'cancelled';

export interface FrameStat {
  status: FrameStatus;
  /** Wall time of the last *completed* boolean, ms; null before the first. */
  lastMs: number | null;
  /** Clock start of the run in flight (only meaningful while computing). */
  startedMs: number;
}

export function formatMs(ms: number): string {
  return ms >= 1000 ? `${(ms / 1000).toFixed(2)} s` : `${ms.toFixed(1)} ms`;
}

/// The "Last Frame" row, or null when there is nothing truthful to say yet
/// (no run has ever started). Exported for tests.
export function frameStatRow(stat: FrameStat, nowMs: number): { label: string; value: string } | null {
  const previous = stat.lastMs === null ? '' : ` (last ${formatMs(stat.lastMs)})`;
  switch (stat.status) {
    case 'computing': {
      const elapsed = Math.max(0, nowMs - stat.startedMs);
      return { label: 'Last Frame', value: `computing… ${formatMs(elapsed)}${previous}` };
    }
    case 'failed':
      return { label: 'Last Frame', value: `failed${previous}` };
    case 'cancelled':
      return { label: 'Last Frame', value: `cancelled${previous}` };
    default:
      return stat.lastMs === null ? null : { label: 'Last Frame', value: formatMs(stat.lastMs) };
  }
}

/** Repaint interval of the ticking clock while a boolean computes. */
const TICK_MS = 100;

export class BooleanReadout {
  private mesh: MeshData | null = null;
  private stat: FrameStat = { status: 'idle', lastMs: null, startedMs: 0 };
  private tick: ReturnType<typeof setInterval> | null = null;

  /** `now` is injectable so tests can drive the clock deterministically. */
  constructor(
    readonly element: HTMLElement,
    private readonly now: () => number = () => performance.now(),
  ) {}

  /** Latest boolean output; the timing row is preserved. */
  setMesh(data: MeshData) {
    this.mesh = data;
    this.render();
  }

  /** Empty the panel (a failed run with nothing valid left to show). */
  clear() {
    this.mesh = null;
    this.render();
  }

  /** Busy transitions from BooleanRunner. Entering busy starts the clock;
   *  leaving it ends the pending state (the result/error hooks below then say
   *  how the run ended). */
  setBusy(busy: boolean) {
    if (busy) {
      this.startClock();
    } else {
      this.stopClock();
      if (this.stat.status === 'computing') this.stat.status = 'idle';
    }
    this.render();
  }

  /** A run completed. If the runner is still busy (a coalesced request took
   *  its place without a busy transition) the clock restarts for that run. */
  setResult(elapsedMs: number) {
    this.stat.lastMs = elapsedMs;
    if (this.stat.status === 'computing') this.stat.startedMs = this.now();
    this.render();
  }

  /** A run failed. Same coalescing rule as setResult: a follow-on run already
   *  in flight keeps the row on its live clock rather than reporting a stale
   *  failure. */
  setFailed() {
    if (this.stat.status === 'computing') this.stat.startedMs = this.now();
    else this.stat.status = 'failed';
    this.render();
  }

  setCancelled() {
    this.stopClock();
    this.stat.status = 'cancelled';
    this.render();
  }

  dispose() {
    this.stopClock();
  }

  private startClock() {
    this.stat.status = 'computing';
    this.stat.startedMs = this.now();
    if (this.tick === null) this.tick = setInterval(() => this.render(), TICK_MS);
  }

  private stopClock() {
    if (this.tick !== null) {
      clearInterval(this.tick);
      this.tick = null;
    }
  }

  private render() {
    const rows: { label: string; value: string }[] = [];
    const data = this.mesh;
    if (data) {
      rows.push(
        { label: 'Vertices', value: String(data.num_vert) },
        { label: 'Triangles', value: String(data.num_tri) },
        { label: 'Volume', value: data.volume.toFixed(4) },
        { label: 'Surface Area', value: data.surface_area.toFixed(4) },
      );
    }
    const frame = frameStatRow(this.stat, this.now());
    if (frame) rows.push(frame);
    updateReadout(this.element, rows);
  }
}
