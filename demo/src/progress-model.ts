// Maps the kernel's per-phase progress notifications onto ONE monotonic
// overall fraction for the busy card's bar.
//
// The kernel (src/progress.rs) reports `(phase, fraction)` where `fraction` is
// the completion of *that phase* — so feeding it straight to the bar makes it
// reset to 0 and refill once per phase, which reads as "stuck in a loop".
// This module keeps a static weight per phase and turns the stream into
//
//   overall = (weight of completed phases + weight[phase] * fraction) / total
//
// clamped so it never moves backwards inside one operation (a run is
// identified by the worker's `seq`; a new seq starts over at 0).
//
// Consumed by boolean-runner.ts, which hands the result to
// BusyIndicator.setPhase(). The phase identifier arrives either as the stable
// display name (what demo/wasm/src/progress.rs forwards today) or as the
// numeric phase id from the FFI phase table — both are accepted so this keeps
// working if the worker switches to ids.

/**
 * Pipeline phases in the order `manifold_rust::progress::Phase` declares them
 * (array index == the stable FFI phase id), with a rough relative duration.
 *
 * THE ONE PLACE weights live. They only need to be right relative to each
 * other; they are a UI pacing hint, never a measurement. Phases the robust
 * engine does not drive with a work total (weight 0) advance the bar to their
 * start and no further.
 */
const PHASES: { name: string; weight: number }[] = [
  { name: 'narrow phase', weight: 5 },
  { name: 'self intersections', weight: 20 },
  { name: 'candidate points', weight: 10 },
  { name: 'registries', weight: 10 },
  { name: 'arrangements', weight: 30 },
  { name: 'cells', weight: 15 },
  { name: 'winding', weight: 5 },
  { name: 'assemble', weight: 5 },
  // The exact engine is a single opaque phase: no weight, so it never claims
  // any share of the bar and stays in the indeterminate sweep.
  { name: 'exact boolean', weight: 0 },
];

const TOTAL_WEIGHT = PHASES.reduce((sum, p) => sum + p.weight, 0);

/** Cumulative weight *before* each phase, i.e. where its span starts. */
const PHASE_START = PHASES.reduce<number[]>((starts, p, i) => {
  starts.push(i === 0 ? 0 : starts[i - 1]! + PHASES[i - 1]!.weight);
  return starts;
}, []);

const BY_NAME = new Map(PHASES.map((p, i) => [p.name, i]));

function phaseIndex(phase: string | number): number {
  if (typeof phase === 'number') return Number.isInteger(phase) ? phase : -1;
  const byName = BY_NAME.get(phase);
  if (byName !== undefined) return byName;
  // Ids may also arrive stringified through postMessage payloads.
  const asNum = Number(phase);
  return Number.isInteger(asNum) && phase.trim() !== '' ? asNum : -1;
}

/** Display name for a phase we know; unknown ids show as-is so a newer kernel
 *  phase still says *something* useful instead of blanking the label. */
export function phaseLabel(phase: string | number): string {
  const i = phaseIndex(phase);
  const known = PHASES[i];
  if (known) return known.name;
  return typeof phase === 'number' ? `phase ${phase}` : String(phase);
}

export interface OverallProgress {
  /** Phase name plus the overall percentage, e.g. `arrangements — 42%`. */
  label: string;
  /** Monotonic overall completion in [0,1], or null while indeterminate. */
  fraction: number | null;
}

/**
 * Accumulates phase updates into a monotonic overall fraction.
 *
 * One instance per consumer; it resets itself whenever a new `seq` (a new
 * operation) shows up, so callers never have to remember to.
 */
export class ProgressModel {
  private seq: number | null = null;
  private best = 0;
  /** Whether any phase has reported a real fraction in this operation. Until
   *  one does, the bar stays in its indeterminate sweep — a determinate 0%
   *  that may never move is worse than an honest "working". */
  private determinate = false;

  /** Forget the current operation's progress (next update starts at 0). */
  reset() {
    this.seq = null;
    this.best = 0;
    this.determinate = false;
  }

  update(seq: number, phase: string | number, fraction: number | null): OverallProgress {
    if (seq !== this.seq) {
      this.reset();
      this.seq = seq;
    }

    const i = phaseIndex(phase);
    const known = PHASES[i];
    if (known && TOTAL_WEIGHT > 0) {
      const start = PHASE_START[i]! / TOTAL_WEIGHT;
      const span = known.weight / TOTAL_WEIGHT;
      let candidate = start;
      if (fraction !== null && Number.isFinite(fraction)) {
        candidate = start + span * Math.max(0, Math.min(1, fraction));
        this.determinate = true;
      }
      // Never step back: phases can report slightly out of order under
      // parallel `advance`, and a bar that retreats is the bug we are fixing.
      this.best = Math.max(this.best, Math.min(1, candidate));
    }
    // Unknown phases (a newer kernel, or the weightless exact boolean) leave
    // the bar exactly where it was.

    return {
      label: this.format(phaseLabel(phase)),
      fraction: this.determinate ? this.best : null,
    };
  }

  private format(name: string): string {
    if (!this.determinate) return name;
    return `${name} — ${Math.round(this.best * 100)}%`;
  }
}
