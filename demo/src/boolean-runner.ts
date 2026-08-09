// Main-thread wrapper around boolean-worker.ts for the Boolean Gallery.
//
// Serializes worker access: at most one boolean run is in flight, and run
// requests made meanwhile coalesce (latest wins) so slider drags and
// animation never queue up stale frames. cancel() implements "stop this
// computation": wasm cannot be preempted, so the worker is terminated and a
// fresh one spawned, with the last imported operands re-sent transparently.

import type { MeshData, BooleanEngine } from './wasm.ts';
import { BusyIndicator } from './busy-indicator.ts';

const ENGINE_NAMES = ['Exact', 'Robust', 'Auto'];
const OP_NAMES = ['union', 'intersection', 'difference'];

function describeRun(params: RunParams): string {
  const engine = ENGINE_NAMES[params.engine] ?? `Engine ${params.engine}`;
  const op = OP_NAMES[params.op] ?? `op ${params.op}`;
  return `Computing… ${engine} ${op}`;
}

export interface OperandArrays {
  positions: Float32Array;
  indices: Uint32Array;
}

export interface ImportInfo {
  ok: boolean;
  status: string;
  is_soup: boolean;
  /** Own triangles genuinely intersect (crossing, overlap, or coincident
   *  surface) — why Auto may pick Robust. null when it was not asked for. */
  self_intersecting: boolean | null;
  num_vert: number;
  num_tri: number;
}

export interface RunParams {
  source: 'thingi' | 'builtin';
  shapeA?: number;
  shapeB?: number;
  op: number;
  engine: BooleanEngine;
  /** Rewind inside-out shells of both operands before the boolean. */
  repair: boolean;
  ox: number; oy: number; oz: number;
  rx: number; ry: number; rz: number;
  /** Opaque payload handed back with the result (debug info, flags). */
  tag?: any;
}

export interface RunResult {
  data: MeshData;
  /** Worker-side wall time of the boolean evaluation itself. */
  elapsedMs: number;
}

export class BooleanRunner {
  private worker: Worker | null = null;
  private seq = 0;
  private importWaiters: ((info: ImportInfo) => void)[] = [];
  private runInFlight = false;
  private currentRun: RunParams | null = null;
  private queuedRun: RunParams | null = null;
  private operandCache: { a: OperandArrays | null; b: OperandArrays | null } = { a: null, b: null };
  private disposed = false;
  /** Shared busy overlay; owned here so every demo using the runner gets
   *  feedback without wiring anything up itself. */
  private indicator: BusyIndicator | null;

  constructor(showBusyIndicator = true) {
    this.indicator = showBusyIndicator ? new BusyIndicator() : null;
    this.indicator?.setCancel(() => this.cancel());
  }

  /** Single funnel for busy transitions: keeps the overlay in sync no matter
   *  what the consumer does with the public onBusyChange hook. */
  private setBusy(busy: boolean, params: RunParams | null) {
    if (busy && params) this.indicator?.show(describeRun(params));
    else this.indicator?.hide();
    this.onBusyChange(busy);
  }

  onResult: (result: RunResult, params: RunParams) => void = () => {};
  onError: (error: string, params: RunParams) => void = () => {};
  onBusyChange: (busy: boolean) => void = () => {};
  onCancelled: () => void = () => {};

  private ensureWorker(): Worker {
    if (!this.worker) {
      // Dev server serves /src/*.ts transpiled per file; the production
      // build (build.ts) bundles the worker as its own entrypoint next to
      // main.js. import.meta.url still ends in .ts exactly in dev.
      const url = import.meta.url.endsWith('.ts')
        ? new URL('./boolean-worker.ts', import.meta.url)
        : new URL('./boolean-worker.js', window.location.href);
      this.worker = new Worker(url, { type: 'module' });
      this.worker.onmessage = (ev: MessageEvent) => this.handleMessage(ev.data);
      // Same page-relative pkg resolution as wasm.ts — correct under both
      // localhost and the GitHub Pages subpath.
      const pkgUrl = new URL('./public/pkg/manifold_wasm.js', window.location.href).href;
      this.worker.postMessage({ type: 'init', pkgUrl });
    }
    return this.worker;
  }

  private handleMessage(msg: any) {
    if (msg.type === 'imported') {
      console.debug('[boolean-runner] imported', msg.slot, msg.status, `${msg.num_tri} tris`);
      const waiter = this.importWaiters.shift();
      waiter?.(msg as ImportInfo);
    } else if (msg.type === 'result') {
      const params = this.currentRun!;
      this.runInFlight = false;
      this.currentRun = null;
      if (this.queuedRun) {
        const next = this.queuedRun;
        this.queuedRun = null;
        // Still busy — retitle in place (the elapsed clock keeps running so a
        // drag that coalesces frames reads as one continuous computation).
        this.indicator?.show(describeRun(next));
        this.send(next);
      } else {
        this.setBusy(false, null);
      }
      if (msg.ok) {
        this.onResult({ data: msg.data as MeshData, elapsedMs: msg.elapsedMs }, params);
      } else {
        this.onError(msg.error, params);
      }
    }
  }

  /** Import one operand into the worker; the arrays are also cached so a
   *  cancel can re-import them into the replacement worker. */
  importOperand(
    slot: 'a' | 'b',
    arrays: OperandArrays,
    wantSelfIntersecting = false,
  ): Promise<ImportInfo> {
    this.operandCache[slot] = arrays;
    const w = this.ensureWorker();
    return new Promise<ImportInfo>(resolve => {
      this.importWaiters.push(resolve);
      // Structured clone (no transfer): the main thread keeps its copy.
      w.postMessage({
        type: 'import', slot, positions: arrays.positions, indices: arrays.indices,
        want_self_intersecting: wantSelfIntersecting,
      });
    });
  }

  /** Request a boolean evaluation. If one is already running the request is
   *  queued, replacing any not-yet-started earlier request (latest wins). */
  requestRun(params: RunParams) {
    if (this.runInFlight) {
      this.queuedRun = params;
    } else {
      this.setBusy(true, params);
      this.send(params);
    }
  }

  private send(params: RunParams) {
    this.runInFlight = true;
    this.currentRun = params;
    const { tag, ...rest } = params;
    this.ensureWorker().postMessage({ type: 'run', seq: ++this.seq, ...rest });
  }

  get busy(): boolean {
    return this.runInFlight;
  }

  /** Hard-stop the current computation: terminate the worker, spawn a fresh
   *  one, and re-import the cached operands so the next run just works. */
  cancel() {
    if (!this.worker) return;
    this.worker.terminate();
    this.worker = null;
    this.runInFlight = false;
    this.currentRun = null;
    this.queuedRun = null;
    // Unblock any import awaiters with a failure status.
    for (const waiter of this.importWaiters.splice(0)) {
      waiter({ ok: false, status: 'cancelled', is_soup: false, self_intersecting: null, num_vert: 0, num_tri: 0 });
    }
    this.setBusy(false, null);
    this.onCancelled();
    if (this.disposed) return;
    // Warm the replacement worker with the operands (fire and forget — the
    // resolved infos were already reported when first imported).
    if (this.operandCache.a) this.importOperand('a', this.operandCache.a);
    if (this.operandCache.b) this.importOperand('b', this.operandCache.b);
  }

  dispose() {
    this.disposed = true;
    this.worker?.terminate();
    this.worker = null;
    this.indicator?.dispose();
    this.indicator = null;
  }
}
