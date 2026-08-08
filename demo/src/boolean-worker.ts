// Boolean Gallery web worker — owns its own wasm instance plus the two
// imported mesh operands, so long boolean evaluations never block the UI
// thread. The main-thread side is boolean-runner.ts.
//
// Worker-bound messages:
//   { type: 'init', pkgUrl }   — must be first; the page-relative wasm pkg
//     URL is computed on the main thread (the worker's own location differs
//     between the dev server and the GitHub Pages subpath).
//   { type: 'import', slot: 'a'|'b', positions, indices }
//   { type: 'run', seq, source: 'thingi'|'builtin', shapeA?, shapeB?,
//     op, engine, repair, ox, oy, oz, rx, ry, rz }
// Replies (exactly one per import/run message, in order):
//   { type: 'imported', slot, ok, status, is_soup, num_vert, num_tri }
//   { type: 'result', seq, ok, data?, error?, elapsedMs }
//
// There is no cancel message: wasm computations cannot be preempted, so the
// main thread cancels by terminating this worker outright and spawning a
// fresh one (the runner re-imports the operands transparently).

let wasmReady: Promise<any> | null = null;
let pkgUrl: string | null = null;
const handles: { a: any; b: any } = { a: null, b: null };

/// Memoized on the PROMISE, not the result: concurrent callers must share
/// one wasm instantiation — calling mod.default() twice re-initializes the
/// module and detaches every handle created against the first instance.
function ensureWasm(): Promise<any> {
  if (!wasmReady) {
    if (!pkgUrl) throw new Error('worker not initialized (missing init message)');
    const url = pkgUrl;
    wasmReady = (async () => {
      const mod = await import(/* @vite-ignore */ url);
      await mod.default();
      return mod;
    })();
  }
  return wasmReady;
}

// Handlers are async (they await the wasm load), and the browser delivers
// the next message as soon as a handler suspends — so chain them: replies
// must go out in message order (the runner's FIFO waiters depend on it),
// and a run must never start between the two imports it follows.
let chain: Promise<void> = Promise.resolve();
self.onmessage = (ev: MessageEvent) => {
  chain = chain.then(() => handleMessage(ev.data));
};

async function handleMessage(msg: any) {
  if (msg.type === 'init') {
    pkgUrl = msg.pkgUrl;
    return;
  }
  try {
    const w = await ensureWasm();
    if (msg.type === 'import') {
      const slot = msg.slot as 'a' | 'b';
      if (handles[slot]) {
        handles[slot].free();
        handles[slot] = null;
      }
      const h = new w.ImportedMesh(msg.positions, msg.indices);
      const reply = {
        type: 'imported', slot,
        ok: h.ok, status: h.status, is_soup: h.is_soup,
        num_vert: h.num_vert, num_tri: h.num_tri,
      };
      if (h.ok) {
        handles[slot] = h;
      } else {
        h.free();
      }
      (self as any).postMessage(reply);
    } else if (msg.type === 'run') {
      const t0 = performance.now();
      let raw: any;
      if (msg.source === 'thingi') {
        if (!handles.a || !handles.b) throw new Error('operands not imported');
        raw = w.imported_boolean(
          handles.a, handles.b, msg.op, msg.engine,
          msg.ox, msg.oy, msg.oz, msg.rx, msg.ry, msg.rz, msg.repair ?? false,
        );
      } else {
        raw = w.boolean_gallery_mesh_repair(
          msg.shapeA, msg.shapeB, msg.op,
          msg.ox, msg.oy, msg.oz, msg.rx, msg.ry, msg.rz, msg.engine,
          msg.repair ?? false,
        );
      }
      // Same extraction as wasm.ts toMeshData: capture the typed-array
      // copies, then free the wasm-side struct.
      const data = {
        positions: raw.positions as Float32Array,
        normals: raw.normals as Float32Array,
        indices: raw.indices as Uint32Array,
        has_colors: (raw.has_colors ?? false) as boolean,
        colors: raw.has_colors ? (raw.colors as Float32Array) : null,
        num_vert: raw.num_vert as number,
        num_tri: raw.num_tri as number,
        volume: raw.volume as number,
        surface_area: raw.surface_area as number,
      };
      if (typeof raw.free === 'function') raw.free();
      const transfers: ArrayBuffer[] = [
        data.positions.buffer as ArrayBuffer,
        data.normals.buffer as ArrayBuffer,
        data.indices.buffer as ArrayBuffer,
      ];
      if (data.colors) transfers.push(data.colors.buffer as ArrayBuffer);
      (self as any).postMessage(
        { type: 'result', seq: msg.seq, ok: true, data, elapsedMs: performance.now() - t0 },
        transfers,
      );
    }
  } catch (e) {
    if (msg.type === 'run') {
      (self as any).postMessage({ type: 'result', seq: msg.seq, ok: false, error: String(e) });
    } else if (msg.type === 'import') {
      (self as any).postMessage({
        type: 'imported', slot: msg.slot,
        ok: false, status: String(e), is_soup: false, num_vert: 0, num_tri: 0,
      });
    }
  }
}
