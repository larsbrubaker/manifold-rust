// Dev script: runs server + watches Rust sources to rebuild WASM automatically

import { spawn } from "bun";
import { watch } from "fs";
import { join, resolve } from "path";

const ROOT = import.meta.dir;
const PROJECT_ROOT = resolve(ROOT, "..");

// Directories to watch for Rust changes
const WATCH_DIRS = [
  join(PROJECT_ROOT, "src"),
  join(ROOT, "wasm", "src"),
];

let building = false;
let pendingBuild = false;

async function buildWasm() {
  if (building) {
    pendingBuild = true;
    return;
  }
  building = true;
  console.log("\n[dev] Rust source changed — rebuilding WASM...");

  const start = performance.now();
  const proc = spawn({
    cmd: ["wasm-pack", "build", join(ROOT, "wasm"), "--target", "web", "--out-dir", join(ROOT, "public", "pkg")],
    cwd: PROJECT_ROOT,
    stdout: "inherit",
    stderr: "inherit",
  });

  const exitCode = await proc.exited;
  const elapsed = ((performance.now() - start) / 1000).toFixed(1);

  if (exitCode === 0) {
    console.log(`[dev] WASM rebuild succeeded (${elapsed}s)`);
  } else {
    console.error(`[dev] WASM rebuild failed (exit code ${exitCode})`);
  }

  building = false;
  if (pendingBuild) {
    pendingBuild = false;
    buildWasm();
  }
}

// Debounce: collect rapid changes into one rebuild
let debounceTimer: ReturnType<typeof setTimeout> | null = null;
const DEBOUNCE_MS = 300;

function onRustChange(filename: string | null) {
  if (filename && !filename.endsWith(".rs")) return;
  if (debounceTimer) clearTimeout(debounceTimer);
  debounceTimer = setTimeout(() => buildWasm(), DEBOUNCE_MS);
}

// Watch Rust source directories
for (const dir of WATCH_DIRS) {
  try {
    watch(dir, { recursive: true }, (_event, filename) => {
      onRustChange(filename as string | null);
    });
    console.log(`[dev] Watching ${dir}`);
  } catch (e: any) {
    console.warn(`[dev] Could not watch ${dir}: ${e.message}`);
  }
}

// Start the dev server
console.log("[dev] Starting server...");
const server = await import("./server.ts");
