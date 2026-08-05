// Thingi10K dataset access for the Boolean Gallery's real-world mesh mode.
//
// Loads the 10K-model metadata index from the Thingi10K browser repo (CDN),
// filters it into manifold / non-manifold pools of reasonably sized closed
// STL models, and fetches + decompresses + parses individual meshes into
// triangle soup ready for the WASM ImportedMesh constructor.
//
// Mesh zips live in three companion repos served via jsDelivr; the URL
// scheme matches the Thingi10K browser app (larsbrubaker/Thingi10K).

import { unzip } from 'fflate';

export interface ThingiModel {
  id: number;
  thing_id: number;
  name: string;
  format: string;
  repo: number;
  closed: boolean;
  edge_manifold: boolean;
  vertex_manifold: boolean;
  faces: number;
  vertices: number;
}

export interface ParsedMesh {
  positions: Float32Array;   // xyz interleaved triangle soup, normalized
  indices: Uint32Array;      // 3 per triangle
  // Normalization applied (original = normalized / scale + center),
  // recorded for debug reproduction.
  center: [number, number, number];
  scale: number;
}

export type PairKind = 'mm' | 'mn' | 'nn';

const METADATA_URL = 'https://cdn.jsdelivr.net/gh/larsbrubaker/Thingi10K@main/docs/data/models.json';
const MESH_CDN = 'https://cdn.jsdelivr.net/gh/larsbrubaker';

// Keep interactive: robust (exact-rational) booleans on big meshes get slow.
const MAX_FACES = 20000;
const MIN_FACES = 4;

export function meshZipUrl(model: ThingiModel): string {
  return `${MESH_CDN}/Thingi10K-meshes-${model.repo}@main/meshes/${model.id}.${model.format}.zip`;
}

let manifoldPool: ThingiModel[] | null = null;
let nonManifoldPool: ThingiModel[] | null = null;
let metadataPromise: Promise<void> | null = null;

/** Fetch and index the metadata once per session. */
async function ensureMetadata(): Promise<void> {
  if (manifoldPool && nonManifoldPool) return;
  if (!metadataPromise) {
    metadataPromise = (async () => {
      const resp = await fetch(METADATA_URL);
      if (!resp.ok) throw new Error(`Thingi10K metadata fetch failed: HTTP ${resp.status}`);
      const all: ThingiModel[] = await resp.json();
      const usable = all.filter(m =>
        m.format === 'stl' && m.closed && m.faces >= MIN_FACES && m.faces <= MAX_FACES);
      manifoldPool = usable.filter(m => m.edge_manifold && m.vertex_manifold);
      nonManifoldPool = usable.filter(m => !(m.edge_manifold && m.vertex_manifold));
    })().catch(e => { metadataPromise = null; throw e; });
  }
  return metadataPromise;
}

function randomFrom(pool: ThingiModel[], exclude?: ThingiModel): ThingiModel {
  let pick = pool[Math.floor(Math.random() * pool.length)];
  // Avoid booleaning a model against itself when both come from one pool.
  while (exclude && pick.id === exclude.id && pool.length > 1) {
    pick = pool[Math.floor(Math.random() * pool.length)];
  }
  return pick;
}

/**
 * Pick one random model of the given manifoldness ('m' = manifold per the
 * dataset flags, 'n' = non-manifold). Callers retry with a fresh pick when a
 * model fails import — some dataset "closed" flags disagree with the robust
 * importer's stricter closed check.
 */
export async function pickRandomModel(kind: 'm' | 'n', exclude?: ThingiModel): Promise<ThingiModel> {
  await ensureMetadata();
  return randomFrom(kind === 'm' ? manifoldPool! : nonManifoldPool!, exclude);
}

/** Operand manifoldness for each pair combination. */
export function pairOperandKinds(kind: PairKind): ['m' | 'n', 'm' | 'n'] {
  return kind === 'mm' ? ['m', 'm'] : kind === 'mn' ? ['m', 'n'] : ['n', 'n'];
}

async function fetchAndUnzip(url: string): Promise<ArrayBuffer> {
  const resp = await fetch(url);
  if (!resp.ok) throw new Error(`HTTP ${resp.status} — ${url}`);
  const compressed = new Uint8Array(await resp.arrayBuffer());
  return new Promise((resolve, reject) => {
    unzip(compressed, (err, files) => {
      if (err) return reject(err);
      const data = Object.values(files)[0];
      if (!data) return reject(new Error(`Empty zip — ${url}`));
      resolve(data.buffer.slice(data.byteOffset, data.byteOffset + data.byteLength) as ArrayBuffer);
    });
  });
}

// Some dataset files declare one more face than the data holds (truncated by
// exactly one 50-byte record); clamp the header so parsing doesn't run off
// the end. Same repair the Thingi10K browser applies.
function fixTruncatedBinaryStl(buffer: ArrayBuffer): ArrayBuffer {
  if (buffer.byteLength < 84) return buffer;
  const view = new DataView(buffer);
  const nFaces = view.getUint32(80, true);
  if (80 + 4 + nFaces * 50 <= buffer.byteLength) return buffer;
  const actualFaces = Math.floor((buffer.byteLength - 84) / 50);
  const fixed = buffer.slice(0);
  new DataView(fixed).setUint32(80, actualFaces, true);
  return fixed;
}

function isAsciiStl(buffer: ArrayBuffer): boolean {
  const head = new TextDecoder().decode(new Uint8Array(buffer, 0, Math.min(512, buffer.byteLength)));
  return head.trimStart().startsWith('solid') && head.includes('facet');
}

function parseBinaryStl(buffer: ArrayBuffer): Float32Array {
  const view = new DataView(buffer);
  const nFaces = view.getUint32(80, true);
  const positions = new Float32Array(nFaces * 9);
  let off = 84;
  for (let f = 0; f < nFaces; f++) {
    off += 12; // skip facet normal
    for (let i = 0; i < 9; i++) {
      positions[f * 9 + i] = view.getFloat32(off, true);
      off += 4;
    }
    off += 2; // attribute byte count
  }
  return positions;
}

function parseAsciiStl(buffer: ArrayBuffer): Float32Array {
  const text = new TextDecoder().decode(buffer);
  const coords: number[] = [];
  const re = /vertex\s+([-+eE\d.]+)\s+([-+eE\d.]+)\s+([-+eE\d.]+)/g;
  let match: RegExpExecArray | null;
  while ((match = re.exec(text)) !== null) {
    coords.push(parseFloat(match[1]), parseFloat(match[2]), parseFloat(match[3]));
  }
  return new Float32Array(coords);
}

/**
 * Fetch, decompress, parse, and normalize a Thingi10K STL model.
 * The mesh is centered at the origin and uniformly scaled so its longest
 * bounding-box side is 2 units — comparable to the gallery's primitives,
 * so the existing offset sliders stay meaningful.
 */
export async function fetchMesh(model: ThingiModel): Promise<ParsedMesh> {
  let buffer = await fetchAndUnzip(meshZipUrl(model));
  let positions: Float32Array;
  if (isAsciiStl(buffer)) {
    positions = parseAsciiStl(buffer);
  } else {
    buffer = fixTruncatedBinaryStl(buffer);
    positions = parseBinaryStl(buffer);
  }
  const nVerts = Math.floor(positions.length / 3);
  if (nVerts < 3) throw new Error(`No triangles parsed from ${model.id}.stl`);

  // Bounding box -> center + uniform scale to a 2-unit box.
  let minX = Infinity, minY = Infinity, minZ = Infinity;
  let maxX = -Infinity, maxY = -Infinity, maxZ = -Infinity;
  for (let i = 0; i < nVerts; i++) {
    const x = positions[i * 3], y = positions[i * 3 + 1], z = positions[i * 3 + 2];
    if (x < minX) minX = x; if (x > maxX) maxX = x;
    if (y < minY) minY = y; if (y > maxY) maxY = y;
    if (z < minZ) minZ = z; if (z > maxZ) maxZ = z;
  }
  const cx = (minX + maxX) / 2, cy = (minY + maxY) / 2, cz = (minZ + maxZ) / 2;
  const maxSide = Math.max(maxX - minX, maxY - minY, maxZ - minZ);
  const scale = maxSide > 0 ? 2 / maxSide : 1;
  for (let i = 0; i < nVerts; i++) {
    positions[i * 3] = (positions[i * 3] - cx) * scale;
    positions[i * 3 + 1] = (positions[i * 3 + 1] - cy) * scale;
    positions[i * 3 + 2] = (positions[i * 3 + 2] - cz) * scale;
  }

  const indices = new Uint32Array(nVerts);
  for (let i = 0; i < nVerts; i++) indices[i] = i;

  return { positions, indices, center: [cx, cy, cz], scale };
}
