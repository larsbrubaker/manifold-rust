// Static tables for the Boolean Gallery: the option lists its sidebar renders
// (boolean-panel.ts), the enum-name arrays the debug-info records use, and the
// curated Thingi10K trouble cases. Split out of boolean-gallery.ts so that
// file holds behavior only.

import type { PairKind, ThingiModel } from '../thingi.ts';
import type { OpChoice } from './gallery-widgets.ts';

export const SHAPES = [
  { value: '0', text: 'Cube' },
  { value: '1', text: 'Sphere' },
  { value: '2', text: 'Cylinder' },
  { value: '3', text: 'Spiky Dodecahedron' },
];

/** Venn-glyph cells of the Operation grid; `value` is the wasm op index. */
export const OP_CHOICES: OpChoice[] = [
  { value: 0, glyph: 'union', text: 'Union' },
  { value: 1, glyph: 'intersection', text: 'Intersect' },
  { value: 2, glyph: 'difference', text: 'Difference' },
];

export const PAIR_KINDS = [
  { value: 'mm', text: 'Manifold + Manifold' },
  { value: 'mn', text: 'Manifold + Non-manifold' },
  { value: 'nn', text: 'Non-manifold + Non-manifold' },
];

/// Pair Type is really a mode selector: the three kinds above roll a random
/// pair, while the two below reveal a sub-control that names an exact pair
/// instead — a curated trouble case, or model ids typed in by hand.
export type PairMode = PairKind | 'trouble' | 'pick';

export const PAIR_MODES = PAIR_KINDS.concat([
  { value: 'trouble', text: 'Trouble Cases…' },
  { value: 'pick', text: 'Pick Models…' },
]);

export const isRandomMode = (m: PairMode): m is PairKind =>
  m === 'mm' || m === 'mn' || m === 'nn';

export const OP_COLORS = [0x4488cc, 0x44aa44, 0xcc4444];
export const OP_NAMES = ['union', 'intersection', 'difference'];
export const ENGINE_NAMES = ['exact', 'robust', 'auto'];
export const SHAPE_NAMES = ['cube', 'sphere', 'cylinder', 'spiky-dodecahedron'];

export interface TroubleCase {
  value: string; text: string; a: number; b: number;
  op: string; engine: string; offset: [number, number, number];
  rot: [number, number, number]; pairKind: PairKind;
}

/// Known Thingi10K trouble cases — the pairs the perf and bug hunts keep
/// returning to. Selecting one loads the exact configuration (models by id,
/// op, engine, offsets, rotation) through the same path as pasted debug info.
export const TROUBLE_CASES: TroubleCase[] = [
  {
    value: '1663774-51334',
    text: '1663774 ∪ 51334 — heavy fins (perf stress)',
    a: 1663774, b: 51334, op: 'union', engine: 'robust',
    offset: [0.3, 0, 0], rot: [231.39999999999753, 124, 273.6000000000049], pairKind: 'nn',
  },
  {
    value: '91946-61459',
    text: '91946 ∪ 61459 — doubled-surface windows (regression watch)',
    a: 91946, b: 61459, op: 'union', engine: 'robust',
    offset: [0.3, 0, 0], rot: [236, 231, 42], pairKind: 'nn',
  },
  {
    value: '939888-93557',
    text: '939888 ∪ 93557 — self-overlapping operand',
    a: 939888, b: 93557, op: 'union', engine: 'robust',
    offset: [0.3, 0, 0], rot: [356, 140, 322], pairKind: 'nn',
  },
  {
    value: '1075458-91115',
    text: '1075458 − 91115 — CDT recovery (bit-identity gate)',
    a: 1075458, b: 91115, op: 'difference', engine: 'robust',
    offset: [0.7, -0.2, 0.4], rot: [311, 55, 345], pairKind: 'mn',
  },
  {
    value: '92068-39926',
    text: '92068 ∪ 39926 — tripled facets (fixed NotClosed)',
    a: 92068, b: 39926, op: 'union', engine: 'robust',
    offset: [0.3, 0, 0], rot: [0, 0, 0], pairKind: 'nn',
  },
];

/// Rebuild a ThingiModel from a Copy-Debug-Info operand record: the URL
/// carries the repo and format, so a pasted report reproduces without the
/// metadata index (whose pools might have changed since the capture).
export function modelFromDebugRecord(rec: any): ThingiModel {
  if (typeof rec?.id !== 'number' || typeof rec?.url !== 'string') {
    throw new Error('operand record is missing id/url');
  }
  const repoMatch = /Thingi10K-meshes-(\d+)@/.exec(rec.url);
  const formatMatch = /\.(\w+)\.zip$/.exec(rec.url);
  if (!repoMatch) throw new Error(`unrecognized mesh URL: ${rec.url}`);
  return {
    id: rec.id,
    thing_id: 0,
    name: rec.name ?? `#${rec.id}`,
    format: formatMatch ? formatMatch[1] : 'stl',
    repo: parseInt(repoMatch[1]),
    closed: true,
    edge_manifold: !!rec.edge_manifold,
    vertex_manifold: !!rec.vertex_manifold,
    faces: rec.faces ?? 0,
    vertices: 0,
  };
}
