// WASM module loader and typed wrappers for Manifold operations

let wasmModule: any = null;

export async function initWasm(): Promise<void> {
  if (wasmModule) return;
  const wasmUrl = new URL('./public/pkg/manifold_wasm.js', window.location.href).href;
  const mod = await import(wasmUrl);
  await mod.default();
  wasmModule = mod;
}

function getWasm(): any {
  if (!wasmModule) throw new Error('WASM not initialized. Call initWasm() first.');
  return wasmModule;
}

export interface MeshData {
  positions: Float32Array;
  normals: Float32Array;
  indices: Uint32Array;
  num_vert: number;
  num_tri: number;
  volume: number;
  surface_area: number;
}

function toMeshData(raw: any): MeshData {
  return {
    positions: raw.positions,
    normals: raw.normals,
    indices: raw.indices,
    num_vert: raw.num_vert,
    num_tri: raw.num_tri,
    volume: raw.volume,
    surface_area: raw.surface_area,
  };
}

export function version(): string {
  return getWasm().version();
}

export function cubeMesh(sx: number, sy: number, sz: number, center: boolean): MeshData {
  return toMeshData(getWasm().cube_mesh(sx, sy, sz, center));
}

export function sphereMesh(radius: number, segments: number): MeshData {
  return toMeshData(getWasm().sphere_mesh(radius, segments));
}

export function cylinderMesh(height: number, radiusLow: number, radiusHigh: number, segments: number): MeshData {
  return toMeshData(getWasm().cylinder_mesh(height, radiusLow, radiusHigh, segments));
}

export function tetrahedronMesh(): MeshData {
  return toMeshData(getWasm().tetrahedron_mesh());
}

export function extrudeMesh(radius: number, segments: number, height: number): MeshData {
  return toMeshData(getWasm().extrude_mesh(radius, segments, height));
}

export function revolveMesh(radius: number, segments: number, degrees: number): MeshData {
  return toMeshData(getWasm().revolve_mesh(radius, segments, degrees));
}

export function hullMesh(points: Float32Array): MeshData {
  return toMeshData(getWasm().hull_mesh(points));
}

export function unionMesh(offsetX: number): MeshData {
  return toMeshData(getWasm().union_mesh(offsetX));
}

export function intersectMesh(offsetX: number): MeshData {
  return toMeshData(getWasm().intersect_mesh(offsetX));
}

export function differenceMesh(offsetX: number): MeshData {
  return toMeshData(getWasm().difference_mesh(offsetX));
}

export function mengerSpongeMesh(depth: number): MeshData {
  return toMeshData(getWasm().menger_sponge_mesh(depth));
}

export function booleanGalleryMesh(shapeA: number, shapeB: number, op: number, ox: number, oy: number, oz: number): MeshData {
  return toMeshData(getWasm().boolean_gallery_mesh(shapeA, shapeB, op, ox, oy, oz));
}

export function refinedShapeMesh(shape: number, refineLevel: number): MeshData {
  return toMeshData(getWasm().refined_shape_mesh(shape, refineLevel));
}

export function extrudeTwistMesh(radius: number, segments: number, height: number, twistDegrees: number, nDivisions: number, scaleTop: number): MeshData {
  return toMeshData(getWasm().extrude_twist_mesh(radius, segments, height, twistDegrees, nDivisions, scaleTop));
}

export function revolvePartialMesh(profile: number, segments: number, degrees: number): MeshData {
  return toMeshData(getWasm().revolve_partial_mesh(profile, segments, degrees));
}
