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
  has_colors: boolean;
  colors: Float32Array | null; // RGBA interleaved (4 floats per vertex)
  num_vert: number;
  num_tri: number;
  volume: number;
  surface_area: number;
}

function toMeshData(raw: any): MeshData {
  // Extract all data from the WASM struct before freeing it.
  // The getters copy typed arrays to JS heap, so we must capture them first.
  const result: MeshData = {
    positions: raw.positions,
    normals: raw.normals,
    indices: raw.indices,
    has_colors: raw.has_colors ?? false,
    colors: raw.has_colors ? raw.colors : null,
    num_vert: raw.num_vert,
    num_tri: raw.num_tri,
    volume: raw.volume,
    surface_area: raw.surface_area,
  };
  // Free the WASM-side MeshData to prevent memory leak.
  // Without this, each frame leaks ~100KB+ of WASM linear memory at 60fps.
  if (typeof raw.free === 'function') raw.free();
  return result;
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

export function spikyDodecahedronMesh(spikeHeight: number): MeshData {
  return toMeshData(getWasm().spiky_dodecahedron_mesh(spikeHeight));
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

export function booleanGalleryMeshRotated(shapeA: number, shapeB: number, op: number, ox: number, oy: number, oz: number, rx: number, ry: number, rz: number): MeshData {
  return toMeshData(getWasm().boolean_gallery_mesh_rotated(shapeA, shapeB, op, ox, oy, oz, rx, ry, rz));
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

// Test Gallery visualizations
export function testMirrorUnionMesh(): MeshData {
  return toMeshData(getWasm().test_mirror_union_mesh());
}

export function testSplitByPlaneMesh(half: number): MeshData {
  return toMeshData(getWasm().test_split_by_plane_mesh(half));
}

export function testVugMesh(): MeshData {
  return toMeshData(getWasm().test_vug_mesh());
}

export function testWarpMesh(): MeshData {
  return toMeshData(getWasm().test_warp_mesh());
}

export function testSpiralMesh(): MeshData {
  return toMeshData(getWasm().test_spiral_mesh());
}

export function testSphereDiffMesh(): MeshData {
  return toMeshData(getWasm().test_sphere_diff_mesh());
}

export function testCubesUnionMesh(): MeshData {
  return toMeshData(getWasm().test_cubes_union_mesh());
}

export function testBatchSubtractMesh(): MeshData {
  return toMeshData(getWasm().test_batch_subtract_mesh());
}
