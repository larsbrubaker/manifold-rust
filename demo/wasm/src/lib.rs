// WASM bindings for manifold-rust
// Exposes a compact browser-facing API for the demo site.

use manifold_rust::cross_section::CrossSection;
use manifold_rust::linalg::{Vec2, Vec3};
use manifold_rust::manifold::Manifold;
use manifold_rust::quickhull;
use wasm_bindgen::prelude::*;
use js_sys::{Float32Array, Uint32Array};

#[wasm_bindgen]
pub struct MeshSummary {
    num_vert: u32,
    num_tri: u32,
    volume: f64,
    surface_area: f64,
}

#[wasm_bindgen]
impl MeshSummary {
    #[wasm_bindgen(getter)]
    pub fn num_vert(&self) -> u32 {
        self.num_vert
    }

    #[wasm_bindgen(getter)]
    pub fn num_tri(&self) -> u32 {
        self.num_tri
    }

    #[wasm_bindgen(getter)]
    pub fn volume(&self) -> f64 {
        self.volume
    }

    #[wasm_bindgen(getter)]
    pub fn surface_area(&self) -> f64 {
        self.surface_area
    }
}

// MeshSummary retained for lightweight queries that don't need geometry.
fn summarize(mesh: &Manifold) -> MeshSummary {
    MeshSummary {
        num_vert: mesh.num_vert() as u32,
        num_tri: mesh.num_tri() as u32,
        volume: mesh.volume(),
        surface_area: mesh.surface_area(),
    }
}

#[wasm_bindgen]
pub fn cube_summary(size_x: f64, size_y: f64, size_z: f64, center: bool) -> MeshSummary {
    summarize(&Manifold::cube(Vec3::new(size_x, size_y, size_z), center))
}

/// MeshData holds geometry suitable for Three.js BufferGeometry.
/// Access positions, normals, and indices as typed arrays.
#[wasm_bindgen]
pub struct MeshData {
    positions: Vec<f32>,
    normals: Vec<f32>,
    indices: Vec<u32>,
    num_vert: u32,
    num_tri: u32,
    volume: f64,
    surface_area: f64,
}

#[wasm_bindgen]
impl MeshData {
    #[wasm_bindgen(getter)]
    pub fn positions(&self) -> Float32Array {
        Float32Array::from(&self.positions[..])
    }

    #[wasm_bindgen(getter)]
    pub fn normals(&self) -> Float32Array {
        Float32Array::from(&self.normals[..])
    }

    #[wasm_bindgen(getter)]
    pub fn indices(&self) -> Uint32Array {
        Uint32Array::from(&self.indices[..])
    }

    #[wasm_bindgen(getter)]
    pub fn num_vert(&self) -> u32 {
        self.num_vert
    }

    #[wasm_bindgen(getter)]
    pub fn num_tri(&self) -> u32 {
        self.num_tri
    }

    #[wasm_bindgen(getter)]
    pub fn volume(&self) -> f64 {
        self.volume
    }

    #[wasm_bindgen(getter)]
    pub fn surface_area(&self) -> f64 {
        self.surface_area
    }
}

fn mesh_data_from(m: &Manifold) -> MeshData {
    let gl = m.get_mesh_gl(0);
    let num_prop = gl.num_prop as usize;
    let vert_count = if num_prop > 0 { gl.vert_properties.len() / num_prop } else { 0 };
    let tri_count = gl.tri_verts.len() / 3;

    // Extract positions
    let mut positions = Vec::with_capacity(vert_count * 3);
    for i in 0..vert_count {
        let base = i * num_prop;
        positions.push(gl.vert_properties[base]);
        positions.push(gl.vert_properties[base + 1]);
        positions.push(gl.vert_properties[base + 2]);
    }

    // Compute face normals and accumulate per-vertex normals
    let mut normals = vec![0.0f32; vert_count * 3];
    for t in 0..tri_count {
        let i0 = gl.tri_verts[3 * t] as usize;
        let i1 = gl.tri_verts[3 * t + 1] as usize;
        let i2 = gl.tri_verts[3 * t + 2] as usize;

        let ax = positions[i1 * 3] - positions[i0 * 3];
        let ay = positions[i1 * 3 + 1] - positions[i0 * 3 + 1];
        let az = positions[i1 * 3 + 2] - positions[i0 * 3 + 2];
        let bx = positions[i2 * 3] - positions[i0 * 3];
        let by = positions[i2 * 3 + 1] - positions[i0 * 3 + 1];
        let bz = positions[i2 * 3 + 2] - positions[i0 * 3 + 2];

        let nx = ay * bz - az * by;
        let ny = az * bx - ax * bz;
        let nz = ax * by - ay * bx;

        for idx in [i0, i1, i2] {
            normals[idx * 3] += nx;
            normals[idx * 3 + 1] += ny;
            normals[idx * 3 + 2] += nz;
        }
    }

    // Normalize
    for i in 0..vert_count {
        let x = normals[i * 3];
        let y = normals[i * 3 + 1];
        let z = normals[i * 3 + 2];
        let len = (x * x + y * y + z * z).sqrt();
        if len > 1e-10 {
            normals[i * 3] = x / len;
            normals[i * 3 + 1] = y / len;
            normals[i * 3 + 2] = z / len;
        }
    }

    MeshData {
        positions,
        normals,
        indices: gl.tri_verts,
        num_vert: vert_count as u32,
        num_tri: tri_count as u32,
        volume: m.volume(),
        surface_area: m.surface_area(),
    }
}

// ---------------------------------------------------------------------------
// Summary-only API (existing)
// ---------------------------------------------------------------------------

#[wasm_bindgen]
pub fn version() -> String {
    env!("CARGO_PKG_VERSION").to_string()
}

// ---------------------------------------------------------------------------
// Mesh-geometry API (for Three.js rendering)
// ---------------------------------------------------------------------------

#[wasm_bindgen]
pub fn cube_mesh(size_x: f64, size_y: f64, size_z: f64, center: bool) -> MeshData {
    let m = Manifold::cube(Vec3::new(size_x, size_y, size_z), center);
    mesh_data_from(&m)
}

#[wasm_bindgen]
pub fn sphere_mesh(radius: f64, circular_segments: i32) -> MeshData {
    let m = Manifold::sphere(radius, circular_segments);
    mesh_data_from(&m)
}

#[wasm_bindgen]
pub fn cylinder_mesh(height: f64, radius_low: f64, radius_high: f64, circular_segments: i32) -> MeshData {
    let m = Manifold::cylinder(height, radius_low, radius_high, circular_segments);
    mesh_data_from(&m)
}

#[wasm_bindgen]
pub fn tetrahedron_mesh() -> MeshData {
    let m = Manifold::tetrahedron();
    mesh_data_from(&m)
}

#[wasm_bindgen]
pub fn extrude_mesh(radius: f64, segments: i32, height: f64) -> MeshData {
    let cs = CrossSection::circle(radius, segments);
    let m = Manifold::extrude(&cs.to_polygons(), height, 0, 0.0, Vec2::new(1.0, 1.0));
    mesh_data_from(&m)
}

#[wasm_bindgen]
pub fn revolve_mesh(radius: f64, segments: i32, revolve_degrees: f64) -> MeshData {
    let cs = CrossSection::circle(radius, segments);
    let m = Manifold::revolve(&cs.to_polygons(), segments, revolve_degrees);
    mesh_data_from(&m)
}

#[wasm_bindgen]
pub fn hull_mesh(points_flat: &[f32]) -> MeshData {
    let points: Vec<Vec3> = points_flat
        .chunks_exact(3)
        .map(|c| Vec3::new(c[0] as f64, c[1] as f64, c[2] as f64))
        .collect();
    let imp = quickhull::convex_hull(&points);
    let m = Manifold::from_impl(imp);
    mesh_data_from(&m)
}

#[wasm_bindgen]
pub fn union_mesh(offset_x: f64) -> MeshData {
    let a = Manifold::cube(Vec3::new(1.0, 1.0, 1.0), false);
    let b = Manifold::cube(Vec3::new(1.0, 1.0, 1.0), false).translate(Vec3::new(offset_x, 0.0, 0.0));
    let result = a.union(&b);
    mesh_data_from(&result)
}

#[wasm_bindgen]
pub fn intersect_mesh(offset_x: f64) -> MeshData {
    let a = Manifold::cube(Vec3::new(1.0, 1.0, 1.0), false);
    let b = Manifold::cube(Vec3::new(1.0, 1.0, 1.0), false).translate(Vec3::new(offset_x, 0.0, 0.0));
    let result = a.intersection(&b);
    mesh_data_from(&result)
}

#[wasm_bindgen]
pub fn difference_mesh(offset_x: f64) -> MeshData {
    let a = Manifold::cube(Vec3::new(1.0, 1.0, 1.0), false);
    let b = Manifold::cube(Vec3::new(1.0, 1.0, 1.0), false).translate(Vec3::new(offset_x, 0.0, 0.0));
    let result = a.difference(&b);
    mesh_data_from(&result)
}

// ---------------------------------------------------------------------------
// Menger Sponge — recursive boolean subtraction
// ---------------------------------------------------------------------------

fn menger_sponge_impl(depth: i32) -> Manifold {
    let mut result = Manifold::cube(Vec3::new(1.0, 1.0, 1.0), true);
    if depth == 0 {
        return result;
    }
    let mut holes = Vec::new();
    let size = 1.0 / 3.0;
    // Cross-shaped holes through all 3 axes
    for axis in 0..3 {
        for i in [-1.0, 0.0, 1.0] {
            for j in [-1.0, 0.0, 1.0] {
                if (i == 0.0) && (j == 0.0) { continue; }
                if (i != 0.0) && (j != 0.0) { continue; }
                let pos = match axis {
                    0 => Vec3::new(0.0, i * size, j * size),
                    1 => Vec3::new(i * size, 0.0, j * size),
                    _ => Vec3::new(i * size, j * size, 0.0),
                };
                let hole_size = match axis {
                    0 => Vec3::new(1.1, size, size),
                    1 => Vec3::new(size, 1.1, size),
                    _ => Vec3::new(size, size, 1.1),
                };
                holes.push(Manifold::cube(hole_size, true).translate(pos));
            }
        }
    }
    for hole in holes {
        result = result.difference(&hole);
    }
    if depth > 1 {
        // Recursive: scale down and place 20 sub-cubes (8 corners + 12 edges, skip face-centers + center)
        let child = menger_sponge_impl(depth - 1);
        let mut pieces = Vec::new();
        for x in [-1.0, 0.0, 1.0] {
            for y in [-1.0, 0.0, 1.0] {
                for z in [-1.0, 0.0, 1.0] {
                    let zeros = (x == 0.0) as i32 + (y == 0.0) as i32 + (z == 0.0) as i32;
                    if zeros >= 2 { continue; } // skip center and face-centers
                    pieces.push(child.clone().scale(Vec3::splat(size)).translate(Vec3::new(x * size, y * size, z * size)));
                }
            }
        }
        result = pieces.iter().skip(1).fold(pieces[0].clone(), |acc, p| acc.union(p));
    }
    result
}

#[wasm_bindgen]
pub fn menger_sponge_mesh(depth: i32) -> MeshData {
    let m = menger_sponge_impl(depth.clamp(0, 2));
    mesh_data_from(&m)
}

// ---------------------------------------------------------------------------
// Boolean Gallery — sphere-sphere, sphere-cube, cylinder-sphere
// ---------------------------------------------------------------------------

#[wasm_bindgen]
pub fn boolean_gallery_mesh(shape_a: i32, shape_b: i32, op: i32, offset_x: f64, offset_y: f64, offset_z: f64) -> MeshData {
    let a = match shape_a {
        0 => Manifold::cube(Vec3::new(1.0, 1.0, 1.0), true),
        1 => Manifold::sphere(0.6, 32),
        2 => Manifold::cylinder(1.0, 0.5, 0.5, 32).translate(Vec3::new(0.0, 0.0, -0.5)),
        _ => Manifold::cube(Vec3::new(1.0, 1.0, 1.0), true),
    };
    let b = match shape_b {
        0 => Manifold::cube(Vec3::new(1.0, 1.0, 1.0), true),
        1 => Manifold::sphere(0.6, 32),
        2 => Manifold::cylinder(1.0, 0.5, 0.5, 32).translate(Vec3::new(0.0, 0.0, -0.5)),
        _ => Manifold::cube(Vec3::new(1.0, 1.0, 1.0), true),
    };
    let b = b.translate(Vec3::new(offset_x, offset_y, offset_z));
    let result = match op {
        0 => a.union(&b),
        1 => a.intersection(&b),
        2 => a.difference(&b),
        _ => a.union(&b),
    };
    mesh_data_from(&result)
}

// ---------------------------------------------------------------------------
// Refined / Smooth shapes
// ---------------------------------------------------------------------------

#[wasm_bindgen]
pub fn refined_shape_mesh(shape: i32, refine_level: i32) -> MeshData {
    let base = match shape {
        0 => Manifold::tetrahedron(),
        1 => Manifold::cube(Vec3::new(1.0, 1.0, 1.0), true),
        2 => Manifold::cylinder(1.0, 0.5, 0.5, 8).translate(Vec3::new(0.0, 0.0, -0.5)),
        _ => Manifold::tetrahedron(),
    };
    let refined = base.refine(refine_level.clamp(1, 5));
    mesh_data_from(&refined)
}

// ---------------------------------------------------------------------------
// Extrude with twist
// ---------------------------------------------------------------------------

#[wasm_bindgen]
pub fn extrude_twist_mesh(radius: f64, segments: i32, height: f64, twist_degrees: f64, n_divisions: i32, scale_top: f64) -> MeshData {
    let cs = CrossSection::circle(radius, segments);
    let m = Manifold::extrude(&cs.to_polygons(), height, n_divisions, twist_degrees, Vec2::new(scale_top, scale_top));
    mesh_data_from(&m)
}

// ---------------------------------------------------------------------------
// Revolve with partial angle
// ---------------------------------------------------------------------------

#[wasm_bindgen]
pub fn revolve_partial_mesh(profile: i32, segments: i32, degrees: f64) -> MeshData {
    let cs = match profile {
        0 => CrossSection::circle(0.3, 16),
        1 => CrossSection::square(0.5),
        _ => CrossSection::circle(0.3, 16),
    };
    // Offset the profile from the axis
    let cs = cs.translate(Vec2::new(0.7, 0.0));
    let m = Manifold::revolve(&cs.to_polygons(), segments, degrees);
    mesh_data_from(&m)
}
