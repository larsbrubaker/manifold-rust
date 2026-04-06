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
