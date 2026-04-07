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
    colors: Option<Vec<f32>>,
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
    pub fn has_colors(&self) -> bool {
        self.colors.is_some()
    }

    #[wasm_bindgen(getter)]
    pub fn colors(&self) -> Option<Float32Array> {
        self.colors.as_ref().map(|c| Float32Array::from(&c[..]))
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

    // Extract per-vertex colors from properties [3..7] if present (RGBA)
    let colors = if num_prop >= 7 {
        let mut cols = Vec::with_capacity(vert_count * 4);
        for i in 0..vert_count {
            let base = i * num_prop;
            cols.push(gl.vert_properties[base + 3]); // R
            cols.push(gl.vert_properties[base + 4]); // G
            cols.push(gl.vert_properties[base + 5]); // B
            cols.push(gl.vert_properties[base + 6]); // A
        }
        Some(cols)
    } else if num_prop >= 6 {
        // RGB only, add alpha=1.0
        let mut cols = Vec::with_capacity(vert_count * 4);
        for i in 0..vert_count {
            let base = i * num_prop;
            cols.push(gl.vert_properties[base + 3]);
            cols.push(gl.vert_properties[base + 4]);
            cols.push(gl.vert_properties[base + 5]);
            cols.push(1.0);
        }
        Some(cols)
    } else {
        None
    };

    MeshData {
        positions,
        normals,
        indices: gl.tri_verts,
        colors,
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
pub fn spiky_dodecahedron_mesh(spike_height: f64) -> MeshData {
    let m = make_spiky_dodecahedron(spike_height);
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

fn make_spiky_dodecahedron(spike_height: f64) -> Manifold {
    use manifold_rust::types::MeshGL;

    // Golden ratio
    let phi: f64 = (1.0 + 5.0_f64.sqrt()) / 2.0;
    let inv_phi = 1.0 / phi;

    // 20 vertices of a regular dodecahedron (scaled to ~unit size)
    let scale = 0.5;
    let raw_verts: [(f64, f64, f64); 20] = [
        // Cube vertices (±1, ±1, ±1)
        ( 1.0,  1.0,  1.0), ( 1.0,  1.0, -1.0), ( 1.0, -1.0,  1.0), ( 1.0, -1.0, -1.0),
        (-1.0,  1.0,  1.0), (-1.0,  1.0, -1.0), (-1.0, -1.0,  1.0), (-1.0, -1.0, -1.0),
        // (0, ±1/φ, ±φ)
        (0.0,  inv_phi,  phi), (0.0,  inv_phi, -phi), (0.0, -inv_phi,  phi), (0.0, -inv_phi, -phi),
        // (±1/φ, ±φ, 0)
        ( inv_phi,  phi, 0.0), (-inv_phi,  phi, 0.0), ( inv_phi, -phi, 0.0), (-inv_phi, -phi, 0.0),
        // (±φ, 0, ±1/φ)
        ( phi, 0.0,  inv_phi), ( phi, 0.0, -inv_phi), (-phi, 0.0,  inv_phi), (-phi, 0.0, -inv_phi),
    ];

    // 12 pentagonal faces (vertex indices, counterclockwise when viewed from outside)
    let faces: [[usize; 5]; 12] = [
        [0, 8, 10, 2, 16],   [0, 16, 17, 1, 12],  [0, 12, 13, 4, 8],
        [1, 17, 3, 11, 9],   [1, 9, 5, 13, 12],   [2, 10, 6, 15, 14],
        [2, 14, 3, 17, 16],  [4, 13, 5, 19, 18],   [4, 18, 6, 10, 8],
        [5, 9, 11, 7, 19],   [6, 18, 19, 7, 15],   [3, 14, 15, 7, 11],
    ];

    let verts: Vec<(f64, f64, f64)> = raw_verts.iter().map(|&(x, y, z)| (x * scale, y * scale, z * scale)).collect();

    // Build mesh: 20 original verts + 12 spike verts = 32 verts, 60 triangles
    let mut positions: Vec<f32> = Vec::with_capacity(32 * 3);
    let mut tri_verts: Vec<u32> = Vec::with_capacity(60 * 3);

    // Add original vertices
    for &(x, y, z) in &verts {
        positions.extend([x as f32, y as f32, z as f32]);
    }

    // For each face, compute center, push spike vertex, create 5 triangles
    for face in &faces {
        let cx: f64 = face.iter().map(|&i| verts[i].0).sum::<f64>() / 5.0;
        let cy: f64 = face.iter().map(|&i| verts[i].1).sum::<f64>() / 5.0;
        let cz: f64 = face.iter().map(|&i| verts[i].2).sum::<f64>() / 5.0;
        let len = (cx * cx + cy * cy + cz * cz).sqrt();
        // Spike outward along face normal (center direction from origin)
        let nx = cx / len;
        let ny = cy / len;
        let nz = cz / len;
        let spike_idx = (positions.len() / 3) as u32;
        positions.extend([
            (cx + nx * spike_height) as f32,
            (cy + ny * spike_height) as f32,
            (cz + nz * spike_height) as f32,
        ]);
        // 5 triangles fanning from spike to face edges
        for j in 0..5 {
            let a = face[j] as u32;
            let b = face[(j + 1) % 5] as u32;
            tri_verts.extend([spike_idx, a, b]);
        }
    }

    let mut mesh = MeshGL::default();
    mesh.num_prop = 3;
    mesh.vert_properties = positions;
    mesh.tri_verts = tri_verts;
    Manifold::from_mesh_gl(&mesh)
}

fn make_shape(id: i32) -> Manifold {
    match id {
        0 => Manifold::cube(Vec3::new(1.0, 1.0, 1.0), true),
        1 => Manifold::sphere(0.6, 32),
        2 => Manifold::cylinder(1.0, 0.5, 0.5, 32).translate(Vec3::new(0.0, 0.0, -0.5)),
        3 => make_spiky_dodecahedron(0.4),
        _ => Manifold::cube(Vec3::new(1.0, 1.0, 1.0), true),
    }
}

fn color_shape(m: &Manifold, r: f64, g: f64, b: f64, a: f64) -> Manifold {
    // Internal num_prop is extra properties beyond xyz.
    // 4 extra = RGBA color channels, stored at property indices 0,1,2,3.
    // get_mesh_gl prepends xyz, so output has 7 props: xyz + rgba.
    m.set_properties(4, move |new_prop, _pos, _old| {
        new_prop[0] = r;
        new_prop[1] = g;
        new_prop[2] = b;
        new_prop[3] = a;
    })
}

#[wasm_bindgen]
pub fn boolean_gallery_mesh(shape_a: i32, shape_b: i32, op: i32, offset_x: f64, offset_y: f64, offset_z: f64) -> MeshData {
    boolean_gallery_mesh_rotated(shape_a, shape_b, op, offset_x, offset_y, offset_z, 0.0, 0.0, 0.0)
}

#[wasm_bindgen]
pub fn boolean_gallery_mesh_rotated(shape_a: i32, shape_b: i32, op: i32, offset_x: f64, offset_y: f64, offset_z: f64, rot_x: f64, rot_y: f64, rot_z: f64) -> MeshData {
    let a = make_shape(shape_a);
    let b = make_shape(shape_b);
    // Assign distinct colors: shape A = blue (opaque), shape B = off-red (translucent)
    let a = color_shape(&a, 0.27, 0.53, 0.80, 1.0);   // #4488cc blue, fully opaque
    let b = color_shape(&b, 0.85, 0.25, 0.25, 0.6);   // off-red, alpha 0.6
    // Rotate shape B about its offset center, then translate
    let b = b.rotate(rot_x, rot_y, rot_z).translate(Vec3::new(offset_x, offset_y, offset_z));
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

// ---------------------------------------------------------------------------
// Test Gallery — visualizations of key test cases
// ---------------------------------------------------------------------------

/// Mirror Union: cube A at origin + cube B translated, then mirrored across (1,1,0) plane
#[wasm_bindgen]
pub fn test_mirror_union_mesh() -> MeshData {
    let a = Manifold::cube(Vec3::new(5.0, 5.0, 5.0), true);
    let b = Manifold::cube(Vec3::new(5.0, 5.0, 5.0), true)
        .translate(Vec3::new(2.5, 2.5, 2.5));
    let b_mirrored = b.mirror(Vec3::new(1.0, 1.0, 0.0));
    let result = a.union(&b).union(&b_mirrored);
    mesh_data_from(&result)
}

/// Split by Plane: cube rotated and split at z=1
#[wasm_bindgen]
pub fn test_split_by_plane_mesh(half: i32) -> MeshData {
    let cube = Manifold::cube(Vec3::new(2.0, 2.0, 2.0), true)
        .translate(Vec3::new(0.0, 1.0, 0.0))
        .rotate(90.0, 0.0, 0.0);
    let (top, bottom) = cube.split_by_plane(Vec3::new(0.0, 0.0, 1.0), 1.0);
    if half == 0 { mesh_data_from(&top) } else { mesh_data_from(&bottom) }
}

/// Vug (cavity): outer cube with inner cube subtracted
#[wasm_bindgen]
pub fn test_vug_mesh() -> MeshData {
    let outer = Manifold::cube(Vec3::new(4.0, 4.0, 4.0), true);
    let inner = Manifold::cube(Vec3::new(1.0, 1.0, 1.0), false);
    let vug = outer.difference(&inner);
    // Split to show the cavity inside
    let (top, _) = vug.split_by_plane(Vec3::new(0.0, 0.0, 1.0), -0.5);
    mesh_data_from(&top)
}

/// Warp: extruded square with parabolic displacement x += z^2
#[wasm_bindgen]
pub fn test_warp_mesh() -> MeshData {
    let square = CrossSection::square(1.0);
    let m = Manifold::extrude(&square.to_polygons(), 2.0, 10, 0.0, Vec2::new(1.0, 1.0))
        .warp(|v| { v.x += v.z * v.z; });
    mesh_data_from(&m)
}

/// Spiral: recursive spiral of cubes
#[wasm_bindgen]
pub fn test_spiral_mesh() -> MeshData {
    fn spiral(rec: i32, r: f64, add: f64, d: f64) -> Manifold {
        let rot = 360.0 / (std::f64::consts::PI * r * 2.0) * d;
        let r_next = r + add / 360.0 * rot;
        let cube = Manifold::cube(Vec3::splat(1.0), true)
            .translate(Vec3::new(0.0, r, 0.0));
        if rec > 0 {
            spiral(rec - 1, r_next, add, d).rotate(0.0, 0.0, rot).union(&cube)
        } else {
            cube
        }
    }
    let result = spiral(10, 25.0, 2.0, 2.0)
        .scale(Vec3::splat(0.1)); // Scale down from ~50 units to ~5 for viewer
    mesh_data_from(&result)
}

/// Sphere Difference: sphere with spherical cap removed
#[wasm_bindgen]
pub fn test_sphere_diff_mesh() -> MeshData {
    let s1 = Manifold::sphere(1.0, 32);
    let s2 = Manifold::sphere(1.0, 32).translate(Vec3::new(0.5, 0.5, 0.5));
    let result = s1.difference(&s2);
    mesh_data_from(&result)
}

/// Overlapping Cubes: three cubes at different offsets, all unioned
#[wasm_bindgen]
pub fn test_cubes_union_mesh() -> MeshData {
    let a = Manifold::cube(Vec3::new(1.0, 1.0, 1.0), true);
    let b = Manifold::cube(Vec3::new(1.0, 1.0, 1.0), true)
        .translate(Vec3::new(0.3, 0.3, 0.0));
    let c = Manifold::cube(Vec3::new(1.0, 1.0, 1.0), true)
        .translate(Vec3::new(-0.3, 0.3, 0.3));
    let result = a.union(&b).union(&c);
    mesh_data_from(&result)
}

/// Batch Boolean: flat slab with three cylinder holes punched through
#[wasm_bindgen]
pub fn test_batch_subtract_mesh() -> MeshData {
    let slab = Manifold::cube(Vec3::new(10.0, 10.0, 1.0), true);
    let c1 = Manifold::cylinder(2.0, 0.5, 0.5, 32)
        .translate(Vec3::new(-3.0, -3.0, 0.0));
    let c2 = Manifold::cylinder(2.0, 0.5, 0.5, 32)
        .translate(Vec3::new(3.0, 3.0, 0.0));
    let c3 = Manifold::cylinder(2.0, 0.5, 0.5, 32)
        .translate(Vec3::new(0.0, 0.0, 0.0));
    let result = slab.difference(&c1).difference(&c2).difference(&c3);
    mesh_data_from(&result)
}
