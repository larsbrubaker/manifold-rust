// WASM bindings for manifold-rust
// Exposes a compact browser-facing API for the demo site.

use manifold_rust::cross_section::CrossSection;
use manifold_rust::linalg::{Vec2, Vec3};
use manifold_rust::manifold::Manifold;
use wasm_bindgen::prelude::*;

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

fn summarize(mesh: Manifold) -> MeshSummary {
    MeshSummary {
        num_vert: mesh.num_vert() as u32,
        num_tri: mesh.num_tri() as u32,
        volume: mesh.volume(),
        surface_area: mesh.surface_area(),
    }
}

#[wasm_bindgen]
pub fn version() -> String {
    env!("CARGO_PKG_VERSION").to_string()
}

#[wasm_bindgen]
pub fn cube(size_x: f64, size_y: f64, size_z: f64, center: bool) -> MeshSummary {
    summarize(Manifold::cube(Vec3::new(size_x, size_y, size_z), center))
}

#[wasm_bindgen]
pub fn sphere(radius: f64, circular_segments: i32) -> MeshSummary {
    summarize(Manifold::sphere(radius, circular_segments))
}

#[wasm_bindgen]
pub fn cylinder(height: f64, radius_low: f64, radius_high: f64, circular_segments: i32) -> MeshSummary {
    summarize(Manifold::cylinder(height, radius_low, radius_high, circular_segments))
}

#[wasm_bindgen]
pub fn union_cubes(offset_x: f64) -> MeshSummary {
    let a = Manifold::cube(Vec3::new(1.0, 1.0, 1.0), false);
    let b = Manifold::cube(Vec3::new(1.0, 1.0, 1.0), false).translate(Vec3::new(offset_x, 0.0, 0.0));
    summarize(a.union(&b))
}

#[wasm_bindgen]
pub fn intersect_cubes(offset_x: f64) -> MeshSummary {
    let a = Manifold::cube(Vec3::new(1.0, 1.0, 1.0), false);
    let b = Manifold::cube(Vec3::new(1.0, 1.0, 1.0), false).translate(Vec3::new(offset_x, 0.0, 0.0));
    summarize(a.intersection(&b))
}

#[wasm_bindgen]
pub fn difference_cubes(offset_x: f64) -> MeshSummary {
    let a = Manifold::cube(Vec3::new(1.0, 1.0, 1.0), false);
    let b = Manifold::cube(Vec3::new(1.0, 1.0, 1.0), false).translate(Vec3::new(offset_x, 0.0, 0.0));
    summarize(a.difference(&b))
}

#[wasm_bindgen]
pub fn extrude_circle(radius: f64, segments: i32, height: f64) -> MeshSummary {
    let cs = CrossSection::circle(radius, segments);
    summarize(Manifold::extrude(&cs.to_polygons(), height, 0, 0.0, Vec2::new(1.0, 1.0)))
}
