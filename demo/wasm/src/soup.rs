// Triangle-soup import and engine-selectable booleans for the demo site.
//
// Backs the Boolean Gallery's Thingi10K mode: meshes arrive from JS as raw
// positions/indices parsed out of STL files, get their coincident vertices
// welded via `MeshGL::merge`, and import through `from_mesh_gl_robust` so
// closed-but-non-manifold geometry survives as triangle soup for the robust
// boolean engine. Manifold input imports exactly as `from_mesh_gl` would.
//
// See `lib.rs` for the shared `MeshData` geometry bridge to Three.js.

use manifold_rust::linalg::Vec3;
use manifold_rust::manifold::Manifold;
use manifold_rust::types::{BooleanEngine, Error, MeshGL};
use wasm_bindgen::prelude::*;

use crate::{mesh_data_from, MeshData};

fn engine_from_i32(engine: i32) -> BooleanEngine {
    match engine {
        1 => BooleanEngine::Robust,
        2 => BooleanEngine::Auto,
        _ => BooleanEngine::Exact,
    }
}

pub(crate) fn op_with_engine(a: &Manifold, b: &Manifold, op: i32, engine: i32) -> Manifold {
    let engine = engine_from_i32(engine);
    match op {
        1 => a.intersection_with_engine(b, engine),
        2 => a.difference_with_engine(b, engine),
        _ => a.union_with_engine(b, engine),
    }
}

/// A mesh imported from arbitrary (possibly non-manifold) triangle soup,
/// held on the WASM side so repeated booleans don't re-import per call.
#[wasm_bindgen]
pub struct ImportedMesh {
    manifold: Manifold,
}

#[wasm_bindgen]
impl ImportedMesh {
    /// Import raw triangle soup: `positions` is xyz-interleaved, `indices`
    /// is 3 per triangle. Coincident vertices are welded first (STL stores
    /// each triangle's corners independently), then the robust import keeps
    /// whatever connectivity results — manifold or soup — as long as the
    /// geometry is closed and orientable.
    #[wasm_bindgen(constructor)]
    pub fn new(positions: &[f32], indices: &[u32]) -> ImportedMesh {
        let mut mesh = MeshGL::default();
        mesh.num_prop = 3;
        mesh.vert_properties = positions.to_vec();
        mesh.tri_verts = indices.to_vec();
        mesh.merge();
        ImportedMesh {
            manifold: Manifold::from_mesh_gl_robust(&mesh),
        }
    }

    /// Import status: "No Error" on success, "Not Closed" when even the
    /// soup path rejected the geometry (result is empty then).
    #[wasm_bindgen(getter)]
    pub fn status(&self) -> String {
        self.manifold.status().to_str().to_string()
    }

    #[wasm_bindgen(getter)]
    pub fn ok(&self) -> bool {
        self.manifold.status() == Error::NoError && !self.manifold.is_empty()
    }

    /// True when the import fell back to triangle soup — i.e. the mesh is
    /// non-manifold and booleans will need the Robust (or Auto) engine.
    #[wasm_bindgen(getter)]
    pub fn is_soup(&self) -> bool {
        self.manifold.as_impl().is_soup
    }

    #[wasm_bindgen(getter)]
    pub fn num_vert(&self) -> u32 {
        self.manifold.num_vert() as u32
    }

    #[wasm_bindgen(getter)]
    pub fn num_tri(&self) -> u32 {
        self.manifold.num_tri() as u32
    }

    /// Geometry for previewing the imported mesh on its own.
    pub fn mesh(&self) -> MeshData {
        mesh_data_from(&self.manifold)
    }
}

/// Boolean between two imported meshes. `b` is rotated (degrees, XYZ order)
/// about its own origin, then translated. op: 0=union, 1=intersection,
/// 2=difference. engine: 0=Exact, 1=Robust, 2=Auto. Throws the result's
/// status string when the operation yields an error (e.g. "Not Manifold"
/// when the Exact engine is handed soup geometry).
#[wasm_bindgen]
#[allow(clippy::too_many_arguments)]
pub fn imported_boolean(
    a: &ImportedMesh,
    b: &ImportedMesh,
    op: i32,
    engine: i32,
    off_x: f64,
    off_y: f64,
    off_z: f64,
    rot_x: f64,
    rot_y: f64,
    rot_z: f64,
) -> Result<MeshData, JsValue> {
    // Color the operands (A opaque blue, B translucent red). Both engines
    // carry properties through; only soup geometry cannot take them
    // (set_properties requires paired halfedges).
    let colorize = !a.is_soup() && !b.is_soup();
    let (a_m, b_m) = if colorize {
        (
            crate::color_shape(&a.manifold, 0.27, 0.53, 0.80, 1.0),
            crate::color_shape(&b.manifold, 0.85, 0.25, 0.25, 0.6),
        )
    } else {
        (a.manifold.clone(), b.manifold.clone())
    };
    let b_m = b_m
        .rotate(rot_x, rot_y, rot_z)
        .translate(Vec3::new(off_x, off_y, off_z));
    let result = op_with_engine(&a_m, &b_m, op, engine);
    if result.status() != Error::NoError {
        return Err(JsValue::from_str(result.status().to_str()));
    }
    Ok(mesh_data_from(&result))
}
