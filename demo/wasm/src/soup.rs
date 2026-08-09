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
use manifold_rust::progress::ProgressReporter;
use manifold_rust::types::{BooleanEngine, Error, MeshGL, OpType};
use wasm_bindgen::prelude::*;

use crate::{mesh_data_from, MeshData};

fn engine_from_i32(engine: i32) -> BooleanEngine {
    match engine {
        1 => BooleanEngine::Robust,
        2 => BooleanEngine::Auto,
        _ => BooleanEngine::Exact,
    }
}

fn op_from_i32(op: i32) -> OpType {
    match op {
        1 => OpType::Intersect,
        2 => OpType::Subtract,
        _ => OpType::Add,
    }
}

/// Boolean between two operands with an optional progress sink. `None` takes
/// the un-instrumented path in the kernel, so the no-callback demo runs are
/// unaffected.
pub(crate) fn op_with_engine_progress(
    a: &Manifold,
    b: &Manifold,
    op: i32,
    engine: i32,
    nonzero: bool,
    progress: Option<&ProgressReporter>,
) -> Manifold {
    a.boolean_with_engine_rule_and_progress(
        b,
        op_from_i32(op),
        engine_from_i32(engine),
        winding_rule(nonzero),
        None,
        progress,
    )
}

/// `nonzero` selects the {w != 0} solid rule, which keeps inside-out
/// geometry as material instead of discarding it; the default {w >= 1} rule
/// is what every other entry point uses.
fn winding_rule(nonzero: bool) -> manifold_rust::types::WindingRule {
    use manifold_rust::types::WindingRule;
    if nonzero {
        WindingRule::Nonzero
    } else {
        WindingRule::Positive
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

    /// True when two of the mesh's own triangles genuinely intersect — they
    /// cross, they overlap, or they are coincident surface — rather than
    /// merely sharing edges and vertices. Manifold connectivity does not
    /// rule this out, and it is why the Auto engine may pick Robust for an
    /// operand that is not soup. The scan runs on demand (the importer does
    /// not pay for it) and its verdict is cached.
    #[wasm_bindgen(getter)]
    pub fn self_intersecting(&self) -> bool {
        self.manifold.has_self_intersections()
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
/// 2=difference. engine: 0=Exact, 1=Robust, 2=Auto. When `repair` is set,
/// `Manifold::repair_orientation` rewinds inside-out shells of both operands
/// before the boolean, so inverted bodies contribute material instead of
/// vanishing under the {winding >= 1} semantics. When `nonzero` is set the
/// robust engine switches to the {winding != 0} solid rule instead, which
/// keeps inside-out geometry even where per-shell repair cannot help (a
/// single shell wound correctly in one region and inverted in another). Throws the result's
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
    repair: bool,
    nonzero: bool,
) -> Result<MeshData, JsValue> {
    imported_boolean_impl(
        a, b, op, engine, off_x, off_y, off_z, rot_x, rot_y, rot_z, repair, nonzero, None,
    )
}

/// [`imported_boolean`] that reports pipeline progress to `on_progress`,
/// called as `on_progress(phaseName, fraction | null)`. Pass `undefined` for
/// the un-instrumented path.
///
/// The call is synchronous and blocking, so this only helps a caller that can
/// observe side effects while it runs — which is exactly the demo's worker,
/// whose callback does `postMessage` (queued for the main thread while the
/// worker stays busy).
#[wasm_bindgen]
#[allow(clippy::too_many_arguments)]
pub fn imported_boolean_progress(
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
    repair: bool,
    nonzero: bool,
    on_progress: Option<js_sys::Function>,
) -> Result<MeshData, JsValue> {
    crate::progress::with_reporter(on_progress, |reporter| {
        imported_boolean_impl(
            a, b, op, engine, off_x, off_y, off_z, rot_x, rot_y, rot_z, repair, nonzero, reporter,
        )
    })
}

#[allow(clippy::too_many_arguments)]
fn imported_boolean_impl(
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
    repair: bool,
    nonzero: bool,
    progress: Option<&ProgressReporter>,
) -> Result<MeshData, JsValue> {
    let (a_in, b_in) = if repair {
        (
            a.manifold.repair_orientation(),
            b.manifold.repair_orientation(),
        )
    } else {
        (a.manifold.clone(), b.manifold.clone())
    };
    // Color the operands (A opaque blue, B translucent red). Both engines
    // carry properties through; only soup geometry cannot take them
    // (set_properties requires paired halfedges).
    let colorize = !a.is_soup() && !b.is_soup();
    let (a_m, b_m) = if colorize {
        (
            crate::color_shape(&a_in, 0.27, 0.53, 0.80, 1.0),
            crate::color_shape(&b_in, 0.85, 0.25, 0.25, 0.6),
        )
    } else {
        (a_in, b_in)
    };
    let b_m = b_m
        .rotate(rot_x, rot_y, rot_z)
        .translate(Vec3::new(off_x, off_y, off_z));
    let result = op_with_engine_progress(&a_m, &b_m, op, engine, nonzero, progress);
    if result.status() != Error::NoError {
        return Err(JsValue::from_str(result.status().to_str()));
    }
    Ok(mesh_data_from(&result))
}
