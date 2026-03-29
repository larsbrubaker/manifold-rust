// WASM bindings for manifold-rust
// Exposes Manifold operations to JavaScript/TypeScript for the demo site.

use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub fn version() -> String {
    env!("CARGO_PKG_VERSION").to_string()
}
