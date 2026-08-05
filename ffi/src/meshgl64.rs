// Copyright 2026 Lars Brubaker
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

//! Double-precision / 64-bit-index mesh path: the `MeshGL64` analogue of the
//! `MeshGL` import and export in `lib.rs`.
//!
//! Element types follow [`manifold_rust::types::MeshGL64`] =
//! `MeshGLP<f64, u64>` exactly: every `I`-typed field is `u64`
//! (`num_prop`, `tri_verts`, `run_index`, `face_id`) and every `P`-typed field
//! is `f64` (`vert_properties`), while `run_original_id` stays `Vec<u32>` in
//! both precisions. Getting this wrong would silently misread the array on the
//! C side, so the accessors below spell each type out rather than sharing the
//! f32 declarations.
//!
//! The kernel computes in f64 and the core's `from_mesh_gl64`/`get_mesh_gl64`
//! are lossless, so coordinates and tolerance survive this path at full double
//! precision — no narrowing anywhere. Indices are the exception: the kernel
//! indexes vertices with 32 bits (matching the C++ reference, whose import
//! casts every index to `uint32_t`), and that narrowing would *wrap*, which is
//! not a precision loss but a correctness one, so `manifold_rs_from_mesh64`
//! rejects out-of-range indices outright rather than passing them through.

use std::ptr;

use manifold_rust::manifold::Manifold;
use manifold_rust::types::MeshGL64;

use crate::error::{guard, set_last_error};
use crate::{into_handle, length_fits, mesh_array_accessor, slice_or_empty, ManifoldRs};

/// Opaque handle wrapping an exported [`MeshGL64`]. The array accessors hand
/// out pointers that borrow from this handle, so it must outlive their use.
pub struct MeshGl64Rs {
    inner: MeshGL64,
}

/// Double-precision counterpart of `manifold_rs_from_mesh`.
///
/// Same validation ladder and the same sentinels: a handle (possibly with a
/// non-zero status) for anything that describes a mesh, NULL for arguments that
/// cannot.
///
/// # Safety
/// `vert_properties` and `tri_verts` must point to at least the stated number
/// of elements, or be NULL when the corresponding length is 0.
#[no_mangle]
pub unsafe extern "C" fn manifold_rs_from_mesh64(
    vert_properties: *const f64,
    vert_properties_len: usize,
    tri_verts: *const u64,
    tri_verts_len: usize,
    num_prop: u32,
) -> *mut ManifoldRs {
    // SAFETY: same caller contract as from_mesh64_impl.
    unsafe {
        from_mesh64_impl(
            vert_properties,
            vert_properties_len,
            tri_verts,
            tri_verts_len,
            num_prop,
            false,
        )
    }
}

/// Double-precision counterpart of `manifold_rs_from_mesh_robust`: the
/// robust (non-manifold-tolerant) import. See that function for semantics.
///
/// # Safety
/// Same contract as [`manifold_rs_from_mesh64`].
#[no_mangle]
pub unsafe extern "C" fn manifold_rs_from_mesh64_robust(
    vert_properties: *const f64,
    vert_properties_len: usize,
    tri_verts: *const u64,
    tri_verts_len: usize,
    num_prop: u32,
) -> *mut ManifoldRs {
    // SAFETY: same caller contract as from_mesh64_impl.
    unsafe {
        from_mesh64_impl(
            vert_properties,
            vert_properties_len,
            tri_verts,
            tri_verts_len,
            num_prop,
            true,
        )
    }
}

unsafe fn from_mesh64_impl(
    vert_properties: *const f64,
    vert_properties_len: usize,
    tri_verts: *const u64,
    tri_verts_len: usize,
    num_prop: u32,
    robust: bool,
) -> *mut ManifoldRs {
    guard(ptr::null_mut(), || {
        if num_prop < 3 {
            set_last_error(format!("manifold_rs_from_mesh64: num_prop {num_prop} < 3"));
            return ptr::null_mut();
        }
        if vert_properties_len % num_prop as usize != 0 {
            set_last_error(format!(
                "manifold_rs_from_mesh64: vert_properties_len {vert_properties_len} is not a multiple of num_prop {num_prop}"
            ));
            return ptr::null_mut();
        }
        if tri_verts_len % 3 != 0 {
            set_last_error(format!(
                "manifold_rs_from_mesh64: tri_verts_len {tri_verts_len} is not a multiple of 3"
            ));
            return ptr::null_mut();
        }
        if (vert_properties.is_null() && vert_properties_len > 0)
            || (tri_verts.is_null() && tri_verts_len > 0)
        {
            set_last_error("manifold_rs_from_mesh64: null array with non-zero length");
            return ptr::null_mut();
        }
        if !length_fits::<f64>(vert_properties_len) || !length_fits::<u64>(tri_verts_len) {
            set_last_error(format!(
                "manifold_rs_from_mesh64: implausible length (vert_properties_len {vert_properties_len}, tri_verts_len {tri_verts_len})"
            ));
            return ptr::null_mut();
        }

        // SAFETY: both pointers are non-null for the non-empty case, and the
        // lengths are within the isize::MAX byte limit `from_raw_parts`
        // requires (both checked above); the caller guarantees the pointers
        // really cover that many elements.
        let (verts, tris) = unsafe {
            (
                slice_or_empty(vert_properties, vert_properties_len),
                slice_or_empty(tri_verts, tri_verts_len),
            )
        };

        // The core's `from_mesh_gl64` truncates every index to u32 (the
        // kernel's index width, matching the C++ import's uint32_t casts),
        // which *wraps* rather than saturating: an index of exactly 2^32
        // would silently become 0 and produce wrong geometry reporting
        // status 0 — the worst possible failure mode. Reject it here instead.
        // This is a real ceiling on the 64-bit path, not a formality, so it
        // gets its own rung on the ladder.
        if let Some((i, &bad)) = tris
            .iter()
            .enumerate()
            .find(|(_, &v)| v > u64::from(u32::MAX))
        {
            set_last_error(format!(
                "manifold_rs_from_mesh64: tri_verts[{i}] = {bad} exceeds u32::MAX; \
                 the kernel indexes vertices with 32 bits"
            ));
            return ptr::null_mut();
        }

        let mesh = MeshGL64 {
            // num_prop is the `I` type of MeshGLP, i.e. u64 here, but the C
            // signature keeps uint32_t: a per-vertex property count that needs
            // more than 32 bits is not a thing, and matching the f32 signature
            // keeps the two entry points interchangeable at the call site.
            num_prop: u64::from(num_prop),
            vert_properties: verts.to_vec(),
            tri_verts: tris.to_vec(),
            ..Default::default()
        };
        into_handle(if robust {
            Manifold::from_mesh_gl64_robust(&mesh)
        } else {
            Manifold::from_mesh_gl64(&mesh)
        })
    })
}

/// Export `m` as a double-precision mesh handle (normal index -1, matching
/// `manifold_rs_get_meshgl`). An error-status manifold exports as an empty mesh
/// rather than failing. NULL on a NULL handle or a panic.
///
/// # Safety
/// `m` must be NULL or a live handle from this library.
#[no_mangle]
pub unsafe extern "C" fn manifold_rs_get_meshgl64(m: *const ManifoldRs) -> *mut MeshGl64Rs {
    guard(ptr::null_mut(), || {
        // SAFETY: caller contract; as_ref() handles the NULL case.
        let Some(handle) = (unsafe { m.as_ref() }) else {
            set_last_error("manifold_rs_get_meshgl64: null manifold");
            return ptr::null_mut();
        };
        Box::into_raw(Box::new(MeshGl64Rs {
            inner: handle.inner.get_mesh_gl64(-1),
        }))
    })
}

/// Properties per vertex (>= 3; the first three are the position). 0 for a
/// NULL handle.
///
/// Narrowed from the field's `u64` to `uint32_t` for signature parity with
/// `manifold_rs_meshgl_num_prop`; the value is a small count that the import
/// side accepts as `uint32_t` in the first place.
///
/// # Safety
/// `g` must be NULL or a live handle from this library.
#[no_mangle]
pub unsafe extern "C" fn manifold_rs_meshgl64_num_prop(g: *const MeshGl64Rs) -> u32 {
    guard(0, || {
        // SAFETY: caller contract; as_ref() handles the NULL case.
        match unsafe { g.as_ref() } {
            Some(handle) => handle.inner.num_prop as u32,
            None => {
                set_last_error("manifold_rs_meshgl64_num_prop: null mesh");
                0
            }
        }
    })
}

mesh_array_accessor!(
    manifold_rs_meshgl64_vert_properties,
    MeshGl64Rs,
    vert_properties,
    f64
);
mesh_array_accessor!(manifold_rs_meshgl64_tri_verts, MeshGl64Rs, tri_verts, u64);
mesh_array_accessor!(manifold_rs_meshgl64_run_index, MeshGl64Rs, run_index, u64);
// u32 even in the 64-bit mesh: `run_original_id` is `Vec<u32>` for both
// instantiations of MeshGLP, because a mesh ID is not an index into the mesh.
mesh_array_accessor!(
    manifold_rs_meshgl64_run_original_id,
    MeshGl64Rs,
    run_original_id,
    u32
);
mesh_array_accessor!(manifold_rs_meshgl64_face_id, MeshGl64Rs, face_id, u64);

/// Release a mesh handle, invalidating every pointer previously returned by the
/// array accessors. NULL is ignored.
///
/// # Safety
/// `g` must be NULL or a handle from this library that has not already been
/// destroyed, and no other thread may be using it.
#[no_mangle]
pub unsafe extern "C" fn manifold_rs_meshgl64_destroy(g: *mut MeshGl64Rs) {
    guard((), || {
        if !g.is_null() {
            // SAFETY: caller guarantees single ownership of a live handle.
            drop(unsafe { Box::from_raw(g) });
        }
    })
}
