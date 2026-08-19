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

//! C ABI for manifold-rust, built as `manifold_rs.dll` / `libmanifold_rs.so`.
//!
//! The surface is deliberately small: import a triangle soup, run an n-ary
//! boolean, export the result. It mirrors the subset of the C++ `manifoldc`
//! API that MatterCAD actually calls, with the same integer op ordering, so a
//! P/Invoke consumer can switch libraries without changing call sites.
//!
//! Three rules hold for every exported function:
//! - it never unwinds (see [`error::guard`]);
//! - it null-checks every pointer argument and returns a documented sentinel;
//! - handles are owned by the caller and freed with the matching destroy call.
//!
//! `manifold_rs.h` is the normative description of the ABI.

// Every raw-pointer dereference here is a caller-supplied pointer, so each one
// gets its own `unsafe` block with a SAFETY note stating which part of the C
// contract makes it sound. Denying the lint keeps that discipline from eroding
// as functions are added.
#![deny(unsafe_op_in_unsafe_fn)]

mod cancel;
mod error;
mod meshgl64;
mod progress;
#[cfg(test)]
mod tests;
#[cfg(test)]
mod tests_cancel;
#[cfg(test)]
mod tests_mesh64;
#[cfg(test)]
mod tests_progress;
#[cfg(test)]
mod tests_robust;

use std::os::raw::c_char;
use std::ptr;

use manifold_rust::cancel::CancelToken;
use manifold_rust::csg_tree::CsgNode;
use manifold_rust::manifold::Manifold;
use manifold_rust::types::{Error, MeshGL, OpType, WindingRule};

use crate::cancel::CancelTokenRs;
use crate::error::{guard, last_error_message, set_last_error};

/// Opaque handle wrapping a [`Manifold`]. Created by the constructor and
/// operation functions, released by [`manifold_rs_destroy`].
pub struct ManifoldRs {
    inner: Manifold,
}

/// Opaque handle wrapping an exported [`MeshGL`]. The array accessors hand out
/// pointers that borrow from this handle, so it must outlive their use.
pub struct MeshGlRs {
    inner: MeshGL,
}

pub(crate) fn into_handle(manifold: Manifold) -> *mut ManifoldRs {
    Box::into_raw(Box::new(ManifoldRs { inner: manifold }))
}

/// Stable numeric codes for [`Error`], matching the declaration order of the
/// enum in the core crate (which in turn matches C++ `ManifoldError`).
pub(crate) fn error_code(status: Error) -> i32 {
    match status {
        Error::NoError => 0,
        Error::NonFiniteVertex => 1,
        Error::NotManifold => 2,
        Error::VertexOutOfBounds => 3,
        Error::PropertiesWrongLength => 4,
        Error::MissingPositionProperties => 5,
        Error::MergeVectorsDifferentLengths => 6,
        Error::MergeIndexOutOfBounds => 7,
        Error::TransformWrongLength => 8,
        Error::RunIndexWrongLength => 9,
        Error::FaceIdWrongLength => 10,
        Error::InvalidConstruction => 11,
        Error::ResultTooLarge => 12,
        Error::InvalidTangents => 13,
        Error::Cancelled => 14,
        Error::NotClosed => 15,
    }
}

/// NUL terminator is part of the literal so the pointer can be handed to C
/// directly without allocating a CString at runtime.
const VERSION: &str = concat!(
    "manifold-ffi ",
    env!("CARGO_PKG_VERSION"),
    " (manifold-rust ",
    env!("MANIFOLD_RUST_VERSION"),
    ")\0"
);

/// Static NUL-terminated version string. Never NULL; must not be freed.
#[no_mangle]
pub extern "C" fn manifold_rs_version() -> *const c_char {
    VERSION.as_ptr() as *const c_char
}

// ---------------------------------------------------------------------------
// Manifold construction and queries
// ---------------------------------------------------------------------------

/// Build a manifold from interleaved vertex properties and triangle indices.
///
/// Returns a handle even when the mesh fails validation — the handle then
/// carries a non-zero [`manifold_rs_status`], which is how the caller learns
/// *why* (matching the C++ behaviour of returning an error-status manifold).
/// NULL is returned only for arguments that cannot describe a mesh at all, or
/// if the import panics.
///
/// # Safety
/// `vert_properties` and `tri_verts` must point to at least the stated number
/// of elements, or be NULL when the corresponding length is 0.
#[no_mangle]
pub unsafe extern "C" fn manifold_rs_from_mesh(
    vert_properties: *const f32,
    vert_properties_len: usize,
    tri_verts: *const u32,
    tri_verts_len: usize,
    num_prop: u32,
) -> *mut ManifoldRs {
    // SAFETY: same caller contract as from_mesh_impl32.
    unsafe {
        from_mesh_impl32(
            vert_properties,
            vert_properties_len,
            tri_verts,
            tri_verts_len,
            num_prop,
            false,
        )
    }
}

/// [`manifold_rs_from_mesh`] via the robust (non-manifold-tolerant) import:
/// manifold input behaves exactly like `manifold_rs_from_mesh`; closed
/// orientable but non-manifold input is retained for the robust boolean
/// engine (status 0) instead of being rejected; geometry that is not even
/// closed reports status 15 (NotClosed).
///
/// # Safety
/// Same contract as [`manifold_rs_from_mesh`].
#[no_mangle]
pub unsafe extern "C" fn manifold_rs_from_mesh_robust(
    vert_properties: *const f32,
    vert_properties_len: usize,
    tri_verts: *const u32,
    tri_verts_len: usize,
    num_prop: u32,
) -> *mut ManifoldRs {
    // SAFETY: same caller contract as from_mesh_impl32.
    unsafe {
        from_mesh_impl32(
            vert_properties,
            vert_properties_len,
            tri_verts,
            tri_verts_len,
            num_prop,
            true,
        )
    }
}

unsafe fn from_mesh_impl32(
    vert_properties: *const f32,
    vert_properties_len: usize,
    tri_verts: *const u32,
    tri_verts_len: usize,
    num_prop: u32,
    robust: bool,
) -> *mut ManifoldRs {
    guard(ptr::null_mut(), || {
        if num_prop < 3 {
            set_last_error(format!("manifold_rs_from_mesh: num_prop {num_prop} < 3"));
            return ptr::null_mut();
        }
        if vert_properties_len % num_prop as usize != 0 {
            set_last_error(format!(
                "manifold_rs_from_mesh: vert_properties_len {vert_properties_len} is not a multiple of num_prop {num_prop}"
            ));
            return ptr::null_mut();
        }
        if tri_verts_len % 3 != 0 {
            set_last_error(format!(
                "manifold_rs_from_mesh: tri_verts_len {tri_verts_len} is not a multiple of 3"
            ));
            return ptr::null_mut();
        }
        if (vert_properties.is_null() && vert_properties_len > 0)
            || (tri_verts.is_null() && tri_verts_len > 0)
        {
            set_last_error("manifold_rs_from_mesh: null array with non-zero length");
            return ptr::null_mut();
        }
        if !length_fits::<f32>(vert_properties_len) || !length_fits::<u32>(tri_verts_len) {
            set_last_error(format!(
                "manifold_rs_from_mesh: implausible length (vert_properties_len {vert_properties_len}, tri_verts_len {tri_verts_len})"
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

        let mesh = MeshGL {
            num_prop,
            vert_properties: verts.to_vec(),
            tri_verts: tris.to_vec(),
            ..Default::default()
        };
        into_handle(if robust {
            Manifold::from_mesh_gl_robust(&mesh)
        } else {
            Manifold::from_mesh_gl(&mesh)
        })
    })
}

// ---------------------------------------------------------------------------
// Boolean engine selection
// ---------------------------------------------------------------------------

/// Set the process-global default boolean engine: 0 = Exact (default),
/// 1 = Robust, 2 = Auto (Exact unless an operand carries non-manifold soup
/// geometry from `manifold_rs_from_mesh_robust`). Applies to every boolean
/// entry point, including `manifold_rs_batch_boolean`. Returns 0 on success,
/// -1 for an unknown engine value.
#[no_mangle]
pub extern "C" fn manifold_rs_set_boolean_engine(engine: i32) -> i32 {
    use manifold_rust::types::{BooleanConfig, BooleanEngine};
    let engine = match engine {
        0 => BooleanEngine::Exact,
        1 => BooleanEngine::Robust,
        2 => BooleanEngine::Auto,
        other => {
            set_last_error(format!(
                "manifold_rs_set_boolean_engine: unknown engine {other}"
            ));
            return -1;
        }
    };
    BooleanConfig::set_default_engine(engine);
    0
}

/// The current process-global default boolean engine (see
/// [`manifold_rs_set_boolean_engine`] for the value mapping).
#[no_mangle]
pub extern "C" fn manifold_rs_get_boolean_engine() -> i32 {
    use manifold_rust::types::{BooleanConfig, BooleanEngine};
    match BooleanConfig::default_engine() {
        BooleanEngine::Exact => 0,
        BooleanEngine::Robust => 1,
        BooleanEngine::Auto => 2,
    }
}

/// Copy of `m` with inside-out shells rewound so every body reads as solid
/// material under the robust engine's {winding >= 1} semantics: outermost
/// shells end up winding +1 and legitimate cavity shells stay inward-wound.
/// Coincident/doubled sheets are left untouched (the robust boolean's
/// winding-stack arithmetic already classifies those). Works on both strict
/// and robust-imported (soup) meshes, independent of any boolean call.
/// Returns NULL only if `m` is NULL or the repair panics; a mesh that needs
/// no repair comes back as a plain copy.
///
/// # Safety
/// `m` must be NULL or a live handle from this library.
#[no_mangle]
pub unsafe extern "C" fn manifold_rs_repair_orientation(m: *const ManifoldRs) -> *mut ManifoldRs {
    guard(ptr::null_mut(), || {
        // SAFETY: caller contract; as_ref() handles the NULL case.
        let Some(handle) = (unsafe { m.as_ref() }) else {
            set_last_error("manifold_rs_repair_orientation: null manifold");
            return ptr::null_mut();
        };
        into_handle(handle.inner.repair_orientation())
    })
}

/// Rebuild `m` into a fresh, properly paired 2-manifold enclosing the same
/// solid region under `winding_rule` (`MANIFOLD_RS_WINDING_*`: 0 = positive
/// `{w >= 1}`, 1 = nonzero `{w != 0}`).
///
/// Where [`manifold_rs_repair_orientation`] only rewinds triangles, this runs
/// the full robust pipeline on the one mesh: self-intersections are cut,
/// doubled and coincident sheets collapse, interior walls with material on
/// both sides dissolve, and the surviving boundary is rewound from the winding
/// numbers. The output is re-triangulated, so triangle and vertex counts do
/// not match the input.
///
/// The input must already be closed and orientable — that is the soup
/// importer's admission test, so an open mesh is an empty handle with status
/// 15 (`NotClosed`) long before it reaches here, and rebuilding it is a no-op.
/// Everything past closedness is fair game.
///
/// Returns NULL only if `m` is NULL, the rule value is unknown, or the rebuild
/// panics; an empty input comes back as a valid empty handle.
///
/// # Safety
/// `m` must be NULL or a live handle from this library.
#[no_mangle]
pub unsafe extern "C" fn manifold_rs_rebuild_solid(
    m: *const ManifoldRs,
    winding_rule: i32,
) -> *mut ManifoldRs {
    // SAFETY: the caller contract of this function is exactly the one
    // `manifold_rs_rebuild_solid_ct` documents; a NULL token means
    // "uncancellable", which is this function's behaviour.
    unsafe { manifold_rs_rebuild_solid_ct(m, winding_rule, ptr::null()) }
}

/// [`manifold_rs_rebuild_solid`] with an optional cancellation token.
///
/// A rebuild costs what a boolean against a partner costs, so anything
/// interactive wants this form. A NULL `token` is the uncancellable path and
/// behaves identically to [`manifold_rs_rebuild_solid`]. When the token is
/// cancelled the call returns a *valid* handle whose status is the cancelled
/// code (14), never NULL — so the caller can tell cancellation apart from an
/// argument error or a caught panic.
///
/// The token's shared flag is cloned out on entry, with the same ownership
/// reasoning as [`manifold_rs_batch_boolean_ct`]: destroying the token handle
/// mid-call is harmless, it just leaves later cancels nowhere to land.
///
/// # Safety
/// `m` must be NULL or a live handle from this library. `token` must be NULL
/// or a live handle at the moment of the call.
#[no_mangle]
pub unsafe extern "C" fn manifold_rs_rebuild_solid_ct(
    m: *const ManifoldRs,
    winding_rule: i32,
    token: *const CancelTokenRs,
) -> *mut ManifoldRs {
    guard(ptr::null_mut(), || {
        let rule = match winding_rule {
            0 => WindingRule::Positive,
            1 => WindingRule::Nonzero,
            other => {
                set_last_error(format!(
                    "manifold_rs_rebuild_solid: unknown winding rule {other}"
                ));
                return ptr::null_mut();
            }
        };
        // SAFETY: caller contract; as_ref() handles the NULL case.
        let Some(handle) = (unsafe { m.as_ref() }) else {
            set_last_error("manifold_rs_rebuild_solid: null manifold");
            return ptr::null_mut();
        };
        // SAFETY: caller contract; as_ref() handles the NULL case.
        let token: Option<CancelToken> = unsafe { token.as_ref() }.map(|t| t.inner.clone());
        into_handle(handle.inner.rebuild_solid_with_token(rule, token.as_ref()))
    })
}

/// 1 when two of `m`'s own triangles genuinely intersect — they cross, they
/// overlap, or they are coincident surface (a doubled or multiply-wound
/// sheet) — 0 when they merely share edges and vertices as every closed mesh
/// does, -1 for a NULL handle or a panic. A mesh with non-finite positions
/// answers 1, that being the safe verdict for geometry no exact predicate
/// can evaluate.
///
/// A mesh can be topologically manifold and still answer 1; those inputs
/// break the exact boolean engine's assumptions, which is why the `Auto`
/// engine routes them to the robust engine. The verdict is cached on the
/// mesh, so repeated calls (and the booleans that consult it) are free after
/// the first.
///
/// # Safety
/// `m` must be NULL or a live handle from this library.
#[no_mangle]
pub unsafe extern "C" fn manifold_rs_has_self_intersections(m: *const ManifoldRs) -> i32 {
    guard(-1, || {
        // SAFETY: caller contract; as_ref() handles the NULL case.
        match unsafe { m.as_ref() } {
            Some(handle) => i32::from(handle.inner.has_self_intersections()),
            None => {
                set_last_error("manifold_rs_has_self_intersections: null manifold");
                -1
            }
        }
    })
}

/// Copy of `m` re-tagged as an original mesh (a fresh mesh ID, no boolean
/// history). NULL if `m` is NULL or the copy panics.
///
/// # Safety
/// `m` must be NULL or a live handle from this library.
#[no_mangle]
pub unsafe extern "C" fn manifold_rs_as_original(m: *const ManifoldRs) -> *mut ManifoldRs {
    guard(ptr::null_mut(), || {
        // SAFETY: caller contract; as_ref() handles the NULL case.
        let Some(handle) = (unsafe { m.as_ref() }) else {
            set_last_error("manifold_rs_as_original: null manifold");
            return ptr::null_mut();
        };
        into_handle(handle.inner.as_original())
    })
}

/// Original mesh ID, or -1 for a NULL handle or a manifold that is not an
/// original (the core crate uses -1 for both).
///
/// # Safety
/// `m` must be NULL or a live handle from this library.
#[no_mangle]
pub unsafe extern "C" fn manifold_rs_original_id(m: *const ManifoldRs) -> i32 {
    guard(-1, || {
        // SAFETY: caller contract; as_ref() handles the NULL case.
        match unsafe { m.as_ref() } {
            Some(handle) => handle.inner.original_id(),
            None => {
                set_last_error("manifold_rs_original_id: null manifold");
                -1
            }
        }
    })
}

/// Error status code (0 = no error, see the table in `manifold_rs.h`).
/// Returns -1 for a NULL handle, which is not a valid status code.
///
/// # Safety
/// `m` must be NULL or a live handle from this library.
#[no_mangle]
pub unsafe extern "C" fn manifold_rs_status(m: *const ManifoldRs) -> i32 {
    guard(-1, || {
        // SAFETY: caller contract; as_ref() handles the NULL case.
        match unsafe { m.as_ref() } {
            Some(handle) => error_code(handle.inner.status()),
            None => {
                set_last_error("manifold_rs_status: null manifold");
                -1
            }
        }
    })
}

/// Release a manifold handle. NULL is ignored.
///
/// # Safety
/// `m` must be NULL or a handle from this library that has not already been
/// destroyed, and no other thread may be using it.
#[no_mangle]
pub unsafe extern "C" fn manifold_rs_destroy(m: *mut ManifoldRs) {
    guard((), || {
        if !m.is_null() {
            // SAFETY: caller guarantees single ownership of a live handle.
            drop(unsafe { Box::from_raw(m) });
        }
    })
}

// ---------------------------------------------------------------------------
// Boolean operations
// ---------------------------------------------------------------------------

/// N-ary boolean over `count` manifolds.
///
/// `op` is 0 = union, 1 = difference, 2 = intersection, matching the C++
/// `ManifoldOpType` ordering. Difference is *first operand minus the union of
/// the rest*, the same as C++ `BatchBoolean`.
///
/// This routes through the CSG tree rather than folding pairwise booleans, so
/// a union of many parts gets the bounding-box partitioning and heap-ordered
/// reduction the kernel was designed around. Returns NULL for an empty list,
/// an unknown op, a NULL entry, or a panic.
///
/// # Safety
/// `manifolds` must point to `count` live handles.
#[no_mangle]
pub unsafe extern "C" fn manifold_rs_batch_boolean(
    manifolds: *const *const ManifoldRs,
    count: usize,
    op: i32,
) -> *mut ManifoldRs {
    // SAFETY: the caller contract of this function is exactly the one
    // `manifold_rs_batch_boolean_ct` documents; a NULL token means
    // "uncancellable", which is this function's behaviour.
    unsafe { manifold_rs_batch_boolean_ct(manifolds, count, op, ptr::null()) }
}

/// [`manifold_rs_batch_boolean`] with an optional cancellation token.
///
/// A NULL `token` is the uncancellable path and behaves identically to
/// [`manifold_rs_batch_boolean`]. When the token is cancelled the call returns a
/// *valid* handle whose status is the cancelled code, never NULL — so the caller
/// can tell cancellation apart from an argument error or a caught panic.
///
/// The token's shared flag is cloned out on entry, so the operation keeps the
/// flag alive by itself: destroying the handle mid-call leaks nothing and
/// corrupts nothing, it just means later cancels have nowhere to land. (Still
/// poor hygiene — see the header.)
///
/// # Safety
/// `manifolds` must point to `count` live handles. `token` must be NULL or a
/// live handle at the moment of the call.
#[no_mangle]
pub unsafe extern "C" fn manifold_rs_batch_boolean_ct(
    manifolds: *const *const ManifoldRs,
    count: usize,
    op: i32,
    token: *const CancelTokenRs,
) -> *mut ManifoldRs {
    guard(ptr::null_mut(), || {
        let op = match op {
            0 => OpType::Add,
            1 => OpType::Subtract,
            2 => OpType::Intersect,
            other => {
                set_last_error(format!("manifold_rs_batch_boolean: unknown op {other}"));
                return ptr::null_mut();
            }
        };
        if manifolds.is_null() || count == 0 {
            set_last_error("manifold_rs_batch_boolean: no input manifolds");
            return ptr::null_mut();
        }
        if !length_fits::<*const ManifoldRs>(count) {
            set_last_error(format!(
                "manifold_rs_batch_boolean: implausible count {count}"
            ));
            return ptr::null_mut();
        }

        // Clone (an Arc bump) rather than borrow: the operation then owns a
        // share of the flag and no longer depends on the *handle* outliving the
        // call. A C# caller that drops the token while an await is in flight
        // gets a harmless no-op instead of a use-after-free. Only the tokened
        // path pays the refcount bump; the NULL path allocates nothing.
        // SAFETY: caller contract; as_ref() handles the NULL case.
        let token: Option<CancelToken> = unsafe { token.as_ref() }.map(|t| t.inner.clone());
        let token: Option<&CancelToken> = token.as_ref();

        // SAFETY: non-null with a count that is > 0 and within the isize::MAX
        // byte limit (both checked above); the caller guarantees the array
        // really holds that many handles.
        let handles = unsafe { std::slice::from_raw_parts(manifolds, count) };
        let mut inputs = Vec::with_capacity(count);
        for (i, &entry) in handles.iter().enumerate() {
            // SAFETY: caller contract; as_ref() handles the NULL case.
            let Some(handle) = (unsafe { entry.as_ref() }) else {
                set_last_error(format!(
                    "manifold_rs_batch_boolean: null manifold at index {i}"
                ));
                return ptr::null_mut();
            };
            inputs.push(&handle.inner);
        }

        // A single operand has nothing to combine with; every op is identity —
        // but an already-cancelled token still wins, so a caller polling only
        // the status never sees a one-operand call report success after cancel.
        if let [only] = inputs[..] {
            return match token {
                Some(t) if t.is_cancelled() => into_handle(Manifold::make_empty(Error::Cancelled)),
                _ => into_handle(only.clone()),
            };
        }

        let leaves = inputs
            .into_iter()
            .map(|m| CsgNode::leaf(m.as_impl().clone()))
            .collect();
        let result = CsgNode::op_n(op, leaves).evaluate_with_token(token);
        into_handle(Manifold::from_impl(result))
    })
}

// ---------------------------------------------------------------------------
// MeshGL export
// ---------------------------------------------------------------------------

/// Export `m` as a MeshGL handle (normal index -1: normals are exported only
/// if the manifold already carries them). An error-status manifold exports as
/// an empty mesh rather than failing. NULL on a NULL handle or a panic.
///
/// # Safety
/// `m` must be NULL or a live handle from this library.
#[no_mangle]
pub unsafe extern "C" fn manifold_rs_get_meshgl(m: *const ManifoldRs) -> *mut MeshGlRs {
    guard(ptr::null_mut(), || {
        // SAFETY: caller contract; as_ref() handles the NULL case.
        let Some(handle) = (unsafe { m.as_ref() }) else {
            set_last_error("manifold_rs_get_meshgl: null manifold");
            return ptr::null_mut();
        };
        Box::into_raw(Box::new(MeshGlRs {
            inner: handle.inner.get_mesh_gl(-1),
        }))
    })
}

/// Properties per vertex (>= 3; the first three are the position). 0 for a
/// NULL handle.
///
/// # Safety
/// `g` must be NULL or a live handle from this library.
#[no_mangle]
pub unsafe extern "C" fn manifold_rs_meshgl_num_prop(g: *const MeshGlRs) -> u32 {
    guard(0, || {
        // SAFETY: caller contract; as_ref() handles the NULL case.
        match unsafe { g.as_ref() } {
            Some(handle) => handle.inner.num_prop,
            None => {
                set_last_error("manifold_rs_meshgl_num_prop: null mesh");
                0
            }
        }
    })
}

/// Defines a borrowing array accessor over one field of a mesh handle
/// (`MeshGlRs` or `MeshGl64Rs`). All ten share the same contract, and writing
/// them out by hand ten times would only add places for the null/length
/// handling to drift apart.
///
/// `merge_from_vert` / `merge_to_vert` are deliberately not exposed: the export
/// path already resolves merges, so a caller reconstructing geometry from these
/// arrays does not need them.
macro_rules! mesh_array_accessor {
    ($name:ident, $handle:ty, $field:ident, $elem:ty) => {
        /// Borrowed pointer to the field's data, valid until the mesh handle is
        /// destroyed. Writes the element count to `out_len` when it is non-NULL.
        /// NULL (with `*out_len` = 0) for a NULL handle; for an empty field the
        /// pointer is non-NULL but must not be dereferenced.
        ///
        /// # Safety
        /// `g` must be NULL or a live handle; `out_len` must be NULL or
        /// writable.
        #[no_mangle]
        pub unsafe extern "C" fn $name(g: *const $handle, out_len: *mut usize) -> *const $elem {
            // Zero the length up front, outside the guard: if the body panics
            // the caller gets NULL and a length of 0 rather than whatever was
            // in their variable before the call.
            // SAFETY: out_len is null-checked inside write_len.
            unsafe { $crate::write_len(out_len, 0) };
            $crate::error::guard(ptr::null(), || {
                // SAFETY: caller contract; as_ref() handles the NULL case.
                let Some(handle) = (unsafe { g.as_ref() }) else {
                    $crate::error::set_last_error(concat!(stringify!($name), ": null mesh"));
                    return ptr::null();
                };
                let data = &handle.inner.$field;
                // SAFETY: out_len is null-checked inside write_len.
                unsafe { $crate::write_len(out_len, data.len()) };
                data.as_ptr()
            })
        }
    };
}

// `macro_rules` items are textually scoped, so re-export the macro for the
// sibling module that also needs it (`meshgl64`).
pub(crate) use mesh_array_accessor;

mesh_array_accessor!(
    manifold_rs_meshgl_vert_properties,
    MeshGlRs,
    vert_properties,
    f32
);
mesh_array_accessor!(manifold_rs_meshgl_tri_verts, MeshGlRs, tri_verts, u32);
mesh_array_accessor!(manifold_rs_meshgl_run_index, MeshGlRs, run_index, u32);
mesh_array_accessor!(
    manifold_rs_meshgl_run_original_id,
    MeshGlRs,
    run_original_id,
    u32
);
mesh_array_accessor!(manifold_rs_meshgl_face_id, MeshGlRs, face_id, u32);

/// Release a mesh handle, invalidating every pointer previously returned by the
/// array accessors. NULL is ignored.
///
/// # Safety
/// `g` must be NULL or a handle from this library that has not already been
/// destroyed, and no other thread may be using it.
#[no_mangle]
pub unsafe extern "C" fn manifold_rs_meshgl_destroy(g: *mut MeshGlRs) {
    guard((), || {
        if !g.is_null() {
            // SAFETY: caller guarantees single ownership of a live handle.
            drop(unsafe { Box::from_raw(g) });
        }
    })
}

// ---------------------------------------------------------------------------
// Error reporting
// ---------------------------------------------------------------------------

/// Copy this thread's last failure message (UTF-8, *not* NUL-terminated) into
/// `buf`, and return the full message length in bytes — so a caller can size a
/// buffer by calling once with `buf` = NULL. The message is set when a function
/// returns a failure sentinel and is not cleared by later successful calls.
///
/// # Safety
/// `buf` must be NULL or writable for `buf_len` bytes.
#[no_mangle]
pub unsafe extern "C" fn manifold_rs_last_error(buf: *mut u8, buf_len: usize) -> usize {
    guard(0, || {
        // Copy the message out before touching the caller's buffer: a panic
        // while the thread-local was borrowed would poison later reads.
        let message = last_error_message();
        let bytes = message.as_bytes();
        if !buf.is_null() && buf_len > 0 {
            let n = bytes.len().min(buf_len);
            // SAFETY: buf is non-null and the caller guarantees buf_len bytes;
            // n is clamped to both lengths and the regions cannot overlap (the
            // source is a freshly allocated String).
            unsafe { ptr::copy_nonoverlapping(bytes.as_ptr(), buf, n) };
        }
        bytes.len()
    })
}

/// Whether `len` elements of `T` can legally form a slice: `from_raw_parts`
/// requires the total size to fit in `isize::MAX` bytes, and violating that is
/// immediate UB. A caller that widens a negative int to `usize` — an easy
/// mistake in a P/Invoke signature — lands here instead of in the kernel.
pub(crate) fn length_fits<T>(len: usize) -> bool {
    match std::mem::size_of::<T>() {
        0 => true,
        size => len <= (isize::MAX as usize) / size,
    }
}

/// `from_raw_parts` requires a non-null aligned pointer even for length 0, so
/// route the empty case around it.
///
/// # Safety
/// When `len > 0`, `ptr` must be valid for `len` elements and `len` must have
/// passed [`length_fits`].
pub(crate) unsafe fn slice_or_empty<'a, T>(ptr: *const T, len: usize) -> &'a [T] {
    if len == 0 {
        &[]
    } else {
        // SAFETY: caller contract, plus the non-null check at the call site.
        unsafe { std::slice::from_raw_parts(ptr, len) }
    }
}

/// # Safety
/// `out_len` must be NULL or writable.
pub(crate) unsafe fn write_len(out_len: *mut usize, len: usize) {
    if !out_len.is_null() {
        // SAFETY: non-null and caller guarantees it is writable.
        unsafe { *out_len = len };
    }
}
