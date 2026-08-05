/*
 * Copyright 2026 Lars Brubaker
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*
 * C ABI for manifold-rust (library name: manifold_rs).
 *
 * This header is hand written and is the normative description of the ABI:
 * ownership, sentinels and lifetimes are documented here, not inferred from
 * the Rust source. It covers the subset of the C++ manifoldc API needed to
 * import a triangle soup, run an n-ary boolean and export the result, using
 * the same integer op ordering as manifoldc.
 *
 * General contract
 * ----------------
 * - Handles are opaque and owned by the caller. Every handle returned
 *   non-NULL must eventually be passed to its destroy function.
 * - No function ever unwinds or aborts: a panic inside the geometry kernel is
 *   caught and turned into the function's failure sentinel (NULL / -1 / 0),
 *   with a message retrievable via manifold_rs_last_error().
 * - Every pointer argument is null-checked. Passing NULL yields the documented
 *   sentinel instead of a crash.
 * - Thread safety: distinct handles may be used concurrently from different
 *   threads. A single handle must not be used concurrently with any call that
 *   mutates or frees it (in particular its destroy function). The last-error
 *   slot is per thread.
 *   The one deliberate exception is CancelTokenRs: see the cancellation
 *   section below, where concurrent use from two threads is the entire point.
 * - Strings and arrays returned by this library are never freed by the caller;
 *   they are either static or owned by a handle.
 *
 * Error codes returned by manifold_rs_status()
 * --------------------------------------------
 *   -1  null handle (not a real status)
 *    0  no error
 *    1  non-finite vertex
 *    2  not manifold
 *    3  vertex out of bounds
 *    4  properties wrong length
 *    5  missing position properties
 *    6  merge vectors different lengths
 *    7  merge index out of bounds
 *    8  transform wrong length
 *    9  run index wrong length
 *   10  face ID wrong length
 *   11  invalid construction
 *   12  result too large
 *   13  invalid tangents
 *   14  cancelled (the operation was interrupted through a CancelTokenRs)
 */

#ifndef MANIFOLD_RS_H
#define MANIFOLD_RS_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Opaque handle to a manifold solid. Free with manifold_rs_destroy(). */
typedef struct ManifoldRs ManifoldRs;

/* Opaque handle to an exported GL-style mesh. Free with
 * manifold_rs_meshgl_destroy(). */
typedef struct MeshGlRs MeshGlRs;

/* Opaque handle to an exported double-precision / 64-bit-index mesh. Free with
 * manifold_rs_meshgl64_destroy(). */
typedef struct MeshGl64Rs MeshGl64Rs;

/* Opaque handle to a cancellation token. Free with
 * manifold_rs_cancel_token_destroy(). */
typedef struct CancelTokenRs CancelTokenRs;

/* Boolean op codes; the ordering matches C++ ManifoldOpType. */
#define MANIFOLD_RS_OP_ADD 0       /* union */
#define MANIFOLD_RS_OP_SUBTRACT 1  /* difference */
#define MANIFOLD_RS_OP_INTERSECT 2 /* intersection */

/*
 * Static NUL-terminated version string, e.g.
 * "manifold-ffi 0.1.0 (manifold-rust 0.9.3)". Never NULL, never freed.
 */
const char* manifold_rs_version(void);

/*
 * Build a manifold from interleaved vertex properties and triangle indices.
 *
 * vert_properties: num_prop floats per vertex; the first three are x, y, z.
 * vert_properties_len: total float count (must be a multiple of num_prop).
 * tri_verts: three vertex indices per triangle, CCW seen from outside.
 * tri_verts_len: total index count (must be a multiple of 3).
 * num_prop: properties per vertex, >= 3.
 *
 * Both arrays are copied; the caller may free them immediately after the call.
 *
 * Returns a handle, including when the mesh fails validation - in that case
 * the handle carries a non-zero manifold_rs_status() describing the failure,
 * which matches the C++ behaviour of returning an error-status manifold.
 * Returns NULL only for arguments that cannot describe a mesh (num_prop < 3,
 * a length that is not evenly divisible, or a NULL array with a non-zero
 * length) or if the import panics.
 */
ManifoldRs* manifold_rs_from_mesh(const float* vert_properties,
                                  size_t vert_properties_len,
                                  const uint32_t* tri_verts,
                                  size_t tri_verts_len,
                                  uint32_t num_prop);

/*
 * Copy of m re-tagged as an original mesh: it is given a fresh mesh ID and its
 * boolean history is dropped, so results derived from it report that ID in
 * their run_original_id. Returns NULL if m is NULL or the copy panics.
 *
 * An empty manifold - which includes every manifold with a non-zero status -
 * is returned as a plain copy instead, keeping original_id == -1. Do not
 * assume a fresh ID was assigned: check manifold_rs_status(m) first, or check
 * manifold_rs_original_id() on the result.
 */
ManifoldRs* manifold_rs_as_original(const ManifoldRs* m);

/*
 * Original mesh ID of m, or -1 if m is NULL or is not an original mesh
 * (meshes imported with manifold_rs_from_mesh are not originals until
 * manifold_rs_as_original is called).
 */
int32_t manifold_rs_original_id(const ManifoldRs* m);

/* Error status of m (see the code table above). -1 if m is NULL. */
int32_t manifold_rs_status(const ManifoldRs* m);

/* Free a manifold handle. NULL is ignored. */
void manifold_rs_destroy(ManifoldRs* m);

/*
 * N-ary boolean over count manifolds.
 *
 * manifolds: array of count non-NULL handles; they are not consumed and the
 *   caller still owns them afterwards.
 * op: one of MANIFOLD_RS_OP_*.
 *
 * MANIFOLD_RS_OP_SUBTRACT is the first operand minus the union of all the
 * others, i.e. (a - b - c ...), matching C++ BatchBoolean.
 *
 * The operation runs through the CSG tree, so a union of many parts gets the
 * bounding-box partitioning and heap-ordered reduction that make large unions
 * tractable - prefer one call with n operands over n-1 pairwise calls.
 *
 * count == 1 returns a copy of the single operand. Returns NULL for count == 0,
 * a NULL array or entry, an unrecognised op, or a panic.
 *
 * IMPORTANT: an operand carrying a non-zero status is absorbed as empty
 * geometry, and the result can still report status 0 - a failed import
 * therefore shows up as a part silently missing from the output, not as an
 * error here. Check manifold_rs_status() on every operand before combining
 * them.
 */
ManifoldRs* manifold_rs_batch_boolean(const ManifoldRs* const* manifolds,
                                      size_t count,
                                      int32_t op);

/*
 * Cancellation
 * ------------
 *
 * A CancelTokenRs is a shared flag that lets one thread interrupt a long
 * running operation running on another. Usage:
 *
 *   CancelTokenRs* t = manifold_rs_cancel_token_new();
 *   // thread A:  r = manifold_rs_batch_boolean_ct(parts, n, op, t);
 *   // thread B:  manifold_rs_cancel_token_cancel(t);   // user hit Cancel
 *   // thread A, after the call returns:
 *   //   manifold_rs_status(r) == 14  ->  cancelled
 *   manifold_rs_cancel_token_destroy(t);   // see LIFETIME below
 *
 * THREAD SAFETY (the exception to the general rule above): calling
 * manifold_rs_cancel_token_cancel() / manifold_rs_cancel_token_is_cancelled()
 * concurrently with an operation that was handed the same token, from any
 * thread, is explicitly supported - that is what the type is for. The flag is
 * atomic. Any number of threads may hold and cancel the same token.
 *
 * LIFETIME: an operation clones the token's shared flag when it starts, so the
 * flag outlives the handle. Calling manifold_rs_cancel_token_destroy() while an
 * operation that was given the token is still running is therefore MEMORY SAFE
 * - not a use after free, and it leaks nothing. It is still a bug: after the
 * destroy there is no handle left to cancel through, so the operation runs to
 * completion. Keep the handle alive (a GC-rooted field, a `using` scope that
 * outlives the await, etc.) for as long as you might want to cancel.
 *
 * The one ordering the caller must still honour is destroy vs.
 * manifold_rs_cancel_token_cancel() / manifold_rs_cancel_token_is_cancelled():
 * those read the handle itself, so they must not race with its destruction.
 *
 * Cancellation is STICKY: once cancelled a token stays cancelled, and every
 * later operation given that token returns immediately as cancelled. Create a
 * new token for work that should be allowed to run.
 */

/* A fresh, uncancelled token. NULL only if the allocation panics. */
CancelTokenRs* manifold_rs_cancel_token_new(void);

/* Request cancellation. Callable from any thread at any time; NULL is
 * ignored. */
void manifold_rs_cancel_token_cancel(const CancelTokenRs* t);

/* 1 if cancellation has been requested, else 0. A NULL token reports 0, which
 * is also its behaviour: NULL means uncancellable. */
int32_t manifold_rs_cancel_token_is_cancelled(const CancelTokenRs* t);

/* Free a token handle. NULL is ignored. Safe against a running operation that
 * holds the token, but see the LIFETIME note above: doing so loses the ability
 * to cancel that operation. */
void manifold_rs_cancel_token_destroy(CancelTokenRs* t);

/*
 * manifold_rs_batch_boolean with an optional cancellation token.
 *
 * token may be NULL, which means "uncancellable" and makes this call exactly
 * equivalent to manifold_rs_batch_boolean (which is implemented by delegating
 * here with a NULL token). Passing NULL costs nothing: the geometry kernel
 * never reads a flag on that path.
 *
 * On cancellation this returns a VALID handle whose manifold_rs_status() is 14
 * (cancelled) and whose geometry is empty - it does NOT return NULL. NULL is
 * reserved for argument errors and caught panics, so a caller can distinguish
 * "the user cancelled" from "something went wrong". The returned handle must
 * still be destroyed.
 *
 * Cancellation is cooperative and checked at phase boundaries inside the
 * kernel, so the call returns shortly after cancel is requested rather than
 * instantly.
 */
ManifoldRs* manifold_rs_batch_boolean_ct(const ManifoldRs* const* manifolds,
                                         size_t count,
                                         int32_t op,
                                         const CancelTokenRs* token);

/*
 * Export m as a mesh handle. A manifold with a non-zero status exports as an
 * empty mesh rather than failing. Returns NULL if m is NULL or the export
 * panics.
 */
MeshGlRs* manifold_rs_get_meshgl(const ManifoldRs* m);

/* Properties per vertex (>= 3; slots 0-2 are the position). 0 if g is NULL. */
uint32_t manifold_rs_meshgl_num_prop(const MeshGlRs* g);

/*
 * The five array accessors below all follow the same contract:
 *
 * - The returned pointer borrows from g and stays valid until
 *   manifold_rs_meshgl_destroy(g). Copy the data out before destroying it.
 * - out_len receives the element count (not the byte count). out_len may be
 *   NULL if the length is already known. It is set to 0 before any work is
 *   done, so it is never left holding a stale value.
 * - If g is NULL the result is NULL and *out_len is set to 0.
 * - If the array is empty the result is non-NULL but must not be
 *   dereferenced; check *out_len first.
 *
 * The merge_from_vert / merge_to_vert arrays are deliberately not exposed: the
 * export path has already resolved merges, so they are not needed to rebuild
 * the geometry.
 *
 * Field meanings:
 *   vert_properties  interleaved, num_prop floats per vertex.
 *   tri_verts        three vertex indices per triangle, CCW from outside.
 *   run_index        start offsets into tri_verts, one per run plus a trailing
 *                    end sentinel (so it has run_original_id_len + 1 entries).
 *   run_original_id  the source mesh ID for each run.
 *   face_id          the source face each triangle came from, one per triangle.
 */
const float* manifold_rs_meshgl_vert_properties(const MeshGlRs* g,
                                                size_t* out_len);
const uint32_t* manifold_rs_meshgl_tri_verts(const MeshGlRs* g, size_t* out_len);
const uint32_t* manifold_rs_meshgl_run_index(const MeshGlRs* g, size_t* out_len);
const uint32_t* manifold_rs_meshgl_run_original_id(const MeshGlRs* g,
                                                   size_t* out_len);
const uint32_t* manifold_rs_meshgl_face_id(const MeshGlRs* g, size_t* out_len);

/* Free a mesh handle, invalidating every pointer the accessors returned for
 * it. NULL is ignored. */
void manifold_rs_meshgl_destroy(MeshGlRs* g);

/*
 * Double-precision mesh path
 * --------------------------
 *
 * The MeshGl64Rs family mirrors the MeshGlRs family with wider element types,
 * for callers whose geometry is already double / uint64 and who would
 * otherwise have to narrow it by hand. Element types follow the core crate's
 * MeshGL64 exactly:
 *
 *   vert_properties  double
 *   tri_verts        uint64_t
 *   run_index        uint64_t
 *   run_original_id  uint32_t   <- a mesh ID, not an index; u32 in both widths
 *   face_id          uint64_t
 *
 * num_prop stays uint32_t on both paths: it is a small per-vertex count, and
 * keeping the signature identical makes the two entry points interchangeable
 * at the call site.
 *
 * NOTE: the kernel currently narrows to single precision and 32 bit indices
 * internally, so this is a wider interface, not more capacity end to end:
 *   - Coordinates round trip through float, losing ~1e-7 of relative
 *     precision. Harmless at any physical tolerance.
 *   - Indices are narrowed to uint32_t. That would WRAP, not saturate, so
 *     manifold_rs_from_mesh64 rejects any tri_verts entry greater than
 *     UINT32_MAX (returns NULL, with the offending index in the last-error
 *     message) rather than silently building wrong geometry that reports
 *     status 0. In practice this caps a mesh at ~4 billion vertices.
 *
 * Handles produced by manifold_rs_from_mesh64 are ordinary ManifoldRs handles
 * and may be mixed freely with ones from manifold_rs_from_mesh.
 */

/*
 * Double-precision counterpart of manifold_rs_from_mesh. Same argument rules,
 * same sentinels: a handle (possibly carrying a non-zero status) for anything
 * that describes a mesh, NULL for arguments that cannot. Both arrays are
 * copied.
 */
ManifoldRs* manifold_rs_from_mesh64(const double* vert_properties,
                                    size_t vert_properties_len,
                                    const uint64_t* tri_verts,
                                    size_t tri_verts_len,
                                    uint32_t num_prop);

/*
 * Export m as a double-precision mesh handle. A manifold with a non-zero
 * status exports as an empty mesh rather than failing. NULL if m is NULL or
 * the export panics.
 */
MeshGl64Rs* manifold_rs_get_meshgl64(const ManifoldRs* m);

/* Properties per vertex (>= 3; slots 0-2 are the position). 0 if g is NULL. */
uint32_t manifold_rs_meshgl64_num_prop(const MeshGl64Rs* g);

/* Same borrowing / out_len / NULL contract as the MeshGlRs accessors above,
 * and the same field meanings. */
const double* manifold_rs_meshgl64_vert_properties(const MeshGl64Rs* g,
                                                   size_t* out_len);
const uint64_t* manifold_rs_meshgl64_tri_verts(const MeshGl64Rs* g,
                                               size_t* out_len);
const uint64_t* manifold_rs_meshgl64_run_index(const MeshGl64Rs* g,
                                               size_t* out_len);
const uint32_t* manifold_rs_meshgl64_run_original_id(const MeshGl64Rs* g,
                                                     size_t* out_len);
const uint64_t* manifold_rs_meshgl64_face_id(const MeshGl64Rs* g,
                                             size_t* out_len);

/* Free a 64-bit mesh handle, invalidating every pointer the accessors returned
 * for it. NULL is ignored. */
void manifold_rs_meshgl64_destroy(MeshGl64Rs* g);

/*
 * Copy the calling thread's last failure message into buf as UTF-8 and return
 * the full message length in bytes.
 *
 * The message is NOT NUL-terminated: use the returned length. If the message
 * is longer than buf_len only buf_len bytes are written, but the full length
 * is still returned, so calling with buf = NULL and buf_len = 0 is a way to
 * size a buffer. Returns 0 when this thread has not seen a failure yet.
 *
 * The slot is written whenever a function returns a failure sentinel
 * (including caught panics) and is not cleared by later successful calls, so
 * read it immediately after observing the failure.
 */
size_t manifold_rs_last_error(uint8_t* buf, size_t buf_len);

#ifdef __cplusplus
}
#endif

#endif /* MANIFOLD_RS_H */
