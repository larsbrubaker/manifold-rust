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
