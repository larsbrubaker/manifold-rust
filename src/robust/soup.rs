// robust/soup.rs — Triangle-soup import support for the robust boolean
// engine.
//
// When the strict halfedge pairing of `Manifold::from_mesh_gl` fails
// (non-manifold connectivity), `from_mesh_gl_robust` falls through to
// `soupify`: the geometry is kept as an unpaired-halfedge triangle soup
// inside `ManifoldImpl` (is_soup = true), provided it passes the one check
// the robust engine genuinely needs — the soup must be geometrically
// **closed and orientable**: after welding vertices by exact position and
// dropping exactly-degenerate triangles (paper §5), every directed edge
// must be balanced by its reverse. That is precisely the condition for the
// soup to bound a solid via winding numbers.

use std::collections::BTreeMap;

use crate::impl_mesh::ManifoldImpl;
use crate::linalg::{IVec3, Vec3};
use crate::types::{Error, Halfedge};

use super::exact::rational::R3;

/// Weld key: exact position identity (with -0.0 normalized so it equals 0.0).
fn pos_key(v: Vec3) -> (u64, u64, u64) {
    let norm = |x: f64| if x == 0.0 { 0.0f64 } else { x }.to_bits();
    (norm(v.x), norm(v.y), norm(v.z))
}

/// Exact zero-area test on f64 positions.
fn is_degenerate(a: Vec3, b: Vec3, c: Vec3) -> bool {
    R3::from_vec3(b)
        .sub(&R3::from_vec3(a))
        .cross(&R3::from_vec3(c).sub(&R3::from_vec3(a)))
        .is_zero()
}

/// Convert `imp` (with `vert_pos` and `mesh_relation.tri_ref` already
/// populated, and `halfedge` full of whatever strict pairing produced) into
/// a validated triangle soup:
///  - drop exactly-degenerate triangles (and their tri_refs),
///  - verify the remainder is closed and orientable (else `Error::NotClosed`),
///  - rebuild halfedges with best-effort pairing (`paired_halfedge == -1`
///    where no partner exists),
///  - recompute per-face normals (no pairing required),
///  - set `is_soup`.
///
/// `tri_prop` / `tri_vert` mirror the `create_halfedges` inputs: `tri_vert`
/// holds the position indices when properties are mapped separately, else
/// `tri_prop` does double duty.
pub fn soupify(
    imp: &mut ManifoldImpl,
    tri_prop: &[IVec3],
    tri_vert: &[IVec3],
) -> Result<(), Error> {
    let position_tris: &[IVec3] = if tri_vert.is_empty() { tri_prop } else { tri_vert };
    debug_assert_eq!(position_tris.len(), tri_prop.len());

    // Weld vertex ids by exact position for the closedness bookkeeping only
    // (vert_pos itself is left as imported).
    let mut weld: BTreeMap<(u64, u64, u64), i32> = BTreeMap::new();
    let mut welded_id = vec![0i32; imp.vert_pos.len()];
    for (i, &p) in imp.vert_pos.iter().enumerate() {
        let id = *weld.entry(pos_key(p)).or_insert(i as i32);
        welded_id[i] = id;
    }

    // Keep non-degenerate triangles; balance directed edges on welded ids.
    let mut keep: Vec<usize> = Vec::with_capacity(position_tris.len());
    let mut balance: BTreeMap<(i32, i32), i64> = BTreeMap::new();
    for (t, tv) in position_tris.iter().enumerate() {
        let (a, b, c) = (tv.x as usize, tv.y as usize, tv.z as usize);
        let (wa, wb, wc) = (welded_id[a], welded_id[b], welded_id[c]);
        if wa == wb || wb == wc || wc == wa
            || is_degenerate(imp.vert_pos[a], imp.vert_pos[b], imp.vert_pos[c])
        {
            continue; // paper §5: degenerate triangles are safe to drop
        }
        keep.push(t);
        for (u, v) in [(wa, wb), (wb, wc), (wc, wa)] {
            let key = (u.min(v), u.max(v));
            *balance.entry(key).or_insert(0) += if u < v { 1 } else { -1 };
        }
    }
    if keep.len() < 4 {
        return Err(Error::NotClosed);
    }
    if balance.values().any(|&n| n != 0) {
        return Err(Error::NotClosed);
    }

    // Rebuild halfedges for the kept triangles with best-effort pairing:
    // LIFO multimap on welded undirected edges, forward (start < end) pairs
    // with reverse.
    let has_props = !tri_vert.is_empty();
    let mut halfedges: Vec<Halfedge> = Vec::with_capacity(3 * keep.len());
    let mut tri_ref = Vec::with_capacity(keep.len());
    for (new_t, &old_t) in keep.iter().enumerate() {
        let tv = position_tris[old_t];
        let tp = tri_prop[old_t];
        for i in 0..3 {
            let j = (i + 1) % 3;
            halfedges.push(Halfedge {
                start_vert: tv[i],
                end_vert: tv[j],
                paired_halfedge: -1,
                prop_vert: if has_props { tp[i] } else { tv[i] },
            });
        }
        if old_t < imp.mesh_relation.tri_ref.len() {
            tri_ref.push(imp.mesh_relation.tri_ref[old_t]);
        }
        let _ = new_t;
    }
    let mut open: BTreeMap<(i32, i32), Vec<usize>> = BTreeMap::new();
    for (idx, he) in halfedges.iter().enumerate() {
        let (u, v) = (welded_id[he.start_vert as usize], welded_id[he.end_vert as usize]);
        open.entry((u.min(v), u.max(v))).or_default().push(idx);
    }
    for (_key, mut idxs) in open {
        // Pair forwards with reverses greedily; leftovers stay -1.
        let mut fwd: Vec<usize> = Vec::new();
        let mut bwd: Vec<usize> = Vec::new();
        for idx in idxs.drain(..) {
            let he = &halfedges[idx];
            if welded_id[he.start_vert as usize] < welded_id[he.end_vert as usize] {
                fwd.push(idx);
            } else {
                bwd.push(idx);
            }
        }
        while let (Some(f), Some(b)) = (fwd.pop(), bwd.pop()) {
            halfedges[f].paired_halfedge = b as i32;
            halfedges[b].paired_halfedge = f as i32;
        }
    }
    imp.halfedge = halfedges;
    if tri_ref.len() == keep.len() {
        imp.mesh_relation.tri_ref = tri_ref;
    } else {
        imp.mesh_relation.tri_ref.clear();
    }

    // Per-face normals need only the triangle itself.
    imp.face_normal = (0..imp.num_tri())
        .map(|t| {
            let a = imp.vert_pos[imp.halfedge[3 * t].start_vert as usize];
            let b = imp.vert_pos[imp.halfedge[3 * t + 1].start_vert as usize];
            let c = imp.vert_pos[imp.halfedge[3 * t + 2].start_vert as usize];
            let n = crate::linalg::cross(b - a, c - a);
            let len = crate::linalg::length(n);
            if len > 0.0 { n / len } else { Vec3::new(0.0, 0.0, 0.0) }
        })
        .collect();
    imp.vert_normal.clear();
    imp.is_soup = true;
    Ok(())
}

#[cfg(test)]
#[path = "soup_tests.rs"]
mod tests;

/// Per-corner vertex properties of an impl, flattened as
/// `props[(3*tri + corner) * num_prop + channel]`, aligned with
/// `impl_to_tris` ordering. Empty when the impl carries no properties.
pub fn impl_to_corner_props(imp: &ManifoldImpl) -> Vec<f64> {
    let np = imp.num_prop;
    if np == 0 {
        return Vec::new();
    }
    let mut out = Vec::with_capacity(3 * imp.num_tri() * np);
    for t in 0..imp.num_tri() {
        for i in 0..3 {
            let pv = imp.halfedge[3 * t + i].prop_vert as usize;
            out.extend_from_slice(&imp.properties[pv * np..(pv + 1) * np]);
        }
    }
    out
}

/// The triangle list of any impl (soup or manifold) as position triples —
/// the robust engine's working form.
pub fn impl_to_tris(imp: &ManifoldImpl) -> Vec<[Vec3; 3]> {
    (0..imp.num_tri())
        .map(|t| {
            [
                imp.vert_pos[imp.halfedge[3 * t].start_vert as usize],
                imp.vert_pos[imp.halfedge[3 * t + 1].start_vert as usize],
                imp.vert_pos[imp.halfedge[3 * t + 2].start_vert as usize],
            ]
        })
        .collect()
}
