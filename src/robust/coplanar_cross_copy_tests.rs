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

// Tests for the robust engine's phase 3, the coplanar cross-copy
// (robust/intersection_graph.rs, `cross_copy_coplanar_regions`). They pin two
// things: that the cross-copy's optimizations (per-region clip setup,
// bounding-box reject, hashed dedupe — robust/coplanar_clip.rs) leave
// the robust result bit-for-bit where the straight per-call version put it,
// and that the step reports its own progress phase rather than running
// silently under "self intersections". Shared 1:1 with manifold-sharp's
// CoplanarCrossCopyTests (its commit c776525), frozen hash included (the
// SHA-256 is computed here, test-only, rather than taken from a crate).

use crate::linalg::Vec3;
use crate::manifold::Manifold;
use crate::progress::{Phase, ProgressReporter};
use crate::types::{BooleanEngine, Error, OpType};

/// SHA-256 of the robust union's MeshGL64 (tri verts, then vertex property
/// bits; see [`hash`]), captured from the per-call cross-copy before it was
/// optimized — in manifold-sharp and, identically, in this crate. Any change
/// here is a change to a computed result.
const FROZEN_UNION_HASH: &str = "8E910C9A57BC34E4421978487DDFFB679CB62180E3B94B5C162877B30EA5CE76";

/// A 4x4 grid of unit cubes composed without a boolean, so neighbours touch
/// face to face (the self-contact that sends Auto to the robust engine), and a
/// slab whose bottom is coplanar with the grid's tops, offset so its
/// triangles straddle many grid triangles: every slab-bottom triangle
/// collects segments from many coplanar overlap regions, which is the shape
/// that made phase 3 expensive.
pub(crate) fn fixture() -> (Manifold, Manifold) {
    let mut cubes = Vec::new();
    for x in 0..4 {
        for y in 0..4 {
            cubes.push(
                Manifold::cube(Vec3::splat(1.0), false)
                    .translate(Vec3::new(x as f64, y as f64, 0.0)),
            );
        }
    }
    let body = Manifold::compose(&cubes);
    let slab =
        Manifold::cube(Vec3::new(3.0, 3.0, 0.5), false).translate(Vec3::new(0.25, 0.75, 1.0));
    (body, slab)
}

/// The bytes manifold-sharp's `Hash` feeds SHA-256 (a .NET `BinaryWriter`,
/// little-endian): the tri-vert count as i32, each tri vert as u64, the
/// vertex-property count as i32, each property's bits as i64. Uppercase hex.
pub(crate) fn hash(m: &Manifold) -> String {
    let mesh = m.get_mesh_gl64(-1);
    let mut bytes = Vec::new();
    bytes.extend((mesh.tri_verts.len() as i32).to_le_bytes());
    for v in &mesh.tri_verts {
        bytes.extend(v.to_le_bytes());
    }
    bytes.extend((mesh.vert_properties.len() as i32).to_le_bytes());
    for p in &mesh.vert_properties {
        bytes.extend(p.to_bits().to_le_bytes());
    }
    sha256(&bytes).iter().map(|b| format!("{b:02X}")).collect()
}

/// FIPS 180-4 SHA-256, test-only, so the frozen hash needs no dependency.
/// (Checked against `shasum -a 256` on the FIPS "abc" and two-block examples
/// when written; the frozen hash matching manifold-sharp's, computed by .NET's
/// SHA256, is the standing check.)
fn sha256(message: &[u8]) -> [u8; 32] {
    const K: [u32; 64] = [
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4,
        0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe,
        0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f,
        0x4a7484aa, 0x5cb0a9dc, 0x76f988da, 0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
        0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc,
        0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b,
        0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070, 0x19a4c116,
        0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7,
        0xc67178f2,
    ];
    let mut h: [u32; 8] = [
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab,
        0x5be0cd19,
    ];
    let mut data = message.to_vec();
    data.push(0x80);
    while data.len() % 64 != 56 {
        data.push(0);
    }
    data.extend(((message.len() as u64) * 8).to_be_bytes());
    for chunk in data.chunks(64) {
        let mut w = [0u32; 64];
        for (i, word) in chunk.chunks(4).enumerate() {
            w[i] = u32::from_be_bytes([word[0], word[1], word[2], word[3]]);
        }
        for i in 16..64 {
            let s0 = w[i - 15].rotate_right(7) ^ w[i - 15].rotate_right(18) ^ (w[i - 15] >> 3);
            let s1 = w[i - 2].rotate_right(17) ^ w[i - 2].rotate_right(19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16]
                .wrapping_add(s0)
                .wrapping_add(w[i - 7])
                .wrapping_add(s1);
        }
        let mut v = h;
        for i in 0..64 {
            let s1 = v[4].rotate_right(6) ^ v[4].rotate_right(11) ^ v[4].rotate_right(25);
            let ch = (v[4] & v[5]) ^ (!v[4] & v[6]);
            let t1 = v[7]
                .wrapping_add(s1)
                .wrapping_add(ch)
                .wrapping_add(K[i])
                .wrapping_add(w[i]);
            let s0 = v[0].rotate_right(2) ^ v[0].rotate_right(13) ^ v[0].rotate_right(22);
            let maj = (v[0] & v[1]) ^ (v[0] & v[2]) ^ (v[1] & v[2]);
            let t2 = s0.wrapping_add(maj);
            v = [
                t1.wrapping_add(t2),
                v[0],
                v[1],
                v[2],
                v[3].wrapping_add(t1),
                v[4],
                v[5],
                v[6],
            ];
        }
        for (hi, vi) in h.iter_mut().zip(v) {
            *hi = hi.wrapping_add(vi);
        }
    }
    let mut out = [0u8; 32];
    for (i, word) in h.iter().enumerate() {
        out[4 * i..4 * i + 4].copy_from_slice(&word.to_be_bytes());
    }
    out
}

/// The union's bits are the ones the unoptimized cross-copy produced (in both
/// builds: run the suite with and without `--features parallel`).
#[test]
fn robust_union_over_coplanar_faces_is_bit_identical_to_the_frozen_result() {
    let (body, slab) = fixture();
    let union = body.boolean_with_engine(&slab, OpType::Add, BooleanEngine::Robust);
    assert_eq!(union.status(), Error::NoError);
    assert_eq!(
        hash(&union),
        FROZEN_UNION_HASH,
        "parallel={}",
        cfg!(feature = "parallel")
    );
}

/// A robust boolean with coplanar overlap regions reports the cross-copy as
/// its own determinate phase, closed at exactly 1.0, between self
/// intersections and candidate points.
#[test]
fn the_cross_copy_reports_its_own_phase() {
    use std::sync::{Arc, Mutex};

    let (body, slab) = fixture();
    let events = Arc::new(Mutex::new(Vec::<(&'static str, Option<f64>)>::new()));
    let sink = Arc::clone(&events);
    let reporter = ProgressReporter::new(move |phase: Phase, fraction| {
        sink.lock()
            .expect("sink poisoned")
            .push((phase.name(), fraction));
    });
    body.boolean_with_engine_and_progress(
        &slab,
        OpType::Add,
        BooleanEngine::Robust,
        None,
        Some(&reporter),
    );
    let events = events.lock().expect("sink poisoned").clone();

    let mut order: Vec<&str> = Vec::new();
    for (name, _) in &events {
        if order.last() != Some(name) {
            order.push(name);
        }
    }
    let cross = order
        .iter()
        .position(|&n| n == Phase::CoplanarOverlaps.name())
        .unwrap_or_else(|| panic!("phases seen: {order:?}"));
    assert!(cross > 0, "phases seen: {order:?}");
    assert_eq!(order[cross - 1], Phase::SelfIntersections.name());
    assert_eq!(order[cross + 1], Phase::CandidatePoints.name());

    let last = events
        .iter()
        .rev()
        .find(|(name, _)| *name == Phase::CoplanarOverlaps.name())
        .and_then(|(_, fraction)| *fraction);
    assert_eq!(last, Some(1.0));
}
