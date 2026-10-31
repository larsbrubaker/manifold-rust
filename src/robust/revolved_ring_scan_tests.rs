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

// Tests for the self-intersection scan (robust/soup.rs
// `scan_triangles`) on a finely revolved ring: a 10 x 10 rectangle at
// radius 20..30 turned through 20,000 sides. Its flat top and bottom are
// annuli of thin coplanar quads whose boxes overlap many of their
// neighbours', and no vertex reaches a flat star's size, so every one of
// those coplanar pairs reaches the narrow phase.

use std::time::Instant;

use crate::linalg::Vec3;

type Tri = [Vec3; 3];

/// The ring's triangles, wound outward: per side, two each on the bottom,
/// the outer wall, the top and the inner wall.
fn revolved_ring(sides: usize) -> Vec<Tri> {
    let at = |i: usize, r: f64, z: f64| {
        let a = std::f64::consts::TAU * (i % sides) as f64 / sides as f64;
        Vec3::new(r * a.cos(), r * a.sin(), z)
    };
    let mut tris = Vec::with_capacity(8 * sides);
    for i in 0..sides {
        let (ib0, ob0, ot0, it0) = (
            at(i, 20.0, 0.0),
            at(i, 30.0, 0.0),
            at(i, 30.0, 10.0),
            at(i, 20.0, 10.0),
        );
        let j = i + 1;
        let (ib1, ob1, ot1, it1) = (
            at(j, 20.0, 0.0),
            at(j, 30.0, 0.0),
            at(j, 30.0, 10.0),
            at(j, 20.0, 10.0),
        );
        tris.extend([
            [ib0, ib1, ob1],
            [ib0, ob1, ob0],
            [ob0, ob1, ot1],
            [ob0, ot1, ot0],
            [it0, ot0, ot1],
            [it0, ot1, it1],
            [ib0, it1, ib1],
            [ib0, it0, it1],
        ]);
    }
    tris
}

/// Measurement only: where the scan's time goes on the 20,000-sided ring.
#[test]
#[ignore]
fn measure_revolved_ring_scan() {
    use super::super::exact::Sign;
    use super::super::intersection_graph::{orient3d_plane, tri_box, SelfCutStats};
    let tris = revolved_ring(20_000);
    let boxes: Vec<crate::types::Box> = tris.iter().map(tri_box).collect();
    let scene = boxes
        .iter()
        .fold(crate::types::Box::new(), |acc, b| acc.union_box(b));
    let mut order: Vec<usize> = (0..tris.len()).collect();
    order.sort_by_key(|&i| crate::sort::morton_code(boxes[i].center(), &scene));
    let collider = crate::collider::Collider::new(
        order.iter().map(|&i| boxes[i]).collect(),
        order
            .iter()
            .map(|&i| crate::sort::morton_code(boxes[i].center(), &scene))
            .collect(),
    );
    let mut pairs = Vec::new();
    let mut cands = Vec::new();
    for (i, b) in boxes.iter().enumerate() {
        cands.clear();
        collider.collisions_one(b, i, |_, leaf| cands.push(order[leaf]));
        cands.sort_unstable();
        pairs.extend(cands.iter().filter(|&&j| j > i).map(|&j| (i, j)));
    }
    let mut by_shared = [0usize; 4];
    let mut coplanar_by_shared = [0usize; 4];
    for &(i, j) in &pairs {
        let shared = tris[i].iter().filter(|v| tris[j].contains(v)).count();
        by_shared[shared] += 1;
        if tris[j]
            .iter()
            .all(|&v| orient3d_plane(&tris[i], v) == Sign::Zero)
        {
            coplanar_by_shared[shared] += 1;
        }
    }
    let mut stats = SelfCutStats::default();
    let t = Instant::now();
    let any = pairs
        .iter()
        .filter(|&&(i, j)| super::genuine_contact(tris[i], tris[j], &mut stats))
        .count();
    let narrow = t.elapsed().as_secs_f64();
    let t = Instant::now();
    let verdict = super::scan_triangles(&tris, None, None);
    eprintln!(
        "ring: {} tris, {} box pairs, by shared corners {by_shared:?}, coplanar {coplanar_by_shared:?}; \
         narrow phase {:.1} ms ({:.2} us/pair), {any} contacts, full tests {} ({:.1} ms); scan {:.1} ms ({verdict:?})",
        tris.len(),
        pairs.len(),
        narrow * 1000.0,
        narrow * 1e6 / pairs.len() as f64,
        stats.full,
        stats.full_secs * 1000.0,
        t.elapsed().as_secs_f64() * 1000.0
    );
}

/// Every box pair a one-tree walk lists on `tris`, lower index first.
fn box_pairs(tris: &[Tri]) -> Vec<(usize, usize)> {
    use super::super::intersection_graph::tri_box;
    let boxes: Vec<crate::types::Box> = tris.iter().map(tri_box).collect();
    let mut pairs = Vec::new();
    for i in 0..tris.len() {
        for j in i + 1..tris.len() {
            if boxes[i].does_overlap_box(&boxes[j]) {
                pairs.push((i, j));
            }
        }
    }
    pairs
}

#[test]
fn a_finely_revolved_ring_keeps_its_verdict_and_skips_the_full_test_on_its_flat_faces() {
    use super::super::intersection_graph::{real_self_contact, SelfCutStats};
    let tris = revolved_ring(1_000);
    let pairs = box_pairs(&tris);
    // The verdict of the full test on every pair, the answer to keep.
    let mut full = SelfCutStats::default();
    let reference = pairs
        .iter()
        .any(|&(i, j)| real_self_contact(tris[i], tris[j], &mut full).is_some());
    assert!(!reference, "a closed ring has no self-contact");
    assert_eq!(
        super::scan_triangles(&tris, None, None),
        Some(false),
        "the scan's verdict"
    );
    let mut stats = SelfCutStats::default();
    let contacts = pairs
        .iter()
        .filter(|&&(i, j)| super::genuine_contact(tris[i], tris[j], &mut stats))
        .count();
    assert_eq!(contacts, 0);
    assert!(
        stats.full <= MAX_FULL_TESTS,
        "{} of {} pairs reached the full test ({} apart in projection; {} without the scan's shortcuts)",
        stats.full,
        pairs.len(),
        stats.apart_2d,
        full.full
    );
}

/// The most box pairs of the 1,000-sided ring that may reach the full
/// exact test. Of its 275,850 box pairs, the full test alone would see
/// 242,850; with the scan's shortcuts 10,467 reach it, 221,383 of the rest
/// proven apart in projection.
const MAX_FULL_TESTS: usize = 12_000;

fn v(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z)
}

#[test]
fn the_projection_shortcut_never_moves_a_pair_verdict() {
    use super::super::graph_self_cut::apart_in_projection;
    use super::super::intersection_graph::{real_self_contact, SelfCutStats};
    let base = [v(0.0, 0.0, 0.0), v(4.0, 0.0, 0.0), v(0.0, 4.0, 0.0)];
    // (name, other triangle, whether the shortcut must prove it apart)
    let cases: [(&str, Tri, Option<bool>); 9] = [
        (
            "apart",
            [v(5.0, 5.0, 0.0), v(9.0, 5.0, 0.0), v(5.0, 9.0, 0.0)],
            Some(true),
        ),
        (
            "overlapping",
            [v(1.0, 1.0, 0.0), v(6.0, 1.0, 0.0), v(1.0, 6.0, 0.0)],
            Some(false),
        ),
        (
            "corner on an edge",
            [v(2.0, 2.0, 0.0), v(6.0, 2.0, 0.0), v(6.0, 6.0, 0.0)],
            Some(false),
        ),
        (
            "corner on a corner, not shared",
            [v(4.0, 0.0, 0.0), v(8.0, 0.0, 0.0), v(8.0, 4.0, 0.0)],
            None,
        ),
        (
            "along an edge, no shared corner",
            [v(1.0, 0.0, 0.0), v(3.0, 0.0, 0.0), v(2.0, -3.0, 0.0)],
            Some(false),
        ),
        (
            "along an edge, one shared corner",
            [v(0.0, 0.0, 0.0), v(2.0, 0.0, 0.0), v(1.0, -3.0, 0.0)],
            None,
        ),
        (
            "zero-area sliver apart",
            [v(5.0, 5.0, 0.0), v(6.0, 6.0, 0.0), v(7.0, 7.0, 0.0)],
            Some(true),
        ),
        (
            "zero-area sliver across",
            [v(-1.0, 1.0, 0.0), v(1.0, 1.0, 0.0), v(3.0, 1.0, 0.0)],
            Some(false),
        ),
        (
            "crossing out of plane",
            [v(1.0, 1.0, -1.0), v(1.5, 1.0, 1.0), v(1.0, 1.5, 1.0)],
            Some(false),
        ),
    ];
    for (name, other, apart) in cases {
        for (t1, t2) in [(base, other), (other, base)] {
            let shortcut = apart_in_projection(t1, t2);
            if let Some(expected) = apart {
                assert_eq!(shortcut, expected, "{name}: proven apart");
            }
            if name.starts_with("zero-area") {
                // The scan never hands a degenerate triangle to the narrow
                // phase (its `live` filter), and the full test requires a
                // real one; only the shortcut's own answer is checked.
                continue;
            }
            let full = real_self_contact(t1, t2, &mut SelfCutStats::default()).is_some();
            if shortcut {
                assert!(
                    !full,
                    "{name}: proven apart yet the full test finds contact"
                );
            }
            assert_eq!(
                super::genuine_contact(t1, t2, &mut SelfCutStats::default()),
                full,
                "{name}: the scan's pair verdict"
            );
        }
    }
}
