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
    revolved_rectangle(sides, 20.0, 30.0)
}

/// A rectangle from radius `inner` to `outer`, height 10, revolved through
/// `sides` sides, wound outward.
fn revolved_rectangle(sides: usize, inner: f64, outer: f64) -> Vec<Tri> {
    let at = |i: usize, r: f64, z: f64| {
        let a = std::f64::consts::TAU * (i % sides) as f64 / sides as f64;
        Vec3::new(r * a.cos(), r * a.sin(), z)
    };
    let mut tris = Vec::with_capacity(8 * sides);
    for i in 0..sides {
        let (ib0, ob0, ot0, it0) = (
            at(i, inner, 0.0),
            at(i, outer, 0.0),
            at(i, outer, 10.0),
            at(i, inner, 10.0),
        );
        let j = i + 1;
        let (ib1, ob1, ot1, it1) = (
            at(j, inner, 0.0),
            at(j, outer, 0.0),
            at(j, outer, 10.0),
            at(j, inner, 10.0),
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

/// Measurement only: the scan on ordinary meshes, where the pre-tests must
/// cost nothing noticeable: a 204,800-triangle sphere, and a 160,000-triangle
/// hollow cylinder with a 0.1 wall. Best of three runs each.
#[test]
#[ignore]
fn measure_ordinary_mesh_scans() {
    let best = |f: &dyn Fn() -> bool| {
        (0..3)
            .map(|_| {
                let t = Instant::now();
                assert!(!f());
                t.elapsed().as_secs_f64() * 1000.0
            })
            .fold(f64::INFINITY, f64::min)
    };
    // `has_self_intersections` caches its answer, so the scan itself is
    // timed, on the sphere's own triangles.
    let mesh = crate::manifold::Manifold::sphere(10.0, 640).get_mesh_gl64(-1);
    let corner = |k: u64| {
        let k = 3 * k as usize;
        let p = &mesh.vert_properties;
        Vec3::new(p[k], p[k + 1], p[k + 2])
    };
    let sphere: Vec<Tri> = mesh
        .tri_verts
        .chunks(3)
        .map(|t| [corner(t[0]), corner(t[1]), corner(t[2])])
        .collect();
    let sphere_ms = best(&|| super::scan_triangles(&sphere, None, None) == Some(true));
    let wall = revolved_rectangle(20_000, 20.0, 20.1);
    let wall_ms = best(&|| super::scan_triangles(&wall, None, None) == Some(true));
    eprintln!("sphere {sphere_ms:.1} ms; thin wall {wall_ms:.1} ms");
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
    let tris = revolved_ring(500);
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
        "{} of {} pairs reached the full test ({} proven apart first; {} without the pre-tests)",
        stats.full,
        pairs.len(),
        stats.apart,
        full.full
    );
}

/// The most box pairs of the 500-sided ring that may reach the full exact
/// test. Of its 82,114 box pairs, the full test alone would see 65,614;
/// with the scan's pre-tests 5,218 reach it, 54,896 of the rest proven
/// apart before it.
const MAX_FULL_TESTS: usize = 6_000;

fn v(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z)
}

#[test]
fn the_projection_shortcut_never_moves_a_pair_verdict() {
    use super::super::intersection_graph::{real_self_contact, SelfCutStats};
    use super::super::tri_tri::{coplanar_separated_2d, f64_drop_axis};
    let base = [v(0.0, 0.0, 0.0), v(4.0, 0.0, 0.0), v(0.0, 4.0, 0.0)];
    // (name, other triangle, whether the shortcut must prove it apart)
    // Square to the base, its shadow on the base's plane a segment beyond
    // the base's long edge: proven apart through that segment.
    let perpendicular = [v(6.0, 0.0, -1.0), v(6.0, 1.0, 1.0), v(6.0, 0.0, 1.0)];
    assert!(coplanar_separated_2d(
        base,
        perpendicular,
        f64_drop_axis(base)
    ));
    let cases: [(&str, Tri, Option<bool>); 10] = [
        ("perpendicular, apart", perpendicular, None),
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
            "one shared corner, apart beyond it",
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
            let shortcut = coplanar_separated_2d(t1, t2, f64_drop_axis(t1));
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

/// Points at subnormal-product scale (coordinates near 2⁻³⁴⁶, so products of
/// three land below the normal range): `b` and `c` with the origin span a
/// plane; halving or eighthing them stays exactly in it.
fn subnormal_plane() -> (Vec3, Vec3) {
    let b = Vec3::new(
        f64::from_bits(0x2a50549ae3dc2d94),
        f64::from_bits(0x2a51ad2109517d2e),
        f64::from_bits(0x2a5a02929b27eae0),
    );
    let c = Vec3::new(
        f64::from_bits(0x2a450fc229ea9000),
        f64::from_bits(0x2a5cf0ce06c472c2),
        f64::from_bits(0xaa3149aee586a2a8),
    );
    (b, c)
}

/// The full exact test's verdict on one pair: the answer the scan must give.
fn full_verdict(t1: Tri, t2: Tri) -> bool {
    use super::super::intersection_graph::{real_self_contact, SelfCutStats};
    real_self_contact(t1, t2, &mut SelfCutStats::default()).is_some()
}

#[test]
fn a_coplanar_overlap_at_subnormal_scale_is_still_found() {
    let (b, c) = subnormal_plane();
    let o = Vec3::new(0.0, 0.0, 0.0);
    // `t2` lies inside `t1`, sharing no corner.
    let t1 = [o, b, c];
    let t2 = [b / 8.0, c / 2.0, c / 8.0];
    assert!(full_verdict(t1, t2), "the full test finds the overlap");
    assert_eq!(super::scan_triangles(&[t1, t2], None, None), Some(true));
    // Sharing the origin, one wedge inside the other: the float filter
    // without its underflow guard certified both of the small triangle's
    // other corners strictly on one side of the large one's plane (a
    // spurious, same sign for both), so the one-shared-vertex shortcut
    // called the overlap a point contact.
    let small = [o, b / 256.0, c / 256.0];
    let large = [o, b * 4.0, c * 8.0];
    assert!(
        full_verdict(small, large),
        "the full test finds the overlap"
    );
    assert_eq!(
        super::scan_triangles(&[small, large], None, None),
        Some(true)
    );
}
