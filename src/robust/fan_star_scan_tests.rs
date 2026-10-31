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

// Tests for the self-intersection scan on triangle fans (robust/soup.rs
// `has_self_intersections`, whose flat-star settling is
// robust/fan_star.rs). A cylinder capped by two fans around a centre vertex
// is the everyday shape here: every pair of cap triangles shares the centre,
// so every pair of their boxes overlaps, and a scan that walks those pairs
// one by one is quadratic in the cap's side count. The verdicts pinned are
// the ones the pair-by-pair scan gives, including the fans whose sectors
// overlap, touch along a ray without sharing an edge, or wind round more
// than once.

use std::time::Instant;

use crate::manifold::Manifold;
use crate::types::MeshGL64;

/// A closed cylinder of `sides` sides about z, radius 10, height 20, each cap
/// a fan of `sides` triangles around its own centre vertex. `cap_turns`
/// winds the top cap's ring that many times round (1 is an ordinary cap;
/// 2 makes the top cap's fan cover its disc twice over, which with the walls
/// following the same ring still closes but genuinely overlaps itself).
fn fan_capped_cylinder(sides: usize, cap_turns: usize) -> MeshGL64 {
    let mut mesh = MeshGL64 {
        num_prop: 3,
        ..Default::default()
    };
    let mut push = |x: f64, y: f64, z: f64| {
        mesh.vert_properties.extend([x, y, z]);
    };
    // 0: bottom centre, 1: top centre, then the bottom ring, then the top ring.
    push(0.0, 0.0, 0.0);
    push(0.0, 0.0, 20.0);
    for z in [0.0, 20.0] {
        for i in 0..sides {
            let a = std::f64::consts::TAU * (cap_turns * i) as f64 / sides as f64;
            push(10.0 * a.cos(), 10.0 * a.sin(), z);
        }
    }
    let bottom = |i: usize| (2 + i % sides) as u64;
    let top = |i: usize| (2 + sides + i % sides) as u64;
    for i in 0..sides {
        mesh.tri_verts.extend([0, bottom(i + 1), bottom(i)]);
        mesh.tri_verts.extend([1, top(i), top(i + 1)]);
        mesh.tri_verts
            .extend([bottom(i), bottom(i + 1), top(i + 1)]);
        mesh.tri_verts.extend([bottom(i), top(i + 1), top(i)]);
    }
    mesh
}

/// Imported strictly, as a manifold mesh is (the app's path for a closed part).
fn import(mesh: &MeshGL64) -> Manifold {
    Manifold::from_mesh_gl64(mesh)
}

/// Measurement only: the scan's time on a heavy fan-capped cylinder.
#[test]
#[ignore]
fn measure_fan_capped_cylinder_scan() {
    for sides in [1_000usize, 4_000, 30_000] {
        let m = import(&fan_capped_cylinder(sides, 1));
        let t = Instant::now();
        let verdict = m.has_self_intersections();
        eprintln!(
            "sides {sides}: tris {} verdict {verdict} scan {:.1} ms",
            4 * sides,
            t.elapsed().as_secs_f64() * 1000.0
        );
        // The broad phase alone: every box pair a morton-ordered tree over
        // the triangles lists, as the scan's rows walk them.
        let (tris, boxes, live) = mesh_soup(&fan_capped_cylinder(sides, 1));
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
        let t = Instant::now();
        let stars = super::super::fan_star::flat_stars(&tris, &live);
        eprintln!(
            "sides {sides}: flat stars {} found in {:.1} ms; scan of the raw triangles {:.1} ms",
            stars.owned.len(),
            t.elapsed().as_secs_f64() * 1000.0,
            {
                let t = Instant::now();
                super::scan_triangles(&tris, None, None);
                t.elapsed().as_secs_f64() * 1000.0
            }
        );
        eprintln!(
            "sides {sides}: the flat-star walk lists {} pairs for {} triangles",
            super::super::fan_star_scan::listed_pairs(&boxes, &live, &stars),
            tris.len()
        );
        if sides > 4_000 {
            // Listing this one's pairs takes about a minute even optimized.
            continue;
        }
        let t = Instant::now();
        let mut pairs = 0usize;
        let mut cands = Vec::new();
        for (i, b) in boxes.iter().enumerate() {
            cands.clear();
            collider.collisions_one(b, i, |_, leaf| cands.push(order[leaf]));
            cands.sort_unstable();
            pairs += cands.iter().filter(|&&j| j > i).count();
        }
        eprintln!(
            "sides {sides}: broad phase lists {pairs} pairs in {:.1} ms (serial)",
            t.elapsed().as_secs_f64() * 1000.0
        );
    }
}

/// Measurement only: what looking for flat stars adds to the scan of a
/// mesh that has none (a fine sphere, every vertex of degree six or so).
#[test]
#[ignore]
fn measure_flat_star_search_on_a_sphere() {
    let sphere = Manifold::sphere(10.0, 640);
    let (tris, _, live) = mesh_soup(&sphere.get_mesh_gl64(-1));
    let t = Instant::now();
    let stars = super::super::fan_star::flat_stars(&tris, &live);
    let search = t.elapsed().as_secs_f64() * 1000.0;
    let t = Instant::now();
    let verdict = super::scan_triangles(&tris, None, None);
    eprintln!(
        "sphere: {} triangles, {} flat stars found in {search:.1} ms; scan {:.1} ms ({verdict:?})",
        tris.len(),
        stars.owned.len(),
        t.elapsed().as_secs_f64() * 1000.0
    );
}

type Tri = [crate::linalg::Vec3; 3];

fn v(x: f64, y: f64, z: f64) -> crate::linalg::Vec3 {
    crate::linalg::Vec3::new(x, y, z)
}

/// The verdict of testing every overlapping pair one by one, with the scan's
/// own narrow phase and argument order: the answer the scan must keep.
fn pair_by_pair(tris: &[Tri]) -> bool {
    use super::super::intersection_graph::{is_degenerate, tri_box, SelfCutStats};
    let boxes: Vec<crate::types::Box> = tris.iter().map(tri_box).collect();
    let live: Vec<bool> = tris.iter().map(|t| !is_degenerate(t)).collect();
    let mut stats = SelfCutStats::default();
    (0..tris.len()).any(|i| {
        live[i]
            && (i + 1..tris.len()).any(|j| {
                live[j]
                    && boxes[i].does_overlap_box(&boxes[j])
                    && super::genuine_contact(tris[i], tris[j], &mut stats)
            })
    })
}

/// The scan's verdict on plain triangles, with no face tree.
fn scan(tris: &[Tri]) -> bool {
    super::scan_triangles(tris, None, None).expect("not cancelled")
}

/// The settled stars of `tris`: each star's owned triangles.
fn settled_stars(tris: &[Tri]) -> Vec<Vec<usize>> {
    use super::super::intersection_graph::is_degenerate;
    let live: Vec<bool> = tris.iter().map(|t| !is_degenerate(t)).collect();
    super::super::fan_star::flat_stars(tris, &live).owned
}

/// Whether the flat-star path really runs on `tris` (anti-vacuity).
fn has_a_flat_star(tris: &[Tri]) -> bool {
    !settled_stars(tris).is_empty()
}

/// The triangles of `mesh` as the scan reads them, with their boxes and
/// which ones are live (not degenerate).
fn mesh_soup(mesh: &MeshGL64) -> (Vec<Tri>, Vec<crate::types::Box>, Vec<bool>) {
    use super::super::intersection_graph::{is_degenerate, tri_box};
    let p = |k: u64| {
        let k = 3 * k as usize;
        v(
            mesh.vert_properties[k],
            mesh.vert_properties[k + 1],
            mesh.vert_properties[k + 2],
        )
    };
    let tris: Vec<Tri> = mesh
        .tri_verts
        .chunks(3)
        .map(|t| [p(t[0]), p(t[1]), p(t[2])])
        .collect();
    let boxes = tris.iter().map(tri_box).collect();
    let live = tris.iter().map(|t| !is_degenerate(t)).collect();
    (tris, boxes, live)
}

/// `n` triangles round the origin in the z = 0 plane, the k-th spanning ring
/// points `k * step` to `k * step + span` of a ring of `ring` points (radius
/// 10, the j-th point at angle j / ring of a turn, `radius_of(j)` scaling it).
fn fan(
    n: usize,
    ring: usize,
    step: usize,
    span: usize,
    radius_of: impl Fn(usize) -> f64,
) -> Vec<Tri> {
    let point = |j: usize| {
        let a = std::f64::consts::TAU * (j % ring) as f64 / ring as f64;
        let r = 10.0 * radius_of(j % ring);
        v(r * a.cos(), r * a.sin(), 0.0)
    };
    (0..n)
        .map(|k| [v(0.0, 0.0, 0.0), point(k * step), point(k * step + span)])
        .collect()
}

/// Asserts the scan keeps the pair-by-pair verdict, and that verdict is `expected`.
fn assert_verdict(name: &str, tris: &[Tri], expected: bool) {
    assert_eq!(
        pair_by_pair(tris),
        expected,
        "{name}: the pair-by-pair verdict"
    );
    assert_eq!(scan(tris), expected, "{name}: the scan's verdict");
}

#[test]
fn a_flat_disjoint_fan_is_settled_and_clean() {
    let tris = fan(300, 300, 1, 1, |_| 1.0);
    assert!(has_a_flat_star(&tris));
    assert_verdict("closed fan", &tris, false);
    // An open fan (a sector missing) is disjoint too.
    let open = fan(250, 300, 1, 1, |_| 1.0);
    assert!(has_a_flat_star(&open));
    assert_verdict("open fan", &open, false);
}

#[test]
fn an_overlapping_flat_fan_keeps_its_contact() {
    // Every triangle spans two ring steps but starts one step on: each
    // overlaps its neighbours over a whole sector.
    let tris = fan(300, 300, 1, 2, |_| 1.0);
    assert!(
        !has_a_flat_star(&tris),
        "an overlapping star is never settled"
    );
    assert_verdict("double-width fan", &tris, true);
    // One wide triangle among an otherwise disjoint fan.
    let mut one_wide = fan(300, 300, 1, 1, |_| 1.0);
    one_wide[7] = [v(0.0, 0.0, 0.0), one_wide[7][1], one_wide[9][2]];
    assert_verdict("one wide sector", &one_wide, true);
}

#[test]
fn a_fan_winding_round_twice_keeps_its_contact() {
    // 301 points walked two at a time: the fan goes round twice and every
    // sector lies over another.
    let tris = fan(301, 301, 2, 2, |_| 1.0);
    assert_verdict("fan wound twice", &tris, true);
}

#[test]
fn sectors_touching_along_a_ray_without_a_shared_corner_keep_their_verdict() {
    // Neighbouring triangles share the ray to every other ring point but not
    // the corner on it: odd triangles end at half the radius there.
    let mut tris = fan(300, 300, 1, 1, |_| 1.0);
    for (k, tri) in tris.iter_mut().enumerate() {
        if k % 2 == 1 {
            tri[1] = tri[1] * 0.5;
        }
    }
    assert!(
        !has_a_flat_star(&tris),
        "a T along a ray is left to the pair scan"
    );
    assert_verdict("T-junction rays", &tris, true);
}

#[test]
fn triangles_outside_a_settled_star_are_still_tested_against_it() {
    let base = fan(300, 300, 1, 1, |_| 1.0);
    // A triangle standing up through the cap, crossing it.
    let mut pierced = base.clone();
    pierced.push([v(1.0, 1.0, -5.0), v(2.0, 1.5, 5.0), v(1.5, 2.5, 5.0)]);
    assert!(has_a_flat_star(&pierced));
    assert_verdict("pierced cap", &pierced, true);
    // One resting a corner on the cap from below: a point contact.
    let mut touching = base.clone();
    touching.push([v(2.0, 1.0, 0.0), v(1.0, 1.0, -5.0), v(3.0, 1.5, -5.0)]);
    assert!(has_a_flat_star(&touching));
    assert_verdict("corner resting on the cap", &touching, false);
    // A second fan over the first, rotated half a step round the same
    // centre: the two make one 600-triangle star whose wedges overlap, so it
    // is never settled and the pair scan finds the overlap.
    let mut doubled = base.clone();
    let half_step = |j: usize| {
        let a = std::f64::consts::TAU * (j as f64 + 0.5) / 300.0;
        v(10.0 * a.cos(), 10.0 * a.sin(), 0.0)
    };
    doubled.extend((0..300).map(|k| [v(0.0, 0.0, 0.0), half_step(k), half_step(k + 1)]));
    assert!(
        !has_a_flat_star(&doubled),
        "one shared centre makes one overlapping star"
    );
    assert_verdict("two fans over each other", &doubled, true);
    // Two disjoint fans side by side in one plane.
    let mut beside = base;
    beside.extend(
        fan(300, 300, 1, 1, |_| 1.0)
            .iter()
            .map(|t| t.map(|p| p + v(25.0, 0.0, 0.0))),
    );
    assert!(has_a_flat_star(&beside));
    assert_verdict("two fans side by side", &beside, false);
}

#[test]
fn a_fan_capped_cylinder_is_clean_and_quick() {
    let small = import(&fan_capped_cylinder(300, 1));
    assert!(!small.has_self_intersections(), "300 sides");
    // Pair-by-pair this lists about 450 million box pairs (tens of seconds
    // even in an optimized parallel build); settled as stars it is a sort.
    let mesh = fan_capped_cylinder(30_000, 1);
    assert!(!import(&mesh).has_self_intersections(), "30,000 sides");
    // The walk's cost counted in pairs, not seconds: both caps settled whole,
    // and every row listing a few partners, not its cap's 30,000.
    let (tris, boxes, live) = mesh_soup(&mesh);
    let stars = super::super::fan_star::flat_stars(&tris, &live);
    let sizes: Vec<usize> = stars.owned.iter().map(Vec::len).collect();
    assert_eq!(sizes, [30_000, 30_000], "both caps settled whole");
    assert!(stars.grouped, "the bucket count lets the caps through");
    let listed = super::super::fan_star_scan::listed_pairs(&boxes, &live, &stars);
    assert!(
        listed <= LISTED_PAIRS_PER_TRIANGLE * tris.len(),
        "the walk listed {listed} pairs for {} triangles",
        tris.len()
    );
}

#[test]
fn a_mesh_without_fans_skips_the_corner_grouping() {
    // A sphere's vertices each touch about six triangles: the bucket count
    // proves there is no star before any sorting.
    let (tris, _, live) = mesh_soup(&Manifold::sphere(10.0, 64).get_mesh_gl64(-1));
    let stars = super::super::fan_star::flat_stars(&tris, &live);
    assert!(!stars.grouped, "no position reaches a star's size");
    assert!(stars.owned.is_empty() && stars.owner.iter().all(Option::is_none));
    assert_verdict("sphere", &tris, false);
    // One flat fan of exactly the threshold's size is still grouped.
    let tris = fan(64, 64, 1, 1, |_| 1.0);
    let live = vec![true; tris.len()];
    let stars = super::super::fan_star::flat_stars(&tris, &live);
    assert!(stars.grouped, "a 64-triangle star passes the bucket count");
    assert_eq!(stars.owned.len(), 1);
}

/// The most stars one row of the fan grid may reach through the star tree,
/// on average. Measured: 21,964 for 16,384 triangles (1.34 each), the row's
/// own star plus the diagonal neighbours its box overlaps.
const STAR_TREE_HITS_PER_TRIANGLE: usize = 2;

/// The most pairs per triangle the flat-star walk may list on a
/// fan-capped cylinder. The 30,000-sided one lists 510,020 pairs for
/// 120,000 triangles (4.25 each, mostly the walls' neighbours); listing a
/// cap pair by pair would be about 3,750 each.
const LISTED_PAIRS_PER_TRIANGLE: usize = 5;

#[test]
fn overlapping_settled_stars_are_tested_against_each_other() {
    // Two disjoint fans, each settled, the second's centre inside the first.
    let mut tris = fan(300, 300, 1, 1, |_| 1.0);
    tris.extend(
        fan(300, 300, 1, 1, |_| 1.0)
            .iter()
            .map(|t| t.map(|p| p + v(5.0, 0.0, 0.0))),
    );
    assert_eq!(settled_stars(&tris).len(), 2, "both fans settled");
    assert_verdict("two settled fans overlapping", &tris, true);
}

#[test]
fn a_triangle_in_two_stars_is_owned_once_and_keeps_its_verdict() {
    // Fan A round the origin with a gap about +x, fan B the same fan turned
    // half a turn about (12.5, 0, 0), so its gap faces A's; one triangle
    // bridges the two centres through both gaps, a member of both stars.
    let point = |j: usize| {
        let a = std::f64::consts::TAU * j as f64 / 300.0;
        v(10.0 * a.cos(), 10.0 * a.sin(), 0.0)
    };
    let fan_a: Vec<Tri> = (5..295)
        .map(|k| [v(0.0, 0.0, 0.0), point(k), point(k + 1)])
        .collect();
    let turn = |p: crate::linalg::Vec3| v(25.0 - p.x, -p.y, p.z);
    let mut tris = fan_a.clone();
    tris.extend(fan_a.iter().map(|t| t.map(turn)));
    tris.push([v(0.0, 0.0, 0.0), v(25.0, 0.0, 0.0), v(12.5, 1.0, 0.0)]);
    let stars = settled_stars(&tris);
    assert_eq!(stars.len(), 2, "both fans settled with the bridge");
    let mut owned: Vec<usize> = stars.concat();
    owned.sort_unstable();
    assert_eq!(
        owned,
        (0..tris.len()).collect::<Vec<_>>(),
        "every triangle owned exactly once"
    );
    assert_verdict("bridge between two fans", &tris, false);
}

/// A 16 by 16 grid of closed 64-triangle fans in staggered rows, placed so
/// diagonal neighbours' boxes overlap while their discs stay apart.
fn staggered_fan_grid() -> Vec<Vec<Tri>> {
    let mut fans = Vec::new();
    for row in 0..16 {
        for col in 0..16 {
            let x = 22.0 * col as f64 + if row % 2 == 1 { 11.0 } else { 0.0 };
            let y = 18.0 * row as f64;
            fans.push(
                fan(64, 64, 1, 1, |_| 1.0)
                    .iter()
                    .map(|t| t.map(|p| p + v(x, y, 0.0)))
                    .collect(),
            );
        }
    }
    fans
}

#[test]
fn many_settled_stars_with_overlapping_boxes_keep_their_verdict() {
    let fans = staggered_fan_grid();
    let tris: Vec<Tri> = fans.concat();
    assert_eq!(settled_stars(&tris).len(), 256, "every fan settled");
    assert_verdict("staggered fan grid", &tris, false);
    // Each row reaches its own star and its few diagonal neighbours, not
    // all 256: a loop over every star would be 16,384 x 256 = 4,194,304.
    let live = vec![true; tris.len()];
    let stars = super::super::fan_star::flat_stars(&tris, &live);
    let boxes: Vec<crate::types::Box> = tris
        .iter()
        .map(super::super::intersection_graph::tri_box)
        .collect();
    let hits = super::super::fan_star_scan::star_tree_hits(&boxes, &live, &stars);
    assert!(
        hits <= STAR_TREE_HITS_PER_TRIANGLE * tris.len(),
        "the rows reached {hits} stars for {} triangles",
        tris.len()
    );
    // One fan moved halfway onto its neighbour.
    let mut moved = fans;
    moved[1] = moved[0]
        .iter()
        .map(|t| t.map(|p| p + v(5.0, 0.0, 0.0)))
        .collect();
    let tris: Vec<Tri> = moved.concat();
    assert_eq!(settled_stars(&tris).len(), 256, "every fan settled");
    assert_verdict("one fan moved onto its neighbour", &tris, true);
}

/// A fan of `n` triangles in the tilted plane z = -(x + y), on integer
/// points so it is exactly coplanar off every axis, each spanning `span`
/// ring steps.
fn tilted_fan(n: usize, span: usize) -> Vec<Tri> {
    let point = |j: usize| {
        let a = std::f64::consts::TAU * (j % n) as f64 / n as f64;
        let (x, y) = ((1000.0 * a.cos()).round(), (1000.0 * a.sin()).round());
        v(x, y, -(x + y))
    };
    (0..n)
        .map(|k| [v(0.0, 0.0, 0.0), point(k), point(k + span)])
        .collect()
}

#[test]
fn a_tilted_flat_fan_is_settled_by_the_exact_plane_test() {
    let tris = tilted_fan(128, 1);
    assert!(has_a_flat_star(&tris), "exactly coplanar off the axes");
    assert_verdict("tilted fan", &tris, false);
    let overlapping = tilted_fan(128, 2);
    assert!(!has_a_flat_star(&overlapping));
    assert_verdict("tilted double-width fan", &overlapping, true);
}
