use super::*;

fn contour(coords: &[(f64, f64)], first_idx: i32) -> Vec<PolyVert> {
    coords
        .iter()
        .enumerate()
        .map(|(offset, &(x, y))| PolyVert {
            pos: Vec2::new(x, y),
            idx: first_idx + offset as i32,
        })
        .collect()
}

#[test]
fn orders_holes_rightmost_first_for_keyholing() {
    let polygons = vec![
        contour(&[(0.0, 0.0), (12.0, 0.0), (12.0, 10.0), (0.0, 10.0)], 0),
        contour(&[(1.0, 2.0), (1.0, 4.0), (3.0, 4.0), (3.0, 2.0)], 4),
        contour(&[(4.0, 2.0), (4.0, 4.0), (6.0, 4.0), (6.0, 2.0)], 8),
        contour(&[(8.0, 2.0), (8.0, 4.0), (10.0, 4.0), (10.0, 2.0)], 12),
    ];

    let ear_clip = EarClip::new(&polygons, 1.0e-10);
    let hole_xs = ear_clip
        .holes
        .iter()
        .map(|&hole| ear_clip.polygon[hole].pos.x)
        .collect::<Vec<_>>();

    assert_eq!(hole_xs, vec![10.0, 6.0, 3.0]);
}

/// End-to-end companion to `orders_holes_rightmost_first_for_keyholing`:
/// with holes keyholed in contour order the bridges cross, producing
/// inverted triangles. Verifies the triangulation itself is valid.
#[test]
fn multi_hole_triangulation_has_no_inverted_triangles() {
    let polygons = vec![
        contour(&[(0.0, 0.0), (12.0, 0.0), (12.0, 10.0), (0.0, 10.0)], 0),
        contour(&[(1.0, 2.0), (1.0, 4.0), (3.0, 4.0), (3.0, 2.0)], 4),
        contour(&[(4.0, 2.0), (4.0, 4.0), (6.0, 4.0), (6.0, 2.0)], 8),
        contour(&[(8.0, 2.0), (8.0, 4.0), (10.0, 4.0), (10.0, 2.0)], 12),
    ];
    let verts: Vec<Vec2> = polygons.iter().flatten().map(|v| v.pos).collect();

    let (triangles, _eps) = EarClip::new(&polygons, -1.0).triangulate();

    // Outer 12x10 rectangle minus three 2x2 holes.
    let expected_area = 12.0 * 10.0 - 3.0 * 4.0;
    let mut total_area = 0.0;
    for tri in &triangles {
        let (a, b, c) = (
            verts[tri.x as usize],
            verts[tri.y as usize],
            verts[tri.z as usize],
        );
        let area = 0.5 * determinant2x2(b - a, c - a);
        assert!(
            area > 0.0,
            "inverted or degenerate triangle {:?} (area {})",
            tri,
            area
        );
        total_area += area;
    }
    assert!(
        (total_area - expected_area).abs() < 1e-9,
        "triangulation area {} != expected {}",
        total_area,
        expected_area
    );
}

/// FNV-1a over a triangle list, pinning the triangles and their order.
fn fnv(tris: &[IVec3Out]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for t in tris {
        for c in [t.x, t.y, t.z] {
            for b in c.to_le_bytes() {
                h ^= u64::from(b);
                h = h.wrapping_mul(0x0100_0000_01b3);
            }
        }
    }
    h
}

fn octagon(cx: f64, cy: f64, r: f64, hole: bool) -> Vec<Vec2> {
    let mut ring: Vec<Vec2> = [
        (1.0, -0.5),
        (0.5, -1.0),
        (-0.5, -1.0),
        (-1.0, -0.5),
        (-1.0, 0.5),
        (-0.5, 1.0),
        (0.5, 1.0),
        (1.0, 0.5),
    ]
    .iter()
    .map(|&(x, y)| Vec2::new(cx + r * x, cy + r * y))
    .collect();
    if !hole {
        ring.reverse();
    }
    ring
}

/// Pins the triangles of a 24x24 grid of octagonal holes, every other row
/// shifted so the bridges run between holes rather than to the outer edge.
#[test]
fn keyholing_many_holes_keeps_its_triangles() {
    let n = 24;
    let size = 3.0 * f64::from(n) + 3.0;
    let mut polys = vec![vec![
        Vec2::new(0.0, 0.0),
        Vec2::new(size, 0.0),
        Vec2::new(size, size),
        Vec2::new(0.0, size),
    ]];
    for i in 0..n {
        for j in 0..n {
            let cx = 3.0 * f64::from(i) + 2.0 + if j % 2 == 1 { 0.75 } else { 0.0 };
            let cy = 3.0 * f64::from(j) + 2.0;
            polys.push(octagon(cx, cy, 1.0, true));
        }
    }
    let tris = crate::polygon::triangulate(&polys, 1e-9, true);
    // Every vert is kept, so v + 2h - 2 triangles.
    let verts = 4 + 8 * n * n;
    assert_eq!(tris.len() as i32, verts + 2 * n * n - 2);
    assert_eq!(fnv(&tris), 0xb444_9c9b_cd61_b83e, "hash {:#x}", fnv(&tris));
}

/// A hole that retraces its outer ring (ring 2 is ring 1, the unit square,
/// drawn the other way) collapses that ring to two verts when it is joined
/// in, and the ring stays in `outers`. The second hole's bridge searches then
/// both walk the collapsed ring: the walk reports it degenerate, and each
/// search skips it and restores the connector it had before the ring.
///
/// Pins the triangles, taken on `main` before the searches walked rings in
/// place, and checks the state that puts this input on that path: after the
/// first cut ring 1 is degenerate, and the walk says so before visiting any
/// vert, so the restore has nothing to undo.
#[test]
fn keyholing_skips_an_outer_ring_collapsed_by_an_earlier_hole() {
    let rings: [&[(f64, f64)]; 4] = [
        &[(-3.0, -1.0), (3.0, -1.0), (3.0, 2.0), (-3.0, 2.0)],
        &[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        &[(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0)],
        &[(-2.0, 0.25), (-2.0, 0.75), (-1.0, 0.75), (-1.0, 0.25)],
    ];
    let polys: Vec<Vec<Vec2>> = rings
        .iter()
        .map(|r| r.iter().map(|&(x, y)| Vec2::new(x, y)).collect())
        .collect();
    let tris: Vec<[i32; 3]> = crate::polygon::triangulate(&polys, 1e-9, true)
        .iter()
        .map(|t| [t.x, t.y, t.z])
        .collect();
    let expected = [
        [9, 10, 6],
        [9, 6, 7],
        [8, 9, 7],
        [8, 7, 4],
        [11, 8, 4],
        [11, 4, 5],
        [10, 11, 5],
        [10, 5, 6],
        [1, 2, 14],
        [1, 14, 15],
        [14, 2, 3],
        [13, 14, 3],
        [13, 3, 0],
        [12, 13, 0],
        [12, 0, 1],
        [12, 1, 15],
    ];
    assert_eq!(tris, expected);

    let polygons: Vec<Vec<PolyVert>> = rings
        .iter()
        .scan(0, |first, r| {
            let c = contour(r, *first);
            *first += r.len() as i32;
            Some(c)
        })
        .collect();
    let mut ear_clip = EarClip::new(&polygons, 1e-9);
    assert_eq!(ear_clip.outers.len(), 2);
    assert_eq!(ear_clip.holes.len(), 2);
    let first_hole = ear_clip.holes[0];
    ear_clip.cut_keyhole(first_hole);
    let collapsed = ear_clip.outers[1];
    let mut visited = Vec::new();
    assert!(!ear_clip.for_each_loop_vert(collapsed, |v| visited.push(v)));
    assert!(visited.is_empty(), "visited {visited:?}");
}
