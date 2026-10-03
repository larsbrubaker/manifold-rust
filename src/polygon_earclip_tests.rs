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

/// Pins the triangles of offset rows of rings with holes and islands, where
/// each hole sees several candidate rings, at a fixed and the automatic epsilon.
#[test]
fn keyholing_many_outer_rings_keeps_its_triangles() {
    let mut polys = Vec::new();
    for row in 0..9 {
        for col in 0..12 {
            let x0 = 5.0 * f64::from(col) + if row % 2 == 1 { 1.25 } else { 0.0 };
            let y0 = 3.75 * f64::from(row) + 0.25 * f64::from(col % 3);
            let (w, h) = (4.0, 4.5);
            polys.push(vec![
                Vec2::new(x0, y0),
                Vec2::new(x0 + w, y0),
                Vec2::new(x0 + w, y0 + h),
                Vec2::new(x0, y0 + h),
            ]);
            if (row + col) % 3 == 0 {
                polys.push(octagon(x0 + 2.0, y0 + 1.25, 0.75, true));
                polys.push(octagon(x0 + 2.0, y0 + 3.25, 0.75, true));
            } else {
                polys.push(octagon(x0 + 2.0, y0 + 2.25, 1.5, true));
                if col % 2 == 0 {
                    polys.push(octagon(x0 + 2.0, y0 + 2.25, 0.5, false));
                }
            }
        }
    }
    // The two epsilons happen to pick the same triangles here.
    for epsilon in [1e-9, -1.0] {
        let tris = crate::polygon::triangulate(&polys, epsilon, true);
        assert_eq!(tris.len(), 1872);
        assert_eq!(fnv(&tris), 0x57ad_354c_2a93_a394, "hash {:#x}", fnv(&tris));
    }
}

fn tri_list(tris: &[IVec3Out]) -> Vec<(i32, i32, i32)> {
    tris.iter().map(|t| (t.x, t.y, t.z)).collect()
}

/// `find_closer_bridge` skips a ring whose bounding box lies outside the
/// line from start to connector, but `ccw` squares the area, and at L =
/// 2^-270 the square of the L-shaped ring's reflex corner's area (-1.26 L^2)
/// underflows to 0, so `ccw` calls it collinear and the tie-break takes it
/// as the bridge. The cull must not apply at that scale. Without the
/// magnitude window the hole bridges to the box corner (L, 2L) instead,
/// giving 11 different triangles. The expected list is the output before
/// the ring culls existed.
#[test]
fn keyhole_cull_keeps_a_bridge_whose_ccw_underflows() {
    let l = 2f64.powi(-270);
    let p = |x: f64, y: f64| Vec2::new(x * l, y * l);
    let polys = vec![
        // Outer box; its right edge, x = L, is crossed by the hole's ray and
        // ends at the first connector (L, 2L).
        vec![
            p(-3000.0, -3000.0),
            p(1.0, -3000.0),
            p(1.0, 2.0),
            p(-3000.0, 2.0),
        ],
        // Hole (clockwise); its rightmost reflex vert, the start, is (0, 0).
        vec![p(0.0, 0.0), p(-1000.0, -1000.0), p(-1000.0, 1.0)],
        // Second outer ring, an L whose reflex inner corner (.12L, 1.5L) lies
        // left of start -> connector. Its long edges keep that corner's own
        // reflex test clear of underflow.
        vec![
            p(-1000.0, 1.2),
            p(0.4, 1.2),
            p(0.4, 1.5),
            p(0.12, 1.5),
            p(0.12, 1000.0),
            p(-1000.0, 1000.0),
        ],
    ];
    let tris = crate::polygon::triangulate(&polys, 0.0, true);
    assert_eq!(
        tri_list(&tris),
        [
            (4, 9, 10),
            (8, 9, 4),
            (3, 0, 1),
            (3, 1, 2),
            (10, 11, 12),
            (10, 12, 7),
            (10, 7, 8),
            (4, 10, 8),
            (6, 4, 5),
        ]
    );
}

/// The overflow side of `keyhole_cull_keeps_a_bridge_whose_ccw_underflows`.
/// At L = 2^254 with epsilon 0.01 L, both sides of `ccw`'s comparison for
/// the Γ-shaped ring's reflex corner (.12L, 1.5L) overflow to infinity, so
/// it is called collinear and taken as the bridge, though its box lies
/// outside the line from start to the connector (1000L, 2000L). Without the
/// magnitude window the hole bridges to that connector instead. The
/// expected list is the output before the ring culls existed.
#[test]
fn keyhole_cull_keeps_a_bridge_whose_ccw_overflows() {
    let l = 2f64.powi(254);
    let p = |x: f64, y: f64| Vec2::new(x * l, y * l);
    let polys = vec![
        vec![
            p(-3000.0, -3000.0),
            p(1000.0, -3000.0),
            p(1000.0, 2000.0),
            p(-3000.0, 2000.0),
        ],
        vec![p(0.0, 0.0), p(-1000.0, -1000.0), p(-1000.0, 1.0)],
        vec![
            p(0.02, 1.2),
            p(0.12, 1.2),
            p(0.12, 1.5),
            p(0.4, 1.5),
            p(0.4, 1.8),
            p(0.02, 1.8),
        ],
    ];
    let tris = crate::polygon::triangulate(&polys, 0.01 * l, true);
    assert_eq!(
        tri_list(&tris),
        [
            (3, 0, 1),
            (1, 2, 3),
            (9, 10, 11),
            (12, 7, 8),
            (8, 9, 4),
            (4, 5, 6),
            (9, 11, 12),
            (8, 4, 6),
            (4, 9, 12),
            (12, 8, 6),
            (6, 4, 12),
        ]
    );
}
