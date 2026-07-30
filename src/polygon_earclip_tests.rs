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
        contour(
            &[(8.0, 2.0), (8.0, 4.0), (10.0, 4.0), (10.0, 2.0)],
            12,
        ),
    ];

    let ear_clip = EarClip::new(&polygons, 1.0e-10);
    let hole_xs = ear_clip
        .holes
        .iter()
        .map(|&hole| ear_clip.polygon[hole].pos.x)
        .collect::<Vec<_>>();

    assert_eq!(hole_xs, vec![10.0, 6.0, 3.0]);
}
