// Tests that `Manifold::decompose` hands each component exactly its own
// share of a body that carries properties. The ports of the C++ Decompose
// tests live in api.rs.

use super::*;

/// Three spheres of different sizes with normals, composed and decomposed:
/// every component is the solid it came from, property vertices included.
#[test]
fn decompose_gives_each_component_its_own_properties() {
    let spheres: Vec<Manifold> = [(1.0, 48, 0.0), (0.5, 12, 5.0), (2.0, 24, 12.0)]
        .iter()
        .map(|&(radius, segments, x)| {
            Manifold::sphere(radius, segments)
                .translate(Vec3::new(x, 0.0, 0.0))
                .calculate_normals(0, 60.0)
        })
        .collect();
    let parts = Manifold::compose(&spheres).decompose();
    assert_eq!(parts.len(), 3);

    for (part, sphere) in parts.iter().zip(&spheres) {
        let (part, sphere) = (
            part.get_mesh_glp::<f64, u32>(-1),
            sphere.get_mesh_glp::<f64, u32>(-1),
        );
        assert_eq!(part.vert_properties, sphere.vert_properties);
        assert_eq!(part.tri_verts, sphere.tri_verts);
        assert_eq!(part.merge_from_vert, sphere.merge_from_vert);
        assert_eq!(part.merge_to_vert, sphere.merge_to_vert);
    }
}
