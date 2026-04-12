use super::*;

/// Helper: square with a square hole, offset along x
fn square_hole(x_offset: f64) -> Polygons {
    vec![
        vec![
            Vec2::new(2.0 + x_offset, 2.0),
            Vec2::new(-2.0 + x_offset, 2.0),
            Vec2::new(-2.0 + x_offset, -2.0),
            Vec2::new(2.0 + x_offset, -2.0),
        ],
        vec![
            Vec2::new(-1.0 + x_offset, 1.0),
            Vec2::new(1.0 + x_offset, 1.0),
            Vec2::new(1.0 + x_offset, -1.0),
            Vec2::new(-1.0 + x_offset, -1.0),
        ],
    ]
}

/// Load an OBJ file from the C++ test models directory.
fn read_test_obj(filename: &str) -> Manifold {
    let path = format!(
        "{}/cpp-reference/manifold/test/models/{}",
        env!("CARGO_MANIFEST_DIR"),
        filename
    );
    let contents = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read OBJ file {}: {}", path, e));
    let mut verts: Vec<f32> = Vec::new();
    let mut tri_verts: Vec<u32> = Vec::new();
    for line in contents.lines() {
        let line = line.trim();
        if line.starts_with("v ") {
            let parts: Vec<f64> = line[2..]
                .split_whitespace()
                .filter_map(|s| s.parse().ok())
                .collect();
            if parts.len() >= 3 {
                verts.push(parts[0] as f32);
                verts.push(parts[1] as f32);
                verts.push(parts[2] as f32);
            }
        } else if line.starts_with("f ") {
            let indices: Vec<u32> = line[2..]
                .split_whitespace()
                .filter_map(|s| s.split('/').next()?.parse::<u32>().ok().map(|i| i - 1))
                .collect();
            for i in 1..indices.len() - 1 {
                tri_verts.push(indices[0]);
                tri_verts.push(indices[i]);
                tri_verts.push(indices[i + 1]);
            }
        }
    }
    let mut mesh = MeshGL::default();
    mesh.num_prop = 3;
    mesh.vert_properties = verts;
    mesh.tri_verts = tri_verts;
    Manifold::from_mesh_gl(&mesh)
}

mod basic;
mod boolean;
mod complex;
mod advanced;
mod hull;
mod api;
mod sdf;
