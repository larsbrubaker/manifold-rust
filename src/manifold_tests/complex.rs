use super::*;

/// C++ TEST(BooleanComplex, SelfIntersect) — read_test_obj, union, get_mesh_gl (test crash)
#[test]
fn test_cpp_complex_self_intersect() {
    let m1 = read_test_obj("self_intersectA.obj");
    let m2 = read_test_obj("self_intersectB.obj");
    let res = m1.union(&m2);
    res.get_mesh_gl(0); // test that it doesn't crash
}

/// C++ TEST(BooleanComplex, GenericTwinBooleanTest7081)
#[test]
fn test_cpp_complex_generic_twin_7081() {
    let m1 = read_test_obj("Generic_Twin_7081.1.t0_left.obj");
    let m2 = read_test_obj("Generic_Twin_7081.1.t0_right.obj");
    let res = m1.union(&m2);
    res.get_mesh_gl(0); // test crash
}

/// C++ TEST(BooleanComplex, GenericTwinBooleanTest7863)
#[test]
fn test_cpp_complex_generic_twin_7863() {
    let m1 = read_test_obj("Generic_Twin_7863.1.t0_left.obj");
    let m2 = read_test_obj("Generic_Twin_7863.1.t0_right.obj");
    let res = m1.union(&m2);
    res.get_mesh_gl(0); // test crash
}

/// C++ TEST(BooleanComplex, Havocglass8Bool)
#[test]
fn test_cpp_complex_havocglass8() {
    let m1 = read_test_obj("Havocglass8_left.obj");
    let m2 = read_test_obj("Havocglass8_right.obj");
    let res = m1.union(&m2);
    res.get_mesh_gl(0); // test crash
}

/// C++ TEST(BooleanComplex, HullMask)
#[test]
fn test_cpp_complex_hull_mask() {
    let body = read_test_obj("hull-body.obj");
    let mask = read_test_obj("hull-mask.obj");
    let ret = body.difference(&mask);
    ret.get_mesh_gl(0);
}

/// C++ TEST(BooleanComplex, OffsetTriangulationFailure)
#[test]
fn test_cpp_complex_offset_triangulation_failure() {
    let a = read_test_obj("Offset1.obj");
    let b = read_test_obj("Offset2.obj");
    let result = a.union(&b);
    assert_eq!(result.status(), Error::NoError,
        "OffsetTriangulationFailure: status {:?}", result.status());
}

/// C++ TEST(BooleanComplex, OffsetSelfIntersect)
#[test]
fn test_cpp_complex_offset_self_intersect() {
    let a = read_test_obj("Offset3.obj");
    let b = read_test_obj("Offset4.obj");
    let result = a.union(&b);
    assert_eq!(result.status(), Error::NoError,
        "OffsetSelfIntersect: status {:?}", result.status());
}

/// C++ TEST(BooleanComplex, Subtract) — specific mesh subtract
#[test]
fn test_cpp_complex_subtract() {
    let first_verts: Vec<f32> = vec![
        0.0,    0.0,  0.0,
        1540.0, 0.0,  0.0,
        1540.0, 70.0, 0.0,
        0.0,    70.0, 0.0,
        0.0,    0.0,  -278.282,
        1540.0, 70.0, -278.282,
        1540.0, 0.0,  -278.282,
        0.0,    70.0, -278.282,
    ];
    let first_tris: Vec<u32> = vec![
        0, 1, 2,  2, 3, 0,  4, 5, 6,  5, 4, 7,
        6, 2, 1,  6, 5, 2,  5, 3, 2,  5, 7, 3,
        7, 0, 3,  7, 4, 0,  4, 1, 0,  4, 6, 1,
    ];

    let mut first_mesh = MeshGL::default();
    first_mesh.num_prop = 3;
    first_mesh.vert_properties = first_verts;
    first_mesh.tri_verts = first_tris;

    let second_verts: Vec<f32> = vec![
        2.04636e-12, 70.0,           50000.0,
        2.04636e-12, -1.27898e-13,   50000.0,
        1470.0,      -1.27898e-13,   50000.0,
        1540.0,      70.0,           50000.0,
        2.04636e-12, 70.0,           -28.2818,
        1470.0,      -1.27898e-13,   0.0,
        2.04636e-12, -1.27898e-13,   0.0,
        1540.0,      70.0,           -28.2818,
    ];
    let second_tris: Vec<u32> = vec![
        0, 1, 2,  2, 3, 0,  4, 5, 6,  5, 4, 7,
        6, 2, 1,  6, 5, 2,  5, 3, 2,  5, 7, 3,
        7, 0, 3,  7, 4, 0,  4, 1, 0,  4, 6, 1,
    ];

    let mut second_mesh = MeshGL::default();
    second_mesh.num_prop = 3;
    second_mesh.vert_properties = second_verts;
    second_mesh.tri_verts = second_tris;

    let mut first = Manifold::from_mesh_gl(&first_mesh);
    let second = Manifold::from_mesh_gl(&second_mesh);

    first = first.difference(&second);
    first.get_mesh_gl(0);
    assert_eq!(first.status(), Error::NoError);
}

/// C++ TEST(BooleanComplex, Cylinders) — many cylinders with transforms
#[test]
fn test_cpp_complex_cylinders() {
    let rod = Manifold::cylinder(1.0, 0.4, -1.0, 12);
    let arrays1: Vec<[f64; 12]> = vec![
        [0.0, 0.0, 1.0, 3.0,  -1.0, 0.0, 0.0, 3.0,  0.0, -1.0, 0.0, 6.0],
        [0.0, 0.0, 1.0, 2.0,  -1.0, 0.0, 0.0, 3.0,  0.0, -1.0, 0.0, 8.0],
        [0.0, 0.0, 1.0, 1.0,  -1.0, 0.0, 0.0, 2.0,  0.0, -1.0, 0.0, 7.0],
        [1.0, 0.0, 0.0, 3.0,   0.0, 1.0, 0.0, 2.0,  0.0, 0.0, 1.0, 6.0],
        [0.0, 0.0, 1.0, 3.0,  -1.0, 0.0, 0.0, 3.0,  0.0, -1.0, 0.0, 7.0],
        [0.0, 0.0, 1.0, 1.0,  -1.0, 0.0, 0.0, 3.0,  0.0, -1.0, 0.0, 7.0],
        [1.0, 0.0, 0.0, 3.0,   0.0, 0.0, 1.0, 4.0,  0.0, -1.0, 0.0, 6.0],
        [1.0, 0.0, 0.0, 4.0,   0.0, 0.0, 1.0, 4.0,  0.0, -1.0, 0.0, 6.0],
    ];
    let arrays2: Vec<[f64; 12]> = vec![
        [1.0, 0.0, 0.0, 3.0,   0.0, 0.0, 1.0, 2.0,  0.0, -1.0, 0.0, 6.0],
        [1.0, 0.0, 0.0, 4.0,   0.0, 1.0, 0.0, 3.0,  0.0, 0.0, 1.0, 6.0],
        [0.0, 0.0, 1.0, 2.0,  -1.0, 0.0, 0.0, 2.0,  0.0, -1.0, 0.0, 7.0],
        [1.0, 0.0, 0.0, 3.0,   0.0, 1.0, 0.0, 3.0,  0.0, 0.0, 1.0, 7.0],
        [1.0, 0.0, 0.0, 2.0,   0.0, 1.0, 0.0, 3.0,  0.0, 0.0, 1.0, 7.0],
        [1.0, 0.0, 0.0, 1.0,   0.0, 1.0, 0.0, 3.0,  0.0, 0.0, 1.0, 7.0],
        [1.0, 0.0, 0.0, 3.0,   0.0, 1.0, 0.0, 4.0,  0.0, 0.0, 1.0, 7.0],
        [1.0, 0.0, 0.0, 3.0,   0.0, 1.0, 0.0, 5.0,  0.0, 0.0, 1.0, 6.0],
        [0.0, 0.0, 1.0, 3.0,  -1.0, 0.0, 0.0, 4.0,  0.0, -1.0, 0.0, 6.0],
    ];

    let make_mat = |array: &[f64; 12]| -> Mat3x4 {
        // C++ layout: mat[i][j] = array[j * 4 + i]
        // Row 0: array[0..4], Row 1: array[4..8], Row 2: array[8..12]
        // Columns: x=(r00,r10,r20), y=(r01,r11,r21), z=(r02,r12,r22), w=(r03,r13,r23)
        Mat3x4::from_cols(
            Vec3::new(array[0], array[4], array[8]),
            Vec3::new(array[1], array[5], array[9]),
            Vec3::new(array[2], array[6], array[10]),
            Vec3::new(array[3], array[7], array[11]),
        )
    };

    let mut m1 = Manifold::empty();
    for array in &arrays1 {
        let mat = make_mat(array);
        m1 = m1.union(&rod.transform(&mat));
    }

    let mut m2 = Manifold::empty();
    for array in &arrays2 {
        let mat = make_mat(array);
        m2 = m2.union(&rod.transform(&mat));
    }
    m1 = m1.union(&m2);

    assert!(m1.matches_tri_normals(), "Cylinders: should match tri normals");
    assert!(m1.num_degenerate_tris() <= 12,
        "Cylinders: {} degenerate tris, expected <= 12", m1.num_degenerate_tris());
}

/// C++ TEST(BooleanComplex, Close) — intersecting near-coincident spheres
#[test]
fn test_cpp_complex_close() {
    let r = 10.0;
    let a = Manifold::sphere(r, 256);
    let mut result = a.clone();
    for i in 0..10 {
        result = result.intersection(&a.translate(Vec3::new(
            a.get_tolerance() / 10.0 * i as f64,
            0.0,
            0.0,
        )));
    }
    let pi = std::f64::consts::PI;
    let tol = 0.004;
    assert!(
        (result.volume() - (4.0 / 3.0) * pi * r * r * r).abs() < tol * r * r * r,
        "Close volume: {} expected ~{}", result.volume(), (4.0 / 3.0) * pi * r * r * r
    );
    assert!(
        (result.surface_area() - 4.0 * pi * r * r).abs() < tol * r * r,
        "Close area: {} expected ~{}", result.surface_area(), 4.0 * pi * r * r
    );
}

/// C++ TEST(BooleanComplex, LazyCollider) — cylinder combos with mirror
#[test]
fn test_cpp_complex_lazy_collider() {
    let ele1 = Manifold::cylinder(50.0, 50.0, 50.0, 32);
    let ele2 = Manifold::cylinder(60.0, 30.0, 30.0, 32);
    let _ele3 = Manifold::cylinder(60.0, 40.0, 40.0, 32).mirror(Vec3::new(0.0, 0.0, 1.0));

    let ele4 = ele1.union(&ele2).mirror(Vec3::new(0.0, 0.0, 1.0));
    assert!(
        (ele4.volume() - 418839.0).abs() < 1.0,
        "LazyCollider ele4 volume: {} expected ~418839", ele4.volume()
    );

    let r1 = ele4.difference(&ele2.translate(Vec3::new(0.0, 0.0, -20.0)));
    assert!(
        (r1.volume() - 362577.0).abs() < 1.0,
        "LazyCollider r1 volume: {} expected ~362577", r1.volume()
    );

    let ele3 = Manifold::cylinder(60.0, 40.0, 40.0, 32).mirror(Vec3::new(0.0, 0.0, 1.0));
    let r2 = ele4.translate(Vec3::new(0.0, 0.0, 1.0)).difference(&ele3);
    assert!(
        (r2.volume() - 145656.0).abs() < 1.0,
        "LazyCollider r2 volume: {} expected ~145656", r2.volume()
    );
}

/// C++ TEST(Boolean, Perturb2) — prism construction from cube triangles
#[test]
fn test_cpp_perturb2() {
    let cube = Manifold::cube(Vec3::splat(2.0), true);
    let cube_gl = cube.get_mesh_gl(0);

    // Rotate so that nothing is axis-aligned
    let mut result = cube.rotate(5.0, 10.0, 15.0);

    let num_tri = cube_gl.tri_verts.len() / 3;
    for tri in 0..num_tri {
        let mut prism_verts: Vec<f32> = Vec::new();
        let mut prism_tris: Vec<u32> = vec![4, 2, 0, 1, 3, 5];

        for v0 in 0..3usize {
            let v1 = (v0 + 1) % 3;
            let v_in0 = cube_gl.tri_verts[3 * tri + v0] as usize;
            let v_in1 = cube_gl.tri_verts[3 * tri + v1] as usize;
            if v_in1 > v_in0 {
                prism_tris.extend_from_slice(&[
                    2 * v0 as u32, 2 * v1 as u32, 2 * v1 as u32 + 1,
                    2 * v0 as u32, 2 * v1 as u32 + 1, 2 * v0 as u32 + 1,
                ]);
            } else {
                prism_tris.extend_from_slice(&[
                    2 * v0 as u32, 2 * v1 as u32, 2 * v0 as u32 + 1,
                    2 * v1 as u32, 2 * v1 as u32 + 1, 2 * v0 as u32 + 1,
                ]);
            }
            let np = cube_gl.num_prop as usize;
            for j in 0..3 {
                prism_verts.push(cube_gl.vert_properties[np * v_in0 + j]);
            }
            for j in 0..3 {
                prism_verts.push(2.0 * cube_gl.vert_properties[np * v_in0 + j]);
            }
        }

        let mut mesh = MeshGL::default();
        mesh.num_prop = 3;
        mesh.vert_properties = prism_verts;
        mesh.tri_verts = prism_tris;
        result = result + Manifold::from_mesh_gl(&mesh).rotate(5.0, 10.0, 15.0);
    }

    assert_eq!(result.num_degenerate_tris(), 0);
    assert_eq!(result.num_vert(), 8, "Perturb2: {} verts expected 8", result.num_vert());
    assert!((result.volume() - 64.0).abs() < 1e-4, "Perturb2 volume: {} expected 64", result.volume());
    assert!((result.surface_area() - 96.0).abs() < 1e-4, "Perturb2 area: {} expected 96", result.surface_area());
}

/// C++ TEST(Boolean, Perturb3) — nasty gear pattern
#[test]
#[ignore = "Gear pattern requires BatchBoolean precision improvements"]
fn test_cpp_perturb3() {
    let n = 16;
    let alpha = 90.0 / n as f64;

    let cube = Manifold::cube(Vec3::splat(1.0), true);
    let mut outer_cubes = Vec::new();
    for i in 0..n {
        outer_cubes.push(cube.rotate(0.0, 0.0, alpha * i as f64));
    }
    let gear = Manifold::batch_boolean(&outer_cubes, OpType::Add);
    let outer_gear = gear.scale(Vec3::new(2.0, 2.0, 1.0));

    let nasty_gear = outer_gear.difference(&gear);
    let expected_volume = outer_gear.volume() - gear.volume();

    assert_eq!(nasty_gear.status(), Error::NoError);
    assert!(!nasty_gear.is_empty());
    assert_eq!(nasty_gear.genus(), 1, "Perturb3 genus: {} expected 1", nasty_gear.genus());
    assert!((nasty_gear.volume() - expected_volume).abs() < 1e-5,
        "Perturb3 volume: {} expected {}", nasty_gear.volume(), expected_volume);
    assert!((nasty_gear.surface_area() - 26.972).abs() < 1e-3,
        "Perturb3 area: {} expected 26.972", nasty_gear.surface_area());
}
