// par_tests.rs — the parallel sites must give the sequential output, on any
// thread count.

use crate::linalg::Vec3;
use crate::manifold::Manifold;

/// FNV-1a (64-bit) over every field of the output `MeshGL64` but
/// `run_original_id`, which comes from the process-wide ID counter. In struct
/// order: `num_prop`; then `vert_properties`, `tri_verts`, `merge_from_vert`,
/// `merge_to_vert`, `run_index`, `run_transform`, `face_id`,
/// `halfedge_tangent` and `run_flags`, each as its length then its elements;
/// then `tolerance`. Lengths and integers are u64 little-endian, floats their
/// IEEE bits as u64 little-endian, `run_flags` one byte each.
/// `examples/boolean_perf.rs` hashes the same way.
fn fingerprint(m: &Manifold) -> u64 {
    let gl = m.get_mesh_gl64(-1);
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    let mut eat = |bytes: &[u8]| {
        for &b in bytes {
            h ^= u64::from(b);
            h = h.wrapping_mul(0x0100_0000_01b3);
        }
    };
    let ints = |eat: &mut dyn FnMut(&[u8]), list: &[u64]| {
        eat(&(list.len() as u64).to_le_bytes());
        for x in list {
            eat(&x.to_le_bytes());
        }
    };
    let floats = |eat: &mut dyn FnMut(&[u8]), list: &[f64]| {
        eat(&(list.len() as u64).to_le_bytes());
        for x in list {
            eat(&x.to_bits().to_le_bytes());
        }
    };
    eat(&gl.num_prop.to_le_bytes());
    floats(&mut eat, &gl.vert_properties);
    ints(&mut eat, &gl.tri_verts);
    ints(&mut eat, &gl.merge_from_vert);
    ints(&mut eat, &gl.merge_to_vert);
    ints(&mut eat, &gl.run_index);
    floats(&mut eat, &gl.run_transform);
    ints(&mut eat, &gl.face_id);
    floats(&mut eat, &gl.halfedge_tangent);
    eat(&(gl.run_flags.len() as u64).to_le_bytes());
    eat(&gl.run_flags);
    eat(&gl.tolerance.to_bits().to_le_bytes());
    h
}

/// `f`'s fingerprint, on `threads` threads when the `parallel` feature is on.
fn run_on(threads: usize, f: &(dyn Fn() -> Manifold + Sync)) -> u64 {
    #[cfg(feature = "parallel")]
    {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        pool.install(|| fingerprint(&f()))
    }
    #[cfg(not(feature = "parallel"))]
    {
        let _ = threads;
        fingerprint(&f())
    }
}

/// Pins a union of 500 cubes touching along edges and corners, one at a time,
/// so every union is below the parallel thresholds of the cleanup. Integer
/// coordinates, so the hash holds on every platform.
#[test]
fn checkerboard_union_keeps_the_sequential_output() {
    let model = || {
        let cube = Manifold::cube(Vec3::splat(1.0), false);
        let mut result = Manifold::empty();
        for x in 0..10 {
            for y in 0..10 {
                for z in 0..10 {
                    if (x + y + z) % 2 == 0 {
                        let at = Vec3::new(f64::from(x), f64::from(y), f64::from(z));
                        result = result.union(&cube.translate(at));
                    }
                }
            }
        }
        result
    };
    for threads in [1, 8] {
        assert_eq!(
            run_on(threads, &model),
            0x63ff_dc74_1011_f536,
            "{threads} threads"
        );
    }
}

/// Checkerboards of 20 cubes a side, over the 100,001-halfedge thresholds:
/// with the feature they take the parallel paths, without it the sequential
/// ones, and both must give `main`'s hashes. Unit cubes only touch; cubes of
/// 1.5 overlap, so the unions cut edges. Each of the four classes of cubes is
/// disjoint, so `compose` builds it and three unions join them.
#[test]
fn large_checkerboard_unions_keep_the_sequential_output() {
    for (size, hash) in [(1.0, 0x7802_34fb_b0d5_3e14), (1.5, 0x4a6e_f237_508a_1f18)] {
        let model = || {
            let cube = Manifold::cube(Vec3::splat(size), false);
            let class = |(ox, oy, oz): (i32, i32, i32)| {
                let mut cubes = Vec::new();
                for x in (ox..20).step_by(2) {
                    for y in (oy..20).step_by(2) {
                        for z in (oz..20).step_by(2) {
                            let at = Vec3::new(f64::from(x), f64::from(y), f64::from(z));
                            cubes.push(cube.translate(at));
                        }
                    }
                }
                Manifold::compose(&cubes)
            };
            let a = class((0, 0, 0)).union(&class((0, 1, 1)));
            let b = class((1, 0, 1)).union(&class((1, 1, 0)));
            let result = a.union(&b);
            assert!(
                3 * result.num_tri() > 100_001,
                "{} triangles",
                result.num_tri()
            );
            result
        };
        for threads in [1, 8] {
            assert_eq!(
                run_on(threads, &model),
                hash,
                "size {size}, {threads} threads"
            );
        }
    }
}

/// A checkerboard of `n` unit cubes a side imported from a `MeshGL64` whose
/// cubes share the verts where they touch, so `split_pinched_verts` has
/// thousands of orbits to split. With `tangents`, every halfedge carries a
/// tangent of small integers, which the import keeps and `sort_geometry`
/// gathers.
fn pinched_checkerboard(n: u32, tangents: bool) -> Manifold {
    // Corners as x | y << 1 | z << 2, wound outward.
    const QUADS: [[u32; 4]; 6] = [
        [0, 4, 6, 2],
        [1, 3, 7, 5],
        [0, 1, 5, 4],
        [2, 6, 7, 3],
        [0, 2, 3, 1],
        [4, 5, 7, 6],
    ];
    let side = n + 1;
    let mut gl = crate::types::MeshGL64 {
        num_prop: 3,
        ..Default::default()
    };
    for z in 0..side {
        for y in 0..side {
            for x in 0..side {
                gl.vert_properties
                    .extend([f64::from(x), f64::from(y), f64::from(z)]);
            }
        }
    }
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                if (x + y + z) % 2 != 0 {
                    continue;
                }
                let corner = |c: u32| {
                    let (cx, cy, cz) = (x + (c & 1), y + (c >> 1 & 1), z + (c >> 2));
                    u64::from(cx + side * (cy + side * cz))
                };
                for q in QUADS {
                    gl.tri_verts
                        .extend([q[0], q[1], q[2], q[0], q[2], q[3]].map(corner));
                }
            }
        }
    }
    if tangents {
        for i in 0..gl.tri_verts.len() as u32 {
            gl.halfedge_tangent.extend([
                f64::from(i % 7) - 3.0,
                f64::from(i % 5) - 2.0,
                f64::from(i % 3),
                f64::from(i % 4) * 0.5,
            ]);
        }
    }
    let result = Manifold::from_mesh_gl64(&gl);
    assert!(
        3 * result.num_tri() > 100_001,
        "{} triangles",
        result.num_tri()
    );
    result
}

/// An imported checkerboard of 20 a side, over the 100,001-halfedge
/// thresholds of the orbit and edge-flag scans. Runs in every build, against
/// `main`'s hash.
#[test]
fn pinched_checkerboard_import_keeps_the_sequential_output() {
    for threads in [1, 8] {
        assert_eq!(
            run_on(threads, &|| pinched_checkerboard(20, false)),
            0x1c3d_e063_12f8_0576,
            "{threads} threads"
        );
    }
}

/// The same import, 26 a side with halfedge tangents: over 100,000
/// triangles, so `sort_geometry` sorts the faces and gathers the tangents in
/// parallel with the feature.
#[test]
fn pinched_checkerboard_import_with_tangents_keeps_the_sequential_output() {
    let model = || {
        let result = pinched_checkerboard(26, true);
        assert!(
            result.num_tri() >= 100_000,
            "{} triangles",
            result.num_tri()
        );
        assert_eq!(
            result.get_mesh_gl64(-1).halfedge_tangent.len(),
            12 * result.num_tri()
        );
        result
    };
    for threads in [1, 8] {
        assert_eq!(
            run_on(threads, &model),
            0x3a43_20b2_3394_3b5d,
            "{threads} threads"
        );
    }
}
