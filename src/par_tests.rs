// par_tests.rs — the parallel sites must give the sequential output, on any
// thread count.

use crate::linalg::Vec3;
use crate::manifold::Manifold;

/// FNV-1a over every field of the output mesh, with the length of every list.
/// Of `run_original_id`, which comes from the process-wide ID counter, only
/// the length is hashed.
fn fingerprint(m: &Manifold) -> u64 {
    let gl = m.get_mesh_gl64(-1);
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    let mut eat = |bytes: &[u8]| {
        for &b in bytes {
            h ^= u64::from(b);
            h = h.wrapping_mul(0x0100_0000_01b3);
        }
    };
    eat(&gl.num_prop.to_le_bytes());
    eat(&gl.tolerance.to_bits().to_le_bytes());
    for list in [&gl.vert_properties, &gl.run_transform, &gl.halfedge_tangent] {
        eat(&(list.len() as u64).to_le_bytes());
        for x in list {
            eat(&x.to_bits().to_le_bytes());
        }
    }
    for list in [
        &gl.tri_verts,
        &gl.merge_from_vert,
        &gl.merge_to_vert,
        &gl.run_index,
        &gl.face_id,
    ] {
        eat(&(list.len() as u64).to_le_bytes());
        for x in list {
            eat(&x.to_le_bytes());
        }
    }
    eat(&(gl.run_original_id.len() as u64).to_le_bytes());
    eat(&(gl.run_flags.len() as u64).to_le_bytes());
    eat(&gl.run_flags);
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
            0x0116_641f_953c_52d0,
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
    for (size, hash) in [(1.0, 0x89fd_3a6c_41b3_a52c), (1.5, 0x230c_a78a_1804_c434)] {
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

/// A checkerboard of unit cubes, `n` a side, as one mesh whose cubes share
/// the verts where they touch, so `split_pinched_verts` has an orbit to split
/// at every shared vert.
fn pinched_checkerboard(n: u32) -> crate::types::MeshGL64 {
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
    gl
}

/// An imported checkerboard of 20 a side: `split_pinched_verts` has
/// thousands of orbits to split, in a mesh over the 100,001-halfedge
/// thresholds of the orbit and edge-flag scans. Runs in every build, against
/// `main`'s hash.
#[test]
fn pinched_checkerboard_import_keeps_the_sequential_output() {
    let model = || {
        let result = Manifold::from_mesh_gl64(&pinched_checkerboard(20));
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
            0xaf78_1fb7_4bc9_643f,
            "{threads} threads"
        );
    }
}

/// The same import with a tangent on every halfedge, so the output's
/// `halfedge_tangent` is not empty and has gone through the cleanup and
/// `sort_geometry`'s gather, which permutes it with the halfedges above its
/// 100,000 threshold. The tangents are multiples of 1/4, exact everywhere.
#[test]
fn pinched_checkerboard_import_with_tangents_keeps_the_sequential_output() {
    let model = || {
        let mut gl = pinched_checkerboard(20);
        let quarters = |k: usize, m: usize| (k % m) as f64 / 4.0 - 1.0;
        gl.halfedge_tangent = (0..gl.tri_verts.len())
            .flat_map(|k| {
                [
                    quarters(k, 9),
                    quarters(k, 7),
                    quarters(k, 5),
                    quarters(k, 4),
                ]
            })
            .collect();
        let result = Manifold::from_mesh_gl64(&gl);
        assert!(
            3 * result.num_tri() > 100_001,
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
            0x5829_215b_5954_d206,
            "{threads} threads"
        );
    }
}
