use super::*;
use crate::face_op::set_normals_and_coplanar;
use crate::impl_mesh::ManifoldImpl;
use crate::linalg::Mat3x4;

#[test]
fn test_pair_up() {
    let mut halfedges = vec![
        Halfedge {
            start_vert: 0,
            end_vert: 1,
            paired_halfedge: -1,
            prop_vert: 0,
        },
        Halfedge {
            start_vert: 1,
            end_vert: 0,
            paired_halfedge: -1,
            prop_vert: 1,
        },
    ];
    pair_up(&mut halfedges, 0, 1);
    assert_eq!(halfedges[0].paired_halfedge, 1);
    assert_eq!(halfedges[1].paired_halfedge, 0);
}

#[test]
fn test_cleanup_topology_noop_on_clean_mesh() {
    let mut m = ManifoldImpl::tetrahedron(&Mat3x4::identity());
    set_normals_and_coplanar(&mut m);
    let before_verts = m.vert_pos.len();
    let before_halfedges = m.halfedge.len();
    cleanup_topology(&mut m);
    // Tetrahedron is already 2-manifold; cleanup should not add verts/halfedges
    assert_eq!(m.vert_pos.len(), before_verts);
    assert_eq!(m.halfedge.len(), before_halfedges);
}

#[test]
fn test_simplify_topology_noop_on_clean_mesh() {
    let mut m = ManifoldImpl::cube(&Mat3x4::identity());
    set_normals_and_coplanar(&mut m);
    simplify_topology(&mut m, 0);
    // After simplify, cube should still be 2-manifold (no degenerate edges)
    // (Some verts/edges may be removed, but topology must be valid)
    // Just check it's still 2-manifold where halfedges are valid
    let valid = m
        .halfedge
        .iter()
        .filter(|h| h.paired_halfedge >= 0)
        .all(|h| h.paired_halfedge < m.halfedge.len() as i32);
    assert!(valid, "invalid paired halfedge after simplify");
}

/// Regression fixture shared with manifold-sharp (DedupeEdgesRegressionTests):
/// 852 triangles cut from a 28060-triangle exact union right before
/// DedupeEdges, within three vertex rings of its sixteen duplicate edges.
/// Format: "numVert numTri", then one "x y z" per vertex, then one line per
/// triangle of three (start, end, pair) halfedges, pair -1 where it fell
/// outside the cut. Positions are round-trip formatted, so bit-exact.
const DEDUPE_STALE_DUPLICATE: &str = include_str!("testdata/dedupe-stale-duplicate.txt");

fn load_halfedge_fixture(text: &str) -> ManifoldImpl {
    let mut lines = text.lines();
    let header: Vec<usize> = lines
        .next()
        .unwrap()
        .split_whitespace()
        .map(|s| s.parse().unwrap())
        .collect();
    let (num_vert, num_tri) = (header[0], header[1]);
    let mut m = ManifoldImpl::new();
    for _ in 0..num_vert {
        let p: Vec<f64> = lines
            .next()
            .unwrap()
            .split_whitespace()
            .map(|s| s.parse().unwrap())
            .collect();
        m.vert_pos.push(Vec3::new(p[0], p[1], p[2]));
    }
    for _ in 0..num_tri {
        let p: Vec<i32> = lines
            .next()
            .unwrap()
            .split_whitespace()
            .map(|s| s.parse().unwrap())
            .collect();
        for k in 0..3 {
            m.halfedge.push(Halfedge {
                start_vert: p[3 * k],
                end_vert: p[3 * k + 1],
                paired_halfedge: p[3 * k + 2],
                prop_vert: p[3 * k],
            });
        }
    }
    m
}

/// Splitting duplicated edges only relabels vertices and adds zero-area
/// triangles, so no pre-existing triangle corner may change position. Before
/// the stale-entry check, a duplicate already resolved by an earlier repair in
/// the same pass was "repaired" again, relabelling an orbit to a copy of the
/// wrong vertex: 16 corners jumped ~0.13 and a union lost 2.3e-5 of volume.
///
/// That check is a deliberate divergence from C++ `DedupeEdges`, which repairs
/// every entry it collected. It moved the import pins in `robust::thingi_tests`:
/// Thingi10K #1147177 from 3201 verts / 6418 tris to 3206 / 6424 (its current
/// counts are pinned in `thingi_1147177_import_counts`) and #939888 from
/// 860 / 1716 to 861 / 1718. manifold-sharp has the same check
/// (`EdgeOp.Dedupe.cs`).
#[test]
fn test_dedupe_edges_never_moves_a_triangle_corner() {
    let mut m = load_halfedge_fixture(DEDUPE_STALE_DUPLICATE);
    let before: Vec<Vec3> = m
        .halfedge
        .iter()
        .map(|h| m.vert_pos[h.start_vert as usize])
        .collect();

    dedupe_edges(&mut m);

    let moved = before
        .iter()
        .enumerate()
        .filter(|(h, p)| {
            let after = m.vert_pos[m.halfedge[*h].start_vert as usize];
            after.x != p.x || after.y != p.y || after.z != p.z
        })
        .count();
    assert_eq!(moved, 0, "a corner that moves changes the solid");
}

/// Two operands cut from a union in a BOSL2 `cubetruss` model, around a
/// concave corner. Format: "3 numVert numTri tolerance 0", then one "x y z"
/// per vertex, then one "v0 v1 v2" per triangle. Positions are round-trip
/// formatted, so bit-exact. Neither mesh carries face IDs.
const UNION_CONCAVE_CORNER_A: &str = include_str!("testdata/union-concave-corner-a.txt");
const UNION_CONCAVE_CORNER_B: &str = include_str!("testdata/union-concave-corner-b.txt");

fn load_tri_fixture(text: &str) -> crate::manifold::Manifold {
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap().split_whitespace().collect();
    let (num_vert, num_tri) = (header[1].parse().unwrap(), header[2].parse().unwrap());
    let mut mesh = crate::types::MeshGL64 {
        num_prop: 3,
        tolerance: header[3].parse().unwrap(),
        ..Default::default()
    };
    for line in lines.by_ref().take(num_vert) {
        mesh.vert_properties
            .extend(line.split_whitespace().map(|s| s.parse::<f64>().unwrap()));
    }
    for line in lines.take(num_tri) {
        mesh.tri_verts
            .extend(line.split_whitespace().map(|s| s.parse::<u64>().unwrap()));
    }
    crate::manifold::Manifold::from_mesh_gl64(&mesh)
}

/// The union's volume must match inclusion-exclusion. Before the stale-entry
/// check in `dedupe_edges` (noted at its repair loop), the cleanup filled the
/// concave corner, as it does in C++ Manifold 3.5.2.
#[test]
fn test_union_keeps_concave_corner() {
    use crate::types::OpType;
    let a = load_tri_fixture(UNION_CONCAVE_CORNER_A);
    let b = load_tri_fixture(UNION_CONCAVE_CORNER_B);
    let expected = a.volume() + b.volume() - a.boolean(&b, OpType::Intersect).volume();
    for (x, y) in [(&a, &b), (&b, &a)] {
        let union = x.boolean(y, OpType::Add);
        assert!(
            (union.volume() - expected).abs() < 1e-9 * expected,
            "union volume {}, expected {expected}",
            union.volume()
        );
    }
}

/// The sequential owner scan that `orbit_owners` replaces: ascending, each
/// unvisited eligible halfedge owns its orbit and marks it visited.
#[cfg(feature = "parallel")]
fn sequential_orbit_owners(
    halfedge: &[Halfedge],
    eligible: &dyn Fn(&Halfedge) -> bool,
) -> Vec<usize> {
    let mut visited = vec![false; halfedge.len()];
    let mut owners = Vec::new();
    for i in 0..halfedge.len() {
        if visited[i] || !eligible(&halfedge[i]) {
            continue;
        }
        owners.push(i);
        let mut current = i;
        loop {
            visited[current] = true;
            current = next_halfedge(halfedge[current].paired_halfedge) as usize;
            if current == i {
                break;
            }
        }
    }
    owners
}

/// `orbit_owners` with a threshold of 0 against the sequential scan, under
/// eligibility rules that make every, most, and few halfedges eligible.
#[cfg(feature = "parallel")]
fn assert_orbit_owners_match(halfedge: &[Halfedge]) {
    let rules: [&(dyn Fn(&Halfedge) -> bool + Sync); 3] = [
        &|h| h.start_vert >= 0,
        &|h| h.start_vert >= 0 && h.end_vert % 3 != 0,
        &|h| h.start_vert >= 0 && h.end_vert % 7 == 3,
    ];
    for eligible in rules {
        let expected = sequential_orbit_owners(halfedge, eligible);
        assert_eq!(orbit_owners(halfedge, 0, eligible), Some(expected));
    }
}

/// `orbit_owners` must give each orbit's smallest eligible halfedge, as the
/// sequential scans do, on cubes touching along edges.
#[cfg(feature = "parallel")]
#[test]
fn test_orbit_owners_match_the_sequential_scan() {
    let cube = crate::manifold::Manifold::cube(Vec3::splat(1.0), false);
    let mut model = crate::manifold::Manifold::empty();
    for x in 0..4 {
        for y in 0..4 {
            for z in 0..4 {
                if (x + y + z) % 2 == 0 {
                    let at = Vec3::new(f64::from(x), f64::from(y), f64::from(z));
                    model = model.union(&cube.translate(at));
                }
            }
        }
    }
    assert_orbit_owners_match(&model.as_impl().halfedge);
}

/// The same on orbits longer than `ORBIT_WALK_CAP`: a bipyramid whose two
/// apexes have valence 300, so their orbits are left to the sequential pass,
/// while the ring vertices' orbits (valence 4) resolve in the parallel walks.
#[cfg(feature = "parallel")]
#[test]
fn test_orbit_owners_match_the_sequential_scan_on_orbits_longer_than_the_cap() {
    let n = 300u64;
    let mut mesh = crate::types::MeshGL64 {
        num_prop: 3,
        ..Default::default()
    };
    for i in 0..n {
        let angle = i as f64 * std::f64::consts::TAU / n as f64;
        mesh.vert_properties
            .extend_from_slice(&[100.0 * angle.cos(), 100.0 * angle.sin(), 0.0]);
    }
    mesh.vert_properties
        .extend_from_slice(&[0.0, 0.0, 50.0, 0.0, 0.0, -50.0]);
    let (top, bottom) = (n, n + 1);
    for i in 0..n {
        let j = (i + 1) % n;
        mesh.tri_verts.extend_from_slice(&[i, j, top, j, i, bottom]);
    }
    let model = crate::manifold::Manifold::from_mesh_gl64(&mesh);
    assert_eq!(model.num_tri(), 2 * n as usize);
    let halfedge = &model.as_impl().halfedge;
    let longest = (0..halfedge.len())
        .map(|i| {
            let mut len = 1;
            let mut current = next_halfedge(halfedge[i].paired_halfedge) as usize;
            while current != i {
                len += 1;
                current = next_halfedge(halfedge[current].paired_halfedge) as usize;
            }
            len
        })
        .max();
    assert_eq!(longest, Some(n as usize));
    assert!(n as usize > super::orbits::ORBIT_WALK_CAP);
    assert_orbit_owners_match(halfedge);
}

/// The helix and core of a three-sided, one-turn thread (as a hole) from
/// MatterCAD, whose roots lie within an f32 step of the core's corners.
/// Every position is an f32 value, so the meshes load into a single-precision
/// `MeshGL` exactly.
const THREAD_CORE_UNION_THREAD: &str = include_str!("testdata/thread-core-union-thread.txt");
const THREAD_CORE_UNION_CORE: &str = include_str!("testdata/thread-core-union-core.txt");

fn load_single_precision_fixture(text: &str) -> crate::manifold::Manifold {
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap().split_whitespace().collect();
    let (num_vert, num_tri) = (header[1].parse().unwrap(), header[2].parse().unwrap());
    let mut mesh = crate::types::MeshGL {
        num_prop: 3,
        ..Default::default()
    };
    for line in lines.by_ref().take(num_vert) {
        mesh.vert_properties
            .extend(line.split_whitespace().map(|s| s.parse::<f32>().unwrap()));
    }
    for line in lines.take(num_tri) {
        mesh.tri_verts
            .extend(line.split_whitespace().map(|s| s.parse::<u32>().unwrap()));
    }
    crate::manifold::Manifold::from_mesh_gl(&mesh)
}

/// A single-precision mesh imports with its tolerance floored at an f32 step
/// of its scale, while epsilon stays at the precision of the f64 arithmetic.
#[test]
fn test_single_precision_import_tolerance_is_an_f32_step() {
    let core = load_single_precision_fixture(THREAD_CORE_UNION_CORE);
    let scale = core.as_impl().bbox.scale();
    assert_eq!(core.get_tolerance(), f64::from(f32::EPSILON) * scale);
    assert!(core.as_impl().epsilon < core.get_tolerance() / 1000.0);
}

/// Every edge of the union is at least the tolerance long, and the union
/// exported as a single-precision `MeshGL` has no triangle with two equal
/// corners. The thread's root cuts the core within a billionth of its corner,
/// leaving three new verts around one old one; one pass of the short-edge
/// collapse merges the first two and leaves the third on a 1e-9 edge it had
/// not flagged, a needle whose corners coincide once narrowed to f32.
#[test]
fn test_union_leaves_no_edge_shorter_than_tolerance() {
    use crate::types::OpType;
    let thread = load_single_precision_fixture(THREAD_CORE_UNION_THREAD);
    let core = load_single_precision_fixture(THREAD_CORE_UNION_CORE);
    let union = thread.boolean(&core, OpType::Add);
    let m = union.as_impl();
    let tolerance = m.tolerance;
    let short: Vec<f64> = m
        .halfedge
        .iter()
        .filter(|h| h.start_vert >= 0)
        .map(|h| {
            let d = m.vert_pos[h.end_vert as usize] - m.vert_pos[h.start_vert as usize];
            dot(d, d).sqrt()
        })
        .filter(|&len| len < tolerance)
        .collect();
    assert!(short.is_empty(), "edges under {tolerance}: {short:?}");

    let gl = union.get_mesh_gl(-1);
    let np = gl.num_prop as usize;
    let pos = |v: u32| &gl.vert_properties[v as usize * np..v as usize * np + 3];
    let pinched = gl
        .tri_verts
        .chunks(3)
        .filter(|t| pos(t[0]) == pos(t[1]) || pos(t[1]) == pos(t[2]) || pos(t[2]) == pos(t[0]))
        .count();
    assert_eq!(pinched, 0, "f32 triangles with coincident corners");

    let expected =
        thread.volume() + core.volume() - thread.boolean(&core, OpType::Intersect).volume();
    assert!((union.volume() - expected).abs() < 1e-9 * expected);
}
