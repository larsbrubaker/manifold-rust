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
/// check in `dedupe_edges` (divergence ledger entry 3), the cleanup filled the
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

/// The owners the sequential scans find: each orbit's smallest eligible
/// halfedge, walked with visited flags in ascending order.
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

/// `orbit_owners` must give each orbit's smallest eligible halfedge, as the
/// sequential scans do, on cubes touching along edges, for two rules.
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
    let halfedge = &model.as_impl().halfedge;
    let rules: [&(dyn Fn(&Halfedge) -> bool + Sync); 2] = [&|h| h.start_vert >= 0, &|h| {
        h.start_vert >= 0 && h.end_vert % 3 != 0
    }];
    for eligible in rules {
        let expected = sequential_orbit_owners(halfedge, eligible);
        assert_eq!(orbit_owners(halfedge, 0, eligible), Some(expected));
    }
}

/// The halfedges of a bicone: two fans of `k` triangles around a ring, so
/// both apexes have valence `k`. The top fan is stored in reverse, so the walk
/// around the top apex meets its halfedges in ascending order.
#[cfg(feature = "parallel")]
fn bicone_halfedges(k: usize) -> Vec<Halfedge> {
    let (top, bottom) = (k as i32, k as i32 + 1);
    let ring = |t: usize| (t % k) as i32;
    let top_tri = |t: usize| k - 1 - t % k;
    let bottom_tri = |t: usize| k + t % k;
    let mut halfedge = vec![
        Halfedge {
            start_vert: -1,
            end_vert: -1,
            paired_halfedge: -1,
            prop_vert: -1,
        };
        6 * k
    ];
    let mut set = |tri: usize, verts: [i32; 3], pairs: [usize; 3]| {
        for j in 0..3 {
            halfedge[3 * tri + j] = Halfedge {
                start_vert: verts[j],
                end_vert: verts[(j + 1) % 3],
                paired_halfedge: pairs[j] as i32,
                prop_vert: verts[j],
            };
        }
    };
    for t in 0..k {
        // Top (apex, r_t, r_t+1), bottom (apex, r_t+1, r_t).
        set(
            top_tri(t),
            [top, ring(t), ring(t + 1)],
            [
                3 * top_tri(t + k - 1) + 2,
                3 * bottom_tri(t) + 1,
                3 * top_tri(t + 1),
            ],
        );
        set(
            bottom_tri(t),
            [bottom, ring(t + 1), ring(t)],
            [
                3 * bottom_tri(t + 1) + 2,
                3 * top_tri(t) + 1,
                3 * bottom_tri(t + k - 1),
            ],
        );
    }
    halfedge
}

/// Two vertices of valence 100,000, in a mesh above `ORBIT_PAR_THRESHOLD`.
/// If every halfedge walked its orbit until it met a smaller one, this would
/// take 1e10 steps; walks stop at `OWNER_WALK_CAP` and each long orbit is
/// walked once.
#[cfg(feature = "parallel")]
#[test]
fn test_orbit_owners_of_a_high_valence_vertex() {
    let halfedge = bicone_halfedges(100_000);
    assert!(halfedge.len() >= ORBIT_PAR_THRESHOLD);
    let rules: [&(dyn Fn(&Halfedge) -> bool + Sync); 2] = [&|h| h.start_vert >= 0, &|h| {
        h.start_vert >= 0 && h.end_vert % 3 != 0
    }];
    for eligible in rules {
        let expected = sequential_orbit_owners(&halfedge, eligible);
        let owners = orbit_owners(&halfedge, ORBIT_PAR_THRESHOLD, eligible);
        assert_eq!(owners, Some(expected));
    }
}
