// Copyright 2026 Lars Brubaker
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// minkowski_union_regression_tests.rs — an IGNORED regression test that keeps
// an open exact-engine item on record: an operand pair on which the exact
// engine turns a folded but manifold, NoError operand into lost solid.
// Listed under "Known unresolved mismatches" in docs/CPP_DIVERGENCES.md, and
// shared 1:1 with manifold-sharp's MinkowskiUnionRegressionTests (its commit
// 1931e87), fixtures included.
//
// Origin: Manifold::minkowski_sum of Thingi10K 641145 (demo import,
// sphere(0.05781898171099809, 12)). Before divergence ledger entry 11,
// triangle 109's QuickHull hull was non-convex
// (quickhull::tests::test_thingi641145_triangle109_swept_hull_is_convex), and
// unioning it in left a zero-thickness fin. The two fixtures are the operands
// of the 19th pairwise boolean in minkowski.rs's reduction back then (287 and
// 286 verts), exported with get_mesh_gl64. With the hulls fixed nothing in the
// crate produces these operands any more; the case stays because the engine
// defect below is still there.
//
// The defect, shared bit for bit with manifold-sharp and with C++
// SetNormalsAndCoplanar (impl.cpp v3.5.2): the boolean is right, and
// simplify_topology's swap_degenerates pass loses the solid (skipping only
// that pass restores the exact volume). The coplanar flood fill gives a sound
// triangle 0.041 tall, next to the fin, the reversed normal of an
// opposite-facing coplanar seed. Projected through that normal it reads as
// inverted, so recursive_edge_swap's normal path swaps its long edge. The
// triangle across that edge lies in a different plane (+z), so the swap cuts
// out a real wedge. Importing fixture b runs the same pass
// (from_mesh_gl64), so b arrives 8.5e-4 short.
//
// Two narrow fixes were tried (in manifold-sharp) and rejected. Refusing
// opposite-facing, non-degenerate neighbours in the flood fill breaks
// test_cpp_simplify (40 tris, not 12/20: its internal double wall must
// merge). Never swapping a triangle taller than tolerance breaks
// test_cpp_nonconvex_convex_minkowski_sum (genus 3, not 5). A fix has to
// tell a fin the swap should resolve from a sound triangle that only looks
// inverted through the normal it inherited.

use crate::manifold::Manifold;
use crate::types::{BooleanEngine, Error, MeshGL64, OpType};

/// Format: "numVert numTri", then one "x y z" line per vertex (round-trip
/// doubles) and one "v0 v1 v2" line per triangle.
const UNION_A: &str = include_str!("testdata/minkowski-641145-union-a.txt");
const UNION_B: &str = include_str!("testdata/minkowski-641145-union-b.txt");

fn load_fixture(text: &str) -> Manifold {
    let mut lines = text.lines();
    let header: Vec<usize> = lines
        .next()
        .expect("fixture header")
        .split_whitespace()
        .map(|s| s.parse().expect("header count"))
        .collect();
    let (num_vert, num_tri) = (header[0], header[1]);
    let mut mesh = MeshGL64 {
        num_prop: 3,
        ..Default::default()
    };
    for line in lines.by_ref().take(num_vert) {
        mesh.vert_properties.extend(
            line.split_whitespace()
                .map(|s| s.parse::<f64>().expect("coordinate")),
        );
    }
    for line in lines.take(num_tri) {
        mesh.tri_verts.extend(
            line.split_whitespace()
                .map(|s| s.parse::<u64>().expect("index")),
        );
    }
    Manifold::from_mesh_gl64(&mesh)
}

/// A union contains both operands: nothing of either may lie outside it.
#[test]
#[ignore = "Open exact-engine item: swap_degenerates cuts solid next to a zero-thickness fin (shared with manifold-sharp and C++); see the file header"]
fn exact_union_of_thingi641145_partial_unions_contains_both_operands() {
    let a = load_fixture(UNION_A);
    let b = load_fixture(UNION_B);

    let union = a.union_with_engine(&b, BooleanEngine::Exact);
    let intersection = a.intersection_with_engine(&b, BooleanEngine::Exact);
    let expected = a.volume() + b.volume() - intersection.volume();

    assert_eq!(a.status(), Error::NoError);
    assert_eq!(b.status(), Error::NoError);
    assert_eq!(union.status(), Error::NoError);

    // Measured: 3.9e-3 relative short, 4.2e-4 of B outside the union.
    let relative = (union.volume() - expected).abs() / expected;
    assert!(
        relative < 1e-9,
        "union volume {} short by {relative:e} relative",
        union.volume()
    );
    let b_outside = b.boolean(&union, OpType::Subtract).volume();
    assert!(
        b_outside < 1e-12,
        "{b_outside:e} of B lies outside the union"
    );
    let a_outside = a.boolean(&union, OpType::Subtract).volume();
    assert!(
        a_outside < 1e-12,
        "{a_outside:e} of A lies outside the union"
    );
}
