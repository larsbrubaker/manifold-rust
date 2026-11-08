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

// csg_tree_nary_contact_tests.rs — eight solid tools that abut each other exactly, kept on record
// for an open exact-engine item (PORTING_PLAN.md). Every operand is clean, and the n-ary union and
// a pairwise left fold reach the same solid, but on the auto engine the n-ary answer carries
// zero-width seams that the self-intersection scan reports.
//
// The tools share whole planes (x = 16, z = 4, z = 2, ...), so a partial union can touch itself
// along a line. Symbolic perturbation gives such contact a consistent answer, but which seams
// survive depends on the order the operands meet in, and a seam is left with sheets lying on one
// line and vertices that do not line up. It is not specific to the n-ary tree: one two-operand
// union of these tools already reports contact, the left fold can be the clean one or not, and
// the robust engine is order-dependent the same way. The shape is right; the seams enclose no
// volume. The fixture is the batch from a bevel's cutters, exported with get_mesh_gl64.

use super::*;
use crate::manifold::Manifold;
use crate::types::{BooleanEngine, Error, MeshGL64};

const ABUTTING_TOOLS: &str = include_str!("testdata/nary-union-abutting-tools.txt");

fn abutting_tools() -> Vec<Manifold> {
    let mut lines = ABUTTING_TOOLS.lines();
    let operands: usize = lines
        .next()
        .expect("operand count")
        .parse()
        .expect("a count");
    (0..operands)
        .map(|_| {
            let header: Vec<&str> = lines.next().expect("operand header").split(' ').collect();
            let num_vert: usize = header[0].parse().expect("vertex count");
            let num_tri: usize = header[1].parse().expect("triangle count");
            let mut mesh = MeshGL64 {
                num_prop: 3,
                tolerance: header[2].parse().expect("tolerance"),
                ..Default::default()
            };
            for line in lines.by_ref().take(num_vert) {
                mesh.vert_properties.extend(
                    line.split(' ')
                        .map(|s| s.parse::<f64>().expect("coordinate")),
                );
            }
            for line in lines.by_ref().take(num_tri) {
                mesh.tri_verts
                    .extend(line.split(' ').map(|s| s.parse::<u64>().expect("index")));
            }
            let operand = Manifold::from_mesh_gl64(&mesh);
            assert_eq!(operand.status(), Error::NoError);
            operand
        })
        .collect()
}

fn n_ary_union(tools: &[Manifold], engine: BooleanEngine) -> Manifold {
    let mut leaves: Vec<CsgLeafNode> = tools
        .iter()
        .map(|tool| CsgLeafNode::new(tool.as_impl().clone()))
        .collect();
    let union = batch_union(&mut leaves, None, Some(engine), None);
    Manifold::from_impl(union.get_impl())
}

/// Every operand is clean, and the two ways of joining them reach one solid.
#[test]
fn the_n_ary_union_of_abutting_tools_is_the_same_solid_as_a_pairwise_fold() {
    let tools = abutting_tools();
    assert_eq!(tools.len(), 8);
    assert!(tools.iter().all(|tool| !tool.has_self_intersections()));

    let n_ary = n_ary_union(&tools, BooleanEngine::Auto);
    let fold = tools[1..].iter().fold(tools[0].clone(), |joined, tool| {
        joined.union_with_engine(tool, BooleanEngine::Auto)
    });

    assert_eq!(n_ary.status(), Error::NoError);
    assert_eq!(fold.status(), Error::NoError);
    assert!((n_ary.volume() - fold.volume()).abs() < 1e-9);
}

/// What the open item asks for: no seam left behind.
#[test]
#[ignore = "Open exact-engine item: unions of abutting operands can leave zero-width seams, depending on operand order; see the file header"]
fn the_n_ary_union_of_abutting_tools_leaves_no_self_contact() {
    let n_ary = n_ary_union(&abutting_tools(), BooleanEngine::Auto);

    assert!(!n_ary.has_self_intersections());
}
