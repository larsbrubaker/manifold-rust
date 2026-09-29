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

// cross_section_ops.rs — the Clipper2-backed operations on CrossSection:
// two-operand booleans, batch booleans and compose, decompose, simplify,
// offset, Minkowski sum, warp and convex hull.
//
// Ports the corresponding members of src/cross_section/cross_section.cpp.
// A child module of cross_section.rs (which owns the struct, constructors,
// transforms and queries) so callers keep the same `CrossSection::...` paths.
// Every operation reads its contours through `paths()` (C++ `GetPaths`), so
// pending transforms are applied first.

use clipper2_rust::{
    boolean_op_d, boolean_op_tree_d, difference_d, inflate_paths_d, intersect_d, minkowski_sum_d,
    simplify_paths, union_d, union_subjects_d, ClipType, EndType, FillRule, JoinType, PathsD,
    PolyTreeD,
};

use std::cmp::Ordering;

use super::{from_paths, path_area, to_paths, CrossSection, PRECISION};
use crate::linalg::Vec2;
use crate::polygon::ccw;
use crate::math;
use crate::types::{OpType, Quality, Rect};

impl CrossSection {
    pub fn union(&self, other: &Self) -> Self {
        Self::from_raw(from_paths(&union_d(
            &to_paths(&self.paths()),
            &to_paths(&other.paths()),
            FillRule::Positive,
            PRECISION,
        )))
    }

    pub fn intersection(&self, other: &Self) -> Self {
        Self::from_raw(from_paths(&intersect_d(
            &to_paths(&self.paths()),
            &to_paths(&other.paths()),
            FillRule::Positive,
            PRECISION,
        )))
    }

    pub fn difference(&self, other: &Self) -> Self {
        Self::from_raw(from_paths(&difference_d(
            &to_paths(&self.paths()),
            &to_paths(&other.paths()),
            FillRule::Positive,
            PRECISION,
        )))
    }
    /// Split into topologically disconnected components, each one outline
    /// with zero or more holes. Mirrors C++ `CrossSection::Decompose`: fewer
    /// than two contours return `self` unchanged; otherwise a Positive union
    /// into a Clipper2 PolyTree, whose containment links decide which holes
    /// belong to which outline, walked as `decompose_outline` /
    /// `decompose_hole` do and emitted in reverse push order.
    pub fn decompose(&self) -> Vec<Self> {
        if self.paths().len() < 2 {
            return vec![self.clone()];
        }
        let mut tree = PolyTreeD::new();
        boolean_op_tree_d(
            ClipType::Union,
            FillRule::Positive,
            &to_paths(&self.paths()),
            &PathsD::new(),
            &mut tree,
            PRECISION,
        );
        let mut comps = Vec::new();
        decompose_outlines(&tree, 0, &mut comps);
        comps
            .iter()
            .rev()
            .map(|poly| Self::from_raw(from_paths(poly)))
            .collect()
    }

    /// Remove vertices closer than `epsilon` to the line through their
    /// neighbours. Mirrors C++ `CrossSection::Simplify`: a Positive union
    /// into a Clipper2 PolyTree, flattened as C++ `flatten` does (each node's
    /// descendants before the node, so holes precede their outline), contours
    /// dropped when `|Area| <= max(box width, box height) * epsilon`, then
    /// `SimplifyPaths` on the closed survivors.
    pub fn simplify(&self, epsilon: f64) -> Self {
        let mut tree = PolyTreeD::new();
        boolean_op_tree_d(
            ClipType::Union,
            FillRule::Positive,
            &to_paths(&self.paths()),
            &PathsD::new(),
            &mut tree,
            PRECISION,
        );
        let mut polys = PathsD::new();
        flatten(&tree, 0, &mut polys);
        let filtered: PathsD = polys
            .into_iter()
            .filter(|poly| {
                let area = path_area(poly);
                let mut bx = Rect::new();
                for vert in poly {
                    bx.union_point(Vec2::new(vert.x, vert.y));
                }
                let size = bx.size();
                area.abs() > size.x.max(size.y) * epsilon
            })
            .collect();
        Self::from_raw(from_paths(&simplify_paths(&filtered, epsilon, true)))
    }

    /// Offset with the C++ `CrossSection::Offset` defaults: Round joins,
    /// miter_limit 2.0, circularSegments 0 (segments from Quality).
    pub fn offset(&self, delta: f64) -> Self {
        self.offset_with_params(delta, 1, 2.0, 0)
    }

    /// Offset with explicit join type and segment count.
    /// join_type: 0=Square, 1=Round, 2=Miter, 3=Bevel (the C++
    /// `CrossSection::JoinType` enumerator order). Other codes fall through to
    /// Square, the value C++ `jt()` starts from before its switch.
    pub fn offset_with_params(
        &self,
        delta: f64,
        join_type: i32,
        miter_limit: f64,
        circular_segments: i32,
    ) -> Self {
        let jt = match join_type {
            1 => JoinType::Round,
            2 => JoinType::Miter,
            3 => JoinType::Bevel,
            _ => JoinType::Square,
        };
        // For round joins, compute arc_tolerance from circular_segments (or,
        // when it is <= 2, Quality's count for radius delta) to get the exact
        // segment count. Matches C++ CrossSection::Offset:
        //   arc_tol = (math::cos(π/n) - 1) * -|delta|
        let arc_tol = if jt == JoinType::Round {
            let n = if circular_segments > 2 {
                circular_segments
            } else {
                Quality::get_circular_segments(delta)
            };
            let abs_delta = delta.abs();
            (math::cos(std::f64::consts::PI / n as f64) - 1.0) * -abs_delta
        } else {
            0.0
        };
        Self::from_raw(from_paths(&inflate_paths_d(
            &to_paths(&self.paths()),
            delta,
            jt,
            EndType::Polygon,
            miter_limit,
            PRECISION,
            arc_tol,
        )))
    }

    pub fn minkowski_sum(&self, other: &Self) -> Self {
        let mut result = Vec::new();
        for a in to_paths(&self.paths()) {
            for b in to_paths(&other.paths()) {
                result.extend(minkowski_sum_d(&a, &b, true, PRECISION));
            }
        }
        Self::from_raw(from_paths(&result))
    }

    /// Move every vertex through `f`, then re-union. Mirrors C++
    /// `CrossSection::Warp` / `WarpBatch`: vertices are visited in contour
    /// order, and the moved contours go through a FillRule::Positive union
    /// at `precision_`, so introduced self-intersections are resolved.
    pub fn warp<F: FnMut(&mut Vec2)>(&self, mut f: F) -> Self {
        let mut paths = to_paths(&self.paths());
        for path in paths.iter_mut() {
            for p in path.iter_mut() {
                let mut v = Vec2::new(p.x, p.y);
                f(&mut v);
                p.x = v.x;
                p.y = v.y;
            }
        }
        Self::from_raw(from_paths(&union_subjects_d(
            &paths,
            FillRule::Positive,
            PRECISION,
        )))
    }

    /// Boolean over a list of sections. Mirrors C++
    /// `CrossSection::BatchBoolean`: no sections give an empty section and
    /// one gives that section back untouched; Intersect folds pairwise
    /// `BooleanOp`s, while Add and Subtract run a single `BooleanOp` with the
    /// first section as subject and every later contour as a clip (so
    /// Subtract removes all of the tail from the head).
    pub fn batch_boolean(sections: &[Self], op: OpType) -> Self {
        match sections.len() {
            0 => return Self::default(),
            1 => return sections[0].clone(),
            _ => {}
        }
        let subjs = to_paths(&sections[0].paths());
        if let OpType::Intersect = op {
            let mut res = subjs;
            for s in &sections[1..] {
                res = boolean_op_d(
                    ClipType::Intersection,
                    FillRule::Positive,
                    &res,
                    &to_paths(&s.paths()),
                    PRECISION,
                );
            }
            return Self::from_raw(from_paths(&res));
        }
        let mut clips = PathsD::new();
        for s in &sections[1..] {
            clips.extend(to_paths(&s.paths()));
        }
        Self::from_raw(from_paths(&boolean_op_d(
            cliptype_of_op(op),
            FillRule::Positive,
            &subjs,
            &clips,
            PRECISION,
        )))
    }

    /// Convex hull of every vertex of `sections`, in section then contour
    /// order. Mirrors C++ `CrossSection::Hull(const std::vector<CrossSection>&)`,
    /// which reads each section through a by-value copy, so the inputs' own
    /// pending transforms stay pending (hence `clone().paths()`).
    pub fn hull_cross_sections(sections: &[Self]) -> Self {
        let mut points: Vec<Vec2> = Vec::new();
        for s in sections {
            for path in s.clone().paths().iter() {
                points.extend_from_slice(path);
            }
        }
        Self::hull_points(&points)
    }

    /// Convex hull of a point set. Mirrors C++ `CrossSection::Hull(SimplePolygon)`
    /// (and `Hull(Polygons)`, which flattens its contours into one list):
    /// the result is always exactly one contour, left degenerate as C++
    /// `HullImpl` leaves it — empty for fewer than three points, two vertices
    /// when every point is collinear.
    pub fn hull_points(points: &[Vec2]) -> Self {
        let mut pts = points.to_vec();
        Self::from_raw(vec![hull_impl(&mut pts)])
    }

    /// Batch union of the sections. Mirrors C++ `CrossSection::Compose`,
    /// which is `BatchBoolean(crossSections, OpType::Add)`.
    pub fn compose(sections: &[Self]) -> Self {
        Self::batch_boolean(sections, OpType::Add)
    }
}

/// C++ `cliptype_of_op`: Add is Union, Subtract Difference, Intersect
/// Intersection.
fn cliptype_of_op(op: OpType) -> ClipType {
    match op {
        OpType::Add => ClipType::Union,
        OpType::Subtract => ClipType::Difference,
        OpType::Intersect => ClipType::Intersection,
    }
}

/// C++ `V2Lesser`: by x, then by y. As an ordering, pairs that are neither
/// lesser are equal (only `+0.0` / `-0.0` ties between distinct values).
fn v2_lesser(a: &Vec2, b: &Vec2) -> Ordering {
    let lesser = |a: &Vec2, b: &Vec2| {
        if a.x == b.x {
            a.y < b.y
        } else {
            a.x < b.x
        }
    };
    if lesser(a, b) {
        Ordering::Less
    } else if lesser(b, a) {
        Ordering::Greater
    } else {
        Ordering::Equal
    }
}

/// C++ `HullBacktrack`: pop while the last two stack points and `pt` do not
/// turn strictly counter-clockwise under `CCW(.., 0.0)`.
fn hull_backtrack(pt: Vec2, stack: &mut Vec<Vec2>) {
    let mut sz = stack.len();
    while sz >= 2 && ccw(stack[sz - 2], stack[sz - 1], pt, 0.0) <= 0 {
        stack.pop();
        sz = stack.len();
    }
}

/// C++ `HullImpl` (cross_section.cpp:183-206), Andrew's monotone chain:
/// sorts `pts` in place with `V2Lesser`, builds the lower chain forwards and
/// the upper chain backwards, drops each chain's last point and returns
/// lower then upper. Fewer than three points give an empty path.
fn hull_impl(pts: &mut [Vec2]) -> Vec<Vec2> {
    if pts.len() < 3 {
        return Vec::new();
    }
    pts.sort_by(v2_lesser);
    let mut lower: Vec<Vec2> = Vec::new();
    for &pt in pts.iter() {
        hull_backtrack(pt, &mut lower);
        lower.push(pt);
    }
    let mut upper: Vec<Vec2> = Vec::new();
    for &pt in pts.iter().rev() {
        hull_backtrack(pt, &mut upper);
        upper.push(pt);
    }
    upper.pop();
    lower.pop();
    lower.extend(upper);
    lower
}

/// C++ `decompose_outline` / `decompose_hole` (cross_section.cpp:126-151):
/// for each outline child of `node`, first recurse into every hole's own
/// outline children (islands), then push `[outline, holes...]`. The C++
/// recurses over sibling indices too; iterating them visits the same nodes
/// in the same order without a stack frame per sibling.
fn decompose_outlines(tree: &PolyTreeD, node: usize, polys: &mut Vec<PathsD>) {
    for &outline in tree.nodes[node].children() {
        let holes = tree.nodes[outline].children();
        let mut poly = PathsD::with_capacity(holes.len() + 1);
        poly.push(tree.nodes[outline].polygon().clone());
        for &hole in holes {
            decompose_outlines(tree, hole, polys);
            poly.push(tree.nodes[hole].polygon().clone());
        }
        polys.push(poly);
    }
}

/// C++ `flatten` (cross_section.cpp:153-164): for each child of `node`,
/// its whole subtree first, then the child's own contour. Iterating the
/// siblings replaces the C++'s index recursion with the same visit order.
fn flatten(tree: &PolyTreeD, node: usize, polys: &mut PathsD) {
    for &child in tree.nodes[node].children() {
        flatten(tree, child, polys);
        polys.push(tree.nodes[child].polygon().clone());
    }
}
