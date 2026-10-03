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

// robust/coplanar_clip.rs — phase 3 of the intersection-graph build (the
// coplanar cross-copy, robust/intersection_graph.rs) and the exact clip and
// containment tests it runs, with the polygon's side of the work done once
// instead of once per call.
//
// `clip_segment_to_polygon` and `point_in_polygon_coplanar` (graph_geom.rs,
// since removed in favour of `CoplanarClipRegion`) recomputed the polygon's
// normal, dominant axis, 2D projection and orientation on every segment they
// clipped. Phase 3 clips every primitive of a triangle against the same
// overlap polygon, and a coplanar triangle can carry hundreds of primitives
// across hundreds of regions (a 12-triangle cube resting on an
// 11,652-triangle self-touching body spent 6.7 s here in manifold-sharp), so
// that setup — all rational — was most of the step's cost. The step also
// deduped by a linear scan of the destination's primitive list.
//
// Why no bit can move:
// - `CoplanarClipRegion::prepare` computes exactly the values the per-call
//   code computed, from the polygon alone, by the same operations in the same
//   order; `clip` and `contains` then run the same loops over them. Exact
//   rationals do not round, so computing a value once or N times gives the
//   same value.
// - The bounding-box reject answers "empty" only for a segment (or point)
//   whose projected box misses the polygon's projected box. When the polygon
//   has positive area, a convex CCW polygon is the intersection of its edge
//   half-planes, so the parametric clip's interval holds only parameters
//   whose points lie in the polygon, hence in its box, hence in the segment's
//   box too — disjoint boxes mean the clip is empty and the containment test
//   fails. A zero-area polygon's half-planes meet in an unbounded line, which
//   a box cannot bound, so the reject is off for it.
// - `PrimIndex` answers the scan's exact question — exact rational equality
//   (`R3Key` compares canonical fields, which is value equality), same
//   provenance, either orientation for a segment — and appends still happen
//   in the scan's order, so every primitive list ends up element for element
//   the same. Its sets are probe-only (never iterated), so their hashing
//   cannot reach a result.
//
// Shared with manifold-sharp (CoplanarClipRegion.cs and
// IntersectionGraphBuild.Types.cs's CrossCopyCoplanarRegions, its commit
// c776525), which also reports the step as its own progress phase,
// `Phase::CoplanarOverlaps`.

use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};

use super::exact::backend::{rat_is_zero, rat_one, rat_zero, Rational};
use super::exact::predicates::{orient2d_r, tri_normal_r};
use super::exact::rational::{R3Key, R2, R3};
use super::exact::Sign;
use super::intersection_graph::TriPrims;
use super::tri_tri::dominant_axis;
use crate::cancel::{is_cancelled, CancelToken};
use crate::progress::{begin_phase, complete_phase, Phase, ProgressReporter};

/// A convex coplanar polygon (at least three vertices) prepared for repeated
/// exact clipping and containment.
pub(super) struct CoplanarClipRegion {
    axis: usize,
    /// The polygon projected along `axis`, wound CCW.
    pts2: Vec<R2>,
    /// Whether the polygon has positive area, which is what makes the box
    /// reject sound (see the file header).
    reject_by_box: bool,
    min: R2,
    max: R2,
}

impl CoplanarClipRegion {
    pub(super) fn prepare(poly: &[R3]) -> Self {
        debug_assert!(poly.len() >= 3);
        let n = tri_normal_r(&poly[0], &poly[1], &poly[2]);
        let axis = dominant_axis(&n);
        let mut pts2: Vec<R2> = poly.iter().map(|p| p.project_drop(axis)).collect();
        if orient2d_r(&pts2[0], &pts2[1], &pts2[2]) == Sign::Neg {
            pts2.reverse();
        }
        let reject_by_box =
            (1..pts2.len() - 1).any(|i| orient2d_r(&pts2[0], &pts2[i], &pts2[i + 1]) != Sign::Zero);
        let (mut min, mut max) = (pts2[0].clone(), pts2[0].clone());
        for p in &pts2 {
            if p.x < min.x {
                min.x = p.x.clone();
            }
            if p.x > max.x {
                max.x = p.x.clone();
            }
            if p.y < min.y {
                min.y = p.y.clone();
            }
            if p.y > max.y {
                max.y = p.y.clone();
            }
        }
        Self {
            axis,
            pts2,
            reject_by_box,
            min,
            max,
        }
    }

    /// Clip segment (a,b) to the region (2D test via projection on the
    /// polygon's own plane). Returns a positive-length sub-segment or None.
    pub(super) fn clip(&self, a: &R3, b: &R3) -> Option<(R3, R3)> {
        let a2 = a.project_drop(self.axis);
        let b2 = b.project_drop(self.axis);
        if self.misses_box(&a2, &b2) {
            return None;
        }
        let dir = b2.sub(&a2);
        let pts2 = &self.pts2;

        // Parametric clip of [0,1] against each CCW edge halfplane.
        let mut t0 = rat_zero();
        let mut t1 = rat_one();
        for i in 0..pts2.len() {
            let e0 = &pts2[i];
            let e1 = &pts2[(i + 1) % pts2.len()];
            let edge = e1.sub(e0);
            // Signed distance numerators of a2 + t*dir against the edge line:
            // f(t) = cross(edge, a2 + t*dir - e0) = fa + t * fd.
            let fa = edge.cross(&a2.sub(e0));
            let fd = edge.cross(&dir);
            if rat_is_zero(&fd) {
                if fa < rat_zero() {
                    return None; // parallel and strictly outside
                }
                continue;
            }
            let t_hit = -&fa / &fd;
            if fd > rat_zero() {
                // entering: f grows with t → require t >= t_hit
                if t_hit > t0 {
                    t0 = t_hit;
                }
            } else if t_hit < t1 {
                t1 = t_hit;
            }
            if t0 >= t1 {
                return None;
            }
        }
        if t0 >= t1 {
            return None;
        }
        let seg = |t: &Rational| a.add(&b.sub(a).scale(t));
        Some((seg(&t0), seg(&t1)))
    }

    /// Exact point-in-region test for a point on the polygon's plane: inside
    /// or on the boundary.
    pub(super) fn contains(&self, p: &R3) -> bool {
        let p2 = p.project_drop(self.axis);
        if self.misses_box(&p2, &p2) {
            return false;
        }
        let pts2 = &self.pts2;
        for i in 0..pts2.len() {
            if orient2d_r(&pts2[i], &pts2[(i + 1) % pts2.len()], &p2) == Sign::Neg {
                return false;
            }
        }
        true
    }

    /// True when the projected segment's box is disjoint from the polygon's
    /// box, which (see the file header) proves the segment misses a
    /// positive-area polygon.
    fn misses_box(&self, a2: &R2, b2: &R2) -> bool {
        if !self.reject_by_box {
            return false;
        }
        (a2.x < self.min.x && b2.x < self.min.x)
            || (a2.x > self.max.x && b2.x > self.max.x)
            || (a2.y < self.min.y && b2.y < self.min.y)
            || (a2.y > self.max.y && b2.y > self.max.y)
    }
}

/// Probe-only membership sets mirroring one `TriPrims`'s lists, kept in step
/// by [`copy_through_region`], which is the only writer while phase 3 runs.
#[derive(Default)]
struct PrimIndex {
    segments: HashSet<(R3Key, R3Key, usize)>,
    points: HashSet<(R3Key, usize)>,
}

impl PrimIndex {
    fn of(prims: &TriPrims) -> Self {
        let mut index = Self::default();
        for (a, b, prov) in &prims.segments {
            index.add_segment(a, b, *prov);
        }
        for (pt, prov) in &prims.points {
            index.add_point(pt, *prov);
        }
        index
    }

    /// Records the unordered segment; false when it was already present.
    fn add_segment(&mut self, a: &R3, b: &R3, prov: usize) -> bool {
        let added = self
            .segments
            .insert((R3Key(a.clone()), R3Key(b.clone()), prov));
        if added {
            self.segments
                .insert((R3Key(b.clone()), R3Key(a.clone()), prov));
        }
        added
    }

    /// Records the point; false when it was already present.
    fn add_point(&mut self, pt: &R3, prov: usize) -> bool {
        self.points.insert((R3Key(pt.clone()), prov))
    }
}

/// Phase 3: cross-copy primitives through coplanar overlap regions so both
/// sides see identical geometry inside the shared area, clipped against the
/// region so unrelated geometry is not dragged across. Region by region in
/// the order phase 1 found them. Reports as [`Phase::CoplanarOverlaps`],
/// counted in regions, and only when there are regions to copy through.
/// Returns false when cancelled.
pub(super) fn cross_copy_coplanar_regions(
    prims: &mut [Vec<TriPrims>; 2],
    coplanar_regions: &[(usize, usize, Vec<R3>)],
    token: Option<&CancelToken>,
    progress: Option<&ProgressReporter>,
) -> bool {
    if coplanar_regions.is_empty() {
        return !is_cancelled(token);
    }
    begin_phase(
        progress,
        Phase::CoplanarOverlaps,
        coplanar_regions.len() as u64,
    );
    // Keyed by (mesh, triangle); probe-only.
    let mut indexes: HashMap<(usize, usize), PrimIndex> = HashMap::default();
    for (pi, qi, poly) in coplanar_regions {
        if is_cancelled(token) {
            return false;
        }
        // Both snapshots are taken BEFORE either copy runs: the second copy
        // must not see what the first one just added.
        let region = CoplanarClipRegion::prepare(poly);
        let from_p = prims[0][*pi].clone();
        let from_q = prims[1][*qi].clone();
        for (m, ti, src) in [(1, *qi, &from_p), (0, *pi, &from_q)] {
            let dst = &mut prims[m][ti];
            let index = indexes.entry((m, ti)).or_insert_with(|| PrimIndex::of(dst));
            copy_through_region(src, dst, index, &region);
        }
        if let Some(p) = progress {
            p.advance(1);
        }
    }
    complete_phase(progress);
    true
}

/// Cross-copy one side's primitives into the other's list, clipped to the
/// shared coplanar overlap region.
fn copy_through_region(
    src: &TriPrims,
    dst: &mut TriPrims,
    present: &mut PrimIndex,
    region: &CoplanarClipRegion,
) {
    for (a, b, prov) in &src.segments {
        if let Some((ca, cb)) = region.clip(a, b) {
            if present.add_segment(&ca, &cb, *prov) {
                dst.segments.push((ca, cb, *prov));
            }
        }
    }
    for (pt, prov) in &src.points {
        // Asked of the clip with a zero-length segment first and the
        // containment test second, as before; both answer "inside or on".
        if (region.clip(pt, pt).is_some() || region.contains(pt)) && present.add_point(pt, *prov) {
            dst.points.push((pt.clone(), *prov));
        }
    }
}
