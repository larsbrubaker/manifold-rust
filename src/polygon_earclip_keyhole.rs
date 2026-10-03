// EarClip keyholing — the hole-to-outer bridge searches of the ear clipper
//
// Ports C++ `EarClip::CutKeyhole`, `FindCloserBridge` and `JoinPolygons`
// (src/polygon.cpp). Before ear clipping, `EarClip::triangulate` (in
// polygon_earclip.rs, which defines the struct, the vert predicates and the
// ear-clipping loop) joins each hole, rightmost first, to an outer ring by a
// zero-width bridge so every polygon left is simple. This file adds a second
// `impl EarClip` block holding just those three steps.

use crate::linalg::Vec2;

use super::super::{ccw, INVALID};
use super::EarClip;

impl EarClip {
    /// Attach a hole to an outer polygon via a keyhole.
    pub(super) fn cut_keyhole(&mut self, start: usize) {
        let bbox = *self.hole2bbox.get(&start).unwrap();
        let start_pos = self.polygon[start].pos;
        let on_top: i32 = if start_pos.y >= bbox.max.y - self.epsilon {
            1
        } else if start_pos.y <= bbox.min.y + self.epsilon {
            -1
        } else {
            0
        };
        let mut connector: usize = INVALID;
        let mut ring: usize = INVALID;
        // Speed only, not in C++: a ring wholly above or below start.y -+ eps
        // has no edge with a finite `vert_interp_y2x`, which needs one end at
        // or below start.y + eps and the other at or above start.y - eps, so
        // it cannot take the connector and is not walked. The margin only
        // widens the band: twice eps, plus 1e-9 relative, far above the
        // rounding of the comparisons. A NaN or overflowed slack culls nothing.
        let slack = 2.0 * self.epsilon.abs() + 1e-9 * (1.0 + start_pos.y.abs());

        // Port of the C++ CheckEdge lambda: take `edge` as the new connector
        // when the horizontal ray from `start` crosses it (finite x), `start`
        // lies inside THAT edge's wedge, and it beats the current connector —
        // either the crossing point is CCW of the connector edge, or (for any
        // non-CCW result) the vertical-ordering InsideEdge tie-break holds.
        // A degenerate ring is skipped whole, as `loop_verts` returning `None`
        // skipped it, so the connector is restored if the walk stops part-way.
        // C++ `Loop` keeps what it saw before the degenerate vert; see
        // docs/CPP_DIVERGENCES.md entry 12.
        for (k, &outer_start) in self.outers.iter().enumerate() {
            let rb = &self.outer_bbox[k];
            if rb.min.y > start_pos.y + slack || rb.max.y < start_pos.y - slack {
                continue;
            }
            let before = (connector, ring);
            let complete = self.for_each_loop_vert(outer_start, |edge| {
                let x = self.vert_interp_y2x(edge, start_pos, on_top);
                if x.is_finite()
                    && self.vert_inside_edge(start, edge, true)
                    && (connector == INVALID
                        || ccw(
                            Vec2::new(x, start_pos.y),
                            self.polygon[connector].pos,
                            self.polygon[self.polygon[connector].right].pos,
                            self.epsilon,
                        ) == 1
                        || (if self.polygon[connector].pos.y < self.polygon[edge].pos.y {
                            self.vert_inside_edge(edge, connector, false)
                        } else {
                            !self.vert_inside_edge(connector, edge, false)
                        }))
                {
                    connector = edge;
                    ring = k;
                }
            });
            if !complete {
                (connector, ring) = before;
            }
        }

        if connector == INVALID {
            self.simples.push(start);
            return;
        }

        let (connector, ring) = self.find_closer_bridge(start, connector, ring);
        self.join_polygons(start, connector);
        // The hole's verts are now part of that ring. The joined verts are
        // copies of verts already in one box or the other, and clipping only
        // removes verts, so the box still holds every vert the ring walks.
        let rb = &mut self.outer_bbox[ring];
        rb.union_point(bbox.min);
        rb.union_point(bbox.max);
    }

    /// Refine keyhole connector: find any reflex vert closer to start.
    /// Also returns the `outers` index of the connector's ring, which starts
    /// as `edge_ring`, the ring `edge` came from.
    fn find_closer_bridge(&self, start: usize, edge: usize, edge_ring: usize) -> (usize, usize) {
        let start_pos = self.polygon[start].pos;
        let edge_right = self.polygon[edge].right;
        let mut connector = if self.polygon[edge].pos.x < start_pos.x {
            edge_right
        } else if self.polygon[edge_right].pos.x < start_pos.x {
            edge
        } else if self.polygon[edge_right].pos.y - start_pos.y
            > start_pos.y - self.polygon[edge].pos.y
        {
            edge
        } else {
            edge_right
        };

        if (self.polygon[connector].pos.y - start_pos.y).abs() <= self.epsilon {
            return (connector, edge_ring);
        }
        let above: f64 = if self.polygon[connector].pos.y > start_pos.y {
            1.0
        } else {
            -1.0
        };

        // Degenerate rings are skipped whole, as in `cut_keyhole`.
        //
        // Speed only, not in C++: a ring whose bounding box shows that no
        // vert in it can pass the test below is not walked, so the bridge,
        // and the triangles, are unchanged. A vert must lie right of
        // start.x - eps, on the `above` side of start.y -+ eps, and not
        // clearly outside start -> connector. The coordinate tests are bounded
        // by the box's edges, widened by twice eps plus 1e-9 relative, far
        // above the rounding of the comparisons; they need no guard, as a
        // NaN or overflowed slack culls nothing.
        let eps = self.epsilon.abs();
        let slack = 2.0 * eps + 1e-9 * (1.0 + start_pos.x.abs() + start_pos.y.abs());
        let mut ring = edge_ring;
        for (k, &outer_start) in self.outers.iter().enumerate() {
            let rb = &self.outer_bbox[k];
            if rb.max.x < start_pos.x - slack
                || (above > 0.0 && rb.max.y < start_pos.y - slack)
                || (above < 0.0 && rb.min.y > start_pos.y + slack)
            {
                continue;
            }
            // The `ccw` test rejects a vert only when ccw(start, vert,
            // connector) is nonzero with sign -above, that is when its area
            // a = v1 x v2 (v1 = vert - start, v2 = connector - start) has
            // above * a < 0 and 4a^2 > max(|v1|^2, |v2|^2) eps^2. above * a
            // is linear in the vert, so its maximum over the box is at a
            // corner (`best`); `dist` bounds |v1| and |v2| over the box. In
            // exact arithmetic best < -(dist eps) rejects every vert in the
            // box. The 1e-9 dist |v2| margin covers rounding: each difference
            // and product is rounded to ~1e-16 relative, so the computed a of
            // a culled vert keeps its sign and |a| > dist eps.
            //
            // That argument needs the squares in `ccw` to be finite normal
            // numbers, so the cull applies only inside a magnitude window.
            // There a culled vert has |a| > 1e-9 dist |v2| >= 1e-9 |v2|^2 >=
            // 1e-129, so 4a^2 >= 4e-258 does not underflow, and |a| <= dist
            // |v2| <= 1e120 and max(|v1|^2, |v2|^2) eps^2 <= 1e240 do not
            // overflow. Outside it, a^2 can underflow to 0 (or both sides
            // overflow to inf) and `ccw` returns 0 for a vert the bound
            // calls outside, which the `inside == 0` tie-break can then take;
            // see `keyhole_cull_keeps_a_bridge_whose_ccw_underflows`.
            let v2 = self.polygon[connector].pos - start_pos;
            let v2_len = (v2.x * v2.x + v2.y * v2.y).sqrt();
            let mut best = f64::NEG_INFINITY;
            let mut dist = v2_len;
            for c in [
                Vec2::new(rb.min.x, rb.min.y),
                Vec2::new(rb.max.x, rb.min.y),
                Vec2::new(rb.min.x, rb.max.y),
                Vec2::new(rb.max.x, rb.max.y),
            ] {
                let v1 = c - start_pos;
                best = best.max(above * (v1.x * v2.y - v1.y * v2.x));
                dist = dist.max((v1.x * v1.x + v1.y * v1.y).sqrt());
            }
            // False on NaN.
            let in_window = v2_len >= 1e-60 && dist <= 1e60 && eps <= 1e60;
            if in_window && best < -(dist * eps + 1e-9 * dist * v2_len) {
                continue;
            }
            let before = (connector, ring);
            let complete = self.for_each_loop_vert(outer_start, |vert| {
                let inside = above
                    * ccw(
                        start_pos,
                        self.polygon[vert].pos,
                        self.polygon[connector].pos,
                        self.epsilon,
                    ) as f64;
                let vp = self.polygon[vert].pos;
                let cp = self.polygon[connector].pos;
                if vp.x > start_pos.x - self.epsilon
                    && vp.y * above > start_pos.y * above - self.epsilon
                    && (inside > 0.0
                        || (inside == 0.0 && vp.x < cp.x && vp.y * above < cp.y * above))
                    && self.vert_inside_edge(vert, edge, true)
                    && self.vert_is_reflex(vert)
                {
                    connector = vert;
                    ring = k;
                }
            });
            if !complete {
                (connector, ring) = before;
            }
        }

        (connector, ring)
    }

    /// Create a keyhole between hole `start` and outer polygon `connector`.
    fn join_polygons(&mut self, start: usize, connector: usize) {
        let new_start = self.polygon.len();
        self.polygon.push(self.polygon[start].clone());
        let new_connector = self.polygon.len();
        self.polygon.push(self.polygon[connector].clone());

        let start_right = self.polygon[start].right;
        self.polygon[start_right].left = new_start;
        let connector_left = self.polygon[connector].left;
        self.polygon[connector_left].right = new_connector;

        self.link(start, connector);
        self.link(new_connector, new_start);

        self.clip_if_degenerate(start);
        self.clip_if_degenerate(new_start);
        self.clip_if_degenerate(connector);
        self.clip_if_degenerate(new_connector);
    }
}
