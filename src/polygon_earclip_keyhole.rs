// polygon_earclip_keyhole.rs — keyholing for the ear clipper: joining each
// hole into the outer ring its bridge reaches (cut_keyhole, find_closer_bridge,
// join_polygons), and the ring walk (loop_verts, for_each_loop_vert) that the
// bridge searches run over every outer ring.
//
// Ports the CutKeyhole / FindCloserBridge / JoinPolygons portion of the
// EarClip class in src/polygon.cpp, plus its `Loop`. Extracted from
// polygon_earclip.rs, which owns the Vert list, the predicates these searches
// use (vert_interp_y2x, vert_inside_edge, vert_is_reflex) and the ear-clipping
// loop. A child module, so it keeps access to EarClip's private fields; the
// rest of polygon_earclip.rs also walks rings through loop_verts.

use crate::linalg::Vec2;

use super::{ccw, EarClip, INVALID};

impl EarClip {
    /// The unclipped verts of the ring starting from `first`, or `None` if
    /// the ring is degenerate.
    pub(super) fn loop_verts(&self, first: usize) -> Option<Vec<usize>> {
        let mut result = Vec::new();
        self.for_each_loop_vert(first, |v| result.push(v))
            .then_some(result)
    }

    /// Apply `f` to each vert `loop_verts` would return, without collecting
    /// them, as C++ `Loop` does. Returns `false` if the ring is degenerate.
    /// A ring degenerates whole, to two verts, so the walk reports that at
    /// the first vert it lands on, before `f` has seen any; the bridge
    /// searches still restore their state if it does stop part-way, so that
    /// they skip the ring whole either way, as `loop_verts` returning `None`
    /// made them.
    pub(super) fn for_each_loop_vert(&self, first: usize, mut f: impl FnMut(usize)) -> bool {
        let mut v = first;
        let mut cur_first = first;
        loop {
            if self.clipped(v) {
                cur_first = self.polygon[self.polygon[v].right].left;
                if !self.clipped(cur_first) {
                    v = cur_first;
                    if self.polygon[v].right == self.polygon[v].left {
                        return false;
                    }
                    f(v);
                }
            } else {
                if self.polygon[v].right == self.polygon[v].left {
                    return false;
                }
                f(v);
            }
            v = self.polygon[v].right;
            if v == cur_first {
                break;
            }
        }
        true
    }

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

        // Port of the C++ CheckEdge lambda: take `edge` as the new connector
        // when the horizontal ray from `start` crosses it (finite x), `start`
        // lies inside THAT edge's wedge, and it beats the current connector —
        // either the crossing point is CCW of the connector edge, or (for any
        // non-CCW result) the vertical-ordering InsideEdge tie-break holds.
        // A degenerate ring is skipped whole, as `loop_verts` returning `None`
        // skipped it, so the connector is restored if the walk stops part-way.
        for &outer_start in &self.outers {
            let before = connector;
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
                }
            });
            if !complete {
                connector = before;
            }
        }

        if connector == INVALID {
            self.simples.push(start);
            return;
        }

        connector = self.find_closer_bridge(start, connector);
        self.join_polygons(start, connector);
    }

    /// Refine keyhole connector: find any reflex vert closer to start.
    fn find_closer_bridge(&self, start: usize, edge: usize) -> usize {
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
            return connector;
        }
        let above: f64 = if self.polygon[connector].pos.y > start_pos.y {
            1.0
        } else {
            -1.0
        };

        // Degenerate rings are skipped whole, as in `cut_keyhole`.
        for &outer_start in &self.outers {
            let before = connector;
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
                }
            });
            if !complete {
                connector = before;
            }
        }

        connector
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
