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

        // Port of the C++ CheckEdge lambda: take `edge` as the new connector
        // when the horizontal ray from `start` crosses it (finite x), `start`
        // lies inside THAT edge's wedge, and it beats the current connector —
        // either the crossing point is CCW of the connector edge, or (for any
        // non-CCW result) the vertical-ordering InsideEdge tie-break holds.
        let outers: Vec<usize> = self.outers.clone();
        for outer_start in &outers {
            let verts = match self.loop_verts(*outer_start) {
                None => continue,
                Some(v) => v,
            };
            for &edge in &verts {
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

        let outers: Vec<usize> = self.outers.clone();
        for outer_start in &outers {
            let verts = match self.loop_verts(*outer_start) {
                None => continue,
                Some(v) => v,
            };
            for &vert in &verts {
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
