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

// Phase 15: Subdivision — ported from C++ subdivision.cpp (811 lines)
//
// Implements the full subdivision system with:
// - Partition class with cached triangulations
// - Triangle and quad subdivision
// - Barycentric interpolation for vertices and properties
// - Edge-division-based refinement

use std::collections::HashMap;
use std::sync::Mutex;

use crate::impl_mesh::ManifoldImpl;
use crate::linalg::{BVec4, IVec3, IVec4, Mat3, Mat3x4, Vec3, Vec4};
use crate::types::{next_halfedge, Barycentric, Halfedge, TmpEdge, TriRef};

// ---------------------------------------------------------------------------
// BaryIndices — maps halfedge to triangle + 4D index pair
// Port of C++ Manifold::Impl::BaryIndices
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug)]
struct BaryIndices {
    tri: i32,
    start4: i32,
    end4: i32,
}

// ---------------------------------------------------------------------------
// Partition — cached topological triangulation
// Port of C++ Partition class (lines 31-380)
// ---------------------------------------------------------------------------

#[derive(Clone, Debug)]
struct Partition {
    idx: IVec4,
    sorted_divisions: IVec4,
    vert_bary: Vec<Vec4>,
    tri_vert: Vec<IVec3>,
}

impl Partition {
    fn new() -> Self {
        Self {
            idx: IVec4::default(),
            sorted_divisions: IVec4::default(),
            vert_bary: Vec::new(),
            tri_vert: Vec::new(),
        }
    }

    fn interior_offset(&self) -> i32 {
        self.sorted_divisions[0] + self.sorted_divisions[1]
            + self.sorted_divisions[2] + self.sorted_divisions[3]
    }

    fn num_interior(&self) -> i32 {
        self.vert_bary.len() as i32 - self.interior_offset()
    }

    /// Port of C++ Partition::GetPartition(ivec4 divisions)
    fn get_partition(divisions: IVec4) -> Partition {
        if divisions[0] == 0 {
            return Partition::new(); // skip wrong side of quad
        }

        let mut sorted_div = divisions;
        let mut tri_idx = IVec4::new(0, 1, 2, 3);

        if divisions[3] == 0 {
            // triangle: sort descending
            if sorted_div[2] > sorted_div[1] {
                swap_i32(&mut sorted_div, &mut tri_idx, 2, 1);
            }
            if sorted_div[1] > sorted_div[0] {
                swap_i32(&mut sorted_div, &mut tri_idx, 1, 0);
                if sorted_div[2] > sorted_div[1] {
                    swap_i32(&mut sorted_div, &mut tri_idx, 2, 1);
                }
            }
        } else {
            // quad: rotate to canonical form
            let mut min_idx = 0;
            let mut min = divisions[0];
            let mut next = divisions[1];
            for i in 1..4 {
                let n = divisions[(i + 1) % 4];
                if divisions[i] < min || (divisions[i] == min && n < next) {
                    min_idx = i;
                    min = divisions[i];
                    next = n;
                }
            }
            let tmp = sorted_div;
            for i in 0..4 {
                tri_idx[i] = ((i + min_idx) % 4) as i32;
                sorted_div[i] = tmp[tri_idx[i] as usize];
            }
        }

        let mut partition = Self::get_cached_partition(sorted_div);
        partition.idx = tri_idx;
        partition
    }

    /// Port of C++ Partition::Reindex()
    fn reindex(
        &self,
        tri_verts: IVec4,
        edge_offsets: IVec4,
        mut edge_fwd: BVec4,
        interior_offset: i32,
    ) -> Vec<IVec3> {
        let mut new_verts: Vec<i32> = Vec::with_capacity(self.vert_bary.len());
        let mut tri_idx = self.idx;
        let mut out_tri = IVec4::new(0, 1, 2, 3);

        if tri_verts[3] < 0 && self.idx[1] != next3(self.idx[0]) {
            tri_idx = IVec4::new(self.idx[2], self.idx[0], self.idx[1], self.idx[3]);
            edge_fwd = !edge_fwd;
            // swap outTri[0] and outTri[1]
            let tmp = out_tri[0];
            out_tri[0] = out_tri[1];
            out_tri[1] = tmp;
        }

        for i in 0..4 {
            if tri_verts[tri_idx[i] as usize] >= 0 {
                new_verts.push(tri_verts[tri_idx[i] as usize]);
            }
        }

        for i in 0..4 {
            let n = self.sorted_divisions[i] - 1;
            let mut offset = edge_offsets[self.idx[i] as usize]
                + if edge_fwd[self.idx[i] as usize] { 0 } else { n - 1 };
            for _ in 0..n {
                new_verts.push(offset);
                offset += if edge_fwd[self.idx[i] as usize] { 1 } else { -1 };
            }
        }

        let offset = interior_offset - new_verts.len() as i32;
        let old = new_verts.len();
        new_verts.resize(self.vert_bary.len(), 0);
        for (j, slot) in new_verts[old..].iter_mut().enumerate() {
            *slot = (old + j) as i32 + offset;
        }

        let num_tri = self.tri_vert.len();
        let mut new_tri_vert = vec![IVec3::default(); num_tri];
        for tri in 0..num_tri {
            for j in 0..3 {
                new_tri_vert[tri][out_tri[j] as usize] =
                    new_verts[self.tri_vert[tri][j] as usize];
            }
        }
        new_tri_vert
    }

    /// Port of C++ Partition::GetCachedPartition(ivec4 n)
    fn get_cached_partition(n: IVec4) -> Partition {
        // Check cache
        {
            let cache = partition_cache().lock().unwrap();
            if let Some(cached) = cache.get(&(n.x, n.y, n.z, n.w)) {
                return cached.clone();
            }
        }

        let mut partition = Partition::new();
        partition.sorted_divisions = n;

        if n[3] > 0 {
            // quad
            partition.vert_bary.push(Vec4::new(1.0, 0.0, 0.0, 0.0));
            partition.vert_bary.push(Vec4::new(0.0, 1.0, 0.0, 0.0));
            partition.vert_bary.push(Vec4::new(0.0, 0.0, 1.0, 0.0));
            partition.vert_bary.push(Vec4::new(0.0, 0.0, 0.0, 1.0));

            let mut edge_offsets = IVec4::default();
            edge_offsets[0] = 4;
            for i in 0..4 {
                if i > 0 {
                    edge_offsets[i] = edge_offsets[i - 1] + n[i - 1] - 1;
                }
                let next_bary = partition.vert_bary[(i + 1) % 4];
                for j in 1..n[i] {
                    partition.vert_bary.push(lerp_vec4(
                        partition.vert_bary[i],
                        next_bary,
                        j as f64 / n[i] as f64,
                    ));
                }
            }
            let n_minus_1 = IVec4::new(n[0] - 1, n[1] - 1, n[2] - 1, n[3] - 1);
            partition_quad(
                &mut partition.tri_vert,
                &mut partition.vert_bary,
                IVec4::new(0, 1, 2, 3),
                edge_offsets,
                n_minus_1,
                BVec4::splat(true),
            );
        } else {
            // triangle
            partition.vert_bary.push(Vec4::new(1.0, 0.0, 0.0, 0.0));
            partition.vert_bary.push(Vec4::new(0.0, 1.0, 0.0, 0.0));
            partition.vert_bary.push(Vec4::new(0.0, 0.0, 1.0, 0.0));

            for i in 0..3 {
                let next_bary = partition.vert_bary[(i + 1) % 3];
                for j in 1..n[i] {
                    partition.vert_bary.push(lerp_vec4(
                        partition.vert_bary[i],
                        next_bary,
                        j as f64 / n[i] as f64,
                    ));
                }
            }
            let edge_offsets = IVec4::new(3, 3 + n[0] - 1, 3 + n[0] - 1 + n[1] - 1, 0);

            let f = (n[2] * n[2] + n[0] * n[0]) as f64;
            if n[1] == 1 {
                if n[0] == 1 {
                    partition.tri_vert.push(IVec3::new(0, 1, 2));
                } else {
                    partition_fan(
                        &mut partition.tri_vert,
                        IVec3::new(0, 1, 2),
                        n[0] - 1,
                        edge_offsets[0],
                    );
                }
            } else if (n[1] * n[1]) as f64 > f - (2.0_f64).sqrt() * n[0] as f64 * n[2] as f64 {
                // acute-ish
                partition.tri_vert.push(IVec3::new(
                    edge_offsets[1] - 1,
                    1,
                    edge_offsets[1],
                ));
                partition_quad(
                    &mut partition.tri_vert,
                    &mut partition.vert_bary,
                    IVec4::new(edge_offsets[1] - 1, edge_offsets[1], 2, 0),
                    IVec4::new(-1, edge_offsets[1] + 1, edge_offsets[2], edge_offsets[0]),
                    IVec4::new(0, n[1] - 2, n[2] - 1, n[0] - 2),
                    BVec4::splat(true),
                );
            } else {
                // obtuse -> split into two acute
                let ns = (n[0] - 2).min(
                    ((f - (n[1] * n[1]) as f64) / (2.0 * n[0] as f64)).round() as i32,
                );
                let nh = 1.0_f64
                    .max(((n[2] * n[2] - ns * ns) as f64).sqrt().round())
                    as i32;

                let h_offset = partition.vert_bary.len() as i32;
                let middle_bary = partition.vert_bary[(edge_offsets[0] + ns - 1) as usize];
                for j in 1..nh {
                    partition.vert_bary.push(lerp_vec4(
                        partition.vert_bary[2],
                        middle_bary,
                        j as f64 / nh as f64,
                    ));
                }

                partition.tri_vert.push(IVec3::new(
                    edge_offsets[1] - 1,
                    1,
                    edge_offsets[1],
                ));
                partition_quad(
                    &mut partition.tri_vert,
                    &mut partition.vert_bary,
                    IVec4::new(
                        edge_offsets[1] - 1,
                        edge_offsets[1],
                        2,
                        edge_offsets[0] + ns - 1,
                    ),
                    IVec4::new(-1, edge_offsets[1] + 1, h_offset, edge_offsets[0] + ns),
                    IVec4::new(0, n[1] - 2, nh - 1, n[0] - ns - 2),
                    BVec4::splat(true),
                );

                if n[2] == 1 {
                    partition_fan(
                        &mut partition.tri_vert,
                        IVec3::new(0, edge_offsets[0] + ns - 1, 2),
                        ns - 1,
                        edge_offsets[0],
                    );
                } else if ns == 1 {
                    partition.tri_vert.push(IVec3::new(
                        h_offset,
                        2,
                        edge_offsets[2],
                    ));
                    partition_quad(
                        &mut partition.tri_vert,
                        &mut partition.vert_bary,
                        IVec4::new(h_offset, edge_offsets[2], 0, edge_offsets[0]),
                        IVec4::new(-1, edge_offsets[2] + 1, -1, h_offset + nh - 2),
                        IVec4::new(0, n[2] - 2, ns - 1, nh - 2),
                        BVec4::new(true, true, true, false),
                    );
                } else {
                    partition.tri_vert.push(IVec3::new(
                        h_offset - 1,
                        0,
                        edge_offsets[0],
                    ));
                    partition_quad(
                        &mut partition.tri_vert,
                        &mut partition.vert_bary,
                        IVec4::new(
                            h_offset - 1,
                            edge_offsets[0],
                            edge_offsets[0] + ns - 1,
                            2,
                        ),
                        IVec4::new(
                            -1,
                            edge_offsets[0] + 1,
                            h_offset + nh - 2,
                            edge_offsets[2],
                        ),
                        IVec4::new(0, ns - 2, nh - 1, n[2] - 2),
                        BVec4::new(true, true, false, true),
                    );
                }
            }
        }

        // Store in cache
        let mut cache = partition_cache().lock().unwrap();
        cache.insert((n.x, n.y, n.z, n.w), partition.clone());
        partition
    }
}

// Global partition cache
use std::sync::OnceLock;
fn partition_cache() -> &'static Mutex<HashMap<(i32, i32, i32, i32), Partition>> {
    static CACHE: OnceLock<Mutex<HashMap<(i32, i32, i32, i32), Partition>>> = OnceLock::new();
    CACHE.get_or_init(|| Mutex::new(HashMap::new()))
}

// ---------------------------------------------------------------------------
// Helper functions
// ---------------------------------------------------------------------------

#[inline]
fn next3(i: i32) -> i32 {
    if i == 2 { 0 } else { i + 1 }
}

#[inline]
fn lerp_vec4(a: Vec4, b: Vec4, t: f64) -> Vec4 {
    Vec4::new(
        a.x + (b.x - a.x) * t,
        a.y + (b.y - a.y) * t,
        a.z + (b.z - a.z) * t,
        a.w + (b.w - a.w) * t,
    )
}

fn swap_i32(sorted_div: &mut IVec4, tri_idx: &mut IVec4, a: usize, b: usize) {
    let tmp = sorted_div[a];
    sorted_div[a] = sorted_div[b];
    sorted_div[b] = tmp;
    let tmp = tri_idx[a];
    tri_idx[a] = tri_idx[b];
    tri_idx[b] = tmp;
}

/// Port of C++ Partition::PartitionFan()
fn partition_fan(tri_vert: &mut Vec<IVec3>, corner_verts: IVec3, added: i32, edge_offset: i32) {
    let mut last = corner_verts[0];
    for i in 0..added {
        let next = edge_offset + i;
        tri_vert.push(IVec3::new(last, next, corner_verts[2]));
        last = next;
    }
    tri_vert.push(IVec3::new(last, corner_verts[1], corner_verts[2]));
}

/// Port of C++ Partition::PartitionQuad()
fn partition_quad(
    tri_vert: &mut Vec<IVec3>,
    vert_bary: &mut Vec<Vec4>,
    corner_verts: IVec4,
    edge_offsets: IVec4,
    edge_added: IVec4,
    edge_fwd: BVec4,
) {
    let get_edge_vert = |edge: usize, idx: i32| -> i32 {
        edge_offsets[edge] + if edge_fwd[edge] { 1 } else { -1 } * idx
    };

    debug_assert!(
        edge_added[0] >= 0 && edge_added[1] >= 0 && edge_added[2] >= 0 && edge_added[3] >= 0,
        "negative divisions!"
    );

    let mut corner: i32 = -1;
    let mut last = 3_usize;
    let mut max_edge: i32 = -1;
    for i in 0..4 {
        if corner == -1 && edge_added[i] == 0 && edge_added[last] == 0 {
            corner = i as i32;
        }
        if edge_added[i] > 0 {
            max_edge = if max_edge == -1 { i as i32 } else { -2 };
        }
        last = i;
    }

    if corner >= 0 {
        // terminate
        if max_edge >= 0 {
            let mut edge = IVec4::default();
            for j in 0..4 {
                edge[j] = ((j as i32 + max_edge) % 4) as i32;
            }
            let middle = edge_added[max_edge as usize] / 2;
            tri_vert.push(IVec3::new(
                corner_verts[edge[2] as usize],
                corner_verts[edge[3] as usize],
                get_edge_vert(max_edge as usize, middle),
            ));
            let mut last_v = corner_verts[edge[0] as usize];
            for i in 0..=middle {
                let next = get_edge_vert(max_edge as usize, i);
                tri_vert.push(IVec3::new(
                    corner_verts[edge[3] as usize],
                    last_v,
                    next,
                ));
                last_v = next;
            }
            last_v = corner_verts[edge[1] as usize];
            for i in (middle..edge_added[max_edge as usize]).rev() {
                let next = get_edge_vert(max_edge as usize, i);
                tri_vert.push(IVec3::new(
                    corner_verts[edge[2] as usize],
                    next,
                    last_v,
                ));
                last_v = next;
            }
        } else {
            let corner_u = corner as usize;
            let mut side_vert = corner_verts[0]; // initial value unused
            for j in 1..=2 {
                let side = (corner_u + j) % 4;
                if j == 2 && edge_added[side] > 0 {
                    tri_vert.push(IVec3::new(
                        corner_verts[side],
                        get_edge_vert(side, 0),
                        side_vert,
                    ));
                } else {
                    side_vert = corner_verts[side];
                }
                for i in 0..edge_added[side] {
                    let next_vert = get_edge_vert(side, i);
                    tri_vert.push(IVec3::new(
                        corner_verts[corner_u],
                        side_vert,
                        next_vert,
                    ));
                    side_vert = next_vert;
                }
                if j == 2 || edge_added[side] == 0 {
                    tri_vert.push(IVec3::new(
                        corner_verts[corner_u],
                        side_vert,
                        corner_verts[(corner_u + j + 1) % 4],
                    ));
                }
            }
        }
        return;
    }

    // recursively partition
    let partitions = 1 + edge_added[1].min(edge_added[3]);
    let mut new_corner_verts = IVec4::new(corner_verts[1], -1, -1, corner_verts[0]);
    let mut new_edge_offsets = IVec4::new(
        edge_offsets[1],
        -1,
        get_edge_vert(3, edge_added[3] + 1),
        edge_offsets[0],
    );
    let mut new_edge_added = IVec4::new(0, -1, 0, edge_added[0]);
    let mut new_edge_fwd = BVec4::new(edge_fwd[1], true, edge_fwd[3], edge_fwd[0]);

    for i in 1..partitions {
        let corner_offset1 = (edge_added[1] * i) / partitions;
        let corner_offset3 = edge_added[3] - 1 - (edge_added[3] * i) / partitions;
        let next_offset1 = get_edge_vert(1, corner_offset1 + 1);
        let next_offset3 = get_edge_vert(3, corner_offset3 + 1);
        let added = ((edge_added[0] as f64)
            + (edge_added[2] as f64 - edge_added[0] as f64) * (i as f64 / partitions as f64))
            .round() as i32;

        new_corner_verts[1] = get_edge_vert(1, corner_offset1);
        new_corner_verts[2] = get_edge_vert(3, corner_offset3);
        new_edge_added[0] = (next_offset1 - new_edge_offsets[0]).abs() - 1;
        new_edge_added[1] = added;
        new_edge_added[2] = (next_offset3 - new_edge_offsets[2]).abs() - 1;
        new_edge_offsets[1] = vert_bary.len() as i32;
        new_edge_offsets[2] = next_offset3;

        for j in 0..added {
            vert_bary.push(lerp_vec4(
                vert_bary[new_corner_verts[1] as usize],
                vert_bary[new_corner_verts[2] as usize],
                (j + 1) as f64 / (added + 1) as f64,
            ));
        }

        partition_quad(
            tri_vert,
            vert_bary,
            new_corner_verts,
            new_edge_offsets,
            new_edge_added,
            new_edge_fwd,
        );

        new_corner_verts[0] = new_corner_verts[1];
        new_corner_verts[3] = new_corner_verts[2];
        new_edge_added[3] = new_edge_added[1];
        new_edge_offsets[0] = next_offset1;
        new_edge_offsets[3] = new_edge_offsets[1] + new_edge_added[1] - 1;
        new_edge_fwd[3] = false;
    }

    new_corner_verts[1] = corner_verts[2];
    new_corner_verts[2] = corner_verts[3];
    new_edge_offsets[1] = edge_offsets[2];
    new_edge_added[0] = edge_added[1] - (new_edge_offsets[0] - edge_offsets[1]).abs();
    new_edge_added[1] = edge_added[2];
    new_edge_added[2] = (new_edge_offsets[2] - edge_offsets[3]).abs() - 1;
    new_edge_offsets[2] = edge_offsets[3];
    new_edge_fwd[1] = edge_fwd[2];

    partition_quad(
        tri_vert,
        vert_bary,
        new_corner_verts,
        new_edge_offsets,
        new_edge_added,
        new_edge_fwd,
    );
}

// ---------------------------------------------------------------------------
// CreateTmpEdges — port of C++ inline CreateTmpEdges()
// ---------------------------------------------------------------------------

fn create_tmp_edges(halfedge: &[Halfedge]) -> Vec<TmpEdge> {
    let mut edges: Vec<TmpEdge> = Vec::with_capacity(halfedge.len());
    for (idx, half) in halfedge.iter().enumerate() {
        if half.is_forward() {
            edges.push(TmpEdge::new(half.start_vert, half.end_vert, idx as i32));
        }
    }
    debug_assert!(
        edges.len() == halfedge.len() / 2,
        "Not oriented! edges={} halfedges={}",
        edges.len(),
        halfedge.len()
    );
    edges
}

// ---------------------------------------------------------------------------
// ManifoldImpl methods for subdivision
// ---------------------------------------------------------------------------

impl ManifoldImpl {
    /// Port of C++ Manifold::Impl::GetNeighbor(int tri)
    fn get_neighbor(&self, tri: i32) -> i32 {
        let mut neighbor: i32 = -1;
        for i in 0..3 {
            if self.is_marked_inside_quad((3 * tri + i) as usize) {
                neighbor = if neighbor == -1 { i } else { -2 };
            }
        }
        neighbor
    }

    /// Port of C++ Manifold::Impl::GetHalfedges(int tri)
    fn get_halfedges_quad(&self, tri: i32) -> IVec4 {
        let mut halfedges = IVec4::new(-1, -1, -1, -1);
        for i in 0..3 {
            halfedges[i] = 3 * tri + i as i32;
        }
        let neighbor = self.get_neighbor(tri);
        if neighbor >= 0 {
            // quad
            let pair = self.halfedge[(3 * tri + neighbor) as usize].paired_halfedge;
            if pair / 3 < tri {
                return IVec4::new(-1, -1, -1, -1); // only process lower tri index
            }
            halfedges[2] = next_halfedge(halfedges[neighbor as usize]);
            halfedges[3] = next_halfedge(halfedges[2]);
            halfedges[0] = next_halfedge(pair);
            halfedges[1] = next_halfedge(halfedges[0]);
        }
        halfedges
    }

    /// Port of C++ Manifold::Impl::GetIndices(int halfedge)
    fn get_indices(&self, halfedge: i32) -> BaryIndices {
        let mut tri = halfedge / 3;
        let mut idx = halfedge % 3;
        let neighbor = self.get_neighbor(tri);
        if idx == neighbor {
            return BaryIndices { tri: -1, start4: -1, end4: -1 };
        }

        if neighbor < 0 {
            // tri
            BaryIndices { tri, start4: idx, end4: next3(idx) }
        } else {
            // quad
            let pair = self.halfedge[(3 * tri + neighbor) as usize].paired_halfedge;
            if pair / 3 < tri {
                tri = pair / 3;
                idx = if next3(neighbor) == idx { 0 } else { 1 };
            } else {
                idx = if next3(neighbor) == idx { 2 } else { 3 };
            }
            BaryIndices { tri, start4: idx, end4: (idx + 1) % 4 }
        }
    }

    /// Port of C++ Manifold::Impl::FillRetainedVerts()
    fn fill_retained_verts(&self, vert_bary: &mut [Barycentric]) {
        let num_tri = self.halfedge.len() / 3;
        for tri in 0..num_tri {
            for i in 0..3 {
                let indices = self.get_indices((3 * tri + i) as i32);
                if indices.start4 < 0 {
                    continue; // skip quad interiors
                }
                let mut uvw = Vec4::splat(0.0);
                uvw[indices.start4 as usize] = 1.0;
                vert_bary[self.halfedge[3 * tri + i].start_vert as usize] = Barycentric {
                    tri: indices.tri,
                    uvw,
                };
            }
        }
    }

    /// Port of C++ Manifold::Impl::Subdivide()
    /// edgeDivisions: takes (edge_vec, tangent0, tangent1) → number of new vertices
    pub fn subdivide(
        &mut self,
        edge_divisions: &dyn Fn(Vec3, Vec4, Vec4) -> i32,
        keep_interior: bool,
    ) -> Vec<Barycentric> {
        let edges = create_tmp_edges(&self.halfedge);
        let num_vert = self.num_vert();
        let num_edge = edges.len();
        let num_tri = self.num_tri();

        // Build half2edge mapping
        let mut half2edge = vec![0i32; 2 * num_edge];
        for (edge, tmp) in edges.iter().enumerate() {
            let idx = tmp.halfedge_idx as usize;
            half2edge[idx] = edge as i32;
            half2edge[self.halfedge[idx].paired_halfedge as usize] = edge as i32;
        }

        // Get face halfedges for each triangle
        let face_halfedges: Vec<IVec4> = (0..num_tri)
            .map(|tri| self.get_halfedges_quad(tri as i32))
            .collect();

        // Compute edge divisions
        let mut edge_added = vec![0i32; num_edge];
        for i in 0..num_edge {
            let edge = &edges[i];
            let h_idx = edge.halfedge_idx as usize;
            if self.is_marked_inside_quad(h_idx) {
                edge_added[i] = 0;
                continue;
            }
            let vec = self.vert_pos[edge.first as usize] - self.vert_pos[edge.second as usize];
            let tangent0 = if self.halfedge_tangent.is_empty() {
                Vec4::splat(0.0)
            } else {
                self.halfedge_tangent[h_idx]
            };
            let tangent1 = if self.halfedge_tangent.is_empty() {
                Vec4::splat(0.0)
            } else {
                self.halfedge_tangent[self.halfedge[h_idx].paired_halfedge as usize]
            };
            edge_added[i] = edge_divisions(vec, tangent0, tangent1);
        }

        // Optional: add extra divisions to short edges for interior thickness
        if keep_interior {
            let orig_edge_added = edge_added.clone();
            for i in 0..num_edge {
                let edge = &edges[i];
                let h_idx = edge.halfedge_idx as usize;
                if self.is_marked_inside_quad(h_idx) {
                    continue;
                }

                let this_added = orig_edge_added[i];
                let added_fn = |mut h: i32| -> i32 {
                    let mut longest = 0;
                    let mut total = 0;
                    for _ in 0..3 {
                        let added = orig_edge_added[half2edge[h as usize] as usize];
                        longest = longest.max(added);
                        total += added;
                        h = next_halfedge(h);
                        if self.is_marked_inside_quad(h as usize) {
                            longest = 0;
                            total = 1;
                            break;
                        }
                    }
                    let min_extra = (longest as f64 * 0.2) as i32 + 1;
                    let extra = 2 * longest + min_extra - total;
                    if longest == 0 {
                        return 0;
                    }
                    if extra > 0 {
                        (extra * (longest - this_added)) / longest
                    } else {
                        0
                    }
                };

                let a1 = added_fn(h_idx as i32);
                let a2 = added_fn(self.halfedge[h_idx].paired_halfedge);
                edge_added[i] = orig_edge_added[i] + a1.max(a2);
            }
        }

        // Compute edge offsets (exclusive scan)
        let mut edge_offset = vec![0i32; num_edge];
        let mut acc = num_vert as i32;
        for i in 0..num_edge {
            edge_offset[i] = acc;
            acc += edge_added[i];
        }

        // Allocate vert_bary
        let total_edge_verts = acc - num_vert as i32;
        let mut vert_bary = vec![
            Barycentric {
                tri: 0,
                uvw: Vec4::splat(0.0)
            };
            acc as usize
        ];
        self.fill_retained_verts(&mut vert_bary);

        // Fill edge vertex barycentric coords
        for i in 0..num_edge {
            let n = edge_added[i];
            let offset = edge_offset[i];
            let indices = self.get_indices(edges[i].halfedge_idx);
            if indices.tri < 0 {
                continue; // inside quad
            }
            let frac = 1.0 / (n as f64 + 1.0);
            for j in 0..n {
                let mut uvw = Vec4::splat(0.0);
                uvw[indices.end4 as usize] = (j + 1) as f64 * frac;
                uvw[indices.start4 as usize] = 1.0 - uvw[indices.end4 as usize];
                vert_bary[(offset + j) as usize] = Barycentric { tri: indices.tri, uvw };
            }
        }

        // Generate partitions for each triangle
        let sub_tris: Vec<Partition> = (0..num_tri)
            .map(|tri| {
                let halfedges = face_halfedges[tri];
                let mut divisions = IVec4::default();
                for i in 0..4 {
                    if halfedges[i] >= 0 {
                        divisions[i] = edge_added[half2edge[halfedges[i] as usize] as usize] + 1;
                    }
                }
                Partition::get_partition(divisions)
            })
            .collect();

        // Compute triangle offsets (exclusive scan)
        let mut tri_offset = vec![0i32; num_tri];
        {
            let mut acc = 0i32;
            for tri in 0..num_tri {
                tri_offset[tri] = acc;
                acc += sub_tris[tri].tri_vert.len() as i32;
            }
        }

        // Compute interior vertex offsets (exclusive scan)
        let mut interior_offset = vec![0i32; num_tri];
        {
            let mut acc = vert_bary.len() as i32;
            for tri in 0..num_tri {
                interior_offset[tri] = acc;
                acc += sub_tris[tri].num_interior();
            }
        }

        // Allocate output arrays
        let total_new_tris = if num_tri > 0 {
            tri_offset[num_tri - 1] + sub_tris[num_tri - 1].tri_vert.len() as i32
        } else {
            0
        };
        let total_new_verts = if num_tri > 0 {
            interior_offset[num_tri - 1] + sub_tris[num_tri - 1].num_interior()
        } else {
            vert_bary.len() as i32
        };

        let mut tri_verts = vec![IVec3::default(); total_new_tris as usize];
        vert_bary.resize(
            total_new_verts as usize,
            Barycentric { tri: 0, uvw: Vec4::splat(0.0) },
        );
        let mut tri_ref_out = vec![TriRef::default(); total_new_tris as usize];
        let mut face_normal_out = vec![Vec3::splat(0.0); total_new_tris as usize];

        // Build new triangles
        for tri in 0..num_tri {
            let halfedges = face_halfedges[tri];
            if halfedges[0] < 0 {
                continue;
            }

            let mut tri3 = IVec4::default();
            let mut edge_offs = IVec4::default();
            let mut edge_fwd = BVec4::splat(false);
            for i in 0..4 {
                if halfedges[i] < 0 {
                    tri3[i] = -1;
                    continue;
                }
                let he = &self.halfedge[halfedges[i] as usize];
                tri3[i] = he.start_vert;
                edge_offs[i] = edge_offset[half2edge[halfedges[i] as usize] as usize];
                edge_fwd[i] = he.is_forward();
            }

            let new_tris = sub_tris[tri].reindex(
                tri3,
                edge_offs,
                edge_fwd,
                interior_offset[tri],
            );

            let start = tri_offset[tri] as usize;
            for (j, t) in new_tris.iter().enumerate() {
                tri_verts[start + j] = *t;
                tri_ref_out[start + j] = self.mesh_relation.tri_ref[tri];
                face_normal_out[start + j] = self.face_normal[tri];
            }

            // Map interior barycentric coordinates
            let idx = sub_tris[tri].idx;
            let v_idx = if halfedges[3] >= 0 || idx[1] == next3(idx[0]) {
                idx
            } else {
                IVec4::new(idx[2], idx[0], idx[1], idx[3])
            };
            let mut r_idx = IVec4::default();
            for i in 0..4 {
                r_idx[v_idx[i] as usize] = i as i32;
            }

            let sub_bary = &sub_tris[tri].vert_bary;
            let int_off = sub_tris[tri].interior_offset() as usize;
            for (j, bary) in sub_bary[int_off..].iter().enumerate() {
                vert_bary[interior_offset[tri] as usize + j] = Barycentric {
                    tri: tri as i32,
                    uvw: Vec4::new(
                        bary[r_idx[0] as usize],
                        bary[r_idx[1] as usize],
                        bary[r_idx[2] as usize],
                        bary[r_idx[3] as usize],
                    ),
                };
            }
        }

        self.mesh_relation.tri_ref = tri_ref_out;
        self.face_normal = face_normal_out;

        // Compute new vertex positions
        let mut new_vert_pos = vec![Vec3::splat(0.0); vert_bary.len()];
        for (vert, bary) in vert_bary.iter().enumerate() {
            let halfedges = face_halfedges[bary.tri as usize];
            if halfedges[3] < 0 {
                // triangle
                let tri_pos = Mat3::from_cols(
                    self.vert_pos[self.halfedge[halfedges[0] as usize].start_vert as usize],
                    self.vert_pos[self.halfedge[halfedges[1] as usize].start_vert as usize],
                    self.vert_pos[self.halfedge[halfedges[2] as usize].start_vert as usize],
                );
                new_vert_pos[vert] = tri_pos * bary.uvw.xyz();
            } else {
                // quad
                let quad_pos = Mat3x4::from_cols(
                    self.vert_pos[self.halfedge[halfedges[0] as usize].start_vert as usize],
                    self.vert_pos[self.halfedge[halfedges[1] as usize].start_vert as usize],
                    self.vert_pos[self.halfedge[halfedges[2] as usize].start_vert as usize],
                    self.vert_pos[self.halfedge[halfedges[3] as usize].start_vert as usize],
                );
                new_vert_pos[vert] = quad_pos * bary.uvw;
            }
        }
        self.vert_pos = new_vert_pos;

        // Handle properties
        if self.num_prop > 0 {
            let num_prop_vert = self.num_prop_vert();
            let added_verts = self.num_vert() - num_vert;
            let prop_offset = num_prop_vert as i32 - num_vert as i32;
            let num_prop = self.num_prop as usize;

            // Allocate new property array
            let mut prop =
                vec![0.0f64; num_prop * (num_prop_vert + added_verts + total_edge_verts as usize)];

            // Copy retained prop verts
            for (i, &v) in self.properties.iter().enumerate() {
                prop[i] = v;
            }

            // Copy interior prop verts and forward edge prop verts
            for i in 0..added_verts {
                let vert = num_prop_vert + i;
                let bary = &vert_bary[num_vert + i];
                let halfedges = face_halfedges[bary.tri as usize];

                for p in 0..num_prop {
                    if halfedges[3] < 0 {
                        // triangle
                        let mut tri_prop = Vec3::splat(0.0);
                        for k in 0..3 {
                            tri_prop[k] = self.properties
                                [self.halfedge[3 * bary.tri as usize + k].prop_vert as usize
                                    * num_prop
                                    + p];
                        }
                        prop[vert * num_prop + p] =
                            tri_prop.x * bary.uvw.x + tri_prop.y * bary.uvw.y + tri_prop.z * bary.uvw.z;
                    } else {
                        // quad
                        let mut quad_prop = Vec4::splat(0.0);
                        for k in 0..4 {
                            quad_prop[k] = self.properties
                                [self.halfedge[halfedges[k] as usize].prop_vert as usize
                                    * num_prop
                                    + p];
                        }
                        prop[vert * num_prop + p] = quad_prop.x * bary.uvw.x
                            + quad_prop.y * bary.uvw.y
                            + quad_prop.z * bary.uvw.z
                            + quad_prop.w * bary.uvw.w;
                    }
                }
            }

            // Copy backward edge prop verts
            for i in 0..num_edge {
                let n = edge_added[i];
                let offset = edge_offset[i] as usize + prop_offset as usize + added_verts;
                let frac = 1.0 / (n as f64 + 1.0);
                let halfedge_idx =
                    self.halfedge[edges[i].halfedge_idx as usize].paired_halfedge as usize;
                let prop0 = self.halfedge[halfedge_idx].prop_vert as usize;
                let prop1 =
                    self.halfedge[next_halfedge(halfedge_idx as i32) as usize].prop_vert as usize;
                for j in 0..n as usize {
                    for p in 0..num_prop {
                        let t = (j + 1) as f64 * frac;
                        prop[(offset + j) * num_prop + p] = self.properties[prop0 * num_prop + p]
                            + (self.properties[prop1 * num_prop + p]
                                - self.properties[prop0 * num_prop + p])
                                * t;
                    }
                }
            }

            // Build property triangles
            let mut tri_prop_out = vec![IVec3::default(); total_new_tris as usize];
            for tri in 0..num_tri {
                let halfedges = face_halfedges[tri];
                if halfedges[0] < 0 {
                    continue;
                }

                let mut tri3 = IVec4::default();
                let mut edge_offs = IVec4::default();
                let mut edge_fwd = BVec4::splat(true);
                for i in 0..4 {
                    if halfedges[i] < 0 {
                        tri3[i] = -1;
                        continue;
                    }
                    let he = &self.halfedge[halfedges[i] as usize];
                    tri3[i] = he.prop_vert;
                    edge_offs[i] = edge_offset[half2edge[halfedges[i] as usize] as usize];
                    if !he.is_forward() {
                        let paired = he.paired_halfedge;
                        if self.halfedge[paired as usize].prop_vert
                            != self.halfedge[next_halfedge(halfedges[i]) as usize].prop_vert
                            || self.halfedge[next_halfedge(paired) as usize].prop_vert
                                != he.prop_vert
                        {
                            // edge doesn't match, point to backward edge propverts
                            edge_offs[i] += added_verts as i32;
                        } else {
                            edge_fwd[i] = false;
                        }
                    }
                }

                // Add prop_offset to edge offsets
                let prop_edge_offs = IVec4::new(
                    edge_offs[0] + prop_offset,
                    edge_offs[1] + prop_offset,
                    edge_offs[2] + prop_offset,
                    edge_offs[3] + prop_offset,
                );

                let new_tris = sub_tris[tri].reindex(
                    tri3,
                    prop_edge_offs,
                    edge_fwd,
                    interior_offset[tri] + prop_offset,
                );

                let start = tri_offset[tri] as usize;
                for (j, t) in new_tris.iter().enumerate() {
                    tri_prop_out[start + j] = *t;
                }
            }

            self.properties = prop;
            self.create_halfedges(&tri_prop_out, &tri_verts);
        } else {
            self.create_halfedges(&tri_verts, &[]);
        }

        vert_bary
    }
}

/// Simple midpoint subdivision: each triangle is split into 4 by inserting
/// edge midpoints. This is a convenience wrapper that uses uniform n=2 divisions.
pub fn subdivide_impl(mesh: &ManifoldImpl, levels: usize) -> ManifoldImpl {
    if levels == 0 || mesh.is_empty() {
        return mesh.clone();
    }

    let mut current = mesh.clone();
    for _ in 0..levels {
        current.subdivide(&|_vec, _t0, _t1| 1, false);
        current.calculate_bbox();
        current.set_epsilon(-1.0, false);
    }

    current
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linalg::Mat3x4;

    #[test]
    fn test_subdivide_cube_once() {
        let cube = ManifoldImpl::cube(&Mat3x4::identity());
        let sub = subdivide_impl(&cube, 1);
        assert_eq!(sub.num_tri(), cube.num_tri() * 4);
        assert!(sub.num_vert() > cube.num_vert());
    }

    #[test]
    fn test_partition_single_triangle() {
        // Simplest case: 1 division per edge = no subdivision
        let part = Partition::get_partition(IVec4::new(1, 1, 1, 0));
        assert_eq!(part.tri_vert.len(), 1);
        assert_eq!(part.vert_bary.len(), 3);
    }

    #[test]
    fn test_partition_two_divisions() {
        // 2 divisions on each edge of a triangle
        let part = Partition::get_partition(IVec4::new(2, 2, 2, 0));
        assert_eq!(part.tri_vert.len(), 4); // 4 sub-triangles
    }

    #[test]
    fn test_subdivide_edge_divisions() {
        // Test that the full Subdivide method works with edge_divisions callback
        let mut cube = ManifoldImpl::cube(&Mat3x4::identity());
        let _vert_bary = cube.subdivide(&|_vec, _t0, _t1| 1, false);
        assert_eq!(cube.num_tri(), 12 * 4); // each of 12 tris becomes 4
    }

    #[test]
    fn test_create_tmp_edges() {
        let cube = ManifoldImpl::cube(&Mat3x4::identity());
        let edges = create_tmp_edges(&cube.halfedge);
        // A cube has 12 triangles = 36 halfedges = 18 edges
        assert_eq!(edges.len(), cube.halfedge.len() / 2);
    }

    #[test]
    fn test_partition_asymmetric() {
        // Asymmetric divisions: 3,2,1
        let part = Partition::get_partition(IVec4::new(3, 2, 1, 0));
        assert!(part.tri_vert.len() > 1);
        // All barycentric coordinates should sum to 1
        for bary in &part.vert_bary {
            let sum = bary.x + bary.y + bary.z + bary.w;
            assert!(
                (sum - 1.0).abs() < 1e-10,
                "Barycentric coords should sum to 1, got {}",
                sum
            );
        }
    }

    #[test]
    fn test_subdivide_preserves_manifold() {
        let cube = ManifoldImpl::cube(&Mat3x4::identity());
        let sub = subdivide_impl(&cube, 1);
        // Every halfedge should have a valid pair
        for (i, he) in sub.halfedge.iter().enumerate() {
            assert!(
                he.paired_halfedge >= 0,
                "Halfedge {} has no pair",
                i
            );
            let pair = he.paired_halfedge as usize;
            assert!(
                pair < sub.halfedge.len(),
                "Halfedge {} pair {} out of range",
                i,
                pair
            );
        }
    }
}
