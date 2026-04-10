// Phase 2: Core Types — ported from include/manifold/common.h, include/manifold/polygon.h,
// include/manifold/manifold.h (MeshGLP/Error), src/shared.h (Halfedge, TriRef, Barycentric, TmpEdge)

use std::collections::HashMap;
use crate::linalg::{Vec2, Vec3, Vec4, Mat3x4};

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

pub const K_PI: f64 = std::f64::consts::PI;
pub const K_TWO_PI: f64 = std::f64::consts::TAU;
pub const K_HALF_PI: f64 = std::f64::consts::FRAC_PI_2;
/// Precision used for epsilon calculations relative to bounding-box scale.
pub const K_PRECISION: f64 = 1e-12;

pub const DEFAULT_SEGMENTS: i32 = 0;
pub const DEFAULT_ANGLE: f64 = 10.0;
pub const DEFAULT_LENGTH: f64 = 1.0;

// ---------------------------------------------------------------------------
// Scalar utilities
// ---------------------------------------------------------------------------

#[inline]
pub fn radians(a: f64) -> f64 {
    a * K_PI / 180.0
}

#[inline]
pub fn degrees(a: f64) -> f64 {
    a * 180.0 / K_PI
}

/// Smooth Hermite interpolation between 0 and 1 when edge0 < x < edge1.
#[inline]
pub fn smoothstep(edge0: f64, edge1: f64, a: f64) -> f64 {
    let x = ((a - edge0) / (edge1 - edge0)).clamp(0.0, 1.0);
    x * x * (3.0 - 2.0 * x)
}

/// Sine function where multiples of 90 degrees come out exact.
pub fn sind(x: f64) -> f64 {
    if !x.is_finite() {
        return x.sin();
    }
    if x < 0.0 {
        return -sind(-x);
    }
    let (remainder, quo) = {
        let q = (x / 90.0).floor() as i64;
        let r = x - q as f64 * 90.0;
        (r, q)
    };
    match ((quo % 4) + 4) % 4 {
        0 => radians(remainder).sin(),
        1 => radians(remainder).cos(),
        2 => -radians(remainder).sin(),
        3 => -radians(remainder).cos(),
        _ => 0.0,
    }
}

/// Cosine function where multiples of 90 degrees come out exact.
#[inline]
pub fn cosd(x: f64) -> f64 {
    sind(x + 90.0)
}

// ---------------------------------------------------------------------------
// Polygon types
// ---------------------------------------------------------------------------

/// Single polygon contour, wound CCW. First and last point are implicitly connected.
pub type SimplePolygon = Vec<Vec2>;

/// Set of polygons with holes (arbitrary nesting).
pub type Polygons = Vec<SimplePolygon>;

/// Polygon vertex with index.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PolyVert {
    pub pos: Vec2,
    pub idx: i32,
}

/// Single indexed polygon contour, wound CCW.
pub type SimplePolygonIdx = Vec<PolyVert>;

/// Set of indexed polygons with holes.
pub type PolygonsIdx = Vec<SimplePolygonIdx>;

// ---------------------------------------------------------------------------
// Box (3D axis-aligned bounding box)
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Box {
    pub min: Vec3,
    pub max: Vec3,
}

impl Default for Box {
    fn default() -> Self {
        Box {
            min: Vec3::splat(f64::INFINITY),
            max: Vec3::splat(f64::NEG_INFINITY),
        }
    }
}

impl Box {
    /// Default is an infinite box containing all space.
    pub fn new() -> Self {
        Self::default()
    }

    /// Box containing the two given points.
    pub fn from_points(p1: Vec3, p2: Vec3) -> Self {
        Box {
            min: Vec3::new(p1.x.min(p2.x), p1.y.min(p2.y), p1.z.min(p2.z)),
            max: Vec3::new(p1.x.max(p2.x), p1.y.max(p2.y), p1.z.max(p2.z)),
        }
    }

    /// A box containing a single point.
    pub fn from_point(p: Vec3) -> Self {
        Box { min: p, max: p }
    }

    /// True when the box has no volume (min > max on any axis).
    pub fn is_empty(&self) -> bool {
        self.min.x > self.max.x || self.min.y > self.max.y || self.min.z > self.max.z
    }

    pub fn size(&self) -> Vec3 {
        self.max - self.min
    }

    pub fn center(&self) -> Vec3 {
        (self.max + self.min) * 0.5
    }

    /// Absolute-largest coordinate value.
    pub fn scale(&self) -> f64 {
        let abs_min = Vec3::new(self.min.x.abs(), self.min.y.abs(), self.min.z.abs());
        let abs_max = Vec3::new(self.max.x.abs(), self.max.y.abs(), self.max.z.abs());
        let m = Vec3::new(
            abs_min.x.max(abs_max.x),
            abs_min.y.max(abs_max.y),
            abs_min.z.max(abs_max.z),
        );
        m.x.max(m.y).max(m.z)
    }

    pub fn contains_point(&self, p: Vec3) -> bool {
        p.x >= self.min.x && p.x <= self.max.x
            && p.y >= self.min.y && p.y <= self.max.y
            && p.z >= self.min.z && p.z <= self.max.z
    }

    pub fn contains_box(&self, other: &Box) -> bool {
        other.min.x >= self.min.x && other.max.x <= self.max.x
            && other.min.y >= self.min.y && other.max.y <= self.max.y
            && other.min.z >= self.min.z && other.max.z <= self.max.z
    }

    /// Expand in-place to include the given point.
    pub fn union_point(&mut self, p: Vec3) {
        self.min.x = self.min.x.min(p.x);
        self.min.y = self.min.y.min(p.y);
        self.min.z = self.min.z.min(p.z);
        self.max.x = self.max.x.max(p.x);
        self.max.y = self.max.y.max(p.y);
        self.max.z = self.max.z.max(p.z);
    }

    /// Return the union of this box with another.
    pub fn union_box(&self, other: &Box) -> Box {
        Box {
            min: Vec3::new(
                self.min.x.min(other.min.x),
                self.min.y.min(other.min.y),
                self.min.z.min(other.min.z),
            ),
            max: Vec3::new(
                self.max.x.max(other.max.x),
                self.max.y.max(other.max.y),
                self.max.z.max(other.max.z),
            ),
        }
    }

    /// Transform by axis-aligned affine transform (Mat3x4 * vec4(pt, 1)).
    pub fn transform(&self, t: &Mat3x4) -> Box {
        use crate::linalg::Vec4 as V4;
        let min_t = *t * V4::new(self.min.x, self.min.y, self.min.z, 1.0);
        let max_t = *t * V4::new(self.max.x, self.max.y, self.max.z, 1.0);
        Box {
            min: Vec3::new(min_t.x.min(max_t.x), min_t.y.min(max_t.y), min_t.z.min(max_t.z)),
            max: Vec3::new(min_t.x.max(max_t.x), min_t.y.max(max_t.y), min_t.z.max(max_t.z)),
        }
    }

    pub fn does_overlap_box(&self, other: &Box) -> bool {
        self.min.x <= other.max.x && self.min.y <= other.max.y && self.min.z <= other.max.z
            && self.max.x >= other.min.x && self.max.y >= other.min.y && self.max.z >= other.min.z
    }

    /// Does the given point project within the XY extent (including equality)?
    pub fn does_overlap_point_xy(&self, p: Vec3) -> bool {
        p.x >= self.min.x && p.x <= self.max.x && p.y >= self.min.y && p.y <= self.max.y
    }

    pub fn is_finite(&self) -> bool {
        self.min.x.is_finite() && self.min.y.is_finite() && self.min.z.is_finite()
            && self.max.x.is_finite() && self.max.y.is_finite() && self.max.z.is_finite()
    }
}

impl std::ops::Add<Vec3> for Box {
    type Output = Box;
    fn add(self, shift: Vec3) -> Box {
        Box { min: self.min + shift, max: self.max + shift }
    }
}
impl std::ops::AddAssign<Vec3> for Box {
    fn add_assign(&mut self, shift: Vec3) {
        self.min = self.min + shift;
        self.max = self.max + shift;
    }
}
impl std::ops::Mul<Vec3> for Box {
    type Output = Box;
    fn mul(self, scale: Vec3) -> Box {
        Box { min: self.min * scale, max: self.max * scale }
    }
}
impl std::ops::MulAssign<Vec3> for Box {
    fn mul_assign(&mut self, scale: Vec3) {
        self.min = self.min * scale;
        self.max = self.max * scale;
    }
}

// ---------------------------------------------------------------------------
// Rect (2D axis-aligned bounding box)
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Rect {
    pub min: Vec2,
    pub max: Vec2,
}

impl Default for Rect {
    fn default() -> Self {
        Rect {
            min: Vec2::splat(f64::INFINITY),
            max: Vec2::splat(f64::NEG_INFINITY),
        }
    }
}

impl Rect {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn from_points(a: Vec2, b: Vec2) -> Self {
        Rect {
            min: Vec2::new(a.x.min(b.x), a.y.min(b.y)),
            max: Vec2::new(a.x.max(b.x), a.y.max(b.y)),
        }
    }

    pub fn size(&self) -> Vec2 {
        self.max - self.min
    }

    pub fn area(&self) -> f64 {
        let sz = self.size();
        sz.x * sz.y
    }

    pub fn scale(&self) -> f64 {
        let abs_min = Vec2::new(self.min.x.abs(), self.min.y.abs());
        let abs_max = Vec2::new(self.max.x.abs(), self.max.y.abs());
        let m = Vec2::new(abs_min.x.max(abs_max.x), abs_min.y.max(abs_max.y));
        m.x.max(m.y)
    }

    pub fn center(&self) -> Vec2 {
        (self.max + self.min) * 0.5
    }

    pub fn contains_point(&self, p: Vec2) -> bool {
        p.x >= self.min.x && p.x <= self.max.x && p.y >= self.min.y && p.y <= self.max.y
    }

    pub fn contains_rect(&self, other: &Rect) -> bool {
        other.min.x >= self.min.x && other.max.x <= self.max.x
            && other.min.y >= self.min.y && other.max.y <= self.max.y
    }

    pub fn does_overlap(&self, other: &Rect) -> bool {
        self.min.x <= other.max.x && self.min.y <= other.max.y
            && self.max.x >= other.min.x && self.max.y >= other.min.y
    }

    pub fn is_empty(&self) -> bool {
        self.max.y <= self.min.y || self.max.x <= self.min.x
    }

    pub fn is_finite(&self) -> bool {
        self.min.x.is_finite() && self.min.y.is_finite()
            && self.max.x.is_finite() && self.max.y.is_finite()
    }

    pub fn union_point(&mut self, p: Vec2) {
        self.min.x = self.min.x.min(p.x);
        self.min.y = self.min.y.min(p.y);
        self.max.x = self.max.x.max(p.x);
        self.max.y = self.max.y.max(p.y);
    }

    pub fn union_rect(&self, other: &Rect) -> Rect {
        Rect {
            min: Vec2::new(self.min.x.min(other.min.x), self.min.y.min(other.min.y)),
            max: Vec2::new(self.max.x.max(other.max.x), self.max.y.max(other.max.y)),
        }
    }
}

impl std::ops::Add<Vec2> for Rect {
    type Output = Rect;
    fn add(self, shift: Vec2) -> Rect {
        Rect { min: self.min + shift, max: self.max + shift }
    }
}
impl std::ops::AddAssign<Vec2> for Rect {
    fn add_assign(&mut self, shift: Vec2) {
        self.min = self.min + shift;
        self.max = self.max + shift;
    }
}
impl std::ops::Mul<Vec2> for Rect {
    type Output = Rect;
    fn mul(self, scale: Vec2) -> Rect {
        Rect { min: self.min * scale, max: self.max * scale }
    }
}
impl std::ops::MulAssign<Vec2> for Rect {
    fn mul_assign(&mut self, scale: Vec2) {
        self.min = self.min * scale;
        self.max = self.max * scale;
    }
}

// ---------------------------------------------------------------------------
// OpType
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum OpType {
    Add,
    Subtract,
    Intersect,
}

// ---------------------------------------------------------------------------
// Error
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Error {
    NoError,
    NonFiniteVertex,
    NotManifold,
    VertexOutOfBounds,
    PropertiesWrongLength,
    MissingPositionProperties,
    MergeVectorsDifferentLengths,
    MergeIndexOutOfBounds,
    TransformWrongLength,
    RunIndexWrongLength,
    FaceIdWrongLength,
    InvalidConstruction,
    ResultTooLarge,
}

impl Error {
    pub fn to_str(self) -> &'static str {
        match self {
            Error::NoError => "No Error",
            Error::NonFiniteVertex => "Non-Finite Vertex",
            Error::NotManifold => "Not Manifold",
            Error::VertexOutOfBounds => "Vertex Out of Bounds",
            Error::PropertiesWrongLength => "Properties Wrong Length",
            Error::MissingPositionProperties => "Missing Position Properties",
            Error::MergeVectorsDifferentLengths => "Merge Vectors Different Lengths",
            Error::MergeIndexOutOfBounds => "Merge Index Out of Bounds",
            Error::TransformWrongLength => "Transform Wrong Length",
            Error::RunIndexWrongLength => "Run Index Wrong Length",
            Error::FaceIdWrongLength => "Face ID Wrong Length",
            Error::InvalidConstruction => "Invalid Construction",
            Error::ResultTooLarge => "Result Too Large",
        }
    }
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.to_str())
    }
}

// ---------------------------------------------------------------------------
// Quality (static global for circle quantization)
// ---------------------------------------------------------------------------

use std::sync::OnceLock;
use std::sync::Mutex;

struct QualityState {
    min_circular_angle: f64,
    min_circular_edge_length: f64,
    circular_segments: i32,
}

static QUALITY_STATE: OnceLock<Mutex<QualityState>> = OnceLock::new();

fn quality_state() -> &'static Mutex<QualityState> {
    QUALITY_STATE.get_or_init(|| {
        Mutex::new(QualityState {
            min_circular_angle: DEFAULT_ANGLE,
            min_circular_edge_length: DEFAULT_LENGTH,
            circular_segments: DEFAULT_SEGMENTS,
        })
    })
}

pub struct Quality;

impl Quality {
    pub fn set_min_circular_angle(angle: f64) {
        quality_state().lock().unwrap().min_circular_angle = angle;
    }

    pub fn set_min_circular_edge_length(length: f64) {
        quality_state().lock().unwrap().min_circular_edge_length = length;
    }

    pub fn set_circular_segments(n: i32) {
        quality_state().lock().unwrap().circular_segments = n;
    }

    pub fn get_circular_segments(radius: f64) -> i32 {
        let q = quality_state().lock().unwrap();
        if q.circular_segments > 0 {
            return q.circular_segments;
        }
        // Match C++ exactly: int truncation (not ceil), fmin (not fmax), round down to multiple of 4
        let n_seg_a = (360.0 / q.min_circular_angle) as i32;
        let n_seg_l = (2.0 * radius.abs() * K_PI / q.min_circular_edge_length) as i32;
        let mut n_seg = n_seg_a.min(n_seg_l) + 3;
        n_seg -= n_seg % 4;
        n_seg.max(4)
    }

    pub fn reset_to_defaults() {
        let mut q = quality_state().lock().unwrap();
        q.min_circular_angle = DEFAULT_ANGLE;
        q.min_circular_edge_length = DEFAULT_LENGTH;
        q.circular_segments = DEFAULT_SEGMENTS;
    }
}

// ---------------------------------------------------------------------------
// ExecutionParams
// ---------------------------------------------------------------------------

#[derive(Clone, Debug)]
pub struct ExecutionParams {
    pub intermediate_checks: bool,
    pub self_intersection_checks: bool,
    pub process_overlaps: bool,
    pub suppress_errors: bool,
    pub cleanup_triangles: bool,
    pub verbose: i32,
}

impl Default for ExecutionParams {
    fn default() -> Self {
        ExecutionParams {
            intermediate_checks: false,
            self_intersection_checks: false,
            process_overlaps: true,
            suppress_errors: false,
            cleanup_triangles: true,
            verbose: 0,
        }
    }
}

// ---------------------------------------------------------------------------
// Smoothness
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Smoothness {
    /// The halfedge index = 3 * tri + i
    pub halfedge: usize,
    /// 0 = sharp, 1 = smooth
    pub smoothness: f64,
}

// ---------------------------------------------------------------------------
// MeshGLP / MeshGL / MeshGL64
// ---------------------------------------------------------------------------

/// GL-style mesh representation. Generic over precision (f32/f64) and index type (u32/u64).
#[derive(Clone, Debug, Default)]
pub struct MeshGLP<P: Copy + Default, I: Copy + Default = u32> {
    /// Number of properties per vertex, always >= 3.
    pub num_prop: I,
    /// Flat interleaved vertex properties: [x, y, z, ...] × num_verts.
    pub vert_properties: Vec<P>,
    /// Triangle vertex indices, 3 per triangle (CCW from outside).
    pub tri_verts: Vec<I>,
    /// Optional: merge-from vertex indices.
    pub merge_from_vert: Vec<I>,
    /// Optional: merge-to vertex indices.
    pub merge_to_vert: Vec<I>,
    /// Optional: run start indices into triVerts.
    pub run_index: Vec<I>,
    /// Optional: original mesh ID per run.
    pub run_original_id: Vec<u32>,
    /// Optional: 3×4 column-major transform per run (12 elements each).
    pub run_transform: Vec<P>,
    /// Optional: source face ID per triangle.
    pub face_id: Vec<I>,
    /// Optional: halfedge tangent vectors (4 per halfedge).
    pub halfedge_tangent: Vec<P>,
    /// Tolerance for mesh simplification.
    pub tolerance: P,
}

impl MeshGLP<f32, u32> {
    pub fn num_vert(&self) -> usize {
        if self.num_prop == 0 { 0 } else { self.vert_properties.len() / self.num_prop as usize }
    }

    pub fn num_tri(&self) -> usize {
        self.tri_verts.len() / 3
    }

    pub fn get_vert_pos(&self, v: usize) -> [f32; 3] {
        let offset = v * self.num_prop as usize;
        [self.vert_properties[offset], self.vert_properties[offset + 1], self.vert_properties[offset + 2]]
    }

    pub fn get_tri_verts(&self, t: usize) -> [u32; 3] {
        let offset = 3 * t;
        [self.tri_verts[offset], self.tri_verts[offset + 1], self.tri_verts[offset + 2]]
    }

    pub fn get_tangent(&self, h: usize) -> [f32; 4] {
        let offset = 4 * h;
        [
            self.halfedge_tangent[offset],
            self.halfedge_tangent[offset + 1],
            self.halfedge_tangent[offset + 2],
            self.halfedge_tangent[offset + 3],
        ]
    }
}

impl MeshGLP<f64, u64> {
    pub fn num_vert(&self) -> usize {
        if self.num_prop == 0 { 0 } else { self.vert_properties.len() / self.num_prop as usize }
    }

    pub fn num_tri(&self) -> usize {
        self.tri_verts.len() / 3
    }

    pub fn get_vert_pos(&self, v: usize) -> [f64; 3] {
        let offset = v * self.num_prop as usize;
        [self.vert_properties[offset], self.vert_properties[offset + 1], self.vert_properties[offset + 2]]
    }

    pub fn get_tri_verts(&self, t: usize) -> [u64; 3] {
        let offset = 3 * t;
        [self.tri_verts[offset], self.tri_verts[offset + 1], self.tri_verts[offset + 2]]
    }

    pub fn get_tangent(&self, h: usize) -> [f64; 4] {
        let offset = 4 * h;
        [
            self.halfedge_tangent[offset],
            self.halfedge_tangent[offset + 1],
            self.halfedge_tangent[offset + 2],
            self.halfedge_tangent[offset + 3],
        ]
    }
}

/// Single-precision mesh (standard for graphics).
pub type MeshGL = MeshGLP<f32, u32>;

/// Double-precision, 64-bit index mesh (for huge meshes).
pub type MeshGL64 = MeshGLP<f64, u64>;

// ---------------------------------------------------------------------------
// Halfedge (from src/shared.h)
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Halfedge {
    pub start_vert: i32,
    pub end_vert: i32,
    pub paired_halfedge: i32,
    pub prop_vert: i32,
}

impl Halfedge {
    pub fn is_forward(&self) -> bool {
        self.start_vert < self.end_vert
    }
}

impl PartialOrd for Halfedge {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Halfedge {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        if self.start_vert == other.start_vert {
            self.end_vert.cmp(&other.end_vert)
        } else {
            self.start_vert.cmp(&other.start_vert)
        }
    }
}

// ---------------------------------------------------------------------------
// Barycentric (from src/shared.h)
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Barycentric {
    pub tri: i32,
    pub uvw: Vec4,
}

// ---------------------------------------------------------------------------
// TriRef (from src/shared.h)
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct TriRef {
    /// Unique ID of the mesh instance of this triangle.
    pub mesh_id: i32,
    /// OriginalID of the mesh this triangle came from.
    pub original_id: i32,
    /// Source face ID.
    pub face_id: i32,
    /// Triangles with same coplanar_id are coplanar.
    pub coplanar_id: i32,
}

impl TriRef {
    pub fn same_face(&self, other: &TriRef) -> bool {
        self.mesh_id == other.mesh_id
            && self.coplanar_id == other.coplanar_id
            && self.face_id == other.face_id
    }
}

// ---------------------------------------------------------------------------
// TmpEdge (from src/shared.h)
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct TmpEdge {
    pub first: i32,
    pub second: i32,
    pub halfedge_idx: i32,
}

impl TmpEdge {
    pub fn new(start: i32, end: i32, idx: i32) -> Self {
        TmpEdge {
            first: start.min(end),
            second: start.max(end),
            halfedge_idx: idx,
        }
    }
}

impl PartialOrd for TmpEdge {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for TmpEdge {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        if self.first == other.first {
            self.second.cmp(&other.second)
        } else {
            self.first.cmp(&other.first)
        }
    }
}

// ---------------------------------------------------------------------------
// MeshRelationD (from src/impl.h)
// ---------------------------------------------------------------------------

/// Transform relation between meshes.
#[derive(Clone, Debug)]
pub struct Relation {
    pub original_id: i32,
    pub transform: Mat3x4,
    pub back_side: bool,
}

impl Default for Relation {
    fn default() -> Self {
        Relation {
            original_id: -1,
            transform: Mat3x4::identity(),
            back_side: false,
        }
    }
}

/// Mesh relation table stored on ManifoldImpl.
#[derive(Clone, Debug, Default)]
pub struct MeshRelationD {
    /// originalID of this Manifold if it is an original; -1 otherwise.
    pub original_id: i32,
    pub mesh_id_transform: HashMap<i32, Relation>,
    pub tri_ref: Vec<TriRef>,
}

impl MeshRelationD {
    pub fn new() -> Self {
        MeshRelationD {
            original_id: -1,
            mesh_id_transform: HashMap::new(),
            tri_ref: Vec::new(),
        }
    }
}

// ---------------------------------------------------------------------------
// Inline utility from shared.h
// ---------------------------------------------------------------------------

/// Return next halfedge index within the same triangle (wraps 0→1→2→0).
#[inline]
pub fn next_halfedge(current: i32) -> i32 {
    let n = current + 1;
    if n % 3 == 0 { n - 3 } else { n }
}

/// Return next index within 0..3 (wraps 0→1→2→0).
#[inline]
pub fn next3(i: i32) -> i32 {
    let n = i + 1;
    if n == 3 { 0 } else { n }
}

// ---------------------------------------------------------------------------
#[cfg(test)]
#[path = "types_tests.rs"]
mod tests;
