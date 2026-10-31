// robust/fan_star.rs — Settles, in one sort, every triangle pair of a flat
// fan that the self-intersection scan (robust/soup.rs
// `compute_self_intersections`) would otherwise test pair by pair.
//
// A fan of n triangles round one vertex (a cylinder's cap, a revolved
// part's end, a CAD export's flat face) puts n triangles' boxes on top of
// each other at the centre, so a box broad phase lists all n²/2 pairs and
// the narrow phase (`graph_self_cut::real_self_contact`) settles each one
// with several exact predicates: 450 million pairs and tens of seconds for a
// 30,000-sided cap. Here a vertex's whole star is settled at once instead:
// when every triangle round it lies exactly in one plane, the triangles are
// sectors (wedges) of that plane at the vertex, and sorting the wedges by
// angle proves in O(n log n) whether any two of them overlap. A star proven
// to be non-overlapping is "flat and disjoint": any two of its triangles
// meet either in the centre alone, or along the one edge they share where
// two neighbouring wedges touch. The narrow phase calls both of those
// ordinary adjacency (an exact coplanar clip down to one point is
// `TriTriIsect::Point`; a coplanar shared edge whose opposite corners lie
// strictly on either side is its edge-neighbour shortcut), so the scan may
// skip those pairs without moving any verdict.
//
// Exact and conservative: coplanarity is the exact `orient3d`, the angular
// order and every overlap test the exact `orient2d`, and a star is settled
// only when everything is proven - any zero sign the proof cannot use (two
// wedges starting on one ray, wedges touching along a ray without sharing
// the far corner, a projection that flattens the plane) leaves the whole
// star to the pair-by-pair scan. Pairs between a settled star's triangle
// and any triangle outside the star are never settled here.

use std::cmp::Ordering;

use crate::linalg::{Vec2, Vec3};

use super::exact::filtered::orient2d;
use super::exact::Sign;
use super::intersection_graph::orient3d_plane;

/// Stars smaller than this are left to the pair-by-pair scan, whose cost on
/// them is small; the threshold only keeps ordinary meshes (vertex degree
/// four to eight) on exactly the path they always took.
pub(super) const FLAT_STAR_MIN_TRIANGLES: usize = 64;

/// The flat, disjoint stars of a mesh: which star (if any) owns each
/// triangle. A triangle in more than one settled star is owned by the first
/// in the stars' (deterministic) key order; every pair of triangles owned by
/// one star is proven benign, which is all the scan reads.
pub(super) struct FlatStars {
    /// Per triangle, the owning star's index, or `None`.
    pub(super) owner: Vec<Option<usize>>,
    /// Per star, the triangles it owns, ascending.
    pub(super) owned: Vec<Vec<usize>>,
    /// Whether the corners were sorted and grouped, rather than skipped by
    /// [`may_hold_a_star`]: how the tests see the skip without a clock.
    #[cfg(test)]
    pub(super) grouped: bool,
}

/// Weld key: exact position identity, with -0.0 read as 0.0 (the same key
/// robust/soup.rs welds by, and the same identity `real_self_contact`'s f64
/// `==` gives).
fn pos_key(v: Vec3) -> (u64, u64, u64) {
    let norm = |x: f64| if x == 0.0 { 0.0f64 } else { x }.to_bits();
    (norm(v.x), norm(v.y), norm(v.z))
}

/// Whether some position could be the corner of at least
/// `FLAT_STAR_MIN_TRIANGLES` live triangles, counted in hashed buckets: a
/// few milliseconds where the sort and grouping below cost about a sixth of
/// the whole scan on an ordinary mesh (48 of 330 ms for a 200k-triangle
/// sphere), which has no star to find.
///
/// Exact in the direction that matters: equal positions have equal
/// `pos_key`s and so share a bucket, and other positions hashing into it
/// only add to its count, so a bucket's count is at least the size of every
/// star in it. When no bucket reaches the threshold, no star does, and
/// skipping the grouping gives exactly the stars (none) it would have found.
/// A collision can only send a mesh to the grouping needlessly, never skip
/// a star; the scan's verdict is the same either way.
pub(super) fn may_hold_a_star(tris: &[[Vec3; 3]], live: &[bool]) -> bool {
    let corners = 3 * live.iter().filter(|&&l| l).count();
    if corners < FLAT_STAR_MIN_TRIANGLES {
        return false;
    }
    // At least one bucket per corner, so collisions stay rare and an
    // ordinary mesh passes; a power of two, so the top hash bits index it.
    // At least 64 corners reach here, so `bits` is at least 6.
    let bits = corners.next_power_of_two().trailing_zeros();
    let mut counts = vec![0u32; 1usize << bits];
    for (t, tri) in tris.iter().enumerate() {
        if !live[t] {
            continue;
        }
        for &v in tri {
            let (x, y, z) = pos_key(v);
            let mut h = 0u64;
            for part in [x, y, z] {
                h = (h.rotate_left(26) ^ part).wrapping_mul(0x9E37_79B9_7F4A_7C15);
            }
            let count = &mut counts[(h >> (64 - bits)) as usize];
            *count += 1;
            // Every triangle with a corner at a position adds at least one to
            // that position's bucket, so reaching the threshold is the only
            // way a star of that size can exist.
            if *count as usize >= FLAT_STAR_MIN_TRIANGLES {
                return true;
            }
        }
    }
    false
}

/// Finds and proves the flat, disjoint stars among the live triangles.
pub(super) fn flat_stars(tris: &[[Vec3; 3]], live: &[bool]) -> FlatStars {
    if !may_hold_a_star(tris, live) {
        return FlatStars {
            owner: vec![None; tris.len()],
            owned: Vec::new(),
            #[cfg(test)]
            grouped: false,
        };
    }
    let mut corners: Vec<((u64, u64, u64), usize)> = Vec::with_capacity(3 * tris.len());
    for (t, tri) in tris.iter().enumerate() {
        if live[t] {
            corners.extend(tri.iter().map(|&v| (pos_key(v), t)));
        }
    }
    corners.sort_unstable();

    let mut owner: Vec<Option<usize>> = vec![None; tris.len()];
    let mut owned: Vec<Vec<usize>> = Vec::new();
    let mut start = 0;
    while start < corners.len() {
        let key = corners[start].0;
        let mut end = start + 1;
        while end < corners.len() && corners[end].0 == key {
            end += 1;
        }
        let star: Vec<usize> = corners[start..end].iter().map(|&(_, t)| t).collect();
        if star.len() >= FLAT_STAR_MIN_TRIANGLES && is_flat_and_disjoint(tris, &star, key) {
            let index = owned.len();
            let mine: Vec<usize> = star
                .iter()
                .copied()
                .filter(|&t| owner[t].is_none())
                .collect();
            for &t in &mine {
                owner[t] = Some(index);
            }
            owned.push(mine);
        }
        start = end;
    }
    FlatStars {
        owner,
        owned,
        #[cfg(test)]
        grouped: true,
    }
}

/// One triangle of a star as a wedge at the centre: its two other corners,
/// counter-clockwise (in the projection) from `start` to `end`.
struct Wedge {
    start: Vec2,
    end: Vec2,
}

/// Whether every triangle of `star` (all live, all with a corner at
/// `centre_key`) lies in one plane and no two of their wedges overlap.
fn is_flat_and_disjoint(tris: &[[Vec3; 3]], star: &[usize], centre_key: (u64, u64, u64)) -> bool {
    let reference = tris[star[0]];
    // A star in a plane of constant x, y or z (every cap a CAD tool or a
    // revolve makes square to an axis) is coplanar by plain f64 equality;
    // exactly coplanar points are the case `orient3d`'s float filter can
    // never decide, so the general test pays the exact fallback per corner.
    let on_axis_plane = (0..3).any(|axis| {
        star.iter()
            .all(|&t| tris[t].iter().all(|v| v[axis] == reference[0][axis]))
    });
    if !on_axis_plane
        && star.iter().any(|&t| {
            tris[t]
                .iter()
                .any(|&v| orient3d_plane(&reference, v) != Sign::Zero)
        })
    {
        return false;
    }
    // The drop axis is only a projection choice: if it flattens the plane,
    // every wedge's orientation below comes back Zero and the star is left
    // alone.
    let axis = dominant_axis(reference);
    let project = |v: Vec3| match axis {
        0 => Vec2::new(v.y, v.z),
        1 => Vec2::new(v.z, v.x),
        _ => Vec2::new(v.x, v.y),
    };
    let Some(&centre3) = reference.iter().find(|&&v| pos_key(v) == centre_key) else {
        return false;
    };
    let centre = project(centre3);

    let mut wedges: Vec<Wedge> = Vec::with_capacity(star.len());
    for &t in star {
        let mut others = tris[t].iter().filter(|&&v| pos_key(v) != centre_key);
        let (Some(&a), Some(&b), None) = (others.next(), others.next(), others.next()) else {
            return false;
        };
        let (a, b) = (project(a), project(b));
        match orient2d(centre, a, b) {
            Sign::Pos => wedges.push(Wedge { start: a, end: b }),
            Sign::Neg => wedges.push(Wedge { start: b, end: a }),
            Sign::Zero => return false,
        }
    }
    wedges.sort_by(|x, y| angular_order(centre, x.start, y.start));

    (0..wedges.len()).all(|k| {
        let this = &wedges[k];
        let next = &wedges[(k + 1) % wedges.len()];
        ends_before(centre, this, next.start)
    })
}

/// Whether wedge `w` ends no later than direction `next` (from the centre
/// `c`), counter-clockwise from `w`'s start, where `next` is the start of
/// the wedge after it in angular order. Touching is allowed only where the
/// two wedges share the far corner itself - a shared edge.
fn ends_before(c: Vec2, w: &Wedge, next: Vec2) -> bool {
    if angular_order(c, w.start, next) == Ordering::Equal {
        // Two wedges starting on one ray overlap.
        return false;
    }
    if orient2d(c, w.start, next) != Sign::Pos {
        // `next` is half a turn or more round from the start, and the wedge
        // (a triangle's corner) spans strictly less than half a turn.
        return true;
    }
    match orient2d(c, w.end, next) {
        Sign::Pos => true,
        Sign::Neg => false,
        // Both lie strictly within half a turn of the start, so collinear
        // means the same ray: fine only when it is the same corner.
        Sign::Zero => w.end == next,
    }
}

/// Counter-clockwise order of directions from `c`, starting along +x:
/// exact, from f64 comparisons and `orient2d` on the original points.
fn angular_order(c: Vec2, x: Vec2, y: Vec2) -> Ordering {
    let upper = |p: Vec2| p.y > c.y || (p.y == c.y && p.x > c.x);
    match (upper(x), upper(y)) {
        (true, false) => Ordering::Less,
        (false, true) => Ordering::Greater,
        _ => match orient2d(c, x, y) {
            Sign::Pos => Ordering::Less,
            Sign::Neg => Ordering::Greater,
            Sign::Zero => Ordering::Equal,
        },
    }
}

/// Axis of the largest |component| of the (f64) triangle normal.
fn dominant_axis(t: [Vec3; 3]) -> usize {
    let n = crate::linalg::cross(t[1] - t[0], t[2] - t[0]);
    let (ax, ay, az) = (n.x.abs(), n.y.abs(), n.z.abs());
    if az >= ax && az >= ay {
        2
    } else if ay >= ax {
        1
    } else {
        0
    }
}
