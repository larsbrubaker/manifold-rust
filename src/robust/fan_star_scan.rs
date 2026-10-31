// robust/fan_star_scan.rs — The self-intersection scan's candidate walk for a
// mesh with flat, disjoint stars (robust/fan_star.rs), called from
// robust/soup.rs `compute_self_intersections` in place of its one-tree walk.
//
// The pairs a one-tree walk would test split four ways, and each is tested
// exactly once, as `genuine_contact(lower index, higher index)` like the
// one-tree walk tests it:
//  - neither triangle owned by a star: the row of the lower one, against a
//    tree over the unowned triangles;
//  - one owned, one not: the owned triangle's row, against that same tree
//    (unowned rows never see owned triangles, so nothing is tested twice);
//  - both owned by different stars: the row of the triangle whose star comes
//    first, against a tree over the later star's triangles;
//  - both owned by one star: proven ordinary adjacency, never listed.
// So the verdict is the one-tree walk's, without the star's n²/2 box pairs
// ever being listed. Rows run through `par::maybe_par_any_ct`, so the serial
// and parallel builds agree and the cancellation contract is the scan's.

use crate::collider::Collider;
use crate::linalg::Vec3;
use crate::types::Box;

use super::fan_star::FlatStars;
use super::intersection_graph::SelfCutStats;
use super::soup::{genuine_contact, SELF_INTERSECT_PAR_THRESHOLD};

/// Some triangles of the mesh, findable by box.
struct TriangleSet {
    /// Leaf index -> triangle index.
    tris: Vec<usize>,
    collider: Collider,
}

impl TriangleSet {
    /// A morton-ordered tree over `members`, as the scan builds for soups.
    fn new(members: Vec<usize>, boxes: &[Box]) -> Self {
        let scene = members
            .iter()
            .fold(Box::new(), |acc, &i| acc.union_box(&boxes[i]));
        let mut order = members;
        order.sort_by_key(|&i| (crate::sort::morton_code(boxes[i].center(), &scene), i));
        let collider = Collider::new(
            order.iter().map(|&i| boxes[i]).collect(),
            order
                .iter()
                .map(|&i| crate::sort::morton_code(boxes[i].center(), &scene))
                .collect(),
        );
        TriangleSet {
            tris: order,
            collider,
        }
    }

    /// Every member whose box overlaps `query`, ascending. A tree of one leaf
    /// has no internal node to report it, so small sets are walked directly.
    fn overlapping(&self, query: &Box, boxes: &[Box], out: &mut Vec<usize>) {
        out.clear();
        if self.tris.len() < 2 {
            out.extend(
                self.tris
                    .iter()
                    .copied()
                    .filter(|&j| query.does_overlap_box(&boxes[j])),
            );
        } else {
            self.collider
                .collisions_one(query, 0, |_, leaf| out.push(self.tris[leaf]));
        }
        out.sort_unstable();
    }
}

/// The scan's verdict for a mesh with at least one flat, disjoint star:
/// `Some(true)` on the first genuine contact, `None` when cancelled first.
pub(super) fn scan_with_flat_stars(
    tris: &[[Vec3; 3]],
    boxes: &[Box],
    live: &[bool],
    stars: &FlatStars,
    token: Option<&crate::cancel::CancelToken>,
) -> Option<bool> {
    let unowned = TriangleSet::new(
        (0..tris.len())
            .filter(|&i| live[i] && stars.owner[i].is_none())
            .collect(),
        boxes,
    );
    let owned: Vec<TriangleSet> = stars
        .owned
        .iter()
        .map(|members| TriangleSet::new(members.clone(), boxes))
        .collect();
    let star_boxes: Vec<Box> = stars
        .owned
        .iter()
        .map(|members| {
            members
                .iter()
                .fold(Box::new(), |acc, &i| acc.union_box(&boxes[i]))
        })
        .collect();

    let contact = |i: usize, j: usize, stats: &mut SelfCutStats| {
        live[j]
            && boxes[i].does_overlap_box(&boxes[j])
            && genuine_contact(tris[i.min(j)], tris[i.max(j)], stats)
    };

    crate::par::maybe_par_any_ct(
        tris.len(),
        SELF_INTERSECT_PAR_THRESHOLD,
        token,
        || (Vec::<usize>::new(), SelfCutStats::default()),
        |i, (cands, stats)| {
            if !live[i] {
                return false;
            }
            unowned.overlapping(&boxes[i], boxes, cands);
            match stars.owner[i] {
                None => {
                    if cands.iter().any(|&j| j > i && contact(i, j, stats)) {
                        return true;
                    }
                }
                Some(star) => {
                    if cands.iter().any(|&j| contact(i, j, stats)) {
                        return true;
                    }
                    for later in star + 1..owned.len() {
                        if !star_boxes[later].does_overlap_box(&boxes[i]) {
                            continue;
                        }
                        owned[later].overlapping(&boxes[i], boxes, cands);
                        if cands.iter().any(|&j| contact(i, j, stats)) {
                            return true;
                        }
                    }
                }
            }
            false
        },
    )
}
