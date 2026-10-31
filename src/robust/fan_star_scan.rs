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
//    first, against a tree over each later star's triangles, the later stars
//    being those a tree over the stars' boxes finds the row's box reaching.
//    One tree over every owned triangle would list a fan triangle's whole
//    own star (every box in a fan overlaps at the centre), the n²/2 pairs
//    this walk exists to skip, so each star keeps its own tree;
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

/// Some members of a box list (triangles, or stars), findable by box.
struct BoxSet {
    /// Leaf index -> member index.
    members: Vec<usize>,
    collider: Collider,
}

impl BoxSet {
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
        BoxSet {
            members: order,
            collider,
        }
    }

    /// Every member whose box overlaps `query`, ascending. A tree of one leaf
    /// has no internal node to report it, so small sets are walked directly.
    fn overlapping(&self, query: &Box, boxes: &[Box], out: &mut Vec<usize>) {
        out.clear();
        if self.members.len() < 2 {
            out.extend(
                self.members
                    .iter()
                    .copied()
                    .filter(|&j| query.does_overlap_box(&boxes[j])),
            );
        } else {
            self.collider
                .collisions_one(query, 0, |_, leaf| out.push(self.members[leaf]));
        }
        out.sort_unstable();
    }
}

/// The broad phase of the flat-star walk: the trees a row is queried against.
pub(super) struct FlatStarWalk<'a> {
    boxes: &'a [Box],
    live: &'a [bool],
    stars: &'a FlatStars,
    /// The live triangles no star owns.
    unowned: BoxSet,
    /// Per star, its owned triangles.
    owned: Vec<BoxSet>,
    /// Per star, the union of its owned triangles' boxes.
    star_boxes: Vec<Box>,
    /// A tree over `star_boxes`, so a row finds the stars its box reaches
    /// without testing every star: with S stars, a loop over them all is S
    /// tests per owned triangle, quadratic in S (a few thousand small fans,
    /// a part with many round holes, would pay billions).
    star_tree: BoxSet,
}

impl<'a> FlatStarWalk<'a> {
    pub(super) fn new(boxes: &'a [Box], live: &'a [bool], stars: &'a FlatStars) -> Self {
        let unowned = BoxSet::new(
            (0..boxes.len())
                .filter(|&i| live[i] && stars.owner[i].is_none())
                .collect(),
            boxes,
        );
        let owned: Vec<BoxSet> = stars
            .owned
            .iter()
            .map(|members| BoxSet::new(members.clone(), boxes))
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
        let star_tree = BoxSet::new((0..star_boxes.len()).collect(), &star_boxes);
        FlatStarWalk {
            boxes,
            live,
            stars,
            unowned,
            owned,
            star_boxes,
            star_tree,
        }
    }

    /// Calls `test(j)` on every partner `j` row `i` is responsible for (the
    /// four-way split in the header), in a fixed order - the unowned
    /// candidates ascending, then each later star's ascending, stars in
    /// ascending order - and answers whether some call answered true,
    /// stopping at the first. `cands` and `hit_stars` are scratch.
    pub(super) fn row_any(
        &self,
        i: usize,
        cands: &mut Vec<usize>,
        hit_stars: &mut Vec<usize>,
        mut test: impl FnMut(usize) -> bool,
    ) -> bool {
        if !self.live[i] {
            return false;
        }
        let query = &self.boxes[i];
        self.unowned.overlapping(query, self.boxes, cands);
        let Some(star) = self.stars.owner[i] else {
            return cands.iter().any(|&j| j > i && test(j));
        };
        if cands.iter().any(|&j| test(j)) {
            return true;
        }
        self.star_tree
            .overlapping(query, &self.star_boxes, hit_stars);
        for &later in hit_stars.iter().filter(|&&s| s > star) {
            self.owned[later].overlapping(query, self.boxes, cands);
            if cands.iter().any(|&j| test(j)) {
                return true;
            }
        }
        false
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
    let walk = FlatStarWalk::new(boxes, live, stars);
    crate::par::maybe_par_any_ct(
        tris.len(),
        SELF_INTERSECT_PAR_THRESHOLD,
        token,
        || {
            (
                Vec::<usize>::new(),
                Vec::<usize>::new(),
                SelfCutStats::default(),
            )
        },
        |i, (cands, hit_stars, stats)| {
            walk.row_any(i, cands, hit_stars, |j| {
                live[j]
                    && boxes[i].does_overlap_box(&boxes[j])
                    && genuine_contact(tris[i.min(j)], tris[i.max(j)], stats)
            })
        },
    )
}

/// How many pairs the flat-star walk hands the narrow phase's filter (every
/// call [`scan_with_flat_stars`]'s rows make, with no early exit): the
/// walk's whole cost in pairs, a count the tests pin where a clock would
/// flake on a loaded machine.
#[cfg(test)]
pub(super) fn listed_pairs(boxes: &[Box], live: &[bool], stars: &FlatStars) -> usize {
    let walk = FlatStarWalk::new(boxes, live, stars);
    let (mut cands, mut hit_stars, mut count) = (Vec::new(), Vec::new(), 0usize);
    for i in 0..boxes.len() {
        walk.row_any(i, &mut cands, &mut hit_stars, |_| {
            count += 1;
            false
        });
    }
    count
}

/// How many stars the rows' boxes reach through the star tree, summed over
/// every row: the walk's cross-star cost, which a loop over every star
/// would make rows times stars.
#[cfg(test)]
pub(super) fn star_tree_hits(boxes: &[Box], live: &[bool], stars: &FlatStars) -> usize {
    let walk = FlatStarWalk::new(boxes, live, stars);
    let (mut cands, mut hit_stars, mut count) = (Vec::new(), Vec::new(), 0usize);
    for i in 0..boxes.len() {
        hit_stars.clear();
        walk.row_any(i, &mut cands, &mut hit_stars, |_| false);
        count += hit_stars.len();
    }
    count
}
