// boolean_stage_progress.rs — NOT A C++ PORT. The exact boolean's optional
// stage sink, mirrored 1:1 from manifold-sharp's BooleanStageProgress.cs (its
// RUST_DIVERGENCES.md entry 6; here docs/CPP_DIVERGENCES.md entry 13): a
// `StageSink` that `boolean3::boolean_with_token_and_stage`,
// `Boolean3::new_with_token_and_stage` and
// `boolean_result::boolean_result_with_token_and_stage` take, and invoke with
// the cumulative marks below at the cancel gates that close each heavy stage.
// Its one caller is convex_dilation.rs's union tree.
//
// A side channel only: the sink is invoked between stages, reads nothing and is
// handed nothing but a constant, so every computed value is the one the
// variants without it produce. `None` is the pre-existing code.
//
// The marks are one boolean's own completed fraction, from stage times
// measured on two warm unions (the edge-face intersection passes 32-54%,
// simplify_topology 22-33%, triangulation 7-16%, the closing sort 2-21%,
// winding and edge assembly the rest). They are estimates, not a contract; a
// mark only has to be monotone and at most 1. The last is reported after
// simplify_topology, so a caller never hears 1.0 before the boolean returns.

/// A boolean's stage sink: hears this boolean's own completed fraction.
pub type StageSink<'a> = Option<&'a (dyn Fn(f64) + Sync)>;

/// After `intersect12` P→Q.
pub const AFTER_INTERSECT_PQ: f64 = 0.22;
/// After `intersect12` Q→P.
pub const AFTER_INTERSECT_QP: f64 = 0.44;
/// After `winding03` on P.
pub const AFTER_WINDING_P: f64 = 0.46;
/// After `winding03` on Q.
pub const AFTER_WINDING_Q: f64 = 0.48;
/// After the edge assembly (`append_whole_edges`).
pub const AFTER_ASSEMBLY: f64 = 0.53;
/// After `face2tri` and `reorder_halfedges`.
pub const AFTER_TRIANGULATION: f64 = 0.63;
/// After `create_properties` and `update_reference`.
pub const AFTER_REFERENCE: f64 = 0.64;
/// After `simplify_topology`; only the closing sort remains.
pub const AFTER_SIMPLIFY: f64 = 0.90;

/// Invoke `stage` with `mark` when a sink is present.
#[inline]
pub(crate) fn report(stage: StageSink<'_>, mark: f64) {
    if let Some(s) = stage {
        s(mark);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::boolean3::{boolean_with_token, boolean_with_token_and_stage};
    use crate::linalg::Vec3;
    use crate::manifold::Manifold;
    use crate::types::OpType;
    use std::sync::Mutex;

    /// The sink hears all eight marks in order, and the boolean computes the
    /// same bits with and without it.
    #[test]
    fn stage_sink_hears_every_mark_and_changes_no_bit() {
        let a = Manifold::cube(Vec3::splat(2.0), true);
        let b = Manifold::sphere(1.2, 24);
        let heard = Mutex::new(Vec::new());
        let sink = |f: f64| heard.lock().expect("sink").push(f);
        let with =
            boolean_with_token_and_stage(a.as_impl(), b.as_impl(), OpType::Add, None, Some(&sink));
        let without = boolean_with_token(a.as_impl(), b.as_impl(), OpType::Add, None);

        assert_eq!(
            *heard.lock().expect("sink"),
            vec![
                AFTER_INTERSECT_PQ,
                AFTER_INTERSECT_QP,
                AFTER_WINDING_P,
                AFTER_WINDING_Q,
                AFTER_ASSEMBLY,
                AFTER_TRIANGULATION,
                AFTER_REFERENCE,
                AFTER_SIMPLIFY,
            ]
        );
        let bits = |m: &crate::impl_mesh::ManifoldImpl| -> Vec<u64> {
            m.vert_pos
                .iter()
                .flat_map(|v| [v.x.to_bits(), v.y.to_bits(), v.z.to_bits()])
                .collect()
        };
        assert_eq!(bits(&with), bits(&without));
        assert_eq!(with.halfedge.len(), without.halfedge.len());
        for (x, y) in with.halfedge.iter().zip(without.halfedge.iter()) {
            assert_eq!(
                (x.start_vert, x.end_vert, x.paired_halfedge),
                (y.start_vert, y.end_vert, y.paired_halfedge)
            );
        }
    }
}
