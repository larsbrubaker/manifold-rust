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

// manifold_boolean.rs — Manifold's boolean operations: the two-operand
// booleans with their engine / winding-rule / cancellation / progress
// variants, the pairwise batch fold, plane splits and trims built on a
// half-space cutter, and the `+` / `-` / `^` operator overloads (C++
// `operator+`, `operator-`, `operator^`). A child module of manifold.rs so
// it can reach the private `imp` field; the kernels are in boolean3.rs and,
// for the robust engine, src/robust.

use super::Manifold;
use crate::boolean3;
use crate::linalg::{normalize, Vec3};
use crate::math;
use crate::types::{Error, OpType};

impl Manifold {
    /// Split this manifold into two using a cutter manifold.
    /// Returns (intersection, difference).
    pub fn split(&self, cutter: &Self) -> (Self, Self) {
        let intersection = self.intersection(cutter);
        let difference = self.difference(cutter);
        (intersection, difference)
    }

    /// Split this manifold by a plane defined by a normal and offset from origin.
    /// Returns (in direction of normal, opposite direction).
    pub fn split_by_plane(&self, normal: Vec3, origin_offset: f64) -> (Self, Self) {
        // Per C++ #1659: errored manifolds are empty, so the is_empty()
        // early-return below would silently drop their status — guard first.
        if self.imp.status != Error::NoError {
            return (self.clone(), self.clone());
        }
        if self.is_empty() {
            return (Self::empty(), Self::empty());
        }
        let halfspace = Self::halfspace(&self.imp.bbox, normal, origin_offset);
        self.split(&halfspace)
    }

    /// Trim this manifold by a half-space, keeping only the part in the direction
    /// of the normal vector.
    pub fn trim_by_plane(&self, normal: Vec3, origin_offset: f64) -> Self {
        if self.is_empty() {
            return Self::empty();
        }
        let halfspace = Self::halfspace(&self.imp.bbox, normal, origin_offset);
        self.intersection(&halfspace)
    }

    /// Apply batch boolean operations on a list of manifolds.
    pub fn batch_boolean(manifolds: &[Self], op: OpType) -> Self {
        if manifolds.is_empty() {
            return Self::empty();
        }
        let mut result = manifolds[0].clone();
        for m in &manifolds[1..] {
            result = result.boolean(m, op);
        }
        result
    }

    /// Internal helper: create a halfspace (large cube) for plane splitting.
    fn halfspace(bbox: &crate::types::Box, normal: Vec3, origin_offset: f64) -> Self {
        let n = normalize(normal);
        let cutter = Self::cube(Vec3::splat(2.0), true).translate(Vec3::new(1.0, 0.0, 0.0));
        let center = bbox.center();
        let size_len = (bbox.size().x * bbox.size().x
            + bbox.size().y * bbox.size().y
            + bbox.size().z * bbox.size().z)
            .sqrt();
        let dist = ((center.x - n.x * origin_offset).powi(2)
            + (center.y - n.y * origin_offset).powi(2)
            + (center.z - n.z * origin_offset).powi(2))
        .sqrt()
            + 0.5 * size_len;
        let cutter = cutter
            .scale(Vec3::splat(dist))
            .translate(Vec3::new(origin_offset, 0.0, 0.0));
        let y_deg = -math::asin(n.z).to_degrees();
        let z_deg = math::atan2(n.y, n.x).to_degrees();
        cutter.rotate(0.0, y_deg, z_deg)
    }

    pub fn boolean(&self, other: &Self, op: OpType) -> Self {
        self.boolean_with_engine(other, op, crate::types::BooleanConfig::default_engine())
    }

    /// [`Manifold::boolean`] with an explicit engine choice, overriding the
    /// process-global default set via
    /// [`crate::types::BooleanConfig::set_default_engine`].
    pub fn boolean_with_engine(
        &self,
        other: &Self,
        op: OpType,
        engine: crate::types::BooleanEngine,
    ) -> Self {
        Self::from_impl(boolean3::boolean_dispatch(
            &self.imp, &other.imp, op, engine, None,
        ))
    }

    /// [`Manifold::boolean_with_engine`] with cooperative cancellation.
    pub fn boolean_with_engine_and_token(
        &self,
        other: &Self,
        op: OpType,
        engine: crate::types::BooleanEngine,
        token: Option<&crate::cancel::CancelToken>,
    ) -> Self {
        Self::from_impl(boolean3::boolean_dispatch(
            &self.imp, &other.imp, op, engine, token,
        ))
    }

    /// [`Manifold::boolean_with_engine_and_token`] that also reports coarse
    /// pipeline progress.
    ///
    /// Cancellation and progress travel together because callers that want one
    /// almost always want the other (a UI showing a progress bar next to a
    /// cancel button); pass `None` for either independently. `None` progress is
    /// byte-for-byte the un-instrumented path — see [`crate::progress`] for the
    /// phases reported and the throttling contract.
    pub fn boolean_with_engine_and_progress(
        &self,
        other: &Self,
        op: OpType,
        engine: crate::types::BooleanEngine,
        token: Option<&crate::cancel::CancelToken>,
        progress: Option<&crate::progress::ProgressReporter>,
    ) -> Self {
        Self::from_impl(boolean3::boolean_dispatch_with_progress(
            &self.imp, &other.imp, op, engine, token, progress,
        ))
    }

    /// [`Manifold::boolean_with_engine`] with an explicit winding rule.
    ///
    /// [`crate::types::WindingRule::Nonzero`] treats inside-out geometry as
    /// solid (`w != 0` rather than `w >= 1`), which keeps the inverted regions
    /// of inconsistently wound scans instead of dropping them. The rule is a
    /// robust-engine semantic: the exact engine ignores it, and `Auto` routes
    /// to the robust engine whenever the rule is `Nonzero` (see
    /// [`crate::boolean3::boolean_dispatch_full`]).
    pub fn boolean_with_engine_and_rule(
        &self,
        other: &Self,
        op: OpType,
        engine: crate::types::BooleanEngine,
        rule: crate::types::WindingRule,
    ) -> Self {
        self.boolean_with_engine_rule_and_progress(other, op, engine, rule, None, None)
    }

    /// The full per-call boolean path: engine, winding rule, cancellation, and
    /// progress. Every other boolean entry point on `Manifold` is this one with
    /// defaults filled in.
    pub fn boolean_with_engine_rule_and_progress(
        &self,
        other: &Self,
        op: OpType,
        engine: crate::types::BooleanEngine,
        rule: crate::types::WindingRule,
        token: Option<&crate::cancel::CancelToken>,
        progress: Option<&crate::progress::ProgressReporter>,
    ) -> Self {
        Self::from_impl(boolean3::boolean_dispatch_full(
            &self.imp, &other.imp, op, engine, rule, token, progress,
        ))
    }

    /// [`Manifold::batch_boolean`] with an explicit engine choice (pairwise
    /// left fold, like `batch_boolean`).
    pub fn batch_boolean_with_engine(
        manifolds: &[Self],
        op: OpType,
        engine: crate::types::BooleanEngine,
    ) -> Self {
        if manifolds.is_empty() {
            return Self::empty();
        }
        let mut result = manifolds[0].clone();
        for m in &manifolds[1..] {
            result = result.boolean_with_engine(m, op, engine);
        }
        result
    }

    pub fn union_with_engine(&self, other: &Self, engine: crate::types::BooleanEngine) -> Self {
        self.boolean_with_engine(other, OpType::Add, engine)
    }

    pub fn difference_with_engine(
        &self,
        other: &Self,
        engine: crate::types::BooleanEngine,
    ) -> Self {
        self.boolean_with_engine(other, OpType::Subtract, engine)
    }

    pub fn intersection_with_engine(
        &self,
        other: &Self,
        engine: crate::types::BooleanEngine,
    ) -> Self {
        self.boolean_with_engine(other, OpType::Intersect, engine)
    }

    /// [`Manifold::boolean`] with cooperative cancellation.
    ///
    /// Pass `None` for the uncancellable behaviour of [`Manifold::boolean`] —
    /// that path is unchanged and touches no atomics. With `Some(token)`, a
    /// cancel requested from any thread (before or during the call) makes this
    /// return an empty manifold whose [`Manifold::status`] is
    /// [`Error::Cancelled`], mirroring the C++ `ExecutionContext` contract.
    pub fn boolean_with_token(
        &self,
        other: &Self,
        op: OpType,
        token: Option<&crate::cancel::CancelToken>,
    ) -> Self {
        Self::from_impl(boolean3::boolean_with_token(
            &self.imp, &other.imp, op, token,
        ))
    }

    pub fn union(&self, other: &Self) -> Self {
        self.boolean(other, OpType::Add)
    }

    pub fn difference(&self, other: &Self) -> Self {
        self.boolean(other, OpType::Subtract)
    }

    pub fn intersection(&self, other: &Self) -> Self {
        self.boolean(other, OpType::Intersect)
    }
}

// Operator overloads: + for union, - for difference, ^ for intersection
// Matches C++ operator+(Manifold), operator-(Manifold), operator^(Manifold)

impl std::ops::Add for Manifold {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        self.union(&rhs)
    }
}

impl std::ops::Add<&Manifold> for Manifold {
    type Output = Self;
    fn add(self, rhs: &Self) -> Self {
        self.union(rhs)
    }
}

impl std::ops::Add<&Manifold> for &Manifold {
    type Output = Manifold;
    fn add(self, rhs: &Manifold) -> Manifold {
        self.union(rhs)
    }
}

impl std::ops::AddAssign for Manifold {
    fn add_assign(&mut self, rhs: Self) {
        *self = self.union(&rhs);
    }
}

impl std::ops::AddAssign<&Manifold> for Manifold {
    fn add_assign(&mut self, rhs: &Self) {
        *self = self.union(rhs);
    }
}

impl std::ops::Sub for Manifold {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        self.difference(&rhs)
    }
}

impl std::ops::Sub<&Manifold> for Manifold {
    type Output = Self;
    fn sub(self, rhs: &Self) -> Self {
        self.difference(rhs)
    }
}

impl std::ops::Sub<&Manifold> for &Manifold {
    type Output = Manifold;
    fn sub(self, rhs: &Manifold) -> Manifold {
        self.difference(rhs)
    }
}

impl std::ops::SubAssign for Manifold {
    fn sub_assign(&mut self, rhs: Self) {
        *self = self.difference(&rhs);
    }
}

impl std::ops::SubAssign<&Manifold> for Manifold {
    fn sub_assign(&mut self, rhs: &Self) {
        *self = self.difference(rhs);
    }
}

impl std::ops::BitXor for Manifold {
    type Output = Self;
    fn bitxor(self, rhs: Self) -> Self {
        self.intersection(&rhs)
    }
}

impl std::ops::BitXor<&Manifold> for Manifold {
    type Output = Self;
    fn bitxor(self, rhs: &Self) -> Self {
        self.intersection(rhs)
    }
}

impl std::ops::BitXor<&Manifold> for &Manifold {
    type Output = Manifold;
    fn bitxor(self, rhs: &Manifold) -> Manifold {
        self.intersection(rhs)
    }
}

impl std::ops::BitXorAssign for Manifold {
    fn bitxor_assign(&mut self, rhs: Self) {
        *self = self.intersection(&rhs);
    }
}

impl std::ops::BitXorAssign<&Manifold> for Manifold {
    fn bitxor_assign(&mut self, rhs: &Self) {
        *self = self.intersection(rhs);
    }
}
