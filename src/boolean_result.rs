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

// Phase 11: Boolean Result — face assembly from intersection data
//
// C++ source: src/boolean_result.cpp (889 lines)
//
// STATUS: NOT YET IMPLEMENTED
// Requires full port of boolean_result.cpp which takes intersection data
// from Boolean3 and assembles the output mesh faces.

use crate::boolean3;
use crate::impl_mesh::ManifoldImpl;
use crate::types::OpType;

/// Thin wrapper that delegates to boolean3::boolean().
/// Will be replaced with proper BooleanResult when boolean3.cpp is fully ported.
#[derive(Clone)]
pub struct BooleanResult {
    pub result: ManifoldImpl,
}

impl BooleanResult {
    pub fn new(a: &ManifoldImpl, b: &ManifoldImpl, op: OpType) -> Self {
        Self {
            result: boolean3::boolean(a, b, op),
        }
    }

    pub fn result(self) -> ManifoldImpl {
        self.result
    }
}
