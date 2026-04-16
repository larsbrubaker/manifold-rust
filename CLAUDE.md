# Claude Code Guidelines — manifold-rust

## Philosophy

**Quality through iterations** - Start with correct implementations, then improve. In a porting project, every function matters.

**No stubs, no shortcuts** - Every function must be complete and production-ready. No `todo!()`, no `unimplemented!()`, no `panic!("not implemented")`, no partial implementations. If dependencies aren't ready, stop and implement them first.

**Exact numerical match** - This port must reproduce the C++ library with perfect numerical exactness. Not "close enough" — identical results on identical inputs.

## Test-First Bug Fixing (Critical Practice)

When a bug is reported, always follow this workflow:

1. **Write a reproducing test first** — create a failing test
2. **Fix the bug** — minimal change to address the root cause
3. **Verify via passing test** — the test must now pass

Do not skip the reproducing test. Even if the fix seems obvious.

## Testing

- Tests MUST test actual production code, not copies
- Tests must verify **exact behavioral match** with the C++ implementation — same floating-point results, same triangle counts, same vertex positions
- To confirm numerical exactness, instrument both the Rust and C++ implementations and compare output byte-for-byte where applicable
- All tests must pass before advancing to the next phase
- When test failures occur, use the fix-test-failures skill — treat all failures as real bugs, resolve through instrumentation and root cause analysis, never by weakening tests

**Running tests:**
```bash
cargo test
cargo test --lib vec_tests
cargo test test_name -- --exact
cargo test -- --nocapture
```

## Coding Standards

### File length
- **Hard limit: 800 lines.** Files that reach this must be refactored by splitting into
  focused modules before adding more code.
- Never reduce a file's line count by removing comments or blank lines to meet the limit —
  that is not refactoring. Split real logic into separate files/modules.
- **Exceptions:** `linalg.rs` (~2260 lines) is exempt — it is dense type definitions and
  operator overloads that are inherently verbose in Rust and would lose cohesion if split.
  The test block has been extracted to `linalg_tests.rs`. `edge_op.rs` and
  `quickhull_algo.rs` are slightly over (~900/~820) with tightly coupled algorithms.

### Documentation
- Every file must begin with a comment block describing its purpose and how it relates
  to other modules in the project.
- Add doc comments to functions when they clarify intent, non-obvious behavior, or
  relationships to other parts of the codebase. Skip them when the function name and
  signature already tell the full story.
- Explain *why*, not *what*. A comment that restates the code is noise.

### Bug workflow — always follow this order
1. **Write a failing test** that reproduces the bug.
2. **Fix the bug.**
3. **Confirm the test passes** (`cargo test`).

Never commit a bug fix that isn't covered by a test.

### General style
- Prefer `Result`/`Option` over `unwrap` in library code; `expect` is acceptable in
  `main` for startup failures with a clear message.
- Keep handler functions focused — if a handler grows complex, extract helpers.
- Avoid unsafe code unless there is no alternative; document every `unsafe` block.

### Names
Follow Rust conventions (`snake_case` for functions/variables, `PascalCase` for types, `SCREAMING_SNAKE_CASE` for constants). Mirror C++ names where they are clear; use better Rust names where they aren't.

### Comments
Explain *why*, not *what*. When porting, note where Rust differs from C++ and why.

### Refactoring
Improve code when it serves a purpose, not for aesthetics.

## C++ to Rust Porting Rules

### Exact Behavioral Matching
- Same algorithms, same mathematical precision (IEEE 754 float/double semantics)
- Same edge case handling
- No "close enough" implementations
- Use `f32`/`f64` to match C++ `float`/`double`; never round or clamp unless the C++ does

### Dependency-Ordered Implementation
Before implementing any function:
1. Read the corresponding C++ source to identify all functions called by the target
2. Verify all dependencies are implemented and tested in Rust
3. If any dependency is incomplete, implement it first

### C++ Patterns → Rust Patterns
| C++ | Rust |
|-----|------|
| Raw pointer linked lists | Arena allocation with `Vec<T>` + index |
| `std::vector<T>` | `Vec<T>` |
| `std::array<T, N>` | `[T; N]` |
| Template specialization | Trait impls |
| `thrust::` parallel algorithms | Rayon (or sequential first, parallel later) |
| `glm::` math | Our own `vec.rs` / `linalg.rs` |
| `static_cast<int>(x)` truncation | `x as i32` (truncates, same as C++) |

### Forbidden Patterns
- `todo!()` or `unimplemented!()` macros
- `panic!()` for missing functionality
- Stub functions or placeholder implementations
- Marking phases complete when any test fails
- Weakening tests to make them pass

## Instrumentation Strategy

To validate numerical exactness, we build and run both implementations:
- **C++ reference**: Build with `cmake` from `cpp-reference/manifold`, capture output
- **Rust port**: Run `cargo test -- --nocapture`, capture output
- Compare intermediate values at each phase boundary

## C++ Reference

The original C++ source is at `cpp-reference/manifold/` (git submodule).
- Public API: `cpp-reference/manifold/include/manifold/`
- Implementation: `cpp-reference/manifold/src/`
- Tests: `cpp-reference/manifold/test/`

## Shell

This project uses **PowerShell** on Windows for build scripts. Use `bash` syntax in Bash tool calls.

## Porting Plan

See `PORTING_PLAN.md` for the current phase status and implementation order.
