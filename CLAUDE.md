# Claude Code Guidelines — manifold-rust

## Philosophy

**Quality through iterations** - Start with correct implementations, then improve. In a porting project, every function matters.

**No stubs, no shortcuts** - Every function must be complete and production-ready. No `todo!()`, no `unimplemented!()`, no `panic!("not implemented")`, no partial implementations. If dependencies aren't ready, stop and implement them first.

**Accuracy and utility over C++ parity** - The port is complete and has its
own test suite; that suite, not the C++ reference, is the bar. Bit-identical
agreement with the C++ is no longer a goal in itself. Judge every change on
accuracy, robustness, performance, and usefulness to callers. Do not pursue
a fix whose only benefit is matching the C++ more closely, and say plainly
when a proposed change is parity-only. The C++ remains the best reference
for *how* an algorithm is meant to work, and the first place to look when
porting something new or debugging.

**Output changes are documented.** A downstream consumer (manifold-sharp, a
C# twin) verifies against this tree bit-for-bit, so any change that moves
shipped output must say so: record it in `docs/CPP_DIVERGENCES.md` (what
changed, why, the evidence, which expected test values moved) with a note of
what manifold-sharp must mirror. Don't churn output for convenience — move it
only when the result is measurably better.

## Test-First Bug Fixing (Critical Practice)

When a bug is reported, always follow this workflow:

1. **Write a reproducing test first** — create a failing test
2. **Fix the bug** — minimal change to address the root cause
3. **Verify via passing test** — the test must now pass

Do not skip the reproducing test. Even if the fix seems obvious.

## Testing

- Tests MUST test actual production code, not copies
- Pin exact expected values (counts, volumes, positions) where the result is deterministic, so unintended output changes are caught; when an expected value moves, explain why the new value is better before updating it
- When instrumenting against the C++ helps diagnose a bug, do it — as a debugging aid, not a requirement
- All tests must pass before a change lands
- When test failures occur, treat all failures as real bugs, resolve through instrumentation and root cause analysis, never by weakening tests

## Coding Standards

### File length
- **Hard limit: 800 lines.** Files that reach this must be refactored by splitting into
  focused modules before adding more code.
- Never reduce a file's line count by removing comments or blank lines to meet the limit —
  that is not refactoring. Split real logic into separate files/modules.
- **Exceptions:** `linalg.rs` (~2260 lines) is exempt — it is dense type definitions and
  operator overloads that are inherently verbose in Rust and would lose cohesion if split.
  The test block has been extracted to `linalg_tests.rs`. `edge_op.rs` and
  `quickhull_algo.rs` are over (~1150/~920, the latter grown by the exact-visibility test, divergence entry 11) with tightly coupled algorithms; `edge_op.rs`'s orbit-owner search and per-orbit duplicate scan already live in `edge_op_orbits.rs`, and further growth should split out more such helpers rather than raise the figure.

### Documentation
- Every file must begin with a comment block describing its purpose and how it relates
  to other modules in the project.
- Add doc comments to functions when they clarify intent, non-obvious behavior, or
  relationships to other parts of the codebase. Skip them when the function name and
  signature already tell the full story.
- Explain *why*, not *what*. A comment that restates the code is noise.

### General style
- Prefer `Result`/`Option` over `unwrap` in library code; `expect` is acceptable in
  `main` for startup failures with a clear message.
- Keep handler functions focused — if a handler grows complex, extract helpers.
- Avoid unsafe code unless there is no alternative; document every `unsafe` block.

## C++ to Rust Porting Rules

### Porting new functionality
- Start from the C++ algorithm and its edge-case handling; deviate where we can do measurably better, and document the deviation as above
- Use `f32`/`f64` where the C++ uses `float`/`double` unless a different precision is more accurate

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

When comparing against the reference helps (diagnosing a bug, porting new code):
- **C++ reference**: Build with `cmake` from `cpp-reference/manifold`, capture output
- **Rust port**: Run `cargo test -- --nocapture`, capture output
- Compare intermediate values at the stage where results diverge

## C++ Reference

The original C++ source is at `cpp-reference/manifold/` (git submodule).

## Shell

This project uses **PowerShell** on Windows for build scripts. Use `bash` syntax in Bash tool calls.

## Porting Plan

See `PORTING_PLAN.md` for the current phase status and implementation order.

## Orchestration pattern

The main session (Fable 5) acts as planner and orchestrator only and should
not write or edit code directly. All implementation is delegated to the
`implementer` subagent (`.claude/agents/implementer.md`), one scoped step at
a time. All post-change review is delegated to the `reviewer` subagent
(`.claude/agents/reviewer.md`). The main session handles only planning,
architecture decisions, and synthesizing subagent results.
