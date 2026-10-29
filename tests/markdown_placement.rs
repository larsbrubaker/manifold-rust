// tests/markdown_placement.rs — the Markdown placement gate (CLAUDE.md,
// "Documentation"). Documentation lives in the code, in the header and doc
// comments of the code it describes, where whoever works on that code reads
// it. A Markdown file is allowed only where it is not documentation of code:
//
//   - `README.md`, `CLAUDE.md` and `PORTING_PLAN.md` (the open work) at the
//     repository root;
//   - license and notice files (`LICENSE.md`, `NOTICE.md`,
//     `THIRD_PARTY_NOTICES.md`, ...) anywhere;
//   - anything under `.claude/` (agent instructions);
//   - the package READMEs `dotnet/README.md` and `ffi/README.md`;
//   - anything under a `testdata/` directory (fixtures);
//   - anything under `cpp-reference/` (the C++ reference submodule, third
//     party).
//
// Every other `.md` or `.markdown` file (any case) fails the gate, named with
// what to do about it. The walk reads the file system, so an untracked file
// counts and a new doc fails before it is committed; it skips `.git`, build
// output (`target/`, `node_modules/`, the demo's `dist/` and `wasm/pkg/`,
// .NET's `bin/` and `obj/`), the reference submodule and agent worktrees.
//
// The crate is published to crates.io, whose package has no `.git`; there the
// gate has no repository to judge and passes without checking.

use std::path::{Path, PathBuf};

/// What the gate makes of one repository path.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Placement {
    /// Not a Markdown file: the gate has nothing to say about it.
    NotMarkdown,
    /// A Markdown file in one of the allowed places.
    Allowed,
    /// A Markdown file anywhere else.
    Refused,
}

/// True for a Markdown file by its extension, whatever its case.
fn is_markdown(path: &str) -> bool {
    let name = path.rsplit('/').next().unwrap_or(path);
    match name.rsplit_once('.') {
        Some((stem, ext)) if !stem.is_empty() => {
            ext.eq_ignore_ascii_case("md") || ext.eq_ignore_ascii_case("markdown")
        }
        _ => false,
    }
}

/// True for a license or notice file by its name, in any directory.
fn is_license_or_notice(path: &str) -> bool {
    let name = path.rsplit('/').next().unwrap_or(path);
    let stem = name.rsplit_once('.').map_or(name, |(stem, _)| stem);
    let stem = stem.to_ascii_uppercase();
    ["LICENSE", "LICENCE", "COPYING"]
        .iter()
        .any(|prefix| stem.starts_with(prefix))
        || stem.ends_with("NOTICE")
        || stem.ends_with("NOTICES")
}

/// Where `path` (relative to the repository root, `/`-separated) stands.
fn placement(path: &str) -> Placement {
    if !is_markdown(path) {
        return Placement::NotMarkdown;
    }
    let under = |tree: &str| path.strip_prefix(tree).is_some_and(|rest| !rest.is_empty());
    let in_testdata = path
        .rsplit_once('/')
        .is_some_and(|(dir, _)| dir.split('/').any(|part| part == "testdata"));
    let allowed = matches!(
        path,
        "README.md" | "CLAUDE.md" | "PORTING_PLAN.md" | "dotnet/README.md" | "ffi/README.md"
    ) || is_license_or_notice(path)
        || under(".claude/")
        || under("cpp-reference/")
        || in_testdata;
    if allowed {
        Placement::Allowed
    } else {
        Placement::Refused
    }
}

/// True for a directory the walk never enters: version control, build output,
/// the third-party reference and agent worktrees. `rel` is the directory's
/// repository path.
fn skipped_dir(rel: &str) -> bool {
    let name = rel.rsplit('/').next().unwrap_or(rel);
    matches!(name, ".git" | "target" | "node_modules")
        || matches!(
            rel,
            "cpp-reference" | "demo/dist" | "demo/wasm/pkg" | ".claude/worktrees"
        )
        || (rel.starts_with("dotnet/") && matches!(name, "bin" | "obj"))
}

/// Every file under `root` the walk keeps, as sorted `/`-separated relative
/// paths.
fn repo_files(root: &Path) -> Vec<String> {
    let mut files = Vec::new();
    let mut pending: Vec<(PathBuf, String)> = vec![(root.to_path_buf(), String::new())];
    while let Some((dir, rel_dir)) = pending.pop() {
        let entries = std::fs::read_dir(&dir)
            .unwrap_or_else(|e| panic!("read the directory {}: {e}", dir.display()));
        for entry in entries {
            let entry = entry.unwrap_or_else(|e| panic!("read an entry of {}: {e}", dir.display()));
            let name = entry.file_name().to_string_lossy().into_owned();
            let rel = if rel_dir.is_empty() {
                name
            } else {
                format!("{rel_dir}/{name}")
            };
            let kind = entry
                .file_type()
                .unwrap_or_else(|e| panic!("read the type of {rel}: {e}"));
            if kind.is_dir() {
                if !skipped_dir(&rel) {
                    pending.push((entry.path(), rel));
                }
            } else {
                files.push(rel);
            }
        }
    }
    files.sort();
    files
}

/// The refused paths among the files the walk keeps under `root`.
fn misplaced(root: &Path) -> Vec<String> {
    repo_files(root)
        .into_iter()
        .filter(|path| placement(path) == Placement::Refused)
        .collect()
}

/// The gate's failure message for the refused paths, empty when there are none.
fn failure_message(refused: &[String]) -> String {
    if refused.is_empty() {
        return String::new();
    }
    let mut message = format!(
        "{} Markdown file(s) outside the places Markdown may live:\n",
        refused.len()
    );
    for path in refused {
        message.push_str(&format!("  {path}\n"));
    }
    message.push_str(
        "Documentation lives in the code. Move what it explains into the header or doc \
         comments of the code it is about, put open work in PORTING_PLAN.md, drop history \
         and status (git keeps them), then delete the file (CLAUDE.md, \"Documentation\").\n",
    );
    message
}

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// A scratch directory removed again on drop, for the walk's own tests.
struct TempTree(PathBuf);

impl TempTree {
    fn new(tag: &str) -> Self {
        let dir = std::env::temp_dir().join(format!(
            "manifold-rust-{tag}-{}-{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("create the scratch tree");
        TempTree(dir)
    }

    fn write(&self, rel: &str, text: &str) {
        let path = self.0.join(rel);
        std::fs::create_dir_all(path.parent().expect("a file has a parent"))
            .expect("create the scratch directories");
        std::fs::write(path, text).expect("write a scratch file");
    }
}

impl Drop for TempTree {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

#[test]
fn markdown_lives_only_in_the_allowed_places() {
    let root = repo_root();
    // A published package has no repository to judge (see the header).
    if !root.join(".git").exists() {
        return;
    }
    let message = failure_message(&misplaced(&root));
    assert!(message.is_empty(), "{message}");
}

#[test]
fn markdown_gate_sees_the_allowed_files() {
    // Guards against a walk that silently checks nothing (wrong root,
    // over-broad skips).
    let root = repo_root();
    if !root.join(".git").exists() {
        return;
    }
    let files = repo_files(&root);
    for expected in [
        "CLAUDE.md",
        "README.md",
        "PORTING_PLAN.md",
        "ffi/README.md",
        "src/lib.rs",
    ] {
        assert!(
            files.iter().any(|p| p == expected),
            "the walk missed {expected}"
        );
    }
}

#[test]
fn markdown_placement_allows_the_listed_places() {
    for path in [
        "README.md",
        "CLAUDE.md",
        "PORTING_PLAN.md",
        "LICENSE.md",
        "THIRD_PARTY_NOTICES.md",
        "dotnet/ManifoldRust/LICENSE.md",
        ".claude/agents/implementer.md",
        "dotnet/README.md",
        "ffi/README.md",
        "src/testdata/notes.md",
        "src/robust/testdata/README.md",
        "cpp-reference/manifold/README.md",
        "cpp-reference/manifold/docs/usage.md",
    ] {
        assert_eq!(placement(path), Placement::Allowed, "{path}");
    }
}

#[test]
fn markdown_placement_refuses_docs_anywhere_else() {
    for path in [
        "docs/ROBUST_ENGINE_STATUS.md",
        "docs/CPP_DIVERGENCES.md",
        "CONTRIBUTING.md",
        "src/README.md",
        "src/robust/NOTES.md",
        "demo/README.md",
        "dotnet/ManifoldRust/README.md",
        "ffi/src/README.md",
        "examples/README.md",
        // Not a testdata directory, only a name that contains it.
        "src/testdata.md",
        "src/my-testdata/notes.md",
        "cpp-reference.md",
    ] {
        assert_eq!(placement(path), Placement::Refused, "{path}");
    }
    // Not Markdown, so not this gate's concern.
    for path in ["src/lib.rs", "docs/booleans2009.pdf", "README", ".md"] {
        assert_eq!(placement(path), Placement::NotMarkdown, "{path}");
    }
}

#[test]
fn markdown_placement_matches_extensions_in_any_case() {
    for path in [
        "notes.MD",
        "notes.Md",
        "notes.markdown",
        "docs/guide.MARKDOWN",
    ] {
        assert_eq!(placement(path), Placement::Refused, "{path}");
    }
    assert_eq!(placement("Readme.MD"), Placement::Refused);
    assert_eq!(placement("license.MD"), Placement::Allowed);
}

#[test]
fn markdown_gate_reads_untracked_files_and_skips_build_trees() {
    let tree = TempTree::new("markdown");
    tree.write("README.md", "readme\n");
    tree.write("docs/old-status.md", "a stray doc\n");
    tree.write("src/robust/notes.markdown", "another\n");
    tree.write("src/testdata/fixture.md", "fixture\n");
    for skipped in [
        ".git/info/notes.md",
        "target/doc/README.md",
        "demo/node_modules/pkg/README.md",
        "demo/dist/README.md",
        "demo/wasm/pkg/README.md",
        "dotnet/ManifoldRust/bin/Release/README.md",
        "dotnet/ManifoldRust/obj/notes.md",
        "cpp-reference/manifold/CONTRIBUTING.md",
        ".claude/worktrees/slot1/notes.md",
    ] {
        tree.write(skipped, "skipped\n");
    }
    assert_eq!(
        misplaced(&tree.0),
        ["docs/old-status.md", "src/robust/notes.markdown"]
    );
}

#[test]
fn markdown_gate_failure_message_says_what_to_do() {
    let message = failure_message(&["docs/old-status.md".to_string()]);
    assert!(message.contains("  docs/old-status.md\n"));
    assert!(message.contains("into the header or doc comments of the code it is about"));
    assert!(message.contains("PORTING_PLAN.md"));
    assert!(message.contains("then delete the file"));
    assert!(failure_message(&[]).is_empty());
}
