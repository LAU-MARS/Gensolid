//! Repository rules for every contributor, human or agent. Changing this file requires
//! review by the verifier owners (see `.github/CODEOWNERS`).

use std::fs;
use std::path::Path;

use gensolid_policy::{
    code_part, dependency_names, files_with_suffix, negative_exponent_literals, repo_root,
};

/// The only source file allowed to contain literal tolerance values.
const TOLERANCE_MODULE: &str = "crates/gensolid-core/src/tolerance.rs";

/// Crates that implement exact B-rep modelling and must never fall back to meshes.
const BREP_CRATES: [&str; 3] = ["gensolid-geom", "gensolid-topo", "gensolid-check"];

/// The only crates allowed to opt out of `unsafe_code = "deny"`.
const UNSAFE_CRATES: [&str; 1] = ["crates/gensolid-ffi"];

fn relative(path: &Path) -> String {
    let root = repo_root().canonicalize().unwrap_or_else(|_| repo_root());
    let path = path.canonicalize().unwrap_or_else(|_| path.to_path_buf());
    path.strip_prefix(&root)
        .unwrap_or(&path)
        .to_string_lossy()
        .replace('\\', "/")
}

fn kernel_sources() -> Vec<std::path::PathBuf> {
    let root = repo_root();
    ["crates", "bindings"]
        .iter()
        .flat_map(|dir| files_with_suffix(&root.join(dir), ".rs"))
        .collect()
}

#[test]
fn tolerances_come_from_precision_not_literals() {
    let mut violations = Vec::new();
    for path in kernel_sources() {
        let rel = relative(&path);
        if rel == TOLERANCE_MODULE {
            continue;
        }
        let text = fs::read_to_string(&path).unwrap();
        for (n, line) in text.lines().enumerate() {
            for literal in negative_exponent_literals(code_part(line)) {
                violations.push(format!("{rel}:{}: `{literal}`", n + 1));
            }
        }
    }
    assert!(
        violations.is_empty(),
        "literal tolerances must be derived from `Precision` or a `Tolerance` \
         (only {TOLERANCE_MODULE} may define them):\n{}",
        violations.join("\n")
    );
}

#[test]
fn brep_crates_do_not_depend_on_meshing() {
    let root = repo_root();
    for krate in BREP_CRATES {
        let manifest =
            fs::read_to_string(root.join("crates").join(krate).join("Cargo.toml")).unwrap();
        for dep in dependency_names(&manifest) {
            let lower = dep.to_lowercase();
            assert!(
                !(lower.contains("mesh") || lower.contains("tessel") || lower.contains("manifold")),
                "{krate} depends on `{dep}`: exact B-rep operations must not route through meshes"
            );
        }
    }
}

#[test]
fn unsafe_code_is_confined_to_the_ffi_crate() {
    for path in kernel_sources() {
        let rel = relative(&path);
        if UNSAFE_CRATES.iter().any(|c| rel.starts_with(c)) {
            continue;
        }
        let text = fs::read_to_string(&path).unwrap();
        assert!(
            !text.contains("allow(unsafe_code"),
            "{rel} opts out of the unsafe_code lint"
        );
    }
}

#[test]
fn every_crate_uses_workspace_lints() {
    let root = repo_root();
    let manifests: Vec<_> = ["crates", "bindings", "verifier"]
        .iter()
        .flat_map(|dir| files_with_suffix(&root.join(dir), "Cargo.toml"))
        .collect();
    assert!(!manifests.is_empty());
    for manifest in manifests {
        let text = fs::read_to_string(&manifest).unwrap().replace("\r\n", "\n");
        assert!(
            text.contains("[lints]\nworkspace = true"),
            "{} must inherit workspace lints",
            relative(&manifest)
        );
    }
}
