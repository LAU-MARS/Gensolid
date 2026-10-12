//! Repository policy checks. The rules live in `tests/policy.rs`; this library only
//! provides the source scanning helpers they use.

use std::fs;
use std::path::{Path, PathBuf};

/// Repository root.
pub fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..").join("..")
}

/// All files under `dir` (recursively) whose name ends with `suffix`.
pub fn files_with_suffix(dir: &Path, suffix: &str) -> Vec<PathBuf> {
    let mut out = Vec::new();
    let mut stack = vec![dir.to_path_buf()];
    while let Some(current) = stack.pop() {
        let Ok(entries) = fs::read_dir(&current) else {
            continue;
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                if path.file_name().is_some_and(|n| n != "target") {
                    stack.push(path);
                }
            } else if path.to_string_lossy().ends_with(suffix) {
                out.push(path);
            }
        }
    }
    out.sort();
    out
}

/// The code part of a line: everything before a `//` comment.
pub fn code_part(line: &str) -> &str {
    line.find("//").map_or(line, |i| &line[..i])
}

/// Float literals with a negative exponent (`1e-8`, `2.5E-12`, `1_0e-3`) in `code`.
pub fn negative_exponent_literals(code: &str) -> Vec<&str> {
    let bytes = code.as_bytes();
    let mut found = Vec::new();
    let mut i = 1;
    while i + 2 < bytes.len() {
        let is_exponent = matches!(bytes[i], b'e' | b'E')
            && bytes[i + 1] == b'-'
            && bytes[i + 2].is_ascii_digit()
            && (bytes[i - 1].is_ascii_digit() || bytes[i - 1] == b'.');
        if is_exponent {
            let mut start = i;
            while start > 0
                && (bytes[start - 1].is_ascii_digit() || matches!(bytes[start - 1], b'.' | b'_'))
            {
                start -= 1;
            }
            let mut end = i + 2;
            while end < bytes.len() && bytes[end].is_ascii_digit() {
                end += 1;
            }
            let preceded_by_ident =
                start > 0 && (bytes[start - 1].is_ascii_alphanumeric() || bytes[start - 1] == b'_');
            if bytes[start].is_ascii_digit() && !preceded_by_ident {
                found.push(&code[start..end]);
            }
            i = end;
        } else {
            i += 1;
        }
    }
    found
}

/// Dependency names listed in the `[dependencies]` table of a manifest.
pub fn dependency_names(manifest: &str) -> Vec<String> {
    let mut in_dependencies = false;
    let mut names = Vec::new();
    for line in manifest.lines().map(str::trim) {
        if line.starts_with('[') {
            in_dependencies = line == "[dependencies]";
            continue;
        }
        if in_dependencies && let Some((name, _)) = line.split_once(['=', '.']) {
            let name = name.trim();
            if !name.is_empty() && !name.starts_with('#') {
                names.push(name.to_owned());
            }
        }
    }
    names
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finds_negative_exponents_only_in_literals() {
        assert_eq!(
            negative_exponent_literals("let a = 1e-8 + 2.5E-12;"),
            vec!["1e-8", "2.5E-12"]
        );
        assert!(negative_exponent_literals("let x = 1e3; let name_e-1 = 0;").is_empty());
        assert!(negative_exponent_literals("let v = vec1e-2;").is_empty());
    }

    #[test]
    fn reads_dependency_names() {
        let manifest = "[package]\nname = \"x\"\n[dependencies]\nfoo.workspace = true\nbar = \"1\"\n[dev-dependencies]\nbaz = \"1\"\n";
        assert_eq!(dependency_names(manifest), vec!["foo", "bar"]);
    }
}
