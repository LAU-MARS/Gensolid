//! Frozen corpus selection and byte-integrity checks, independent of body validity.

use std::collections::BTreeMap;
use std::fs;
use std::path::{Component, Path, PathBuf};

use gensolid_core::{Error, Result};
use serde::Deserialize;
use sha2::{Digest, Sha256};

#[derive(Deserialize)]
struct Manifest {
    files: Vec<Entry>,
}

#[derive(Deserialize)]
struct Entry {
    file: String,
    sha256: String,
    bytes: Option<usize>,
}

#[derive(Default)]
pub(crate) struct FrozenInputs {
    pub(crate) expected: BTreeMap<PathBuf, (String, Option<usize>)>,
    pub(crate) errors: Vec<String>,
}

pub(crate) fn digest(bytes: &[u8]) -> String {
    Sha256::digest(bytes)
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect()
}

impl FrozenInputs {
    pub(crate) fn read(root: &Path, files: &mut Vec<PathBuf>) -> Result<Self> {
        let mut result = Self::default();
        result.visit(root)?;
        // Expected CAD inputs stay in the workload even if absent from disk.
        files.extend(
            result
                .expected
                .keys()
                .filter(|path| crate::format(path).is_some())
                .cloned(),
        );
        Ok(result)
    }

    fn visit(&mut self, directory: &Path) -> Result<()> {
        for entry in fs::read_dir(directory).map_err(|error| Error::parse(error.to_string()))? {
            let entry = entry.map_err(|error| Error::parse(error.to_string()))?;
            let kind = entry
                .file_type()
                .map_err(|error| Error::parse(error.to_string()))?;
            if kind.is_dir() {
                self.visit(&entry.path())?;
            } else if kind.is_file() && entry.file_name() == "manifest.json" {
                match self.load(&entry.path()) {
                    Ok(()) => {}
                    Err(error) => self
                        .errors
                        .push(format!("{}: {error}", entry.path().display())),
                }
            }
        }
        Ok(())
    }

    fn load(&mut self, path: &Path) -> Result<()> {
        let parent = path
            .parent()
            .ok_or_else(|| Error::parse("manifest has no parent directory"))?;
        let text = fs::read_to_string(path).map_err(|error| Error::parse(error.to_string()))?;
        let manifest: Manifest =
            serde_json::from_str(&text).map_err(|error| Error::parse(error.to_string()))?;
        if manifest.files.is_empty() {
            return Err(Error::parse("manifest contains no expected files"));
        }
        for entry in manifest.files {
            let relative = Path::new(&entry.file);
            if relative.as_os_str().is_empty()
                || relative
                    .components()
                    .any(|component| !matches!(component, Component::Normal(_)))
            {
                return Err(Error::parse(
                    "manifest paths must stay relative to their directory",
                ));
            }
            if entry.sha256.len() != 64
                || !entry.sha256.bytes().all(|byte| byte.is_ascii_hexdigit())
            {
                return Err(Error::parse(format!(
                    "{} has an invalid SHA-256 digest",
                    entry.file
                )));
            }
            let target = parent.join(relative);
            if target.is_file() {
                let resolved =
                    fs::canonicalize(&target).map_err(|error| Error::parse(error.to_string()))?;
                let allowed =
                    fs::canonicalize(parent).map_err(|error| Error::parse(error.to_string()))?;
                if !resolved.starts_with(allowed) {
                    return Err(Error::parse("manifest file link escapes its directory"));
                }
            }
            if self
                .expected
                .insert(target, (entry.sha256, entry.bytes))
                .is_some()
            {
                return Err(Error::parse(format!(
                    "duplicate expected file {}",
                    entry.file
                )));
            }
        }
        Ok(())
    }

    pub(crate) fn violations(&self) -> BTreeMap<PathBuf, String> {
        self.expected
            .iter()
            .filter_map(|(path, (expected, bytes))| {
                let issue = match fs::read(path) {
                    Ok(content) if !digest(&content).eq_ignore_ascii_case(expected) => {
                        Some("SHA-256 differs from the frozen manifest".to_owned())
                    }
                    Ok(content) if bytes.is_some_and(|expected| expected != content.len()) => {
                        Some("byte count differs from the frozen manifest".to_owned())
                    }
                    Ok(_) => None,
                    Err(error) => Some(format!("expected file is unreadable: {error}")),
                };
                issue.map(|issue| (path.clone(), issue))
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deleting_or_altering_a_frozen_input_cannot_improve_coverage() -> Result<()> {
        let suffix = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_err(|error| Error::internal(error.to_string()))?
            .as_nanos();
        let directory =
            std::env::temp_dir().join(format!("gensolid-manifest-{}-{suffix}", std::process::id()));
        fs::create_dir(&directory).map_err(|error| Error::internal(error.to_string()))?;
        let result = (|| {
            let original = b"original frozen input";
            let file = directory.join("a.step");
            fs::write(&file, original).map_err(|error| Error::internal(error.to_string()))?;
            let manifest = serde_json::json!({"files":[{"file":"a.step","sha256":digest(original),"bytes":original.len()}]});
            fs::write(directory.join("manifest.json"), manifest.to_string())
                .map_err(|error| Error::internal(error.to_string()))?;
            let mut files = Vec::new();
            let frozen = FrozenInputs::read(&directory, &mut files)?;
            assert_eq!(files, vec![file.clone()]);
            assert!(frozen.violations().is_empty());
            fs::write(&file, b"altered").map_err(|error| Error::internal(error.to_string()))?;
            assert_eq!(frozen.violations().len(), 1);
            fs::remove_file(&file).map_err(|error| Error::internal(error.to_string()))?;
            assert_eq!(frozen.violations().len(), 1);
            let mut absent = Vec::new();
            FrozenInputs::read(&directory, &mut absent)?;
            assert_eq!(absent, vec![file]);
            Ok(())
        })();
        fs::remove_dir_all(&directory).map_err(|error| Error::internal(error.to_string()))?;
        result
    }
}
