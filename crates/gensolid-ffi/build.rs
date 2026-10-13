//! Optional native SMLib build driver; default builds require no C++ toolchain.
use std::{env, error::Error, path::PathBuf, process::Command};

fn main() -> Result<(), Box<dyn Error>> {
    println!("cargo:rerun-if-env-changed=SMLIB_SOURCE_DIR");
    println!("cargo:rerun-if-env-changed=CXX");
    println!("cargo:rerun-if-env-changed=SMLIB_BUILD_DIR");
    println!("cargo:rerun-if-changed=native");
    println!("cargo:rerun-if-changed=../../tools/build-smlib.py");
    if env::var_os("CARGO_FEATURE_SMLIB").is_none() {
        return Ok(());
    }
    // Browser SMLib is packaged by tools/build-smlib.py --wasm as a companion module.
    if env::var("TARGET")? == "wasm32-unknown-unknown" {
        return Ok(());
    }
    if env::var("HOST")? != env::var("TARGET")? {
        return Err("SMLib cross-compilation is not supported by this builder".into());
    }
    let source = env::var_os("SMLIB_SOURCE_DIR").ok_or(
        "smlib feature requires SMLIB_SOURCE_DIR pointing to the pinned usd-brep checkout",
    )?;
    println!(
        "cargo:rerun-if-changed={}",
        PathBuf::from(&source).join("source/SMLib").display()
    );
    println!(
        "cargo:rerun-if-changed={}",
        PathBuf::from(&source).join("source/SM_API").display()
    );
    let out = PathBuf::from(
        env::var_os("SMLIB_BUILD_DIR")
            .or_else(|| env::var_os("OUT_DIR"))
            .ok_or("OUT_DIR missing")?,
    );
    if !out.is_absolute() {
        return Err("SMLIB_BUILD_DIR must be an absolute path".into());
    }
    let status = Command::new("python3")
        .args(["../../tools/build-smlib.py", "--source"])
        .arg(source)
        .arg("--output")
        .arg(&out)
        .status()?;
    if !status.success() {
        return Err("SMLib native compilation failed; see compiler output above".into());
    }
    println!("cargo:rustc-link-search=native={}", out.display());
    println!("cargo:rustc-link-lib=static=gensolid_smlib");
    println!(
        "cargo:rustc-link-lib={}",
        if env::var("CARGO_CFG_TARGET_OS")? == "macos" {
            "c++"
        } else {
            "stdc++"
        }
    );
    Ok(())
}
