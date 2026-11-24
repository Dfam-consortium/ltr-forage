// Sets LTRFORAGE_VERSION to "<package version> (<commit>)", or "<package version> (<commit>-dirty)" when a
// file that builds the binary has uncommitted changes, so `--version` says which code a binary was built
// from (a stale target/release once ran in place of the current code without a sign).  Outside a git
// checkout (a source tarball) it is the package version alone.  After rmblast's build-support/git_build_id.rs.

use std::path::{Path, PathBuf};
use std::process::Command;

/// Paths, relative to the crate root, whose content goes into the binary.
const BUILD_INPUTS: &[&str] = &["Cargo.toml", "Cargo.lock", "build.rs", "src"];

fn git(root: &Path, args: &[&str]) -> Option<String> {
    let out = Command::new("git").arg("-C").arg(root).args(args).output().ok()?;
    out.status.success().then(|| String::from_utf8_lossy(&out.stdout).trim().to_string())
}

fn main() {
    let root = PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").expect("CARGO_MANIFEST_DIR"));
    let pkg = std::env::var("CARGO_PKG_VERSION").expect("CARGO_PKG_VERSION");
    for p in BUILD_INPUTS {
        println!("cargo:rerun-if-changed={}", root.join(p).display());
    }
    let version = match git(&root, &["rev-parse", "--short=7", "HEAD"]) {
        None => pkg,
        Some(commit) => {
            // rerun when HEAD moves: a commit, checkout or reset
            if let Some(dir) = git(&root, &["rev-parse", "--absolute-git-dir"]) {
                let dir = PathBuf::from(dir);
                println!("cargo:rerun-if-changed={}", dir.join("HEAD").display());
                println!("cargo:rerun-if-changed={}", dir.join("packed-refs").display());
                if let Some(r) = git(&root, &["symbolic-ref", "-q", "HEAD"]) {
                    println!("cargo:rerun-if-changed={}", dir.join(r).display());
                }
            }
            let mut args = vec!["status", "--porcelain", "--"];
            args.extend_from_slice(BUILD_INPUTS);
            let dirty = git(&root, &args).map(|s| !s.is_empty()).unwrap_or(false);
            format!("{pkg} ({commit}{})", if dirty { "-dirty" } else { "" })
        }
    };
    println!("cargo:rustc-env=LTRFORAGE_VERSION={version}");
}
