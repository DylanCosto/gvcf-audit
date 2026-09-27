#![allow(dead_code)]

use serde_json::Value;
use std::{
    fs,
    path::{Path, PathBuf},
    process::{Command, Output},
    sync::atomic::{AtomicU64, Ordering},
};

static NEXT: AtomicU64 = AtomicU64::new(0);

pub struct Workspace(pub PathBuf);

impl Workspace {
    pub fn new() -> Self {
        loop {
            let path = std::env::temp_dir().join(format!(
                "gvcf-audit-test-{}-{}",
                std::process::id(),
                NEXT.fetch_add(1, Ordering::Relaxed)
            ));
            match fs::create_dir(&path) {
                Ok(()) => return Self(path),
                Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(e) => panic!("create test directory: {e}"),
            }
        }
    }

    pub fn run(&self, args: &[&str], code: i32) -> Output {
        let output = Command::new(env!("CARGO_BIN_EXE_gvcf-audit"))
            .current_dir(&self.0)
            .args(args)
            .output()
            .expect("run gvcf-audit");
        assert_eq!(
            output.status.code(),
            Some(code),
            "{args:?}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        for entry in fs::read_dir(&self.0).unwrap() {
            assert!(
                !entry
                    .unwrap()
                    .file_name()
                    .to_string_lossy()
                    .starts_with(".gvcf-"),
                "temporary output was not cleaned up"
            );
        }
        output
    }

    pub fn copy(&self, from: &Path, name: &str) {
        fs::copy(from, self.0.join(name)).unwrap();
    }

    pub fn write(&self, name: &str, text: &str) {
        fs::write(self.0.join(name), text).unwrap();
    }

    pub fn read(&self, name: &str) -> String {
        fs::read_to_string(self.0.join(name)).unwrap()
    }

    pub fn json(&self, name: &str) -> Value {
        serde_json::from_str(&self.read(name)).unwrap()
    }

    pub fn core() -> Self {
        let temp = Self::new();
        for name in [
            "reference.fa",
            "reference.fa.fai",
            "sample.g.vcf",
            "targets.tsv",
        ] {
            temp.copy(&fixture(name), name);
        }
        temp
    }

    pub fn audit(&self, out: &str, flags: &[&str], code: i32) -> Output {
        let mut args = vec![
            "--gvcf",
            "sample.g.vcf",
            "--reference",
            "reference.fa",
            "--out",
            out,
        ];
        args.extend_from_slice(flags);
        self.run(&args, code)
    }
}

impl Drop for Workspace {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

pub fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/core")
        .join(name)
}
