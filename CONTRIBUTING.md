# Contributing

Use GitHub issues for bugs and questions, and pull requests for changes. A useful
bug report includes the version, command, error and a small input that reproduces
it. Small synthetic examples are especially helpful.

Build with Rust 1.88 or later. Before sending a change, run:

```sh
cargo fmt --check
cargo clippy --locked --all-targets -- -D warnings
cargo test --locked
cargo build --locked
python3 tests/oracles/run.py target/debug/gvcf-audit
python3 examples/giab/test_overlap.py
```

The Python checks use only the standard library. See [tests/README.md](tests/README.md)
for optional MultiQC checks and for replaying a randomized failure. Full public
genomes are examples, not required development downloads.

Keep changes focused and explain the behavior that changes. For a bug fix, include
a small test with an independently calculated expectation. Changes to callability
rules need an explanation of the evidence being accepted or excluded, plus updates
to the policy documentation and affected outputs. Discuss those changes in an issue
first so their scientific implications are clear.

The maintainer reviews issues and pull requests in the repository. There is no
formal support schedule. Contributions are covered by the project's MIT license.
