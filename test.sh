#!/bin/sh
set -euox pipefail

# Perform the Rust tests and benchmarks.
cargo test -r
cargo bench

# Check whether a version bump satisfies the semantic versioning checks.
cargo semver-checks

# Test the Python bindings.
cd bindings/python
source venv/bin/activate
pip install .
python tests/verification.py
