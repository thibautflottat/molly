#!/bin/bash
set -x
set -e

trajectory=$1

# Build all examples to prevent running stale executables.
cargo build --release --examples

# Verify that the output is equal between chemfiles and ours.
time target/release/examples/compare $trajectory

hyperfine --shell=none --warmup 10 "target/release/examples/cfreader $trajectory" "target/release/examples/reader $trajectory"
