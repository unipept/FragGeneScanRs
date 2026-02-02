#!/bin/bash
set -euo pipefail

# Benchmark script for FragGeneScanRs using hyperfine
# Runs benchmarks on all example files with at least 10 iterations

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
BINARY="$PROJECT_ROOT/target/release/FragGeneScanRs"

# Check for hyperfine
if ! command -v hyperfine &> /dev/null; then
    echo "Error: hyperfine is not installed."
    echo "Install with: brew install hyperfine or similar"
    exit 1
fi

# Build release binary
echo "Building release binary..."
cargo build --release --manifest-path "$PROJECT_ROOT/Cargo.toml"

# Create temp directory for outputs
TEMP_DIR=$(mktemp -d)
trap 'rm -rf "$TEMP_DIR"' EXIT

echo ""
echo "Running benchmarks (minimum 10 runs, 3 warmup runs each)..."
echo "============================================================"

# Benchmark 1: Short reads (NC_000913-454.fna with 454_10 training)
echo ""
echo "Benchmark: Short reads (NC_000913-454.fna)"
hyperfine \
    --warmup 3 \
    --min-runs 20 \
    "$BINARY -s $PROJECT_ROOT/example/NC_000913-454.fna -t 454_10 -w 0 -o $TEMP_DIR/NC_000913-454"

# Benchmark 2: Complete genome (NC_000913.fna with complete training)
echo ""
echo "Benchmark: Complete genome (NC_000913.fna)"
hyperfine \
    --warmup 3 \
    --min-runs 20 \
    "$BINARY -s $PROJECT_ROOT/example/NC_000913.fna -t complete -w 1 -o $TEMP_DIR/NC_000913"

# Benchmark 3: Long reads (contigs.fna with complete training)
echo ""
echo "Benchmark: Long reads (contigs.fna)"
hyperfine \
    --warmup 3 \
    --min-runs 10 \
    "$BINARY -s $PROJECT_ROOT/example/contigs.fna -t complete -w 1 -o $TEMP_DIR/contigs"

echo ""
echo "Benchmarks complete!"
