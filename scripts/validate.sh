#!/bin/bash
set -euo pipefail

# Validation script for FragGeneScanRs
# Usage:
#   ./scripts/validate.sh --baseline   Generate baseline output files
#   ./scripts/validate.sh --check      Compare current output against baseline (default)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
BASELINE_DIR="$PROJECT_ROOT/baseline"
BINARY="$PROJECT_ROOT/target/release/FragGeneScanRs"

# Example files and their configurations
# Format: "input_file:training_file:whole_genome_flag:output_name"
EXAMPLES=(
    "example/NC_000913-454.fna:454_10:0:NC_000913-454"
    "example/NC_000913.fna:complete:1:NC_000913"
    "example/contigs.fna:complete:1:contigs"
)

usage() {
    echo "Usage: $0 [--baseline|--check]"
    echo "  --baseline  Generate baseline output files"
    echo "  --check     Compare current output against baseline (default)"
    exit 1
}

build_release() {
    echo "Building release binary..."
    cargo build --release --manifest-path "$PROJECT_ROOT/Cargo.toml"
}

run_example() {
    local input="$1"
    local train="$2"
    local whole="$3"
    local output_prefix="$4"

    "$BINARY" \
        -s "$PROJECT_ROOT/$input" \
        -t "$train" \
        -w "$whole" \
        -o "$output_prefix"
}

generate_baseline() {
    echo "Generating baseline outputs..."
    mkdir -p "$BASELINE_DIR"

    for example in "${EXAMPLES[@]}"; do
        IFS=':' read -r input train whole name <<< "$example"
        echo "  Processing $name..."
        run_example "$input" "$train" "$whole" "$BASELINE_DIR/$name"
    done

    echo "Baseline generated in $BASELINE_DIR"
}

check_against_baseline() {
    if [[ ! -d "$BASELINE_DIR" ]]; then
        echo "Error: Baseline directory not found. Run with --baseline first."
        exit 1
    fi

    local temp_dir
    temp_dir=$(mktemp -d)
    trap 'rm -rf "$temp_dir"' EXIT

    echo "Running current version and comparing against baseline..."
    local failed=0

    for example in "${EXAMPLES[@]}"; do
        IFS=':' read -r input train whole name <<< "$example"
        echo "  Processing $name..."
        run_example "$input" "$train" "$whole" "$temp_dir/$name"

        for ext in out faa ffn; do
            local baseline_file="$BASELINE_DIR/$name.$ext"
            local current_file="$temp_dir/$name.$ext"

            if [[ ! -f "$baseline_file" ]]; then
                echo "    Warning: Baseline file $baseline_file not found"
                continue
            fi

            if diff -q "$baseline_file" "$current_file" > /dev/null 2>&1; then
                echo "    ✓ $name.$ext matches"
            else
                echo "    ✗ $name.$ext DIFFERS"
                failed=1
            fi
        done
    done

    if [[ $failed -eq 0 ]]; then
        echo "All outputs match baseline!"
        exit 0
    else
        echo "Some outputs differ from baseline!"
        exit 1
    fi
}

# Parse arguments
MODE="check"
if [[ $# -gt 0 ]]; then
    case "$1" in
        --baseline)
            MODE="baseline"
            ;;
        --check)
            MODE="check"
            ;;
        *)
            usage
            ;;
    esac
fi

# Main
build_release

case "$MODE" in
    baseline)
        generate_baseline
        ;;
    check)
        check_against_baseline
        ;;
esac
