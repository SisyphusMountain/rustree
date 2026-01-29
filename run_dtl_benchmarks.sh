#!/bin/bash

# Run comprehensive DTL benchmarks
# This script runs all major benchmark tests and saves output

echo "========================================"
echo "DTL Simulator Comprehensive Benchmarks"
echo "========================================"
echo ""

# Quick test (always runs)
echo "==> Running quick benchmark..."
cargo test --test dtl_benchmark benchmark_dtl_quick_test -- --nocapture

# Varying rates (fast, ~1 second)
echo ""
echo "==> Running varying rates benchmark..."
cargo test --test dtl_benchmark benchmark_dtl_varying_rates -- --ignored --nocapture

# Loss impact (fast, ~0.2 seconds)
echo ""
echo "==> Running loss impact benchmark..."
cargo test --test dtl_benchmark benchmark_dtl_loss_impact -- --ignored --nocapture

echo ""
echo "========================================"
echo "Benchmarks complete!"
echo "See DTL_BENCHMARK_RESULTS.md for full analysis"
echo "========================================"
