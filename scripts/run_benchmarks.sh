#!/usr/bin/env bash
#
# Master launcher for dynamic_multiplex benchmarks.
#
# Usage:
#   bash scripts/run_benchmarks.sh          # run both R and Python
#   bash scripts/run_benchmarks.sh python   # Python only
#   bash scripts/run_benchmarks.sh r        # R only
#
set -euo pipefail
cd "$(dirname "$0")/.."

TARGET="${1:-all}"

# ---- Python benchmark ----
if [[ "$TARGET" == "all" || "$TARGET" == "python" ]]; then
    echo "=== Installing Python package ==="
    pip install -e ./python_code[louvain,dev] --quiet

    # scikit-learn is optional but preferred for NMI/ARI
    pip install scikit-learn --quiet 2>/dev/null || true

    echo "=== Running Python tests ==="
    python -m pytest python_code/tests -q

    echo ""
    echo "=== Running Python benchmark ==="
    python scripts/benchmark_python.py
fi

# ---- R benchmark ----
if [[ "$TARGET" == "all" || "$TARGET" == "r" ]]; then
    echo ""
    echo "=== Installing R package ==="
    Rscript -e 'if (!requireNamespace("remotes", quietly=TRUE)) install.packages("remotes", repos="https://cloud.r-project.org"); remotes::install_local("r_code", force=TRUE, quiet=TRUE)'

    echo "=== Running R benchmark ==="
    Rscript scripts/benchmark_r.R
fi

echo ""
echo "=== Done ==="
echo "Results:"
ls -la scripts/benchmark_*_results.csv 2>/dev/null || echo "  (no result files yet)"
