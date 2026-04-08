#!/bin/bash
#SBATCH --job-name=dynmux_py_bench
#SBATCH --account=open
#SBATCH --partition=open
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=slurm/logs/python_benchmarks_%j.out
#SBATCH --error=slurm/logs/python_benchmarks_%j.err

# ============================================================
# Python benchmarks for dynamic_multiplex
# ============================================================
set -euo pipefail

PROJECT_DIR="/scratch/jfe4/dynamic_multiplex"
cd "$PROJECT_DIR"
mkdir -p slurm/logs

# Load modules
module purge
module load python/3.13.0

# Activate virtual environment
source venv/bin/activate

# Install scikit-learn (needed for NMI/ARI)
pip install scikit-learn --quiet 2>/dev/null || true

# Run Python tests first
echo "=== Running Python tests ==="
python -m pytest python_code/tests -q

# Run Python benchmark
echo ""
echo "=== Running Python benchmark ==="
python scripts/benchmark_python.py

echo ""
echo "=== Python benchmarks complete ==="
ls -la scripts/benchmark_*_results.csv 2>/dev/null || echo "  (no result files yet)"
