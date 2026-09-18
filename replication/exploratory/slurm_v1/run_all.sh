#!/bin/bash
# ============================================================
# Master launcher: submit all publication workflow jobs to SLURM
#
# Usage:
#   bash slurm/run_all.sh             # submit all jobs
#   bash slurm/run_all.sh python      # Python benchmarks only
#   bash slurm/run_all.sh synthetic   # synthetic experiments only
#   bash slurm/run_all.sh empirical   # empirical application only
#   bash slurm/run_all.sh manuscript  # synthetic + empirical (no Python)
#
# Jobs are submitted with dependencies so that empirical waits
# for synthetic to finish (they share R package installation).
# ============================================================
set -euo pipefail

PROJECT_DIR="/scratch/jfe4/dynamic_multiplex"
cd "$PROJECT_DIR"
mkdir -p slurm/logs

TARGET="${1:-all}"

echo "=== dynamic_multiplex SLURM pipeline ==="
echo "Target: $TARGET"
echo ""

# Python benchmarks
if [[ "$TARGET" == "all" || "$TARGET" == "python" ]]; then
    PY_JOB=$(sbatch --parsable slurm/01_python_benchmarks.sh)
    echo "Submitted Python benchmarks: job $PY_JOB"
fi

# Synthetic experiments
if [[ "$TARGET" == "all" || "$TARGET" == "synthetic" || "$TARGET" == "manuscript" ]]; then
    SYN_JOB=$(sbatch --parsable slurm/02_synthetic_experiments.sh)
    echo "Submitted synthetic experiments: job $SYN_JOB"
fi

# Empirical application (depends on synthetic finishing to avoid R lib conflicts)
if [[ "$TARGET" == "all" || "$TARGET" == "empirical" || "$TARGET" == "manuscript" ]]; then
    if [[ -n "${SYN_JOB:-}" ]]; then
        EMP_JOB=$(sbatch --parsable --dependency=afterok:$SYN_JOB slurm/03_empirical_application.sh)
        echo "Submitted empirical application: job $EMP_JOB (after $SYN_JOB)"
    else
        EMP_JOB=$(sbatch --parsable slurm/03_empirical_application.sh)
        echo "Submitted empirical application: job $EMP_JOB"
    fi
fi

echo ""
echo "Monitor with: squeue -u jfe4"
echo "Outputs will go to: manuscript/output/ and scripts/"
