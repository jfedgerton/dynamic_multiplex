#!/bin/bash
#SBATCH --job-name=dm_cov_grid
#SBATCH --account=jfe4_cr_default
#SBATCH --partition=basic
#SBATCH --array=1-990%64
#SBATCH --time=10:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --output=replication/slurm/logs/02_cov_grid_%A_%a.out
# Binary coverage grid: 3564 tasks (594 cells x 6 specs), 250 sims each
# (sim/02). The array cap is 990 per submission, so run_all.sh submits this
# file four times with COV_OFFSET = 0, 990, 1980, 2970 (the last with
# --array=1-594). Completed tasks self-skip.
set -euo pipefail
cd "${DM_ROOT:?set DM_ROOT to the project root}"
module load r/4.5.0
mkdir -p replication/slurm/logs
COV_CORES=8 COV_OFFSET=${COV_OFFSET:-0} Rscript replication/sim/02_coverage_grid.R
