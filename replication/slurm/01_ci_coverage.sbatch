#!/bin/bash
# =============================================================================
# replication/slurm/01_ci_coverage.sbatch
#
# CI-coverage validation, all three arms in one array. The global array index
# maps to (arm, local task) inside 01_ci_coverage.R:
#     1 .. 3564         -> binary  (planted-partition SBM)
#     3565 .. 4284      -> valued  (weighted networks)
#     4285 .. 4716      -> misspec (degree-corrected SBM stress test)
# Total 4716 tasks. Each task = 250 sims x 100 bootstrap x 100 fresh-draw.
#
# NOTE on the queue cap: some Slurm associations limit queued array elements
# (AssocMaxSubmitJobLimit, often ~4000). If the single 4716-element submit is
# rejected, submit the three arms separately instead:
#     sbatch --export=ALL,COV_ARM=binary  --array=1-3564%64 01_ci_coverage.sbatch
#     sbatch --export=ALL,COV_ARM=valued  --array=1-720%64  01_ci_coverage.sbatch
#     sbatch --export=ALL,COV_ARM=misspec --array=1-432%64  01_ci_coverage.sbatch
# (COV_ARM makes the script use the local index directly.)
#
# Set --account / --partition to your allocation before submitting.
# =============================================================================
#SBATCH --job-name=repl_ci
#SBATCH --output=replication/slurm/logs/ci_%A_%a.out
#SBATCH --array=1-4716%64
#SBATCH --time=8:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --account=jfe4_cr_default
#SBATCH --partition=basic
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
module load r/4.5.0
export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/libs}"
mkdir -p replication/slurm/logs replication/output
COV_CORES=8 Rscript replication/01_ci_coverage.R
