#!/bin/bash
#SBATCH --job-name=dm_cov_misspec
#SBATCH --account=jfe4_cr_default
#SBATCH --partition=basic
#SBATCH --array=1-432%64
#SBATCH --time=10:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --output=replication/slurm/logs/04_cov_misspec_%A_%a.out
# Degree-corrected (DC-SBM) coverage arm: 432 tasks, 250 sims each (sim/04).
set -euo pipefail
cd "${DM_ROOT:?set DM_ROOT to the project root}"
module load r/4.5.0
mkdir -p replication/slurm/logs
COV_CORES=8 Rscript replication/sim/04_coverage_misspec.R
