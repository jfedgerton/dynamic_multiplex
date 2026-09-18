#!/bin/bash
#SBATCH --job-name=dm_cov_valued
#SBATCH --account=jfe4_cr_default
#SBATCH --partition=basic
#SBATCH --array=1-720%64
#SBATCH --time=10:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --output=replication/slurm/logs/03_cov_valued_%A_%a.out
# Weighted-network coverage arm: 720 tasks, 250 sims each (sim/03).
set -euo pipefail
cd "${DM_ROOT:?set DM_ROOT to the project root}"
module load r/4.5.0
mkdir -p replication/slurm/logs
COV_CORES=8 Rscript replication/sim/03_coverage_valued.R
