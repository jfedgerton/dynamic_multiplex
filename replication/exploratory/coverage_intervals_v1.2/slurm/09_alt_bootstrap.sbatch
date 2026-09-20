#!/bin/bash
#SBATCH --job-name=dm_alt_boot
#SBATCH --account=jfe4_cr_default
#SBATCH --partition=basic
#SBATCH --array=1-594%64
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --output=replication/slurm/logs/09_alt_boot_%A_%a.out
# EXPLORATORY: three bootstrap variants + per-pair feature sample on the
# 594-cell grid (sim/05). One cell per task, 50 sims. Self-skips if done.
set -euo pipefail
cd "${DM_ROOT:?set DM_ROOT to the project root}"
module load r/4.5.0
mkdir -p replication/slurm/logs
COV_CORES=8 Rscript replication/sim/05_alt_bootstrap.R
