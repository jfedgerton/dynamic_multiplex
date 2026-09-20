#!/bin/bash
#SBATCH --job-name=dm_stab
#SBATCH --account=jfe4_cr_default
#SBATCH --partition=basic
#SBATCH --array=1-594%64
#SBATCH --time=8:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --output=replication/slurm/logs/11_stab_%A_%a.out
# EXPLORATORY: does bootstrap stability predict accuracy? (sim/06), 594 cells x 50 sims.
set -euo pipefail
cd "${DM_ROOT:?set DM_ROOT to the project root}"
module load r/4.5.0
mkdir -p replication/slurm/logs
COV_CORES=8 Rscript replication/sim/06_stability.R
