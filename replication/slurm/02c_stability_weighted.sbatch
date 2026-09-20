#!/bin/bash
#SBATCH --job-name=dm_stab_wt
#SBATCH --account=jfe4_cr_default
#SBATCH --partition=basic
#SBATCH --array=1-72%48
#SBATCH --time=8:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --output=replication/slurm/logs/02c_stab_weighted_%A_%a.out
# Stability vs accuracy, weighted-network arm (sim/02, 72 cells).
set -euo pipefail
cd "${DM_ROOT:?set DM_ROOT to the project root}"
module load r/4.5.0
mkdir -p replication/slurm/logs
STAB_ARM=weighted COV_CORES=8 Rscript replication/sim/02_stability.R
