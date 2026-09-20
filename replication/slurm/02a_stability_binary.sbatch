#!/bin/bash
#SBATCH --job-name=dm_stab_bin
#SBATCH --account=jfe4_cr_default
#SBATCH --partition=basic
#SBATCH --array=1-594%64
#SBATCH --time=8:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --output=replication/slurm/logs/02a_stab_binary_%A_%a.out
# Stability vs accuracy, binary planted-partition arm (sim/02, 594 cells x 50 sims x B=100): the calibration arm.
set -euo pipefail
cd "${DM_ROOT:?set DM_ROOT to the project root}"
module load r/4.5.0
mkdir -p replication/slurm/logs
STAB_ARM=binary COV_CORES=8 Rscript replication/sim/02_stability.R
