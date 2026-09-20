#!/bin/bash
#SBATCH --job-name=dm_coupling
#SBATCH --account=jfe4_cr_default
#SBATCH --partition=basic
#SBATCH --array=1-48
#SBATCH --time=6:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --output=replication/slurm/logs/13_coupling_%A_%a.out
# Jaccard vs overlap coupling regimes (sim/07): 48 configs x 30 reps x 5 couplings.
set -euo pipefail
cd "${DM_ROOT:?set DM_ROOT to the project root}"
module load r/4.5.0
mkdir -p replication/slurm/logs
Rscript replication/sim/07_coupling_regimes.R
