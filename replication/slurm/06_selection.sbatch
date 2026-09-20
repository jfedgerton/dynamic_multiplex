#!/bin/bash
#SBATCH --job-name=dm_selection
#SBATCH --account=jfe4_cr_default
#SBATCH --partition=basic
#SBATCH --array=1-72
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --output=replication/slurm/logs/06_selection_%A_%a.out
# Stability as a selection rule (sim/06): 72 configs x 10 reps, DynMux and multislice each bootstrapped (B=50).
set -euo pipefail
cd "${DM_ROOT:?set DM_ROOT to the project root}"
module load r/4.5.0
mkdir -p replication/slurm/logs
Rscript replication/sim/06_selection_rule.R
