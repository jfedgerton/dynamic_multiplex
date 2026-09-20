#!/bin/bash
#SBATCH --job-name=dm_regime
#SBATCH --account=jfe4_cr_default
#SBATCH --partition=basic
#SBATCH --array=1-72
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=replication/slurm/logs/01_regime_%A_%a.out
# Regime comparison: 72 configs x 30 reps x every method (sim/01). Finished tasks self-skip on resubmission.
set -euo pipefail
cd "${DM_ROOT:?set DM_ROOT to the project root}"
module load r/4.5.0
mkdir -p replication/slurm/logs
Rscript replication/sim/01_regime_comparison.R
