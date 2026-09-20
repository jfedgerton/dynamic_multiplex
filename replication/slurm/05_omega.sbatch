#!/bin/bash
#SBATCH --job-name=dm_omega
#SBATCH --account=jfe4_cr_default
#SBATCH --partition=basic
#SBATCH --array=1-72
#SBATCH --time=12:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --output=replication/slurm/logs/05_omega_%A_%a.out
# Multislice omega sweep on the regime grid (sim/05): 72 configs x 10 reps x 11 variants.
set -euo pipefail
cd "${DM_ROOT:?set DM_ROOT to the project root}"
module load r/4.5.0
mkdir -p replication/slurm/logs
Rscript replication/sim/05_omega_sweep.R
