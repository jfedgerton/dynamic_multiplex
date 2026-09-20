#!/bin/bash
#SBATCH --job-name=dm_mech
#SBATCH --account=jfe4_cr_default
#SBATCH --partition=basic
#SBATCH --array=1-50
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=replication/slurm/logs/03_mechanism_%A_%a.out
# Mechanism tests (sim/03): turnover / coreperiphery / longT, 50 configs x 10 reps x 3 methods.
set -euo pipefail
cd "${DM_ROOT:?set DM_ROOT to the project root}"
module load r/4.5.0
mkdir -p replication/slurm/logs
Rscript replication/sim/03_mechanism_tests.R
