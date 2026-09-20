#!/bin/bash
#SBATCH --job-name=dm_emp_fit
#SBATCH --account=jfe4_cr_default
#SBATCH --partition=basic
#SBATCH --array=1-4
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --output=replication/slurm/logs/08_emp_fit_%A_%a.out
# Fit every method on one empirical network (empirical/08); array index -> atop, dca, igo, trade. multinet on ATOP/IGO needs 64G.
set -euo pipefail
cd "${DM_ROOT:?set DM_ROOT to the project root}"
module load r/4.5.0
mkdir -p replication/slurm/logs
NETS=(atop dca igo trade); Rscript replication/empirical/08_fit_networks.R ${NETS[$((SLURM_ARRAY_TASK_ID-1))]}
