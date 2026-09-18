#!/bin/bash
#SBATCH --job-name=dm_emp_fit
#SBATCH --account=jfe4_cr_default
#SBATCH --partition=basic
#SBATCH --array=1-4
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --output=replication/slurm/logs/06_emp_fit_%A_%a.out
# Fit the 12 methods on each empirical network (empirical/06). One array
# task per network; multinet on ATOP/IGO is the memory hog (64G).
set -euo pipefail
cd "${DM_ROOT:?set DM_ROOT to the project root}"
module load r/4.5.0
mkdir -p replication/slurm/logs
NETS=(atop dca igo trade)
NET=${NETS[$((SLURM_ARRAY_TASK_ID - 1))]}
echo "network: $NET"
Rscript replication/empirical/06_fit_networks.R "$NET"
