#!/bin/bash
#SBATCH --job-name=dm_emp_score
#SBATCH --account=jfe4_cr_default
#SBATCH --partition=basic
#SBATCH --time=2:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --output=replication/slurm/logs/07_emp_score_%j.out
# Score recovery of the Braumoeller-coded orders on all four networks
# (empirical/07) -> output/empirical/order_recovery.csv
set -euo pipefail
cd "${DM_ROOT:?set DM_ROOT to the project root}"
module load r/4.5.0
mkdir -p replication/slurm/logs
Rscript replication/empirical/07_score_orders.R
