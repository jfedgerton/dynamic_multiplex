#!/bin/bash
#SBATCH --job-name=dm_emp_stab
#SBATCH --account=jfe4_cr_default
#SBATCH --partition=basic
#SBATCH --array=1-4
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --output=replication/slurm/logs/10_emp_stab_%A_%a.out
# Bootstrap stability of DynMux and multislice on one empirical network (empirical/10, B=100).
set -euo pipefail
cd "${DM_ROOT:?set DM_ROOT to the project root}"
module load r/4.5.0
mkdir -p replication/slurm/logs
NETS=(atop dca igo trade); Rscript replication/empirical/10_stability_networks.R ${NETS[$((SLURM_ARRAY_TASK_ID-1))]}
