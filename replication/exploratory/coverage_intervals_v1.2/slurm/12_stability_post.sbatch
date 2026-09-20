#!/bin/bash
#SBATCH --job-name=dm_stab_post
#SBATCH --account=jfe4_cr_default
#SBATCH --partition=basic
#SBATCH --time=1:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --output=replication/slurm/logs/12_stab_post_%j.out
# EXPLORATORY: applies the pre-registered decision rule (post/16).
set -euo pipefail
cd "${DM_ROOT:?set DM_ROOT to the project root}"
module load r/4.5.0
mkdir -p replication/slurm/logs
Rscript replication/post/16_stability.R
