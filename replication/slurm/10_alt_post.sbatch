#!/bin/bash
#SBATCH --job-name=dm_alt_post
#SBATCH --account=jfe4_cr_default
#SBATCH --partition=basic
#SBATCH --time=4:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=2
#SBATCH --output=replication/slurm/logs/10_alt_post_%j.out
# EXPLORATORY: compares interval alternatives (post/15). Reads sim/02 and
# sim/05 output; writes output/alternatives/*.csv and prints the comparison.
set -euo pipefail
cd "${DM_ROOT:?set DM_ROOT to the project root}"
module load r/4.5.0
mkdir -p replication/slurm/logs
Rscript replication/post/15_alternatives.R
