#!/bin/bash
#SBATCH --job-name=dm_coupling_post
#SBATCH --account=jfe4_cr_default
#SBATCH --partition=basic
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=replication/slurm/logs/14_coupling_post_%j.out
# Tables and figure for the coupling appendix (post/17).
set -euo pipefail
cd "${DM_ROOT:?set DM_ROOT to the project root}"
module load r/4.5.0
mkdir -p replication/slurm/logs
Rscript replication/post/17_coupling.R
