#!/bin/bash
#SBATCH --job-name=dm_post
#SBATCH --account=jfe4_cr_default
#SBATCH --partition=basic
#SBATCH --time=4:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=2
#SBATCH --output=replication/slurm/logs/08_postprocess_%j.out
# Every table and figure in the manuscript, plus the calibrated-interval
# lookup table, from the outputs of steps 01-07. Also runnable by hand:
#   bash replication/run_all.sh post
set -euo pipefail
cd "${DM_ROOT:?set DM_ROOT to the project root}"
module load r/4.5.0
mkdir -p replication/slurm/logs
bash replication/run_all.sh post
