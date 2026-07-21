#!/bin/bash
#SBATCH --job-name=dm_cov3grid
#SBATCH --account=open
#SBATCH --partition=basic
#SBATCH --output=slurm/logs/16_cov3grid_%A_%a.out
#SBATCH --error=slurm/logs/16_cov3grid_%A_%a.err
#SBATCH --time=10:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-990%64
# Submit three times with COV_OFFSET=0, 990, 1980 to cover tasks 1..2970:
#   COV_OFFSET=0    sbatch --export=ALL,COV_OFFSET=0    slurm/13_coverage2_main.sh
#   COV_OFFSET=990  sbatch --export=ALL,COV_OFFSET=990  slurm/13_coverage2_main.sh
#   COV_OFFSET=1980 sbatch --export=ALL,COV_OFFSET=1980 slurm/13_coverage2_main.sh
set -euo pipefail
cd /storage/group/LiberalArts/default/jfe4_collab/dynamic_multiplex
module load r/4.5.0
export R_LIBS_USER=/storage/home/jfe4/R/libs.old.20260602
mkdir -p slurm/logs manuscript/output
COV_MODE=main COV_CORES=8 COV_OFFSET=${COV_OFFSET:-0} Rscript manuscript/16_coverage3_grid.R
