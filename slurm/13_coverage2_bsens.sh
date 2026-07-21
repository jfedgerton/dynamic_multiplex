#!/bin/bash
#SBATCH --job-name=dm_cov2_bsens
#SBATCH --account=open
#SBATCH --partition=basic
#SBATCH --output=slurm/logs/13_cov2_bsens_%A_%a.out
#SBATCH --error=slurm/logs/13_cov2_bsens_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-6
set -euo pipefail
cd /storage/group/LiberalArts/default/jfe4_collab/dynamic_multiplex
module load r/4.5.0
export R_LIBS_USER=/storage/home/jfe4/R/libs.old.20260602
mkdir -p slurm/logs manuscript/output
COV_MODE=bsens COV_CORES=8 Rscript manuscript/13_coassign_ci_coverage.R
