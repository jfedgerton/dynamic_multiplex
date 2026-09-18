#!/bin/bash
#SBATCH --job-name=dm_cov_valued
#SBATCH --account=open
#SBATCH --partition=basic
#SBATCH --output=slurm/logs/12_cov_valued_%A_%a.out
#SBATCH --error=slurm/logs/12_cov_valued_%A_%a.err
#SBATCH --time=10:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-160%32
set -euo pipefail
cd "${DM_ROOT:?set DM_ROOT to the project root}"
module load r/4.5.0
# R package library: left unset so R resolves its own default user library
# (~/R/<platform>-library/<version>). Export R_LIBS_USER in your environment
# before submitting if your packages live elsewhere.
mkdir -p slurm/logs manuscript/output
COV_MODE=valued COV_CORES=8 Rscript manuscript/12_bootstrap_coverage.R
