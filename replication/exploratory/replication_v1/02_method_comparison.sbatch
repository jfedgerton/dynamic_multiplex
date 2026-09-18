#!/bin/bash
# =============================================================================
# replication/slurm/02_method_comparison.sbatch
#
# Head-to-head method comparison, one array task per config (1..22). Each task
# runs all 100 reps x all 11 methods for its config. Dynamic SBM is ~1000x
# slower than the other methods, so the wall clock is set generously (the
# heaviest T=40 configs take a few hours on 8 cores).
#
# Requires the dynsbm package for the Dynamic SBM column:
#     Rscript -e '.libPaths(Sys.getenv("R_LIBS_USER"));
#                 install.packages("dynsbm",
#                   repos="https://cloud.r-project.org",
#                   lib=Sys.getenv("R_LIBS_USER"))'
# If dynsbm is absent the script runs the other 10 methods and skips it.
#
# Set --account / --partition to your allocation before submitting.
# =============================================================================
#SBATCH --job-name=repl_cmp
#SBATCH --output=replication/slurm/logs/cmp_%A_%a.out
#SBATCH --array=1-22%12
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --account=jfe4_cr_default
#SBATCH --partition=basic
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
module load r/4.5.0
export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/libs}"
mkdir -p replication/slurm/logs replication/output/comparison
CMP_CORES=8 Rscript replication/02_method_comparison.R
