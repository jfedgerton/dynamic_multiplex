#!/bin/bash
# =============================================================================
# replication/extended/slurm/06_regime_comparison.sbatch
#
# Head-to-head method comparison across the four dynamic-community regimes,
# one array task per config (1..72). Each task runs all reps x all available
# methods for its config on the tracked (meta) partition.
#
# RANDOMIZED CONDITION ORDER: the 72-row config grid
# (regime x N{50,100,200} x r{1.5,3,6} x intensity{low,high}) is shuffled once
# inside 06_regime_comparison.R with `set.seed(123); ord <- sample(72)` and the
# array index selects cfgs[ord[TASK], ]. Because the order is randomized, even
# a PARTIAL run (only the first handful of array indices) already samples every
# regime, every N, every separation r and both intensities -- no condition is
# systematically left for last.
#
# Dynamic SBM is ~1000x slower than the other methods and runs only on the
# first 50 reps (fast methods run 100); the wall clock below is generous.
#
# Optional packages for the full method set (install once into R_LIBS_USER):
#   Rscript -e '.libPaths(Sys.getenv("R_LIBS_USER"));
#               install.packages(c("dynsbm","multinet","clue"),
#                 repos="https://cloud.r-project.org",
#                 lib=Sys.getenv("R_LIBS_USER"))'
#   - dynsbm   : "Dynamic SBM" method (skipped if absent)
#   - multinet : "multinet GLouvain" method (recorded as NA if absent)
#   - clue     : optimal Hungarian label matching (a greedy fallback is used
#                if absent, but clue is preferred)
#
# Set --account / --partition to your allocation before submitting.
# =============================================================================
#SBATCH --job-name=repl_regime
#SBATCH --output=replication/extended/slurm/logs/regime_%A_%a.out
#SBATCH --array=1-72%16
#SBATCH --time=8:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --account=jfe4_cr_default
#SBATCH --partition=basic
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
module load r/4.5.0
export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/libs}"
mkdir -p replication/extended/slurm/logs replication/extended/output/comparison
Rscript replication/extended/06_regime_comparison.R
