#!/bin/bash
# =============================================================================
# replication/extended/slurm/06_regime_comparison.sbatch
#
# Head-to-head method comparison across four dynamic-community regimes
# (seasonality, churnswitch, birthdeath, regimeshift), one array task per
# config (1..72 = 4 regimes x N{50,100,200} x r{1.5,3,6} x intensity{low,high}).
# Each task runs 30 reps x all available methods for its config on the tracked
# (meta) partition. Config order is shuffled (set.seed(123)) so a partial run
# still samples every regime / N / r / intensity.
#
# Outputs: replication/extended/output/dynamic/dyn_cfg%02d.csv (one per config;
# resume-safe -- an existing file for a task is skipped).
#
# R_LIBS_USER is pinned to the user library that actually holds igraph>=2.0.0,
# dynsbm, multinet and clue. Do NOT fall back to $HOME/R/libs: that shadows the
# real user library and makes pkgload::load_all() fail the igraph version check.
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
export R_LIBS_USER=/storage/home/jfe4/R/x86_64-pc-linux-gnu-library/4.5
mkdir -p replication/extended/slurm/logs replication/extended/output/dynamic
Rscript replication/extended/06_regime_comparison.R
