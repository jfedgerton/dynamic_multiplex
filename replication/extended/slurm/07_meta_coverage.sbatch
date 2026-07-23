#!/bin/bash
# =============================================================================
# replication/extended/slurm/07_meta_coverage.sbatch
#
# Meta co-membership interval coverage re-validation, one array task per
# (deployment x config) cell (1..150). Each task runs all sims for its cell:
# simulate -> bootstrap_multilayer + co_assignment_ci (both on the tracked META
# communities) -> coverage vs fresh-draw ground-truth co-assignment propensity.
#
# RANDOMIZED CONDITION ORDER: the 150-row task grid
# (10 deployments {jaccard,overlap,weighted_jaccard,weighted_overlap,identity}
#  x {leiden,louvain}, each over 15 configs = 9 balanced + 3 skewed + 3 misspec)
# is shuffled once inside 07_meta_coverage.R with
# `set.seed(123); ord <- sample(150)` and the array index selects tasks[ord[TASK],].
# Because the order is randomized, even a PARTIAL run (only the first handful of
# array indices) already samples every deployment, every block (balanced /
# skewed / misspec), every N and every separation -- nothing is left for last.
#
# THE GATE is width_P_mean < 0.05: coverage is trusted (nominal ~95%) inside the
# gate and refused outside it. Nominal cov_P INSIDE the gate is the ship
# criterion, and the SKEWED and MISSPEC (DC-SBM, severe degree heterogeneity)
# blocks are the ship-or-die check that the META communities keep nominal
# coverage under community-size skew and degree misspecification -- not just on
# the balanced planted-partition grid.
#
# Set --account / --partition to your allocation before submitting.
# =============================================================================
#SBATCH --job-name=repl_metacov
#SBATCH --output=replication/extended/slurm/logs/metacov_%A_%a.out
#SBATCH --array=1-150%16
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --account=jfe4_cr_default
#SBATCH --partition=basic
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
module load r/4.5.0
export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/libs}"
mkdir -p replication/extended/slurm/logs replication/extended/output/coverage
Rscript replication/extended/07_meta_coverage.R
