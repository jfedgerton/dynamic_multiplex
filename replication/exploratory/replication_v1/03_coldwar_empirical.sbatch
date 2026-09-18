#!/bin/bash
# =============================================================================
# replication/slurm/03_coldwar_empirical.sbatch
#
# Cold-War empirical application (alliance, IGO, trade). A single job that
# fits all three networks over 1965-2000, computes community persistence, the
# Cold-War break diagnostic, cross-specification agreement, and the figures.
#
# Requires the peacesciencer package (downloads / caches COW-ATOP-IGO data):
#     Rscript -e '.libPaths(Sys.getenv("R_LIBS_USER"));
#                 install.packages("peacesciencer",
#                   repos="https://cloud.r-project.org",
#                   lib=Sys.getenv("R_LIBS_USER"))'
# Run once on a node with network access so peacesciencer can fetch/cache its
# data; subsequent runs use the cache.
#
# Set --account / --partition to your allocation before submitting.
# =============================================================================
#SBATCH --job-name=repl_emp
#SBATCH --output=replication/slurm/logs/emp_%j.out
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --account=jfe4_cr_default
#SBATCH --partition=basic
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
module load r/4.5.0
export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/libs}"
mkdir -p replication/slurm/logs replication/output/empirical
Rscript replication/03_coldwar_empirical.R
