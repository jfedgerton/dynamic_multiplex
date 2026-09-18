#!/usr/bin/env bash
# =============================================================================
# replication/run_all.sh  -  MASTER REPLICATION SCRIPT
#
# One entry point for the full replication of the paper's computational
# results. It (1) installs the packages, then (2) submits the three analysis
# pipelines to SLURM, then (3) points you at the aggregation/figure steps.
#
# The three pipelines are independent and may be submitted in any order or in
# parallel; there is no cross-pipeline dependency. Every script uses fixed
# seeds, so a rerun reproduces the numbers in the paper exactly.
#
#   Pipeline (a)  01_ci_coverage.R        CI-coverage validation + gate
#                 slurm/01_ci_coverage.sbatch      (binary + valued + misspec)
#   Pipeline (b)  02_method_comparison.R  head-to-head, 11 methods + Dynamic SBM
#                 slurm/02_method_comparison.sbatch
#   Pipeline (c)  03_coldwar_empirical.R  alliance / IGO / trade, Cold-War break
#                 slurm/03_coldwar_empirical.sbatch
#
# Usage:
#   bash replication/run_all.sh deps          # install R + Python packages only
#   bash replication/run_all.sh submit        # submit all three pipelines
#   bash replication/run_all.sh submit ci      # submit pipeline (a) only
#   bash replication/run_all.sh submit compare # submit pipeline (b) only
#   bash replication/run_all.sh submit empirical  # submit pipeline (c) only
#   bash replication/run_all.sh                # deps, then submit all
#
# EDIT before running: set --account / --partition in slurm/*.sbatch to your
# allocation, and R_LIBS_USER below to a writable library path.
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")/.."                       # repo root
export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/libs}"
mkdir -p "$R_LIBS_USER" replication/output replication/slurm/logs

ACTION="${1:-full}"
WHICH="${2:-all}"

install_deps() {
  echo "--- Installing the dynamicmultiplex R package (local source) ---"
  Rscript -e '.libPaths(Sys.getenv("R_LIBS_USER"));
    if (!requireNamespace("remotes", quietly=TRUE))
      install.packages("remotes", repos="https://cloud.r-project.org",
                       lib=Sys.getenv("R_LIBS_USER"));
    remotes::install_local("r_code", force=TRUE, upgrade="never",
                           lib=Sys.getenv("R_LIBS_USER"))'
  echo "--- Installing analysis dependencies (CRAN) ---"
  Rscript -e '.libPaths(Sys.getenv("R_LIBS_USER"));
    pkgs <- c("igraph","ggplot2","dplyr","tidyr","parallel",
              "dynsbm","peacesciencer","ggalluvial");
    miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)];
    if (length(miss)) install.packages(miss,
        repos="https://cloud.r-project.org", lib=Sys.getenv("R_LIBS_USER"));
    cat("R dependencies OK\n")'
  echo "--- Installing the Python package (for cross-language checks) ---"
  pip install --user -e python_code || \
    echo "(python install skipped; not required for the R pipelines)"
}

submit_ci()       { echo "submitting pipeline (a) CI coverage";
                    sbatch replication/slurm/01_ci_coverage.sbatch; }
submit_compare()  { echo "submitting pipeline (b) method comparison";
                    sbatch replication/slurm/02_method_comparison.sbatch; }
submit_empirical(){ echo "submitting pipeline (c) Cold-War empirical";
                    sbatch replication/slurm/03_coldwar_empirical.sbatch; }

case "$ACTION" in
  deps)   install_deps ;;
  submit)
    case "$WHICH" in
      ci)        submit_ci ;;
      compare)   submit_compare ;;
      empirical) submit_empirical ;;
      all)       submit_ci; submit_compare; submit_empirical ;;
      *) echo "unknown pipeline '$WHICH'"; exit 1 ;;
    esac ;;
  full)   install_deps; submit_ci; submit_compare; submit_empirical ;;
  *) echo "usage: run_all.sh [deps|submit [ci|compare|empirical|all]]"; exit 1 ;;
esac

echo ""
echo "Submitted. When jobs finish, aggregate with the analysis notes in"
echo "replication/README.md (section 'Aggregating the results')."
