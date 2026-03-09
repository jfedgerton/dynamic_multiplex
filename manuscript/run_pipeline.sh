#!/usr/bin/env bash
#
# Manuscript pipeline launcher for dynamic_multiplex.
#
# Runs synthetic experiments and empirical application scripts,
# producing tables and figures in manuscript/output/.
#
# Usage:
#   bash manuscript/run_pipeline.sh             # run everything
#   bash manuscript/run_pipeline.sh synthetic    # synthetic experiments only
#   bash manuscript/run_pipeline.sh empirical    # empirical application only
#
set -euo pipefail
cd "$(dirname "$0")/.."

TARGET="${1:-all}"

echo "=== dynamic_multiplex manuscript pipeline ==="

# ---- Install R package ----
echo ""
echo "--- Installing R package ---"
Rscript -e '
  if (!requireNamespace("remotes", quietly=TRUE))
    install.packages("remotes", repos="https://cloud.r-project.org")
  remotes::install_local("r_code", force=TRUE, quiet=TRUE)
'

# ---- Install manuscript dependencies ----
echo ""
echo "--- Checking R dependencies ---"
Rscript -e '
  pkgs <- c("ggplot2", "dplyr", "tidyr")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
  if (length(missing) > 0) {
    cat("Installing:", paste(missing, collapse=", "), "\n")
    install.packages(missing, repos="https://cloud.r-project.org")
  }
  cat("Core dependencies OK.\n")

  # Optional
  opt <- c("ggalluvial", "RColorBrewer", "peacesciencer")
  missing_opt <- opt[!vapply(opt, requireNamespace, logical(1), quietly=TRUE)]
  if (length(missing_opt) > 0) {
    cat("Optional packages not installed:", paste(missing_opt, collapse=", "), "\n")
    cat("Install with: install.packages(c(", paste0("\"", missing_opt, "\"", collapse=", "), "))\n")
  }
'

# ---- Synthetic experiments ----
if [[ "$TARGET" == "all" || "$TARGET" == "synthetic" ]]; then
    echo ""
    echo "--- Running synthetic experiments ---"
    Rscript manuscript/01_synthetic_experiments.R
fi

# ---- Empirical application ----
if [[ "$TARGET" == "all" || "$TARGET" == "empirical" ]]; then
    echo ""
    echo "--- Running empirical application ---"
    Rscript manuscript/02_empirical_application.R
fi

echo ""
echo "=== Pipeline complete ==="
echo "Outputs:"
ls -la manuscript/output/ 2>/dev/null || echo "  (no outputs yet)"
