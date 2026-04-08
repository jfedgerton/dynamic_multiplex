#!/bin/bash
#SBATCH --job-name=dynmux_empirical
#SBATCH --account=open
#SBATCH --partition=open
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=slurm/logs/empirical_%j.out
#SBATCH --error=slurm/logs/empirical_%j.err

# ============================================================
# Empirical application (R) for manuscript
# Uses peacesciencer for IR network data (trade, IGO, alliance)
# ============================================================
set -euo pipefail

PROJECT_DIR="/scratch/jfe4/dynamic_multiplex"
cd "$PROJECT_DIR"
mkdir -p slurm/logs manuscript/output

# Load modules
module purge
module load r/4.5.0

# Install R package and dependencies
echo "=== Installing R package ==="
Rscript -e '
  if (!requireNamespace("remotes", quietly=TRUE))
    install.packages("remotes", repos="https://cloud.r-project.org")
  remotes::install_local("r_code", force=TRUE, quiet=TRUE)
'

echo "=== Checking R dependencies ==="
Rscript -e '
  pkgs <- c("ggplot2", "dplyr", "tidyr", "ggalluvial", "RColorBrewer", "peacesciencer")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
  if (length(missing) > 0) {
    cat("Installing:", paste(missing, collapse=", "), "\n")
    install.packages(missing, repos="https://cloud.r-project.org")
  }
  cat("All dependencies OK.\n")
'

# Run empirical application
echo ""
echo "=== Running empirical application ==="
Rscript manuscript/02_empirical_application.R

echo ""
echo "=== Empirical application complete ==="
ls -la manuscript/output/ 2>/dev/null || echo "  (no outputs yet)"
