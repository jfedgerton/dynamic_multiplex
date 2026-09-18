#!/bin/bash
#SBATCH --job-name=dm_cov_misspec
#SBATCH --account=jfe4_cr_default
#SBATCH --partition=basic
#SBATCH --array=1-54%32
#SBATCH --time=72:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --output=replication/slurm/logs/04_cov_misspec_%A_%a.out
# Degree-corrected (DC-SBM) coverage arm: 432 tasks, 250 sims each (sim/04),
# PACKED 8 tasks per array element (see 03_coverage_valued.sbatch).
set -euo pipefail
cd "${DM_ROOT:?set DM_ROOT to the project root}"
module load r/4.5.0
mkdir -p replication/slurm/logs
PACK=8; NTASK=432
for j in $(seq 1 $PACK); do
  T=$(( (SLURM_ARRAY_TASK_ID - 1) * PACK + j ))
  [ "$T" -le "$NTASK" ] || break
  echo "=== element $SLURM_ARRAY_TASK_ID -> task $T ($(date)) ==="
  SLURM_ARRAY_TASK_ID=$T COV_CORES=8 Rscript replication/sim/04_coverage_misspec.R
done
