#!/bin/bash
#SBATCH --job-name=dm_cov_valued
#SBATCH --account=jfe4_cr_default
#SBATCH --partition=basic
#SBATCH --array=1-90%32
#SBATCH --time=72:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --output=replication/slurm/logs/03_cov_valued_%A_%a.out
# Weighted-network coverage arm: 720 tasks, 250 sims each (sim/03), PACKED
# 8 tasks per array element (element i runs tasks 8(i-1)+1 .. 8i) to stay
# under the 4,000 submitted-job cap on the account. Each task still writes
# its own cov/calib/joint files and self-skips if they exist.
set -euo pipefail
cd "${DM_ROOT:?set DM_ROOT to the project root}"
module load r/4.5.0
mkdir -p replication/slurm/logs
PACK=8; NTASK=720
for j in $(seq 1 $PACK); do
  T=$(( (SLURM_ARRAY_TASK_ID - 1) * PACK + j ))
  [ "$T" -le "$NTASK" ] || break
  echo "=== element $SLURM_ARRAY_TASK_ID -> task $T ($(date)) ==="
  SLURM_ARRAY_TASK_ID=$T COV_CORES=8 Rscript replication/sim/03_coverage_valued.R
done
