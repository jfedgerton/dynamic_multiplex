#!/bin/bash
#SBATCH --job-name=repl_atop
#SBATCH --output=replication/extended/slurm/logs/atop_%j.out
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --account=jfe4_cr_default
#SBATCH --partition=basic

cd "$SLURM_SUBMIT_DIR"
module load r/4.5.0
export R_LIBS_USER=/storage/home/jfe4/R/x86_64-pc-linux-gnu-library/4.5
export EMP_DATA="$SLURM_SUBMIT_DIR/replication/extended/output/empirical_data"
export EMP_OUT="$SLURM_SUBMIT_DIR/replication/extended/output/empirical"
mkdir -p "$EMP_DATA" "$EMP_OUT"
rm -f "$EMP_OUT/atop_partitions.rds"          # no stale cache
echo "=== build ATOP 1816:2018 ==="
Rscript replication/extended/00_build_atop.R
echo "=== run all specs + baselines on ATOP ==="
Rscript replication/extended/08_alliance_dca_empirical.R atop
echo "=== JOB DONE ==="
