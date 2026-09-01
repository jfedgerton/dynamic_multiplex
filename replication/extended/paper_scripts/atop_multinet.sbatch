#!/bin/bash
#SBATCH --job-name=atop_multinet
#SBATCH --output=replication/extended/paper_scripts/atop_multinet_%j.out
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --account=jfe4_cr_default
#SBATCH --partition=basic
cd "$SLURM_SUBMIT_DIR"
module load r/4.5.0
export R_LIBS_USER=/storage/home/jfe4/R/x86_64-pc-linux-gnu-library/4.5
export EMP_DATA="$SLURM_SUBMIT_DIR/replication/extended/output/empirical_data"
export EMP_OUT="$SLURM_SUBMIT_DIR/replication/extended/output/empirical"
echo "=== refill ATOP multinet only (11 methods cached, multinet computed) ==="
Rscript replication/extended/paper_scripts/emp_fix.R atop
echo "=== JOB DONE ==="
