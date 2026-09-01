#!/bin/bash
#SBATCH --job-name=atop_mn_pres
#SBATCH --output=paper_scripts/logs/atop_mn_pres_%j.out
#SBATCH --time=1:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
#SBATCH --account=open
#SBATCH --partition=basic
cd "$SLURM_SUBMIT_DIR"
module load r/4.5.0
export R_LIBS_USER=/storage/home/jfe4/R/x86_64-pc-linux-gnu-library/4.5
export EMP_DATA="$SLURM_SUBMIT_DIR/replication/extended/output/empirical_data"
export EMP_OUT="$SLURM_SUBMIT_DIR/replication/extended/output/empirical"
Rscript paper_scripts/multinet_fix.R atop
echo "ATOP_MN_PRESENT_DONE"
