#!/bin/bash
#SBATCH --job-name=bdfix
#SBATCH --output=paper_scripts/logs/bdfix_%a.out
#SBATCH --array=1-72
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --account=open
#SBATCH --partition=basic
cd "$SLURM_SUBMIT_DIR"
module load r/4.5.0
export R_LIBS_USER=/storage/home/jfe4/R/x86_64-pc-linux-gnu-library/4.5
Rscript paper_scripts/06_bd_fix.R
