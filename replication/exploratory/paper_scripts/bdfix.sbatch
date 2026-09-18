#!/bin/bash
#SBATCH --job-name=bdfix
#SBATCH --output=replication/extended/paper_scripts/logs/bdfix_%a.out
#SBATCH --array=1-72
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --account=open
#SBATCH --partition=basic
cd "$SLURM_SUBMIT_DIR"
module load r/4.5.0
# R package library: left unset so R resolves its own default user library
# (~/R/<platform>-library/<version>). Export R_LIBS_USER in your environment
# before submitting if your packages live elsewhere.
Rscript replication/extended/paper_scripts/06_bd_fix.R
