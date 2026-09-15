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
# R package library: left unset so R resolves its own default user library
# (~/R/<platform>-library/<version>). Export R_LIBS_USER in your environment
# before submitting if your packages live elsewhere.
export EMP_DATA="$SLURM_SUBMIT_DIR/replication/extended/output/empirical_data"
export EMP_OUT="$SLURM_SUBMIT_DIR/replication/extended/output/empirical"
echo "=== refill ATOP multinet only (11 methods cached, multinet computed) ==="
Rscript replication/extended/paper_scripts/emp_fix.R atop
echo "=== JOB DONE ==="
