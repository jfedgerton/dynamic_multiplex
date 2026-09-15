#!/bin/bash
#SBATCH --job-name=repl_dca
#SBATCH --output=replication/extended/slurm/logs/dca_%j.out
#SBATCH --time=2:00:00
#SBATCH --mem=24G
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
rm -f "$EMP_OUT/dca_partitions.rds"          # only the DCA cache; ATOP results untouched
echo "=== build DCA (dcaAnyV1) ==="
Rscript replication/extended/00_build_dca.R
echo "=== run all specs + baselines on DCA ==="
Rscript replication/extended/08_alliance_dca_empirical.R dca
echo "=== JOB DONE ==="
