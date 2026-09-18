#!/bin/bash
#SBATCH --job-name=igo_parts
#SBATCH --account=open
#SBATCH --partition=basic
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --output=replication/extended/paper_scripts/igo_parts.out
cd "${DM_ROOT:?set DM_ROOT to the project root}"
module load r/4.5.0
# R package library: left unset so R resolves its own default user library
# (~/R/<platform>-library/<version>). Export R_LIBS_USER in your environment
# before submitting if your packages live elsewhere.
export EMP_DATA=${DM_ROOT}/replication/extended/output/empirical_data
export EMP_OUT=${DM_ROOT}/replication/extended/output/empirical
Rscript replication/extended/08_alliance_dca_empirical.R igo
