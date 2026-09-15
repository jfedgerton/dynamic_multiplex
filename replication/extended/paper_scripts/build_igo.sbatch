#!/bin/bash
#SBATCH --job-name=build_igo
#SBATCH --account=open
#SBATCH --partition=basic
#SBATCH --time=1:00:00
#SBATCH --mem=48G
#SBATCH --output=replication/extended/paper_scripts/build_igo.out
cd "${DM_ROOT:?set DM_ROOT to the project root}"
module load r/4.5.0
# R package library: left unset so R resolves its own default user library
# (~/R/<platform>-library/<version>). Export R_LIBS_USER in your environment
# before submitting if your packages live elsewhere.
export EMP_DATA=${DM_ROOT}/replication/extended/output/empirical_data
export IGO_THRESH=1
Rscript replication/extended/paper_scripts/build_igo.R
