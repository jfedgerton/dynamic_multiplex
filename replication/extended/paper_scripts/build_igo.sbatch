#!/bin/bash
#SBATCH --job-name=build_igo
#SBATCH --account=open
#SBATCH --partition=basic
#SBATCH --time=1:00:00
#SBATCH --mem=48G
#SBATCH --output=replication/extended/paper_scripts/build_igo.out
cd /storage/work/jfe4/dynamic_multiplex
module load r/4.5.0
export R_LIBS_USER=/storage/home/jfe4/R/x86_64-pc-linux-gnu-library/4.5
export EMP_DATA=/storage/work/jfe4/dynamic_multiplex/replication/extended/output/empirical_data
export IGO_THRESH=1
Rscript replication/extended/paper_scripts/build_igo.R
