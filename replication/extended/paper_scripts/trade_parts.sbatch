#!/bin/bash
#SBATCH --job-name=trade_parts
#SBATCH --account=open
#SBATCH --partition=basic
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --output=replication/extended/paper_scripts/trade_parts.out
cd /storage/work/jfe4/dynamic_multiplex
module load r/4.5.0
export R_LIBS_USER=/storage/home/jfe4/R/x86_64-pc-linux-gnu-library/4.5
export EMP_DATA=/storage/work/jfe4/dynamic_multiplex/replication/extended/output/empirical_data
export EMP_OUT=/storage/work/jfe4/dynamic_multiplex/replication/extended/output/empirical
Rscript replication/extended/08_alliance_dca_empirical.R trade
