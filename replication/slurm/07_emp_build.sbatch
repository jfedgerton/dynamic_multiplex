#!/bin/bash
#SBATCH --job-name=dm_emp_build
#SBATCH --account=jfe4_cr_default
#SBATCH --partition=basic
#SBATCH --time=2:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=2
#SBATCH --output=replication/slurm/logs/07_emp_build_%j.out
# Build the four empirical networks into output/empirical_data/ (empirical/07). Needs peacesciencer extdata and data/DCAD-v1.0-dyadic.csv.
set -euo pipefail
cd "${DM_ROOT:?set DM_ROOT to the project root}"
module load r/4.5.0
mkdir -p replication/slurm/logs
Rscript replication/empirical/07_build_networks.R
