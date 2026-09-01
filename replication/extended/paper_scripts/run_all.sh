#!/bin/bash
module load r/4.5.0 2>/dev/null
cd /storage/work/jfe4/dynamic_multiplex || exit 1
mkdir -p manuscript/figures manuscript/tables
: > replication/extended/paper_scripts/gen.log
echo "START $(date)" >> replication/extended/paper_scripts/gen.log
Rscript replication/extended/paper_scripts/10_make_tables.R   >> replication/extended/paper_scripts/gen.log 2>&1
Rscript replication/extended/paper_scripts/13_bloc_validity.R >> replication/extended/paper_scripts/gen.log 2>&1
Rscript replication/extended/paper_scripts/11_make_figures.R  >> replication/extended/paper_scripts/gen.log 2>&1
for f in fig_regime fig_empirical_persistence fig_empirical_quality fig_bloc_validity fig_omega; do
  [ -f manuscript/figures/$f.pdf ] && gs -sDEVICE=png16m -r120 -dNOPAUSE -dBATCH -o manuscript/figures/$f.png manuscript/figures/$f.pdf >> replication/extended/paper_scripts/gen.log 2>&1
done
echo "ALL_DONE $(date)" >> replication/extended/paper_scripts/gen.log
