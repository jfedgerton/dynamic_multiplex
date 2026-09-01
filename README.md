# dynamic_multiplex

DynMux: temporal community detection for multiplex networks with customizable
interlayer coupling, plus the replication materials for the accompanying paper.

## Layout: packages vs. paper

| Directory | Contents | Released? |
|---|---|---|
| `r_code/` | R package `dynamicmultiplex` (CRAN) | yes |
| `python_code/` | Python package `dynamic_multiplex` (PyPI) | yes |
| `replication/` | Replication pipelines for the paper (Studies I-III; `extended/` for the regime comparison, coverage study, empirical application, and `extended/paper_scripts/` for figure/table generation) | paper only |
| `manuscript/` | Legacy synthetic pipeline; generated figures/tables are gitignored | paper only |
| `scripts/`, `slurm/` | Benchmarks and cluster job files | paper only |

Only `r_code/` and `python_code/` are built into the released packages.
Everything else is paper-side code and is not part of either package.
See `replication/README.md` for how to reproduce the paper's results.
