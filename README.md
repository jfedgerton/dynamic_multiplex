# dynamic_multiplex

DynMux tracks communities in changing temporal networks, where nodes enter and
leave and groups form, split, merge and dissolve. This repository holds the R and
Python packages and the replication materials for the accompanying paper.

## Layout: packages vs. paper

| Directory | Contents | Released? |
|---|---|---|
| `r_code/` | R package `dynamicmultiplex` (CRAN) | yes |
| `python_code/` | Python package `dynamic_multiplex` (PyPI) | yes |
| `replication/` | Replication package for the paper: `sim/` (regime comparison, coverage study), `empirical/` (alliance, DCA, IGO, trade networks and order recovery), `post/` (every table and figure), `slurm/` (job files), `run_all.sh` (one entry point), `exploratory/` (archived, not run) | paper only |
| `manuscript/` | Generated `tables/` and `figures/` land here (gitignored until publication) | paper only |
| `scripts/` | Benchmarks | paper only |

Only `r_code/` and `python_code/` are built into the released packages.
Everything else is paper-side code and is not part of either package.
See `replication/README.md` for how to reproduce the paper's results.
