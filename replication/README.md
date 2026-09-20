# Replication package: DynMux (dynamic multiplex community detection)

Everything the manuscript reports is produced by the scripts here from the
R package in `../r_code` (Python mirror in `../python_code`). One entry point:

```
source replication/env.sh                 # sets DM_ROOT, loads R on Roar
bash replication/run_all.sh deps          # install the package into the user library
bash replication/run_all.sh test          # package regression tests
bash replication/run_all.sh submit        # every job, with dependencies (Roar / SLURM)
bash replication/run_all.sh post          # tables and figures, once the jobs are done
```

## Pipeline

| Step | Script | Tasks | Produces | Used by |
|---|---|---|---|---|
| 1 | `sim/01_regime_comparison.R` | 72 | `output/regime/` | Table 2, appendix A2 (`post/10`, `post/11`) |
| 2 | `sim/02_stability.R` (`STAB_ARM` = binary / dcsbm / weighted) | 594 / 216 / 72 | `output/stability/` | Section 4, Figure 2, appendix A4 (`post/12`) |
| 4 | `sim/04_coupling_regimes.R` | 48 | `output/coupling/` | appendix A3 (`post/14`) |
| 5 | `sim/05_omega_sweep.R` | 72 | `output/omega/` | appendix omega sweep (`post/15`) |
| 6 | `sim/06_selection_rule.R` | 72 | `output/selection/` | appendix selection rule (`post/15`) |
| 7-9 | `empirical/07_build_networks.R`, `08_fit_networks.R`, `09_score_orders.R` | 1 / 4 / 1 | `output/empirical_data/`, `output/empirical/` | Figure 3, appendix empirical tables (`post/10`, `post/13`) |
| 10 | `empirical/10_stability_networks.R` | 4 | `output/empirical/<net>_stability.csv` | appendix empirical stability (`post/15`) |
| post | `post/10`-`15` | | `manuscript/tables/`, `manuscript/figures/`, calibration table copied into both packages | |

`sim/lib_regimes.R` holds the generators, method wrappers and metrics shared
by sims 1, 5 and 6. Every sim task self-skips when its output exists, so
`submit` can be rerun after a partial failure. All seeds are fixed
(`set.seed(123)` for grids; per-task seeds documented in each script header).

## What each result is

* **Table 2** (`tab_metrics_wide.tex`): joint NMI and K MAE by regime for
  DynMux (Jaccard coupling), multislice with adjacent links and with the same
  links DynMux receives, Hungarian matching, dynamic SBM, pooled Leiden and
  multinet. The overlap coupling appears in the appendix paired tables only.
* **Section 4** (`tab_stability_floor.tex`, `fig_stability_floor.pdf`): the
  bootstrap stability score and its calibrated accuracy floor. Calibrated on
  odd-indexed cells of the binary arm, validated on even-indexed cells and on
  the degree-corrected and weighted arms; `post/12` applies the pre-registered
  pass/fail rule and exits non-zero on FAIL.
* **Appendix**: paired differences (A2), Jaccard vs overlap regimes (A3),
  stability calibration details, leave-one-level-out, decided pairs, arms
  (A4), omega sweep, selection rule, empirical stability, order recovery.

## Requirements

R >= 4.1 with igraph >= 2.0, ggplot2, pkgload, clue; optional multinet and
dynsbm (Table 2 rows are skipped when absent); peacesciencer with its extdata
downloaded and `data/DCAD-v1.0-dyadic.csv` for the empirical networks.
Python 3.10+ with numpy, pandas, networkx, python-igraph, leidenalg for the
Python package tests. Slurm headers assume Roar (`--account=jfe4_cr_default
--partition=basic`).

## History

`exploratory/` keeps superseded material: `coverage_intervals_v1.2/` (the
retired node co-membership interval study and why it was retired),
`replication_v1/`, `paper_scripts/`, `manuscript_sims/`, `extended/`.
