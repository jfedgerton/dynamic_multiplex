# Replication package -- DynMux (dynamic multiplex community detection)

Reproduces every table and figure in the manuscript from scratch. All
randomness is seeded (123 for design shuffles; per-task seeds derived from the
task id), so a rerun reproduces the reported numbers exactly, up to
differences in the installed versions of `igraph`/`leidenalg`.

```
replication/
  env.sh                 source it: sets DM_ROOT, loads R on Roar
  run_all.sh             deps | submit | post | status  (see below)
  sim/
    01_regime_comparison.R    72 configs x 30 reps x 16 methods       -> output/regime/
    02_coverage_grid.R        binary SBM coverage, 3564 tasks x 250   -> output/coverage_grid/
    03_coverage_valued.R      weighted arm, 720 tasks                 -> output/coverage_valued/
    04_coverage_misspec.R     degree-corrected arm, 432 tasks         -> output/coverage_misspec/
  empirical/
    05_build_networks.R       ATOP, DCA, IGO, trade year-layers       -> output/empirical_data/
    06_fit_networks.R <net>   12 methods per network                  -> output/empirical/<net>_partitions.rds
    07_score_orders.R         Braumoeller-order recovery              -> output/empirical/order_recovery.csv
  post/
    10_main_text.R            Table 2, Figure 2, Figure 3
    11_appendix_regimes.R     paired-difference figures, CI and Wilcoxon tables
    12_appendix_coverage.R    gate ladder, gate grid, per-spec, width bands, arms
    13_appendix_empirical.R   order-recovery tables
    14_calibration_table.R    calibrated interval lookup + validation tables
  slurm/                 one sbatch per step (01-08), --account=jfe4_cr_default --partition=basic
  exploratory/           archived scripts, not run (see its README)
```

Outputs: simulation and empirical results under `$DM_ROOT/output/`, tables
(bare booktabs tabulars, no captions) under `$DM_ROOT/manuscript/tables/`,
figures (PDF + PNG) under `$DM_ROOT/manuscript/figures/`. `DM_ROOT` is the
repo root; every script reads it and nothing else is configured.

## Running on Roar

```bash
cd /storage/group/LiberalArts/default/jfe4_collab/dynamic_multiplex
source replication/env.sh
cp /path/to/DCAD-v1.0-dyadic.csv data/          # Kinne DCAD v1.0, not redistributed
bash replication/run_all.sh deps                 # R CMD INSTALL r_code (+ pip -e python_code)
bash replication/run_all.sh submit               # everything, with dependencies
bash replication/run_all.sh status               # squeue + task-file counts
```

`submit` queues 4,860 array tasks (01: 72, 02: 3,564 in four 990-task chunks,
03: 720, 04: 432) plus the empirical chain 05 -> 06 (x4) -> 07, and a final
`08_postprocess` job that waits for all of them and runs `run_all.sh post`.
Every sim task writes one file per task and exits immediately if that file
exists, so after a partial failure rerun `submit` and only missing tasks run.
Post-processing can also be run by hand (`bash replication/run_all.sh post`,
about 20 minutes, needs ~30 GB for the joint files).

Roar-specific assumptions: `module load r/4.5.0`; R packages in the user
library (`dynamicmultiplex` is installed by `deps`; `igraph >= 2.0`,
`ggplot2`, `pkgload`, `clue`, `multinet`, `dynsbm`, `peacesciencer`, `dplyr`
must already be there; `peacesciencer::download_extdata()` must have been run
once). The sims never touch the network.

## What the manuscript uses

| Manuscript | File | Generator |
|---|---|---|
| Table 2 | `tables/tab_metrics_wide.tex` (full float with caption) | `post/10` |
| Figure 2 | `figures/fig_coverage_by_config.pdf` | `post/10` |
| Figure 3 | `figures/fig_order_recovery.pdf` | `post/10` |
| App. regimes | `figures/fig_regime_paired_{nmi,kmae}.pdf`, `tables/tab_app_paired_ci_{nmi,kmae}.tex`, `tables/tab_paired_wilcoxon.tex` | `post/11` |
| App. coverage | `tables/tab_coverage.tex`, `tab_coverage_gate_grid.tex`, `tab_coverage_spec.tex`, `tab_coverage_width_bands.tex`, `tab_coverage_valued.tex`, `tab_coverage_misspec{,_hetero,_balance}.tex`; `figures/fig_coverage_curve.pdf`, `fig_misspec_curve.pdf` | `post/12` |
| App. empirical | `tables/tab_order_recovery_summary.tex`, `tab_order_recovery_{atop,dca,igo,trade}.tex` | `post/13` |
| App. calibrated interval | `tables/tab_calibrated_interval.tex`, `tab_calibrated_by_design.tex`, `tab_calibrated_arms.tex`; `figures/fig_calibrated_interval.pdf` | `post/14` |

## The calibrated co-assignment interval

The coverage study evaluates the Wilson interval that treats the B = 100
bootstrap replicates as binomial draws. That interval's width scales as
1/sqrt(B), so it measures Monte Carlo error in the bootstrap rather than
uncertainty about the estimand (the fresh-draw co-assignment propensity p*),
and its coverage depends on where the bootstrap share p-hat falls.

`sim/02`-`04` therefore also record, for every node pair, the joint
distribution of (p-hat bin, p* bin) on a 50 x 50 grid (`joint_task*.csv`).
`post/14` pools those counts over the calibration half of the design, takes
the conditional 2.5% and 97.5% quantiles of p* within each p-hat bin, and
writes `output/calibration/coassign_calibration_table.csv`. That table is the
interval: its width does not depend on B, its conditional coverage given
p-hat is >= 0.95 on the calibration split by construction, and the validation
split, the weighted arm and the degree-corrected arm are the out-of-sample
checks (`tab_calibrated_*.tex`). `run_all.sh post` copies the table into
`r_code/inst/extdata/` and `python_code/src/dynamic_multiplex/data/`, where
`co_assignment_ci(method = "calibrated")` (the new default in both packages)
reads it; `method = "wilson"` keeps the old interval. Reinstall the packages
after post-processing so the bundled table is current, then commit it.

## Package changes made with this compression (r_code 1.2.0, python 1.2.0)

- `co_assignment_ci()` gains `method = c("calibrated", "wilson")` and
  `calibration_table`; the calibrated table is bundled as package data.
- `fit_multilayer_{jaccard,overlap,weighted_jaccard,weighted_overlap}()` gain
  `allow_unequal_nodes` (already present on `fit_multilayer_identity_ties`);
  with it, layers with different node sets are matched by vertex name.
  Unequal node sets without names are rejected. `sim/01` (birth-death regime)
  and `empirical/06` (states entering and exiting) rely on this; the earlier
  `assignInNamespace()` workaround is gone.
- `layer_node_strengths()` keys strengths by vertex name when graphs are
  named (previously by position, which mis-keyed weighted couplings on any
  named graph with unequal node sets).

## Design notes

- Coverage gate: `width_P_mean < 0.05 & n >= 100`, chosen on the calibration
  split (odd-indexed sorted configurations), reported on the validation split.
- `sim/01` folds the birth-death fix (`06_bd_fix.R`) into the main script:
  coupling and tracking methods see only the nodes alive in each layer.
- Order codings in `empirical/07` are Braumoeller/Goodhart as transcribed in
  the archived scripts, including the Warsaw Pact ccode-100 entry, which is
  flagged in a NOTE and left unchanged so results match the sources.
