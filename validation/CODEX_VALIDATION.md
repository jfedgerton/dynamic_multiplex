# Codex: independent validation of the DynMux replication pipeline (v1.3.1)

Purpose: a second, independent check of the numbers in the *Political Analysis*
manuscript. Do NOT rerun the simulations. Do NOT modify anything under
`replication/`, `r_code/`, or `python_code/`. Write everything you produce under
`validation/` (create it) and report findings in `validation/REPORT.md`.

## Setup (Jared does this before you start)

1. Fresh clone of `jfedgerton/dynamic_multiplex` at tag/commit for v1.3.1.
2. `dynmux_output_v131.tgz` unpacked at the repo root (creates `output/`,
   `manuscript/tables`, `manuscript/figures`, `replication/slurm/logs`).
3. R >= 4.5 with igraph, ggplot2, pkgload, clue, dplyr; Python 3.12 with build.

## Tasks, in priority order

### 1. Re-derive every manuscript table from the raw CSVs (the real check)

For each file in `manuscript/tables/*.tex` dated 2026-09-21, recompute every
cell from `output/<arm>/*.csv` using your own aggregation code (pandas or
base R; do not call or copy `replication/post/*.R`). Sources:

| Table | Raw source | Aggregation |
|---|---|---|
| tab_metrics_wide | output/regime/dyn_cfg*.csv | mean joint NMI and K MAE by regime x method (18 cfg x 30 reps) |
| tab_app_paired_ci_nmi / _kmae, tab_paired_wilcoxon | output/regime | paired differences DynMux Jaccard minus each baseline, per series; 95% CI; Wilcoxon |
| tab_decision_tree | output/regime, output/mechanism, output/coupling | joint NMI by scenario x method; best and margin |
| tab_app_mechanism, tab_app_mechanism_cells | output/mechanism/mech_cfg*.csv | means by regime x method; paired DynMux minus multislice; win share |
| tab_app_coupling_regimes, tab_app_coupling_methods | output/coupling | Jaccard minus overlap by regime, break and continue conventions, paired CI |
| tab_app_omega_sweep | output/omega/omega_cfg*.csv | mean joint NMI by regime x variant |
| tab_app_selection_rule | output/selection | share of series where higher stability = higher NMI; mean NMI of rule / always-DynMux / always-multislice / oracle |
| tab_stability_floor, tab_app_stability_levels | output/stability (binary arm) | per stability bin: n, median accuracy, 5th percentile, validation share; odd-indexed cfg = calibration, even = validation |
| tab_app_stability_arms, tab_app_stability_lolo, tab_app_stability_pairs | output/stability | robustness arms; leave-one-level-out; decided/undetermined pairs |
| tab_app_empirical_stability | output/empirical/*_stability*.csv | s, floor lookup from stability_calibration_table.csv, share of pairs |
| tab_order_recovery_* | output/empirical | precision/recall/IoU by network x method |

Report any cell that differs from the .tex by more than rounding (0.005 for
NMI-type quantities, 0.01 for K MAE, 1 for counts). A clean result is
"all cells reproduce"; list the cells you could NOT reproduce and why.

### 2. Design and seed conventions

Check `replication/sim/*.R` and `replication/README.md` agree on: task counts
(regime 72, stability 594/216/72, mechanism 50, coupling 48, omega 72,
selection 72), replicates per task, and the seed formula for each script
(e.g. sim/03 uses 15000 + TASK*1000 + rep; sim/05 uses 9000 + TASK*1000 + rep).
Confirm `set.seed(123)` is used for any grid shuffle. Flag any script whose
config grid does not match the appendix design table.

### 3. Package and test checks

- `bash replication/run_all.sh test` (R regression tests) — must print
  ALL PACKAGE REGRESSION TESTS PASSED.
- `R CMD build r_code && R CMD check --as-cran dynamicmultiplex_1.3.1.tar.gz`
- `cd python_code && python -m build` and `pytest` if tests exist.
- `diff r_code/inst/extdata/stability_calibration_table.csv output/stability/stability_calibration_table.csv`
  and the same for `python_code/src/dynamic_multiplex/data/`.

### 4. Known items, do not re-report as new

- Multislice (`genlouvain_multislice`) is seed-sensitive on abrupt-rewiring
  series at n <= 100: two runs of the same series with different RNG state
  differ by up to 0.3 joint NMI per config. This is why the omega sweep
  (sim/05, 10 reps) shows 0.41 on abrupt rewiring while Table 2 (sim/01,
  30 reps) shows 0.49. Code paths were verified identical.
- IGO multislice at resolution 2 needs ~40 min per fit; the per-method cap in
  empirical/08 is 7200 s.
- Python tests are skipped on Roar (no pip in the module); run them locally.

### 5. After the prose is final (Wednesday)

- Extract every number in `manuscript/main.tex` and `appendix.tex` and match
  it to a cell in `manuscript/tables/*.tex` or a value in
  `replication/slurm/logs/11_postprocess_*.out`. List mismatches.
- List every sentence that claims DynMux outperforms multislice in general
  (the supported claim is asymmetric: DynMux when the node set or community
  set changes; multislice when nodes persist).
- List every `\includegraphics` and `\input` whose file has no generating
  script under `replication/`.

## Output

`validation/REPORT.md` with one section per task, a table of reproduced vs
non-reproduced cells, and the exact commands you ran. Keep code under
`validation/` so it can be committed as an independent check.
