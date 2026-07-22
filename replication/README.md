# Replication package — Dynamic Multiplex Community Detection

This folder reproduces every computational result in the paper. It is the same
code and the same seeds used in the study, reorganized into three self-contained
pipelines plus a master script. Nothing here needs to be run to read the paper;
run it only to regenerate the results (e.g., for a revision).

All randomness is seeded. A rerun reproduces the reported numbers exactly.

## Paper outline → code map

| Paper section | What it shows | Script | SLURM |
|---|---|---|---|
| Study I — Method comparison | DynMux vs. 10 competitors + Dynamic SBM on synthetic dynamic networks; headline result is tracking a *changing* number of communities | `02_method_comparison.R` | `slurm/02_method_comparison.sbatch` |
| Study II — Uncertainty quantification | Bootstrap CI coverage for community count and node co-assignment; the observable reliability **gate**; the degree-heterogeneity stress test | `01_ci_coverage.R` | `slurm/01_ci_coverage.sbatch` |
| Study III — Empirical (Cold War) | Alliance / IGO / trade community persistence 1965–2000 and the Cold-War break | `03_coldwar_empirical.R` | `slurm/03_coldwar_empirical.sbatch` |

The R and Python packages themselves live in `../r_code` and `../python_code`.
The method registry used by Study I (base fitters and competitors) is defined
in `../manuscript/01_synthetic_experiments.R`, which `02_method_comparison.R`
sources.

## Directory layout

```
replication/
  README.md                     this file
  run_all.sh                    master: install deps + submit all pipelines
  01_ci_coverage.R              Study II (binary + valued + misspec arms)
  02_method_comparison.R        Study I (11 methods incl. Dynamic SBM)
  03_coldwar_empirical.R        Study III (alliance / IGO / trade)
  slurm/
    01_ci_coverage.sbatch
    02_method_comparison.sbatch
    03_coldwar_empirical.sbatch
  output/                       all results land here (created on run)
```

## How to reproduce

1. Edit `--account` and `--partition` in `slurm/*.sbatch` to your allocation,
   and set `R_LIBS_USER` to a writable library path.
2. Install dependencies:
   ```
   bash replication/run_all.sh deps
   ```
3. Submit the pipelines (independent; any order, or all at once):
   ```
   bash replication/run_all.sh submit all
   ```
   or individually: `submit ci`, `submit compare`, `submit empirical`.

Each `.sbatch` documents its own resource needs and any per-arm submission
fallback (see the CI note about the queued-array cap).

## The three pipelines in detail

### (a) `01_ci_coverage.R` — CI coverage + gate

One arm-dispatched program sharing a single edge-resampling bootstrap and
estimand core. The global SLURM array index selects the arm:

| Arm | Networks | Array range | Seed band |
|---|---|---|---|
| `binary` | planted-partition SBM, full factorial grid (3564 tasks) | 1–3564 | 15000 |
| `valued` | weighted networks, aligned/orthogonal weights (720 tasks) | 3565–4284 | 20000 |
| `misspec` | degree-corrected SBM, degree heterogeneity + size skew (432 tasks) | 4285–4716 | 12000 |

Each task runs 250 simulations × 100 bootstrap replicates × 100 fresh-draw
ground-truth replicates. Per-sim seed = `(band + task) * 1e5 + sim_id`.
Outputs: `output/coverage_<arm>/cov_task#####.csv` (per-sim coverage) and
`calib_task#####.csv` (20-bin calibration).

**The gate.** Accept an interval set iff
`width_K_mean <= 1  AND  width_P_mean < 0.05`.
Coverage is nominal inside this observable region and refused outside it.

### (b) `02_method_comparison.R` — head-to-head, 11 methods

One array task per config (22 configs). Each task runs 100 reps × all methods.
Methods: 4 DynMux specs, DynMux (filtered, leakage-free), Independent Leiden
(± greedy matching), Aggregated Leiden, Multislice (Adjacent), multinet
GLouvain, and Dynamic SBM (dynsbm, ICL-selected over Q ∈ 2..10). Two arms:
`borrow` (persistent structure below the per-layer detectability threshold) and
`mergesplit` (the true K changes 6→3 or 3→6 mid-series). Metrics per method:
per-layer NMI, joint (node × layer) NMI, and K MAE. Per-rep seed =
`9000 + cfg*1000 + rep`. Output: `output/comparison/comparison_cfg##.csv`.

Dynamic SBM is ~1000× slower than the other methods; it is given a generous
K range (Qmax = 10) so its results cannot be attributed to too small a
candidate set. If `dynsbm` is not installed the script runs the other 10
methods and skips that column.

### (c) `03_coldwar_empirical.R` — Cold-War empirical

The untouched empirical analysis, with the output directory redirected here.
Fits alliance, IGO, and trade networks (peacesciencer, 1965–2000) under all
five interlayer weightings, compares against multislice Louvain, and measures
community persistence — the diagnostic in which the Cold-War break appears
(sharp in alliances, corroborated in IGOs, absent/gradual in trade).
Confidence intervals are intentionally not drawn on the empirical figures;
the CI machinery is validated separately in pipeline (a).

## Aggregating the results

After the jobs finish, the headline numbers are simple aggregations of the
per-task CSVs.

**Gate coverage (Study II).** Pool `output/coverage_*/cov_task*.csv`, then:
```r
g <- d$width_K_mean <= 1 & d$width_P_mean < 0.05
c(pass = mean(g),
  cov_K = mean(d$cov_K_mean[g]),
  cov_P = mean(d$cov_P_mean[g]))
```
Out-of-sample validation: calibrate the thresholds on an early chunk of tasks
and evaluate on the rest. Break out `misspec` by `hetero` (none/moderate/
severe) to read the ship-or-die stress test.

**Method comparison (Study I).** Pool `output/comparison/comparison_cfg*.csv`,
then aggregate the three metrics by `arm` and `method`; 95% CIs are
`mean ± 1.96 * sd/sqrt(reps)` across the 100 reps.

## Dependencies

R (≥ 4.5): `dynamicmultiplex` (from `../r_code`), `igraph`, `ggplot2`, `dplyr`,
`tidyr`, `parallel`, `dynsbm`, `peacesciencer`, `ggalluvial`.
Python (optional, cross-language checks): the package in `../python_code`.

## Determinism / provenance

Every simulation seed is fixed and documented above. Iteration order within
each CI arm is a fixed shuffle (`set.seed(123)`), which only affects which
tasks finish first, not their contents. Given the same code, seeds, and package
version, results are bit-for-bit reproducible.
