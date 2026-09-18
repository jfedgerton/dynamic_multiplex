# exploratory/ -- archived scripts (not part of the replication)

Everything here was superseded by `replication/{sim,empirical,post}` in the
September 2026 code compression. Nothing in this folder is run by
`run_all.sh`, nothing in the manuscript depends on it, and paths inside these
files (`/storage/work/jfe4/...`, `setwd()`, `manuscript/output/...`) are the
ones that were in use at the time and are not maintained. Kept for provenance
only: each current script names the archived file(s) it consolidates in its
header.

| Folder | What it was | Superseded by |
|---|---|---|
| `manuscript_sims/` | the sim scripts that lived in Dropbox `manuscript/` and were never in git: `16_coverage3_grid.R`, `17_coverage3_valued.R`, `19_coverage_misspec.R` (the coverage study), plus earlier coverage designs `12`-`15`, `18_comparison_extensions.R`, `20_dynsbm_comparison.R` | `sim/02`-`04` |
| `extended/` | `04_regime_comparison.R` (regime comparison, birth-death methods fed the full node set), `00_build_atop.R`, `00_build_dca.R`, `06_alliance_dca_empirical.R`, `09`-`13` (earlier empirical analyses, omega sweep, bloc validity) | `sim/01`, `empirical/05`-`06` |
| `paper_scripts/` | the 60-odd per-figure/per-table scripts: `06_bd_fix.R` + `splice_bdfix.R` (birth-death rerun, now folded into `sim/01`), `15b`, `21`, `33` (paired figures/tests/CI tables -> `post/11`), `19` (Table 2 -> `post/10`), `20`, `34` (coverage gate -> `post/12`), `23` (Figure 2 -> `post/10`), `29`, `32` (order codings and set-level scoring -> `empirical/07`), `emp_fix.R`, `multinet_fix.R`, `build_igo.R`, `build_trade.R` (-> `empirical/05`-`06`), and diagnostics that never reached the paper (`dca_check*.R`, `leakage_*.R`, `atop_*.R`, `frac_analysis.R`, `net_analysis.R`, `temporal_plot.R`, `24`-`28`, `30`-`31`) | see column 2 |
| `replication_v1/` | the first replication layout (`01_ci_coverage.R`, `02_method_comparison.R`, `03_coldwar_empirical.R`, their sbatch files, README) that ran 3 of the 60 scripts | `run_all.sh` |
| `slurm_v1/` | sbatch files for the Dropbox sims (`--account=open`, task 2971-3564 gap, hard-coded paths) | `slurm/` |

Known defects that were fixed on the way out, listed so nobody re-imports them:

- `dca_check.R`, `dca_check2.R`: literal `"${DM_ROOT}"` string instead of `Sys.getenv`.
- `18_comparison_extensions.R`: three independent breakages; never produced the table it was written for.
- `16_coverage3_grid.sh`: array covered tasks 1-2970 of 3564.
- `23_coverage_alt_figs.R`: hard-coded `/storage/group/.../jfe4_collab` path.
- `tab_coverage.tex` had three writers, `tab_metrics_wide.tex` two, `fig_temporal.png` two.
- `06_bd_fix.R` relied on a package version whose fitters accepted unequal node sets silently; the current package requires `allow_unequal_nodes = TRUE`, which `sim/01` passes.
- `32_setlevel_allorders.R`: Warsaw Pact coding lists ccode 100 (Colombia) where Bulgaria (355) was meant; `empirical/07` reproduces the coding as-is and flags it in a NOTE.
