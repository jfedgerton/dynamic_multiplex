#!/usr/bin/env bash
# =============================================================================
# replication/run_all.sh  --  one entry point for the full replication
#
#   bash replication/run_all.sh deps      install the R package (+ Python) locally
#   bash replication/run_all.sh submit    submit every SLURM job with dependencies
#   bash replication/run_all.sh submit sims       only sim/01-04
#   bash replication/run_all.sh submit empirical  only empirical/05-07
#   bash replication/run_all.sh submit rest J1:J2  sim/03-04 + empirical + post, with the
#                                             post job also waiting on job ids J1:J2
#                                             (use after a partial submit hit the
#                                             4,000 submitted-job cap)
#   bash replication/run_all.sh submit alt   EXPLORATORY: sim/05 (594 tasks) + post/15
#   bash replication/run_all.sh submit stab  EXPLORATORY: sim/06 (594 tasks) + post/16
#   bash replication/run_all.sh submit refit  sim/01 + empirical/05-07 + post, after a
#                                             package fix (archives output/regime first
#                                             so the 72 regime tasks do not self-skip)
#   bash replication/run_all.sh submit coupling  sim/07 (48 tasks) + post/17 (appendix:
#                                             when Jaccard vs overlap coupling is right)
#   bash replication/run_all.sh test      package regression tests (R + Python)
#   bash replication/run_all.sh post      run post/10-14 in the current shell
#   bash replication/run_all.sh post alt  run post/15 (interval alternatives) only
#   bash replication/run_all.sh post coupling  run post/17 only
#   bash replication/run_all.sh status    squeue for this user's dm_* jobs
#   bash replication/run_all.sh clean-tables   remove generated tables/figures
#
# Pipeline (every script reads DM_ROOT; defaults to this repo's root):
#   sim/01_regime_comparison.R    -> output/regime/dyn_cfgNN.csv          (72 tasks)
#   sim/02_coverage_grid.R        -> output/coverage_grid/*_taskNNNNN.csv (3564 tasks)
#   sim/03_coverage_valued.R      -> output/coverage_valued/               (720 tasks)
#   sim/04_coverage_misspec.R     -> output/coverage_misspec/              (432 tasks)
#   empirical/05_build_networks.R -> output/empirical_data/<net>_{series,union}.rds
#   empirical/06_fit_networks.R   -> output/empirical/<net>_partitions.rds (4 tasks)
#   empirical/07_score_orders.R   -> output/empirical/order_recovery.csv
#   post/10_main_text.R           -> Table 2, Figure 2, Figure 3
#   post/11_appendix_regimes.R    -> paired-difference figures and CI tables
#   post/12_appendix_coverage.R   -> gate ladder, gate grid, spec, width bands, arms
#   post/13_appendix_empirical.R  -> order-recovery tables
#   post/14_calibration_table.R   -> calibrated interval lookup + validation tables
#                                    (copied into r_code/inst/extdata and
#                                     python_code/src/dynamic_multiplex/data)
# Tables land in manuscript/tables/, figures in manuscript/figures/.
#
# Every sim task self-skips when its output file exists, so `submit` can be
# rerun after a partial failure and only the missing tasks execute.
# Roar caps submitted jobs at 4,000 per account and every array task counts:
# 01 (72) + 02 (3,564) = 3,636, so 03 and 04 are packed 8 tasks per array
# element (90 + 54 jobs) and the empirical chain adds 7.
#
# Data you must supply: $DM_ROOT/data/DCAD-v1.0-dyadic.csv (Kinne DCAD v1.0).
# peacesciencer's extdata (ATOP, COW trade, IGO) must already be downloaded
# into the R user library: Rscript -e 'peacesciencer::download_extdata()'.
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export DM_ROOT="${DM_ROOT:-$(dirname "$HERE")}"
cd "$DM_ROOT"
ACTION="${1:-help}"
WHICH="${2:-all}"                              # all | sims | empirical | rest
SB="replication/slurm"
mkdir -p "$SB/logs" output manuscript/tables manuscript/figures

have_module() { command -v module >/dev/null 2>&1; }
load_r() { if have_module; then module load r/4.5.0; fi; }

case "$ACTION" in

  deps)
    load_r
    echo "--- R: installing dynamicmultiplex from r_code/ into the user library ---"
    LIB="$(Rscript -e 'cat(Sys.getenv("R_LIBS_USER"))')"; mkdir -p "$LIB"
    R CMD INSTALL -l "$LIB" r_code
    echo "--- R: checking optional packages used by the sims ---"
    Rscript -e 'for (p in c("igraph","ggplot2","pkgload","clue","multinet","dynsbm","peacesciencer","dplyr"))
                  cat(sprintf("%-14s %s\n", p, if (requireNamespace(p, quietly=TRUE)) "ok" else "MISSING"))'
    if command -v pip >/dev/null 2>&1; then
      echo "--- Python: installing dynamic_multiplex from python_code/ ---"
      pip install -e python_code --quiet || echo "(pip install failed; Python package is optional for the replication)"
    fi
    ;;

  submit)
    command -v sbatch >/dev/null 2>&1 || { echo "sbatch not found: run this on the cluster login node"; exit 1; }
    DEPS=()
    EXTRA="${3:-}"                                   # job ids already queued (submit rest J1:J2)
    if [[ -n "$EXTRA" ]]; then IFS=':' read -r -a ex <<< "$EXTRA"; DEPS+=("${ex[@]}"); fi
    if [[ "$WHICH" == "refit" ]]; then
      if [[ -d output/regime ]] && ls output/regime/dyn_cfg*.csv >/dev/null 2>&1; then
        stamp=$(date +%Y%m%d_%H%M); mkdir -p output/_archive
        mv output/regime "output/_archive/regime_$stamp"; echo "archived old regime output to output/_archive/regime_$stamp"
      fi
    fi
    if [[ "$WHICH" == "all" || "$WHICH" == "sims" || "$WHICH" == "refit" ]]; then
      j=$(sbatch --parsable --export=ALL "$SB/01_regime.sbatch");                            echo "01 regime            $j"; DEPS+=("$j")
      j=$(sbatch --parsable --export=ALL,COV_OFFSET=0    "$SB/02_coverage_grid.sbatch");     echo "02 grid   tasks 1-990     $j"; DEPS+=("$j")
      j=$(sbatch --parsable --export=ALL,COV_OFFSET=990  "$SB/02_coverage_grid.sbatch");     echo "02 grid   tasks 991-1980  $j"; DEPS+=("$j")
      j=$(sbatch --parsable --export=ALL,COV_OFFSET=1980 "$SB/02_coverage_grid.sbatch");     echo "02 grid   tasks 1981-2970 $j"; DEPS+=("$j")
      j=$(sbatch --parsable --export=ALL,COV_OFFSET=2970 --array=1-594%64 "$SB/02_coverage_grid.sbatch"); echo "02 grid   tasks 2971-3564 $j"; DEPS+=("$j")
    fi
    if [[ "$WHICH" == "all" || "$WHICH" == "sims" || "$WHICH" == "rest" ]]; then
      j=$(sbatch --parsable --export=ALL "$SB/03_coverage_valued.sbatch");                   echo "03 valued  (90 packed)  $j"; DEPS+=("$j")
      j=$(sbatch --parsable --export=ALL "$SB/04_coverage_misspec.sbatch");                  echo "04 misspec (54 packed)  $j"; DEPS+=("$j")
    fi
    if [[ "$WHICH" == "all" || "$WHICH" == "empirical" || "$WHICH" == "rest" || "$WHICH" == "refit" ]]; then
      if [[ ! -f data/DCAD-v1.0-dyadic.csv ]]; then
        echo "WARNING: data/DCAD-v1.0-dyadic.csv not found; empirical/05 will stop at the DCA network."
      fi
      b=$(sbatch --parsable --export=ALL "$SB/05_emp_build.sbatch");                          echo "05 emp build         $b"
      f=$(sbatch --parsable --export=ALL --dependency=afterok:$b "$SB/06_emp_fit.sbatch");    echo "06 emp fit  (after 05) $f"
      s=$(sbatch --parsable --export=ALL --dependency=afterok:$f "$SB/07_emp_score.sbatch");  echo "07 emp score (after 06) $s"
      DEPS+=("$s")
    fi
    if [[ "$WHICH" == "alt" ]]; then
      a=$(sbatch --parsable --export=ALL "$SB/09_alt_bootstrap.sbatch");                        echo "09 alt bootstrap (594 tasks)  $a"
      p=$(sbatch --parsable --export=ALL --dependency=afterany:$a "$SB/10_alt_post.sbatch");    echo "10 alt post (after 09)        $p"
      echo "Results: replication/slurm/logs/10_alt_post_${p}.out and output/alternatives/"
    fi
    if [[ "$WHICH" == "coupling" ]]; then
      a=$(sbatch --parsable --export=ALL "$SB/13_coupling.sbatch");                             echo "13 coupling (48 tasks)        $a"
      p=$(sbatch --parsable --export=ALL --dependency=afterany:$a "$SB/14_coupling_post.sbatch"); echo "14 coupling post (after 13)   $p"
      echo "Results: replication/slurm/logs/14_coupling_post_${p}.out, manuscript/tables/tab_app_coupling_*.tex"
    fi
    if [[ "$WHICH" == "stab" ]]; then
      a=$(sbatch --parsable --export=ALL "$SB/11_stability.sbatch");                            echo "11 stability (594 tasks)      $a"
      p=$(sbatch --parsable --export=ALL --dependency=afterany:$a "$SB/12_stability_post.sbatch"); echo "12 stability post (after 11)  $p"
      echo "Decision: replication/slurm/logs/12_stab_post_${p}.out"
    fi
    if [[ "$WHICH" == "all" || "$WHICH" == "rest" || "$WHICH" == "refit" ]]; then
      dep=$(IFS=:; echo "${DEPS[*]}")
      p=$(sbatch --parsable --export=ALL --dependency=afterok:$dep "$SB/08_postprocess.sbatch")
      echo "08 postprocess (after all of the above) $p"
      echo
      echo "Post-processing runs automatically once every job above ends OK."
      echo "If any array task fails, rerun 'bash replication/run_all.sh submit' --"
      echo "finished tasks skip themselves -- or run 'bash replication/run_all.sh post' by hand."
    fi
    ;;

  post)
    load_r
    if [[ "$WHICH" == "alt" ]]; then Rscript replication/post/15_alternatives.R; exit 0; fi
    if [[ "$WHICH" == "coupling" ]]; then Rscript replication/post/17_coupling.R; exit 0; fi
    for s in 10_main_text 11_appendix_regimes 12_appendix_coverage 13_appendix_empirical 14_calibration_table; do
      echo "=================== post/$s.R ==================="
      Rscript "replication/post/$s.R"
    done
    if [[ -f output/calibration/coassign_calibration_table.csv ]]; then
      mkdir -p r_code/inst/extdata python_code/src/dynamic_multiplex/data
      cp output/calibration/coassign_calibration_table.csv r_code/inst/extdata/
      cp output/calibration/coassign_calibration_table.csv python_code/src/dynamic_multiplex/data/
      echo "calibration table copied into r_code/inst/extdata and python_code/src/dynamic_multiplex/data"
      echo "(reinstall the packages: bash replication/run_all.sh deps)"
    fi
    echo "tables:  $(ls manuscript/tables  | wc -l) files in manuscript/tables"
    echo "figures: $(ls manuscript/figures | wc -l) files in manuscript/figures"
    ;;

  test)
    load_r
    Rscript replication/tests/test_package_regressions.R
    if command -v python >/dev/null 2>&1; then python replication/tests/test_package_regressions.py || echo "(Python tests failed or dynamic_multiplex not installed)"; fi
    ;;

  status)
    squeue -u "$USER" -o "%.10i %.16j %.8T %.10M %.6D %R" | { head -1; grep dm_ || true; }
    for d in regime coverage_grid coverage_valued coverage_misspec; do
      n=$(ls output/$d 2>/dev/null | grep -c '^cov_task\|^dyn_cfg' || true)
      echo "output/$d: $n task files"
    done
    ls output/empirical 2>/dev/null || true
    ;;

  clean-tables)
    rm -f manuscript/tables/*.tex manuscript/figures/*.pdf manuscript/figures/*.png
    echo "removed generated tables and figures (simulation output untouched)"
    ;;

  *)
    sed -n '2,40p' "$0"
    ;;
esac
