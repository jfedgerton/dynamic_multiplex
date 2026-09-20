#!/usr/bin/env bash
# =============================================================================
# replication/run_all.sh  --  one entry point for the full replication
#
#   bash replication/run_all.sh deps            install the R package (+ Python) into the user library
#   bash replication/run_all.sh test            package regression tests (R + Python)
#   bash replication/run_all.sh submit          submit EVERYTHING with dependencies (front-to-end rerun)
#   bash replication/run_all.sh submit sims     sim/01, 02 (three arms), 03, 04, 05, 06 only
#   bash replication/run_all.sh submit empirical   empirical/07-10 only
#   bash replication/run_all.sh submit post     post-processing only (after the above finished)
#   bash replication/run_all.sh post            run post/10-16 in the current shell
#   bash replication/run_all.sh status          queue + output counts
#   bash replication/run_all.sh archive         move existing output/ and manuscript tables/figures
#                                               to output/_archive/<stamp>/ before a clean rerun
#   bash replication/run_all.sh clean-tables    remove generated tables/figures
#
# Pipeline (every script reads DM_ROOT; defaults to this repo's root):
#   sim/01_regime_comparison.R      -> output/regime/dyn_cfgNN.csv        72 tasks   Table 2, App. A2
#   sim/02_stability.R (3 arms)     -> output/stability/<arm>_*_taskNNNNN 594+216+72 Section 4, App. A4
#   sim/03_mechanism_tests.R        -> output/mechanism/mech_cfgNN.csv    50 tasks   decision table, App. mechanism
#   sim/04_coupling_regimes.R       -> output/coupling/coup_cfgNN.csv     48 tasks   App. A3
#   sim/05_omega_sweep.R            -> output/omega/omega_cfgNN.csv       72 tasks   App. omega sweep
#   sim/06_selection_rule.R         -> output/selection/sel_cfgNN.csv     72 tasks   App. selection rule
#   empirical/07_build_networks.R   -> output/empirical_data/<net>_{series,union}.rds
#   empirical/08_fit_networks.R     -> output/empirical/<net>_partitions.rds  4 tasks
#   empirical/09_score_orders.R     -> output/empirical/order_recovery.csv
#   empirical/10_stability_networks.R -> output/empirical/<net>_stability.csv  4 tasks
#   post/10_main_text.R             -> Table 2, Figure 3
#   post/11_appendix_regimes.R      -> paired-difference figures and CI tables (Jaccard + overlap columns)
#   post/12_stability.R             -> calibration table, Figure 2, stability tables (main + appendix)
#   post/13_appendix_empirical.R    -> order-recovery tables
#   post/14_coupling.R              -> Jaccard vs overlap appendix table + figure
#   post/15_omega_selection.R       -> omega sweep, selection rule, empirical stability tables
#   post/16_mechanism_decision.R    -> main-text decision table (tab_decision_tree), mechanism appendix tables
# Tables land in manuscript/tables/, figures in manuscript/figures/. The
# calibration table is copied into r_code/inst/extdata and
# python_code/src/dynamic_multiplex/data by `post`.
#
# Every sim task self-skips when its output exists, so `submit` can be rerun
# after a partial failure and only the missing tasks execute. Total array
# tasks: 72 + 882 + 50 + 48 + 72 + 72 + 9 = 1,205 (Roar caps submitted jobs at
# 4,000 per account).
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
WHICH="${2:-all}"
SB="replication/slurm"
mkdir -p "$SB/logs" output manuscript/tables manuscript/figures

have_module() { command -v module >/dev/null 2>&1; }
load_r() { if have_module; then module load r/4.5.0; fi; }
sb() { sbatch --parsable --export=ALL "$@"; }

case "$ACTION" in

  deps)
    load_r
    echo "--- R: installing dynamicmultiplex from r_code/ into the user library ---"
    LIB="$(Rscript -e 'cat(Sys.getenv("R_LIBS_USER"))')"; mkdir -p "$LIB"
    R CMD INSTALL -l "$LIB" r_code
    Rscript -e 'for (p in c("igraph","ggplot2","pkgload","clue","multinet","dynsbm","peacesciencer","dplyr"))
                  cat(sprintf("%-14s %s\n", p, if (requireNamespace(p, quietly=TRUE)) "ok" else "MISSING"))'
    if command -v pip >/dev/null 2>&1; then
      echo "--- Python: installing dynamic_multiplex from python_code/ ---"
      pip install -e python_code --quiet || echo "(pip install failed; the Python package is optional for the replication)"
    fi
    ;;

  test)
    load_r
    Rscript replication/tests/test_package_regressions.R
    if command -v python >/dev/null 2>&1; then python replication/tests/test_package_regressions.py || echo "(Python tests failed or dynamic_multiplex not installed)"; fi
    ;;

  archive)
    stamp=$(date +%Y%m%d_%H%M); mkdir -p "output/_archive/$stamp"
    for d in regime stability mechanism coupling omega selection empirical empirical_data; do
      [[ -d output/$d ]] && mv "output/$d" "output/_archive/$stamp/" && echo "archived output/$d"
    done
    mkdir -p "output/_archive/$stamp/manuscript"
    for d in tables figures; do [[ -d manuscript/$d ]] && cp -r "manuscript/$d" "output/_archive/$stamp/manuscript/" ; done
    echo "archive: output/_archive/$stamp (empirical_data archived too; 07 rebuilds it)"
    ;;

  submit)
    command -v sbatch >/dev/null 2>&1 || { echo "sbatch not found: run this on the cluster login node"; exit 1; }
    DEPS=()
    if [[ "$WHICH" == "all" || "$WHICH" == "sims" ]]; then
      j=$(sb "$SB/01_regime.sbatch");             echo "01 regime            (72)   $j"; DEPS+=("$j")
      j=$(sb "$SB/02a_stability_binary.sbatch");  echo "02a stability binary (594)  $j"; DEPS+=("$j")
      j=$(sb "$SB/02b_stability_dcsbm.sbatch");   echo "02b stability dcsbm  (216)  $j"; DEPS+=("$j")
      j=$(sb "$SB/02c_stability_weighted.sbatch");echo "02c stability weighted (72) $j"; DEPS+=("$j")
      j=$(sb "$SB/03_mechanism.sbatch");          echo "03 mechanism tests   (50)   $j"; DEPS+=("$j")
      j=$(sb "$SB/04_coupling.sbatch");           echo "04 coupling          (48)   $j"; DEPS+=("$j")
      j=$(sb "$SB/05_omega.sbatch");              echo "05 omega sweep       (72)   $j"; DEPS+=("$j")
      j=$(sb "$SB/06_selection.sbatch");          echo "06 selection rule    (72)   $j"; DEPS+=("$j")
    fi
    if [[ "$WHICH" == "all" || "$WHICH" == "empirical" ]]; then
      [[ -f data/DCAD-v1.0-dyadic.csv ]] || echo "WARNING: data/DCAD-v1.0-dyadic.csv not found; empirical/07 will stop at the DCA network."
      b=$(sb "$SB/07_emp_build.sbatch");                               echo "07 emp build              $b"
      f=$(sb --dependency=afterok:$b "$SB/08_emp_fit.sbatch");         echo "08 emp fit   (after 07)   $f"
      s=$(sb --dependency=afterok:$f "$SB/09_emp_score.sbatch");       echo "09 emp score (after 08)   $s"
      t=$(sb --dependency=afterok:$b "$SB/10_emp_stability.sbatch");   echo "10 emp stability (after 07) $t"
      DEPS+=("$s" "$t")
    fi
    if [[ "$WHICH" == "all" || "$WHICH" == "post" ]]; then
      if [[ ${#DEPS[@]} -gt 0 ]]; then dep=$(IFS=:; echo "${DEPS[*]}"); p=$(sb --dependency=afterok:$dep "$SB/11_postprocess.sbatch")
      else p=$(sb "$SB/11_postprocess.sbatch"); fi
      echo "11 postprocess            $p"
      echo "Post-processing runs once every job above ends OK; on a partial failure rerun 'submit' (finished tasks skip) or 'post' by hand."
    fi
    ;;

  post)
    load_r
    for s in 10_main_text 11_appendix_regimes 12_stability 13_appendix_empirical 14_coupling 15_omega_selection 16_mechanism_decision; do
      echo "=================== post/$s.R ==================="
      Rscript "replication/post/$s.R"
    done
    if [[ -f output/stability/stability_calibration_table.csv ]]; then
      mkdir -p r_code/inst/extdata python_code/src/dynamic_multiplex/data
      cp output/stability/stability_calibration_table.csv r_code/inst/extdata/
      cp output/stability/stability_calibration_table.csv python_code/src/dynamic_multiplex/data/
      echo "calibration table copied into r_code/inst/extdata and python_code/src/dynamic_multiplex/data (reinstall: run_all.sh deps)"
    fi
    echo "tables:  $(ls manuscript/tables  | wc -l) files in manuscript/tables"
    echo "figures: $(ls manuscript/figures | wc -l) files in manuscript/figures"
    ;;

  status)
    squeue -u "$USER" -h -r -o "%j %T" | sort | uniq -c || true
    echo "regime     $(ls output/regime 2>/dev/null | grep -c dyn_cfg)/72"
    echo "stability  binary $(ls output/stability 2>/dev/null | grep -c '^binary_stab')/594  dcsbm $(ls output/stability 2>/dev/null | grep -c '^dcsbm_stab')/216  weighted $(ls output/stability 2>/dev/null | grep -c '^weighted_stab')/72"
    echo "mechanism  $(ls output/mechanism 2>/dev/null | grep -c mech_cfg)/50"
    echo "coupling   $(ls output/coupling 2>/dev/null | grep -c coup_cfg)/48   omega $(ls output/omega 2>/dev/null | grep -c omega_cfg)/72   selection $(ls output/selection 2>/dev/null | grep -c sel_cfg)/72"
    ls output/empirical 2>/dev/null || true
    ;;

  clean-tables)
    rm -f manuscript/tables/*.tex manuscript/figures/*.pdf manuscript/figures/*.png
    echo "removed generated tables and figures (simulation output untouched)"
    ;;

  *)
    sed -n '2,45p' "$0"
    ;;
esac
