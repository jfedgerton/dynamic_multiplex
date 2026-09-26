#!/usr/bin/env bash
# =============================================================================
# replication/run_local.sh  --  run the full replication on one machine (no SLURM)
#
# Usage (from a fresh clone, so the rerun cannot touch your working copy):
#   JOBS=16 REF_TABLES=/path/to/overleaf/tables bash replication/run_local.sh all
#
# Stages (run one at a time, or `all` for check -> smoke -> packages -> sims ->
# empirical -> post -> compare):
#   check      R version, R and Python packages, input data, cores, disk
#   smoke      one mini task of every simulation in a scratch folder (2 reps);
#              catches setup errors in minutes without touching output/
#   packages   R CMD check --as-cran, R and Python unit tests, regression tests,
#              R vs Python multislice cross-check
#   sims       sim/01, 03, 04, 05, 06 (314 tasks) in parallel with JOBS workers
#   empirical  empirical/07 -> 08 (4 networks in parallel) -> 09 -> 11
#   post       post/10, 11, 13, 14, 15, 16, 17, 18 (tables and figures)
#   compare    regenerated manuscript/tables vs REF_TABLES (the tables in the paper)
#   status     print local_status.txt
#
# Not run (not used in the paper): sim/02 stability arms, empirical/10 network
# stability, and post/12. Run those on Roar with run_all.sh if ever needed.
#
# Every sim task skips itself when its output file exists, so rerunning a
# stage after an interruption only runs what is missing.
# Progress: replication/local_logs/local_status.txt (updated every 60 s).
# Inputs you must supply: data/DCAD-v1.0-dyadic.csv (Kinne DCAD v1.0) and the
# peacesciencer extdata (Rscript -e 'peacesciencer::download_extdata()').
# =============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export DM_ROOT="${DM_ROOT:-$(dirname "$HERE")}"
cd "$DM_ROOT"
STAGE="${1:-help}"
NCPU="$(getconf _NPROCESSORS_ONLN)"
JOBS="${JOBS:-$(( NCPU > 2 ? NCPU - 2 : 1 ))}"
REF_TABLES="${REF_TABLES:-}"
LOGDIR="replication/local_logs"
STATUS="$LOGDIR/local_status.txt"
mkdir -p "$LOGDIR" replication/slurm/logs output manuscript/tables manuscript/figures

# one thread per R process; parallelism comes from JOBS
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

# simulation tasks: script | number of tasks | output folder | output prefix
SIMS="01_regime_comparison.R 72 regime dyn_cfg
03_mechanism_tests.R 50 mechanism mech_cfg
04_coupling_regimes.R 48 coupling coup_cfg
05_omega_sweep.R 72 omega omega_cfg
06_selection_rule.R 72 selection sel_cfg"

say() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOGDIR/run_local.log"; }

write_status() {
  {
    echo "updated: $(date '+%Y-%m-%d %H:%M:%S')   JOBS=$JOBS   DM_ROOT=$DM_ROOT"
    total_done=0; total_all=0
    while read -r script n dir prefix; do
      done_n=$(find "output/$dir" -name "${prefix}*.csv" 2>/dev/null | wc -l | tr -d ' ')
      fail_n=$(grep -L '^0$' "$LOGDIR"/sim_"${script%.R}"_*.exit 2>/dev/null | wc -l | tr -d ' ')
      printf "  %-24s %3s / %3s done   %s failed\n" "$script" "$done_n" "$n" "$fail_n"
      total_done=$(( total_done + done_n )); total_all=$(( total_all + n ))
    done <<< "$SIMS"
    echo "  sims total: $total_done / $total_all"
    if [[ -f "$LOGDIR/sims_started" && $total_done -gt 0 ]]; then
      t0=$(cat "$LOGDIR/sims_started"); now=$(date +%s); el=$(( now - t0 ))
      new=$(( total_done - $(cat "$LOGDIR/sims_done_at_start") ))
      if [[ $new -gt 0 ]]; then
        eta=$(( el * (total_all - total_done) / new ))
        printf "  elapsed %dh%02dm   rough ETA %dh%02dm (tasks differ a lot in size)\n" $((el/3600)) $((el%3600/60)) $((eta/3600)) $((eta%3600/60))
      fi
    fi
    for f in "$LOGDIR"/stage_*; do [[ -f "$f" ]] && echo "  $(basename "$f"): $(cat "$f")"; done
  } > "$STATUS.tmp" && mv "$STATUS.tmp" "$STATUS"
}

stage_mark() { echo "$2  $(date '+%Y-%m-%d %H:%M:%S')" > "$LOGDIR/stage_$1"; write_status; }

# ---------------------------------------------------------------------------
do_check() {
  stage_mark check running
  ok=1
  say "cores: $NCPU   JOBS: $JOBS   free disk: $(df -h . | tail -1 | awk '{print $4}')"
  command -v Rscript >/dev/null || { say "MISSING: Rscript"; ok=0; }
  rv=$(Rscript -e 'cat(paste(R.version$major, R.version$minor, sep="."))' 2>/dev/null)
  say "R version: $rv (paper runs used 4.5.0; a different version can change random draws)"
  Rscript -e 'p <- c("pkgload","igraph","clue","multinet","dynsbm","peacesciencer","dplyr","ggplot2","testthat","leiden")
               for (x in p) cat(sprintf("  %-14s %s\n", x, if (requireNamespace(x, quietly=TRUE)) as.character(packageVersion(x)) else "MISSING"))' | tee -a "$LOGDIR/run_local.log"
  Rscript -e 'p <- c("pkgload","igraph","clue","multinet","dynsbm","peacesciencer","dplyr","ggplot2")
               m <- p[!vapply(p, requireNamespace, logical(1), quietly=TRUE)]; if (length(m)) quit(status=1)' || { say "MISSING R packages (see list above; install before running)"; ok=0; }
  if [[ -f data/DCAD-v1.0-dyadic.csv ]]; then say "data/DCAD-v1.0-dyadic.csv: ok"; else say "MISSING: data/DCAD-v1.0-dyadic.csv (copy it from Roar)"; ok=0; fi
  if command -v python3 >/dev/null; then say "python3: $(python3 --version 2>&1)"; else say "python3 not found (needed for the packages stage only)"; fi
  if [[ $ok -eq 1 ]]; then stage_mark check ok; else stage_mark check FAILED; say "check FAILED: fix the items above"; return 1; fi
}

# ---------------------------------------------------------------------------
do_smoke() {
  stage_mark smoke running
  SMOKE="$(mktemp -d "${TMPDIR:-/tmp}/dm_smoke.XXXX")"
  ln -s "$DM_ROOT/r_code" "$SMOKE/r_code"; ln -s "$DM_ROOT/replication" "$SMOKE/replication"
  [[ -d "$DM_ROOT/data" ]] && ln -s "$DM_ROOT/data" "$SMOKE/data"
  fails=0
  while read -r script n dir prefix; do
    t0=$(date +%s)
    if DM_ROOT="$SMOKE" CMP_MINI=1 CMP_CFG=1 Rscript "replication/sim/$script" < /dev/null > "$LOGDIR/smoke_${script%.R}.log" 2>&1; then
      say "smoke $script ok ($(( $(date +%s) - t0 )) s, 2 reps of task 1)"
    else
      say "smoke $script FAILED (see $LOGDIR/smoke_${script%.R}.log)"; fails=$((fails+1))
    fi
  done <<< "$SIMS"
  rm -rf "$SMOKE"
  if [[ $fails -eq 0 ]]; then stage_mark smoke ok; else stage_mark smoke "FAILED ($fails)"; return 1; fi
}

# ---------------------------------------------------------------------------
do_packages() {
  stage_mark packages running
  P="$LOGDIR/packages"; mkdir -p "$P"; res=""
  run() { local name="$1"; shift; if "$@" > "$P/$name.log" 2>&1; then res+="  $name: ok\n"; else res+="  $name: FAILED (see $P/$name.log)\n"; fi; }

  # R: build, CRAN check, unit tests, regression tests
  rm -f dynamicmultiplex_*.tar.gz
  run r_build R CMD build --no-build-vignettes r_code
  TARBALL=$(ls dynamicmultiplex_*.tar.gz 2>/dev/null | head -1)
  if [[ -n "$TARBALL" ]]; then
    run r_cmd_check env _R_CHECK_FORCE_SUGGESTS_=false R CMD check --as-cran --no-manual -o "$P" "$TARBALL"
  fi
  run r_testthat Rscript -e 'pkgload::load_all("r_code", quiet=TRUE); testthat::test_dir("r_code/tests/testthat", stop_on_failure=TRUE)'
  run r_regressions Rscript replication/tests/test_package_regressions.R

  # Python: fresh virtual environment, editable install, tests
  if command -v python3 >/dev/null; then
    VENV="$DM_ROOT/.venv_dm"
    [[ -d "$VENV" ]] || python3 -m venv "$VENV"
    run py_install "$VENV/bin/pip" install --quiet -e python_code pytest
    run py_pytest "$VENV/bin/python" -m pytest -q python_code/tests
    run py_regressions "$VENV/bin/python" replication/tests/test_package_regressions.py
    run crosscheck_r_vs_python env PATH="$VENV/bin:$PATH" Rscript replication/tests/crosscheck_multislice_r_vs_python.R
  else
    res+="  python: skipped (python3 not found)\n"
  fi
  printf "packages stage results:\n$res" | tee "$P/SUMMARY.txt" | tee -a "$LOGDIR/run_local.log"
  if grep -q FAILED "$P/SUMMARY.txt"; then stage_mark packages "FAILED (see $P/SUMMARY.txt)"; else stage_mark packages ok; fi
}

# ---------------------------------------------------------------------------
run_one_sim() {   # called by xargs: $1 = script, $2 = task
  local script="$1" task="$2" tag="${1%.R}_$(printf %02d "$2")"
  SLURM_ARRAY_TASK_ID="$task" Rscript "replication/sim/$script" > "$LOGDIR/sim_${tag}.log" 2>&1
  echo $? > "$LOGDIR/sim_${tag}.exit"
}
export -f run_one_sim
export LOGDIR

do_sims() {
  stage_mark sims running
  date +%s > "$LOGDIR/sims_started"
  n0=0; while read -r script n dir prefix; do n0=$(( n0 + $(find "output/$dir" -name "${prefix}*.csv" 2>/dev/null | wc -l) )); done <<< "$SIMS"
  echo "$n0" > "$LOGDIR/sims_done_at_start"
  ( while true; do write_status; sleep 60; done ) & MON=$!
  # one queue across all five scripts. Longest tasks first: on Roar the selection
  # tasks took up to 7.3 h each and regime tasks up to 3.1 h, so starting them
  # early keeps the last few workers from finishing alone.
  QUEUE="$LOGDIR/sim_queue.txt"; : > "$QUEUE"
  for s in 06_selection_rule.R 01_regime_comparison.R 05_omega_sweep.R 03_mechanism_tests.R 04_coupling_regimes.R; do
    n=$(echo "$SIMS" | awk -v s="$s" '$1==s {print $2}')
    for i in $(seq 1 "$n"); do echo "$s $i" >> "$QUEUE"; done
  done
  say "sims: $(wc -l < "$QUEUE" | tr -d ' ') tasks queued, $JOBS at a time"
  xargs -P "$JOBS" -L 1 bash -c 'run_one_sim "$0" "$1"' < "$QUEUE"
  kill "$MON" 2>/dev/null; write_status
  nfail=$(grep -L '^0$' "$LOGDIR"/sim_*.exit 2>/dev/null | wc -l | tr -d ' ')
  if [[ "$nfail" -eq 0 ]]; then stage_mark sims ok; else stage_mark sims "FAILED ($nfail tasks; rerun the stage to retry only those)"; return 1; fi
}

# ---------------------------------------------------------------------------
do_empirical() {
  stage_mark empirical running
  say "empirical/07 build networks"
  Rscript replication/empirical/07_build_networks.R > "$LOGDIR/emp_07.log" 2>&1 || { stage_mark empirical "FAILED at 07"; return 1; }
  say "empirical/08 fit four networks in parallel"
  # log names match what post/17_runtime.R reads for the empirical timings
  i=0; pids=""
  for net in atop dca igo trade; do
    i=$((i+1)); Rscript replication/empirical/08_fit_networks.R "$net" > "replication/slurm/logs/08_emp_fit_0_${i}.out" 2>&1 & pids="$pids $!"
  done
  f=0; for p in $pids; do wait "$p" || f=$((f+1)); done
  [[ $f -eq 0 ]] || { stage_mark empirical "FAILED at 08 ($f networks)"; return 1; }
  say "empirical/09 score orders"
  Rscript replication/empirical/09_score_orders.R > "$LOGDIR/emp_09.log" 2>&1 || { stage_mark empirical "FAILED at 09"; return 1; }
  say "empirical/11 aligned and rival dyads"
  Rscript replication/empirical/11_alignment_dyads.R > "$LOGDIR/emp_11.log" 2>&1 || { stage_mark empirical "FAILED at 11"; return 1; }
  stage_mark empirical ok
}

# ---------------------------------------------------------------------------
do_post() {
  stage_mark post running
  f=0
  for s in 10_main_text 11_appendix_regimes 13_appendix_empirical 14_coupling 15_omega_selection 16_mechanism_decision 17_runtime 18_mechanism_paired; do
    if Rscript "replication/post/$s.R" > "$LOGDIR/post_$s.log" 2>&1; then say "post/$s ok"; else say "post/$s FAILED"; f=$((f+1)); fi
  done
  if [[ $f -eq 0 ]]; then stage_mark post ok; else stage_mark post "FAILED ($f scripts)"; return 1; fi
}

# ---------------------------------------------------------------------------
do_compare() {
  if [[ -z "$REF_TABLES" ]]; then say "compare: set REF_TABLES to the tables/ folder from the Overleaf export"; return 1; fi
  stage_mark compare running
  Rscript replication/compare_tables.R manuscript/tables "$REF_TABLES" "$LOGDIR/compare_report.md" | tee -a "$LOGDIR/run_local.log"
  stage_mark compare "done (see $LOGDIR/compare_report.md)"
}

case "$STAGE" in
  check)     do_check ;;
  smoke)     do_smoke ;;
  packages)  do_packages ;;
  sims)      do_sims ;;
  empirical) do_empirical ;;
  post)      do_post ;;
  compare)   do_compare ;;
  status)    write_status; cat "$STATUS" ;;
  all)
    do_check && do_smoke || exit 1
    do_packages            # failures here are reported but do not stop the run
    do_sims && do_empirical && do_post && do_compare
    write_status; cat "$STATUS" ;;
  *) sed -n '2,30p' "$0" ;;
esac
