# =============================================================================
# replication/sim/01_regime_comparison.R  --  Table 2 and the paired-CI appendix
#
# Every method on the SAME simulated series, four dynamic-community regimes:
#   seasonality  a bank of latent partitions recurs on a period (period-lagged
#                interlayer links supplied by the generator)
#   churnswitch  edge churn plus Markov switching of node memberships
#   birthdeath   open population: whole communities are born and die
#   regimeshift  change-point: era-1 communities replaced by disjoint era-2 ones
# Grid: regime x N{50,100,200} x r{1.5,3,6} x intensity{low,high} = 72 configs,
# 30 replicates each, shuffled with set.seed(123) so a partial run samples every
# condition. Per-rep seed = 9000 + TASK*1000 + rep.
#
# Methods (sim/lib_regimes.R, METHODS): DynMux Jaccard / Overlap / weighted
# variants (Leiden and Louvain), multislice identity ties with adjacent links
# and with the generator's links, pooled Leiden/Louvain, cross-sectional +
# Hungarian, multinet GLouvain and dynsbm (both skipped when not installed).
# Metrics: per-layer and joint NMI, transition-layer NMI, K MAE, co-membership
# accuracy, runtime, community counts (active nodes only in open populations).
#
# Usage (smoke): CMP_MINI=1 CMP_CFG=1 DM_ROOT=. Rscript replication/sim/01_regime_comparison.R
# Array: one SLURM_ARRAY_TASK_ID per config (1..72); slurm/01_regime.sbatch.
# Output: $DM_ROOT/output/regime/dyn_cfg%02d.csv (a task with an existing file
# exits at once, so the array can be resubmitted after a partial failure).
# =============================================================================
ROOT <- Sys.getenv("DM_ROOT", unset = getwd())
suppressMessages(pkgload::load_all(file.path(ROOT, "r_code"), quiet = TRUE))
suppressPackageStartupMessages(library(igraph))
TASK <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", Sys.getenv("CMP_CFG", "1")))
MINI <- identical(Sys.getenv("CMP_MINI", "0"), "1")
REPS_FAST   <- if (MINI) 2L else 30L
REPS_DYNSBM <- if (MINI) 2L else 30L
QMIN <- 2L; QMAX <- 8L; NSTART <- 3L; RHO <- 0.10
source(file.path(ROOT, "replication", "sim", "lib_regimes.R"))

outdir <- file.path(ROOT, "output", "regime"); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
cfgs <- regime_grid()
stopifnot(TASK >= 1L, TASK <= nrow(cfgs))
cfg <- cfgs[TASK, ]
outfile <- file.path(outdir, sprintf("dyn_cfg%02d.csv", TASK))
if (nzchar(Sys.getenv("SLURM_ARRAY_TASK_ID")) && file.exists(outfile)) {
  cat("[dyn skip] task", TASK, "already complete:", outfile, "\n"); quit(save = "no")
}
if (!HAVE_MULTINET) message("[note] multinet not installed; 'multinet GLouvain' recorded as NA.")
if (!HAVE_DYNSBM)   message("[note] dynsbm not installed; 'Dynamic SBM' skipped.")
cat(sprintf("[dyn] task=%d/%d -> regime=%s n=%d r=%.1f intensity=%s | reps=%d | clue=%s dynsbm=%s multinet=%s\n",
            TASK, nrow(cfgs), cfg$regime, cfg$n, cfg$r, cfg$intensity, REPS_FAST, HAVE_CLUE, HAVE_DYNSBM, HAVE_MULTINET))

# ---- run: one rep at a time, every method on the same simulated series ----
metric_cols <- c("nmi_layer", "nmi_joint", "nmi_change", "k_mae", "comembership_acc", "mean_n_comm", "total_n_comm")
run_rep <- function(rep) {
  seed <- 9000L + TASK * 1000L + rep
  sim  <- simulate_regime(cfg, seed)
  names_this <- names(METHODS)
  if (rep > REPS_DYNSBM || !HAVE_DYNSBM) names_this <- setdiff(names_this, "Dynamic SBM")
  rows <- vector("list", length(names_this)); ri <- 0L
  for (mname in names_this) {
    t0  <- proc.time()[["elapsed"]]
    det <- tryCatch(METHODS[[mname]](sim), error = function(e) {
      message(sprintf("  [rep %d] method '%s' errored: %s", rep, mname, conditionMessage(e))); NULL })
    el  <- proc.time()[["elapsed"]] - t0
    m   <- if (is.null(det)) setNames(rep(NA_real_, length(metric_cols)), metric_cols) else eval_method(det, sim)
    ri  <- ri + 1L
    rows[[ri]] <- data.frame(regime = cfg$regime, n = cfg$n, r = cfg$r, intensity = cfg$intensity, rep = rep, method = mname,
      nmi_layer = round(m[["nmi_layer"]], 4), nmi_joint = round(m[["nmi_joint"]], 4), nmi_change = round(m[["nmi_change"]], 4),
      k_mae = round(m[["k_mae"]], 4), comembership_acc = round(m[["comembership_acc"]], 4), runtime_s = round(el, 3),
      mean_n_comm = round(m[["mean_n_comm"]], 4), total_n_comm = m[["total_n_comm"]], stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}
t0 <- proc.time()[["elapsed"]]; all_rows <- vector("list", REPS_FAST)
for (rep in seq_len(REPS_FAST)) {
  all_rows[[rep]] <- tryCatch(run_rep(rep), error = function(e) { message(sprintf("rep %d failed: %s", rep, conditionMessage(e))); NULL })
  if (rep %% 10L == 0L || MINI) cat(sprintf("  ... rep %d/%d done (%.1f min)\n", rep, REPS_FAST, (proc.time()[["elapsed"]] - t0) / 60))
}
res <- do.call(rbind, Filter(Negate(is.null), all_rows))
stopifnot(nrow(res) > 0, all(metric_cols %in% names(res)))
write.csv(res, outfile, row.names = FALSE)
agg <- aggregate(cbind(nmi_layer, nmi_joint, comembership_acc) ~ method, data = res, FUN = function(x) mean(x, na.rm = TRUE))
cat(sprintf("[dyn] task %d done: %d rows in %.1f min -> %s\n", TASK, nrow(res), (proc.time()[["elapsed"]] - t0) / 60, outfile))
print(agg[order(-agg$nmi_joint), ], row.names = FALSE)
