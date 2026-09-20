# =============================================================================
# replication/sim/05_omega_sweep.R  --  multislice coupling sensitivity (appendix)
#
# Mucha multislice (identity ties, correct per-slice null model) on the 72-cell
# regime grid of sim/01 for omega in {0.25, 0.5, 1, 2, 4}, with adjacent links
# and with the generator's links (period-lagged for recurring structure), plus
# DynMux Jaccard on the same series for reference. 10 replicates per cell,
# seeds 9000 + TASK*1000 + rep (the first 10 sim/01 replicates, so the DynMux
# reference reproduces Table 2's first ten reps).
# Usage (mini): CMP_MINI=1 CMP_CFG=1 DM_ROOT=. Rscript replication/sim/05_omega_sweep.R
# Array: 1..72 (slurm/05_omega.sbatch). Output: output/omega/omega_cfg%02d.csv
# =============================================================================
ROOT <- Sys.getenv("DM_ROOT", unset = getwd())
suppressMessages(pkgload::load_all(file.path(ROOT, "r_code"), quiet = TRUE))
suppressPackageStartupMessages(library(igraph))
TASK <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", Sys.getenv("CMP_CFG", "1")))
MINI <- identical(Sys.getenv("CMP_MINI", "0"), "1")
REPS <- if (MINI) 2L else 10L
OMEGAS <- c(0.25, 0.5, 1, 2, 4)
HAVE_DYNSBM <- FALSE; HAVE_MULTINET <- FALSE                 # not needed here
source(file.path(ROOT, "replication", "sim", "lib_regimes.R"))

outdir <- file.path(ROOT, "output", "omega"); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
cfgs <- regime_grid(); stopifnot(TASK >= 1L, TASK <= nrow(cfgs)); cfg <- cfgs[TASK, ]
outfile <- file.path(outdir, sprintf("omega_cfg%02d.csv", TASK))
if (nzchar(Sys.getenv("SLURM_ARRAY_TASK_ID")) && file.exists(outfile)) { cat("[skip] task", TASK, "\n"); quit(save = "no") }
cat(sprintf("[omega] task=%d regime=%s n=%d r=%.1f intensity=%s reps=%d\n", TASK, cfg$regime, cfg$n, cfg$r, cfg$intensity, REPS))

ms <- function(omega, links) function(sim) {
  n <- nrow(sim$layers[[1]])
  f <- if (links) fit_multilayer_identity_ties(AL(sim), algorithm = "leiden", layer_links = sim$links, omega = omega, allow_unequal_nodes = TRUE)
       else       fit_multilayer_identity_ties(AL(sim), algorithm = "leiden", omega = omega, allow_unequal_nodes = TRUE)
  EXPn(extract_meta_membership(f), sim, n)
}
VARIANTS <- list("DynMux Jaccard" = METHODS[["DynMux Jaccard"]])
for (om in OMEGAS) { VARIANTS[[sprintf("Multislice adjacent omega=%g", om)]] <- ms(om, FALSE)
                     VARIANTS[[sprintf("Multislice custom links omega=%g", om)]] <- ms(om, TRUE) }

rows <- list(); t0 <- proc.time()[["elapsed"]]
for (rep in seq_len(REPS)) {
  sim <- simulate_regime(cfg, 9000L + TASK * 1000L + rep)
  for (v in names(VARIANTS)) {
    t1 <- proc.time()[["elapsed"]]
    det <- tryCatch(VARIANTS[[v]](sim), error = function(e) { message(sprintf("rep %d %s: %s", rep, v, conditionMessage(e))); NULL })
    if (is.null(det)) next
    m <- eval_method(det, sim)
    rows[[length(rows) + 1]] <- data.frame(regime = cfg$regime, n = cfg$n, r = cfg$r, intensity = cfg$intensity, rep = rep, variant = v,
      nmi_joint = m[["nmi_joint"]], nmi_layer = m[["nmi_layer"]], k_mae = m[["k_mae"]], total_n_comm = m[["total_n_comm"]],
      runtime_s = proc.time()[["elapsed"]] - t1, stringsAsFactors = FALSE)
  }
  cat(sprintf("  rep %d/%d (%.1f min)\n", rep, REPS, (proc.time()[["elapsed"]] - t0) / 60))
}
out <- do.call(rbind, rows); stopifnot(nrow(out) >= length(VARIANTS))
write.csv(out, outfile, row.names = FALSE)
print(aggregate(cbind(nmi_joint, k_mae) ~ variant, out, function(x) round(mean(x), 3)))
