# =============================================================================
# replication/sim/06_selection_rule.R  --  stability as a model-selection rule
#
# On the 72-cell regime grid of sim/01, fit DynMux (Jaccard) and multislice
# (identity ties, the generator's links) on the same series, bootstrap each
# (B = 50 block-density redraws from its own fitted partition, refit with the
# same method) and record stability (mean NMI replicate vs point) and accuracy
# (joint NMI vs truth) for both. post/15 reports how often the method with the
# higher stability is also the more accurate one, by regime, and the paired
# accuracy gain from following the rule. 10 replicates per cell, seeds
# 9000 + TASK*1000 + rep (same series as sim/05 and the first ten of sim/01).
# Usage (mini): CMP_MINI=1 CMP_CFG=1 DM_ROOT=. Rscript replication/sim/06_selection_rule.R
# Array: 1..72 (slurm/06_selection.sbatch). Output: output/selection/sel_cfg%02d.csv
# =============================================================================
ROOT <- Sys.getenv("DM_ROOT", unset = getwd())
suppressMessages(pkgload::load_all(file.path(ROOT, "r_code"), quiet = TRUE))
suppressPackageStartupMessages(library(igraph))
TASK <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", Sys.getenv("CMP_CFG", "1")))
MINI <- identical(Sys.getenv("CMP_MINI", "0"), "1")
REPS <- if (MINI) 2L else 10L
B_BOOT <- if (MINI) 4L else 50L
HAVE_DYNSBM <- FALSE; HAVE_MULTINET <- FALSE
source(file.path(ROOT, "replication", "sim", "lib_regimes.R"))

outdir <- file.path(ROOT, "output", "selection"); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
cfgs <- regime_grid(); stopifnot(TASK >= 1L, TASK <= nrow(cfgs)); cfg <- cfgs[TASK, ]
outfile <- file.path(outdir, sprintf("sel_cfg%02d.csv", TASK))
if (nzchar(Sys.getenv("SLURM_ARRAY_TASK_ID")) && file.exists(outfile)) { cat("[skip] task", TASK, "\n"); quit(save = "no") }
cat(sprintf("[selection] task=%d regime=%s n=%d r=%.1f intensity=%s reps=%d B=%d\n", TASK, cfg$regime, cfg$n, cfg$r, cfg$intensity, REPS, B_BOOT))

FITS <- list("DynMux Jaccard" = function(L, links) fit_multilayer_jaccard(L, algorithm = "leiden", layer_links = links, allow_unequal_nodes = TRUE),
             "Multislice (custom links)" = function(L, links) fit_multilayer_identity_ties(L, algorithm = "leiden", layer_links = links, allow_unequal_nodes = TRUE))

# block-density redraw on the alive nodes of each layer (dimnames kept) -------
redraw_layers <- function(AL_layers, mems) lapply(seq_along(AL_layers), function(t) {
  A <- AL_layers[[t]]; mem <- as.integer(mems[[t]]); n <- nrow(A)
  if (n < 2) return(A)
  same <- outer(mem, mem, "=="); up <- upper.tri(A); E <- A > 0
  p_all <- mean(E[up]); in_d <- up & same; out_d <- up & !same
  p_in <- if (any(in_d)) mean(E[in_d]) else p_all; p_out <- if (any(out_d)) mean(E[out_d]) else p_all
  M <- matrix(0, n, n, dimnames = dimnames(A)); M[up] <- rbinom(sum(up), 1, ifelse(same[up], p_in, p_out)); M + t(M) })

rows <- list(); t0 <- proc.time()[["elapsed"]]
for (rep in seq_len(REPS)) {
  sim <- simulate_regime(cfg, 9000L + TASK * 1000L + rep); n <- nrow(sim$layers[[1]]); L <- AL(sim)
  for (mname in names(FITS)) {
    fit0 <- tryCatch(FITS[[mname]](L, sim$links), error = function(e) NULL); if (is.null(fit0)) next
    mem0 <- extract_meta_membership(fit0)
    acc <- eval_method(EXPn(mem0, sim, n), sim)
    s <- numeric(0)
    for (b in seq_len(B_BOOT)) {
      bfit <- tryCatch(FITS[[mname]](redraw_layers(L, mem0), sim$links), error = function(e) NULL); if (is.null(bfit)) next
      bmem <- extract_meta_membership(bfit)
      s <- c(s, mean(vapply(seq_along(L), function(t) if (length(mem0[[t]]) < 2) NA_real_ else
        igraph::compare(as.integer(bmem[[t]]), as.integer(mem0[[t]]), method = "nmi"), numeric(1)), na.rm = TRUE))
    }
    rows[[length(rows) + 1]] <- data.frame(regime = cfg$regime, n = cfg$n, r = cfg$r, intensity = cfg$intensity, rep = rep, method = mname,
      stability = mean(s), stability_sd = sd(s), B = length(s), acc_joint = acc[["nmi_joint"]], acc_layer = acc[["nmi_layer"]],
      k_mae = acc[["k_mae"]], stringsAsFactors = FALSE)
  }
  cat(sprintf("  rep %d/%d (%.1f min)\n", rep, REPS, (proc.time()[["elapsed"]] - t0) / 60))
}
out <- do.call(rbind, rows); stopifnot(nrow(out) >= 2)
write.csv(out, outfile, row.names = FALSE)
print(aggregate(cbind(stability, acc_joint) ~ method, out, function(x) round(mean(x), 3)))
