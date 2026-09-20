# =============================================================================
# replication/empirical/10_stability_networks.R  --  stability scores for the
# four empirical networks (Section 4, appendix table tab_app_empirical_stability)
#
# For one network (atop | dca | igo | trade) and each of DynMux (Jaccard) and
# multislice (identity ties, adjacent links): fit on the native per-year graphs
# (non-present states dropped, allow_unequal_nodes = TRUE, as in 08), then
# B bootstrap redraws from the block densities of the fitted meta-partition on
# each year's present states, refit, and record
#   stability   mean NMI(replicate, point) over replicates and years
#   floor       calibrated 5th-percentile accuracy (partition_stability lookup)
#   node        share of present node-years with Jaccard stability >= 0.9
#   pairs       decided-together / decided-apart / undetermined shares
# The bootstrap is the same block-density redraw as bootstrap_multilayer();
# it is inlined here because the package function assumes a fixed node set.
#
# Usage: DM_ROOT=... Rscript replication/empirical/10_stability_networks.R <net>
# Env: STAB_B (default 100), FORCE=1 to recompute.
# Input:  output/empirical_data/<net>_{series,union}.rds
# Output: output/empirical/<net>_stability.csv (one row per method) and
#         output/empirical/<net>_stability_years.csv (per-year stability)
# =============================================================================
ROOT <- Sys.getenv("DM_ROOT", unset = getwd())
suppressMessages(pkgload::load_all(file.path(ROOT, "r_code"), quiet = TRUE))
suppressPackageStartupMessages(library(igraph))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript 10_stability_networks.R <atop|dca|igo|trade>", call. = FALSE)
net <- args[1]; stopifnot(net %in% c("atop", "dca", "igo", "trade"))
B_BOOT <- as.integer(Sys.getenv("STAB_B", "100"))
EMP_DATA <- file.path(ROOT, "output", "empirical_data"); EMP_OUT <- file.path(ROOT, "output", "empirical")
outf <- file.path(EMP_OUT, sprintf("%s_stability.csv", net)); outy <- file.path(EMP_OUT, sprintf("%s_stability_years.csv", net))
if (file.exists(outf) && !identical(Sys.getenv("FORCE", "0"), "1")) { cat("[skip]", outf, "exists\n"); quit(save = "no") }
set.seed(123)

S <- readRDS(file.path(EMP_DATA, sprintf("%s_series.rds", net))); U <- readRDS(file.path(EMP_DATA, sprintf("%s_union.rds", net)))
if (is.null(U$present)) U$present <- U$active
yrs <- S$years; T_ <- length(yrs)
GLc <- lapply(seq_len(T_), function(k) igraph::delete_vertices(S$graph_layers[[k]], which(!U$present[[k]])))
AL <- lapply(GLc, function(g) { A <- as.matrix(igraph::as_adjacency_matrix(g, sparse = FALSE)); A[A > 0] <- 1; A })
stopifnot(all(sapply(AL, function(A) !is.null(rownames(A)))))
cat(sprintf("[stab-emp] %s: %d years, %d-%d present states\n", net, T_, min(sapply(AL, nrow)), max(sapply(AL, nrow))))

FITS <- list("DynMux Jaccard" = function(L) fit_multilayer_jaccard(L, algorithm = "leiden", allow_unequal_nodes = TRUE),
             "Multislice adjacent" = function(L) fit_multilayer_identity_ties(L, algorithm = "leiden", allow_unequal_nodes = TRUE))

redraw_layers <- function(L, mems) lapply(seq_along(L), function(t) {
  A <- L[[t]]; mem <- as.integer(mems[[t]]); n <- nrow(A); if (n < 2) return(A)
  same <- outer(mem, mem, "=="); up <- upper.tri(A); E <- A > 0
  p_all <- mean(E[up]); in_d <- up & same; out_d <- up & !same
  p_in <- if (any(in_d)) mean(E[in_d]) else p_all; p_out <- if (any(out_d)) mean(E[out_d]) else p_all
  M <- matrix(0, n, n, dimnames = dimnames(A)); M[up] <- rbinom(sum(up), 1, ifelse(same[up], p_in, p_out)); M + t(M) })
node_jaccard <- function(a, b) { sa <- outer(a, a, "=="); sb <- outer(b, b, "=="); rowSums(sa & sb) / rowSums(sa | sb) }
tab <- dynamicmultiplex:::.load_stability_table(NULL)

rows <- list(); yrows <- list()
for (mname in names(FITS)) {
  t0 <- proc.time()[["elapsed"]]
  fit0 <- FITS[[mname]](AL); mem0 <- lapply(extract_meta_membership(fit0), as.integer)
  S_rep <- matrix(NA_real_, 0, T_); node_acc <- lapply(seq_len(T_), function(t) numeric(nrow(AL[[t]])))
  co <- lapply(seq_len(T_), function(t) matrix(0, nrow(AL[[t]]), nrow(AL[[t]]))); b_ok <- 0L
  for (b in seq_len(B_BOOT)) {
    bfit <- tryCatch(FITS[[mname]](redraw_layers(AL, mem0)), error = function(e) NULL); if (is.null(bfit)) next
    bmem <- lapply(extract_meta_membership(bfit), as.integer); b_ok <- b_ok + 1L
    S_rep <- rbind(S_rep, vapply(seq_len(T_), function(t) if (length(mem0[[t]]) < 2) NA_real_ else igraph::compare(bmem[[t]], mem0[[t]], method = "nmi"), numeric(1)))
    for (t in seq_len(T_)) if (length(mem0[[t]]) >= 2) { node_acc[[t]] <- node_acc[[t]] + node_jaccard(bmem[[t]], mem0[[t]]); co[[t]] <- co[[t]] + outer(bmem[[t]], bmem[[t]], "==") }
    if (b %% 10 == 0) cat(sprintf("  %s: %d/%d replicates (%.1f min)\n", mname, b, B_BOOT, (proc.time()[["elapsed"]] - t0) / 60))
  }
  stopifnot(b_ok >= 10)
  by_year <- colMeans(S_rep, na.rm = TRUE); per_rep <- rowMeans(S_rep, na.rm = TRUE); s <- mean(per_rep)
  fl <- dynamicmultiplex:::.floor_lookup(s, tab[tab$level == "partition_nmi", ])
  node_st <- unlist(lapply(seq_len(T_), function(t) node_acc[[t]] / b_ok))
  up <- unlist(lapply(seq_len(T_), function(t) { P <- co[[t]] / b_ok; P[upper.tri(P)] }))
  rows[[mname]] <- data.frame(net = net, method = mname, years = T_, B = b_ok, stability = s, stability_mc_se = sd(per_rep) / sqrt(length(per_rep)),
    floor = fl$floor, floor_median = fl$median, bin = fl$bin, node_share_ge_0.9 = mean(node_st >= 0.9),
    pairs_together = mean(up >= 0.9), pairs_apart = mean(up <= 0.1), pairs_undetermined = mean(up > 0.1 & up < 0.9),
    mean_K = mean(vapply(mem0, function(m) length(unique(m)), integer(1))), runtime_min = (proc.time()[["elapsed"]] - t0) / 60,
    stringsAsFactors = FALSE)
  yrows[[mname]] <- data.frame(net = net, method = mname, year = yrs, stability = by_year, present = sapply(AL, nrow),
    K = vapply(mem0, function(m) length(unique(m)), integer(1)), stringsAsFactors = FALSE)
  cat(sprintf("[stab-emp] %s %s: stability %.3f floor %.2f (%s) decided %.0f%%  [%.1f min]\n", net, mname, s, fl$floor, fl$bin,
              100 * (1 - rows[[mname]]$pairs_undetermined), rows[[mname]]$runtime_min))
}
write.csv(do.call(rbind, rows), outf, row.names = FALSE); write.csv(do.call(rbind, yrows), outy, row.names = FALSE)
cat("[stab-emp] wrote", outf, "\n")
