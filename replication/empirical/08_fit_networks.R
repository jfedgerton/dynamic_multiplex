# =============================================================================
# 08_fit_networks.R
# Empirical case study: all DynMux specs + baselines on one real dynamic
# international network (atop | dca | igo | trade), as built by
# 07_build_networks.R. DynSBM excluded (does not scale to the 203-layer open
# population). DynMux runs on the native per-year graphs with non-present
# (isolate) states dropped; fixed-node baselines run on the union node set and
# are then scored on the per-layer `present` mask. multinet GLouvain is fit
# active-per-layer (actors added only in the years they are present), so no
# second pass (multinet_fix.R) is needed.
#
# Consolidates: replication/extended/paper_scripts/emp_fix.R (authoritative
#               fitter logic; supersedes 06_alliance_dca_empirical.R) and
#               replication/extended/paper_scripts/multinet_fix.R
#
# Usage:  Rscript 08_fit_networks.R <atop|dca|igo|trade>
# Env:    DM_ROOT  project root (default getwd()); FORCE=1 recomputes everything
# Input:  $DM_ROOT/output/empirical_data/<net>_series.rds, <net>_union.rds
# Output: $DM_ROOT/output/empirical/<net>_partitions.rds
#         list(net, years, partitions, timing_s); partitions[[method]][[t]] is a
#         named integer vector (names = ccodes present in layer t).
# Idempotent: if the output already holds every method, the script exits
# without refitting; a partial file (earlier crash) resumes the missing methods.
# =============================================================================
set.seed(123)
suppressMessages({ library(igraph) })

ROOT     <- Sys.getenv("DM_ROOT", unset = getwd())
EMP_DATA <- file.path(ROOT, "output", "empirical_data")
EMP_OUT  <- file.path(ROOT, "output", "empirical")
dir.create(EMP_OUT, recursive = TRUE, showWarnings = FALSE)
stopifnot(dir.exists(EMP_DATA), dir.exists(EMP_OUT))
suppressMessages(pkgload::load_all(file.path(ROOT, "r_code"), quiet = TRUE))
stopifnot(exists("fit_multilayer_identity_ties"), exists("fit_multilayer_jaccard"),
          exists("fit_multilayer_overlap"), exists("extract_meta_membership"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript 08_fit_networks.R <atop|dca|igo|trade>", call. = FALSE)
net <- args[1]
stopifnot(net %in% c("atop", "dca", "igo", "trade"))
FORCE <- identical(Sys.getenv("FORCE", "0"), "1")

HAVE_CLUE     <- requireNamespace("clue",     quietly = TRUE)
HAVE_MULTINET <- requireNamespace("multinet", quietly = TRUE)
if (!HAVE_MULTINET) message("package 'multinet' not installed: multinet GLouvain will be skipped")

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
series_f <- file.path(EMP_DATA, sprintf("%s_series.rds", net))
union_f  <- file.path(EMP_DATA, sprintf("%s_union.rds",  net))
if (!file.exists(series_f) || !file.exists(union_f))
  stop("Missing ", series_f, " and/or ", union_f, " -- run 07_build_networks.R ", net, " first.", call. = FALSE)
S <- readRDS(series_f); yrs <- S$years
U <- readRDS(union_f)
if (is.null(U$present)) U$present <- U$active   # older union files saved only `active`
stopifnot(!is.null(U$layers), !is.null(U$present), !is.null(U$names),
          length(S$graph_layers) == length(yrs),
          length(U$layers) == length(yrs), length(U$present) == length(yrs),
          all(sapply(U$layers, nrow) == length(U$names)),
          all(sapply(U$present, length) == length(U$names)),
          all(sapply(S$graph_layers, function(g) identical(igraph::V(g)$name, U$names))))
if (!is.null(U$active)) stopifnot(identical(U$active, U$present))

# DynMux input: native per-year graphs with non-present states dropped
GLc  <- lapply(seq_along(S$graph_layers), function(k) {
  g <- S$graph_layers[[k]]; igraph::delete_vertices(g, which(!U$present[[k]])) })
# baseline input: union node set + per-layer mask (present == active, see 05)
simU <- list(layers = U$layers, active = U$present, links = NULL)
stopifnot(all(sapply(seq_along(GLc), function(k) igraph::vcount(GLc[[k]]) == sum(U$present[[k]]))))

# The per-year graphs have different node sets (states enter and exit). The
# package fitters accept that with allow_unequal_nodes = TRUE (matching nodes
# by vertex name; only nodes present in both layers of a linked pair
# contribute to the interlayer similarity). Padding every layer to the union
# set instead would add thousands of isolate singleton communities to the
# second-stage community graph and change the results.
UNEQ <- list(allow_unequal_nodes = TRUE)

# ---------------------------------------------------------------------------
# Methods (emp_fix.R, with multinet from multinet_fix.R)
# ---------------------------------------------------------------------------
leiden_layer <- function(mat) {
  g <- igraph::graph_from_adjacency_matrix(mat, mode = "undirected", weighted = TRUE, diag = FALSE)
  if (igraph::ecount(g) == 0L) return(rep(1L, nrow(mat)))
  as.integer(igraph::membership(igraph::cluster_leiden(g, objective_function = "modularity", weights = igraph::E(g)$weight)))
}
method_pooled <- function(sim) {
  L <- sim$layers; agg <- Reduce("+", L)
  g <- igraph::graph_from_adjacency_matrix(agg, mode = "undirected", weighted = TRUE, diag = FALSE)
  mem <- igraph::membership(igraph::cluster_leiden(g, objective_function = "modularity", weights = igraph::E(g)$weight))
  replicate(length(L), as.integer(mem), simplify = FALSE)
}
.greedy <- function(M) {
  nr <- nrow(M); a <- rep(NA_integer_, nr); u <- integer(0)
  for (i in order(-apply(M, 1, max))) { c <- order(-M[i, ]); c <- c[!(c %in% u)]; a[i] <- c[1]; u <- c(u, c[1]) }
  a
}
match_hung <- function(mems) {
  o <- vector("list", length(mems)); o[[1]] <- as.integer(mems[[1]]); nf <- max(o[[1]]) + 1L
  for (t in 2:length(mems)) {
    pv <- o[[t - 1L]]; cu <- as.integer(mems[[t]]); cl <- sort(unique(cu)); pl <- sort(unique(pv))
    M <- matrix(0, length(cl), length(pl))
    for (i in seq_along(cl)) for (j in seq_along(pl)) M[i, j] <- sum(cu == cl[i] & pv == pl[j])
    d <- max(nrow(M), ncol(M)); Ms <- matrix(0, d, d); Ms[seq_len(nrow(M)), seq_len(ncol(M))] <- M
    as <- if (HAVE_CLUE) as.integer(clue::solve_LSAP(max(Ms) - Ms)) else .greedy(Ms)
    mp <- rep(NA_integer_, length(cl))
    for (i in seq_along(cl)) { co <- as[i]; if (co <= length(pl) && M[i, co] > 0) mp[i] <- pl[co] }
    for (i in seq_along(cl)) if (is.na(mp[i])) { mp[i] <- nf; nf <- nf + 1L }
    nf <- max(nf, max(mp) + 1L); names(mp) <- as.character(cl); o[[t]] <- unname(mp[as.character(cu)])
  }
  lapply(o, as.integer)
}
# multinet GLouvain, active-per-layer (multinet_fix.R): actors are added to a
# layer only in the years they are present; present-but-isolate states get
# singletons post hoc.
method_multinet <- function(sim) {
  if (!HAVE_MULTINET) return(NULL)
  L <- sim$layers; ac <- sim$active; n <- nrow(L[[1]]); T_ <- length(L); nm <- as.character(seq_len(n))
  nn <- multinet::ml_empty()
  for (t in seq_len(T_)) {
    keep <- which(ac[[t]]); if (length(keep) == 0) next
    A <- L[[t]][keep, keep, drop = FALSE]
    g <- igraph::graph_from_adjacency_matrix(A, mode = "undirected", diag = FALSE)
    igraph::V(g)$name <- nm[keep]
    multinet::add_igraph_layer_ml(nn, g, paste0("layer", t))
  }
  cat("  multinet actor-layer entries:", sum(sapply(ac, sum)), " (vs padded", n * T_, ")\n")
  cm <- multinet::glouvain_ml(nn, gamma = 1, omega = 1)
  lapply(seq_len(T_), function(t) {
    s <- cm[cm$layer == paste0("layer", t), , drop = FALSE]
    v <- rep(NA_integer_, n)
    if (nrow(s) > 0) v[as.integer(s$actor)] <- as.integer(as.factor(s$cid))
    if (anyNA(v)) v[is.na(v)] <- max(0L, v, na.rm = TRUE) + seq_len(sum(is.na(v)))
    v
  })
}
dm <- function(ff, rp, extra = list()) function() {
  f <- do.call(ff, c(list(GLc, algorithm = "leiden", resolution_parameter = rp), extra))
  m <- extract_meta_membership(f)
  lapply(seq_along(GLc), function(t) setNames(as.integer(m[[t]]), igraph::V(GLc[[t]])$name))
}
comp <- function(fn) function() {
  d <- fn(simU); if (is.null(d)) return(NULL)
  lapply(seq_along(d), function(t) { a <- U$present[[t]]; setNames(as.integer(d[[t]][a]), U$names[a]) })
}
METHODS <- list(
  "DynMux multislice r1" = dm(fit_multilayer_identity_ties, 1, UNEQ),
  "DynMux multislice r2" = dm(fit_multilayer_identity_ties, 2, UNEQ),
  "DynMux multislice r4" = dm(fit_multilayer_identity_ties, 4, UNEQ),
  "DynMux Jaccard r1"    = dm(fit_multilayer_jaccard, 1, UNEQ),
  "DynMux Jaccard r2"    = dm(fit_multilayer_jaccard, 2, UNEQ),
  "DynMux Jaccard r4"    = dm(fit_multilayer_jaccard, 4, UNEQ),
  "DynMux Overlap r1"    = dm(fit_multilayer_overlap, 1, UNEQ),
  "DynMux Overlap r2"    = dm(fit_multilayer_overlap, 2, UNEQ),
  "DynMux Overlap r4"    = dm(fit_multilayer_overlap, 4, UNEQ),
  "Pooled Leiden"        = comp(method_pooled),
  "Cross-sectional + Hungarian" = comp(function(s) match_hung(lapply(s$layers, leiden_layer))),
  "multinet GLouvain"    = comp(method_multinet))
# NOTE: "DynMux multislice" here is the identity-tie (Mucha multislice)
# specification with the package's default adjacent-layer coupling
# (layer_links = NULL). The sources fit no all-to-all ("full") multislice
# variant on the empirical networks, so none is fit here.

# ---------------------------------------------------------------------------
# Fit (per-method timing; cached results reused unless FORCE=1)
# ---------------------------------------------------------------------------
outf <- file.path(EMP_OUT, sprintf("%s_partitions.rds", net))
results <- list(); timing <- c()
# multinet GLouvain is optional (as in the sources): required only when installed
required <- setdiff(names(METHODS), if (HAVE_MULTINET) character(0) else "multinet GLouvain")
if (file.exists(outf) && !FORCE) {
  prev <- readRDS(outf); results <- prev$partitions
  if (!is.null(prev$timing_s)) timing <- prev$timing_s
  stopifnot(identical(prev$net, net), identical(prev$years, yrs))
  if (all(required %in% names(results))) {
    cat(sprintf("[%s] %s already holds all %d required methods; skipping (set FORCE=1 to refit)\n",
                net, basename(outf), length(required)))
    quit(save = "no", status = 0)
  }
  cat(sprintf("[%s] resuming: %d/%d required methods cached\n", net, sum(required %in% names(results)), length(required)))
}
if (FORCE) cat(sprintf("[%s] FORCE=1: refitting all methods\n", net))

for (mn in names(METHODS)) {
  if (!is.null(results[[mn]])) { cat(sprintf("[%s] %-28s cached\n", net, mn)); next }
  gc(); t0 <- proc.time()["elapsed"]
  det <- tryCatch({
    setTimeLimit(elapsed = 7200, transient = TRUE); r <- METHODS[[mn]](); setTimeLimit(); r
  }, error = function(e) { setTimeLimit(); message(sprintf("  %s FAILED: %s", mn, conditionMessage(e))); NULL })
  el <- as.numeric(proc.time()["elapsed"] - t0)
  if (!is.null(det)) {
    stopifnot(length(det) == length(yrs),
              all(sapply(seq_along(det), function(t) identical(names(det[[t]]), U$names[U$present[[t]]]))),
              !anyNA(unlist(det)))
  }
  results[[mn]] <- det; timing[mn] <- el
  ok <- !is.null(det); nc <- if (ok) length(unique(unlist(det))) else NA
  saveRDS(list(net = net, years = yrs, partitions = results, timing_s = timing), outf)
  cat(sprintf("[%s] %-28s %s ncomm=%s (%.1f min)\n", net, mn, if (ok) "ok" else "NA", nc, el / 60)); flush.console()
}

# ---------------------------------------------------------------------------
# Summary + assertions
# ---------------------------------------------------------------------------
cat(sprintf("\n[%s] timing summary\n", net))
for (mn in names(METHODS)) {
  ok <- !is.null(results[[mn]])
  status <- if (ok) "ok" else if (mn == "multinet GLouvain" && !HAVE_MULTINET) "skipped (no multinet)" else "FAILED"
  cat(sprintf("  %-28s %-22s %8.1f s  ncomm=%s\n", mn, status,
              if (is.na(timing[mn])) NA_real_ else timing[mn],
              if (ok) length(unique(unlist(results[[mn]]))) else "NA"))
}
missing  <- required[!required %in% names(results)]
if (length(missing))
  stop(sprintf("[%s] methods failed or timed out: %s (partial results saved to %s)",
               net, paste(missing, collapse = "; "), outf), call. = FALSE)
stopifnot(file.exists(outf))
cat(sprintf("ALL_DONE %s -> %s (%.1f min total)\n", net, outf, sum(timing, na.rm = TRUE) / 60))
