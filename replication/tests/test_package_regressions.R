# =============================================================================
# replication/tests/test_package_regressions.R
# Regression tests for two package bugs found 2026-09-20 (fixed in 1.2.1):
#   1. weighted_jaccard_similarity counted layer-wide node strengths for nodes
#      that were NOT members of the community, so disjoint communities scored 1.
#   2. fit_multilayer_identity_ties optimised plain single-graph modularity on
#      the stacked supra-graph instead of Mucha multislice modularity, and
#      returned each layer as one community at omega <= 1.
# Usage: DM_ROOT=. Rscript replication/tests/test_package_regressions.R
# =============================================================================
ROOT <- Sys.getenv("DM_ROOT", unset = getwd())
suppressMessages(pkgload::load_all(file.path(ROOT, "r_code"), quiet = TRUE))
suppressPackageStartupMessages(library(igraph))

# --- 1. weighted Jaccard --------------------------------------------------
wj <- dynamicmultiplex:::weighted_jaccard_similarity
w  <- setNames(rep(5, 10), as.character(1:10))
stopifnot(wj(c("1","2","3","4","5"), c("6","7","8","9","10"), w, w) == 0)
stopifnot(abs(wj(c("1","2","3","4","5"), c("4","5","6","7","8"), w, w) - 0.25) < 1e-12)
stopifnot(abs(wj(c("1","2"), c("1","2"), w, w) - 1) < 1e-12)
cat("weighted Jaccard: ok\n")

# --- 2. multislice on a static planted partition --------------------------
set.seed(123); n <- 100; K <- 4; mem <- sample(seq_len(K), n, replace = TRUE)
build <- function(m, p_in, p_out) {
  Pr <- ifelse(outer(m, m, "=="), p_in, p_out); diag(Pr) <- 0
  A <- matrix(0, n, n); up <- upper.tri(Pr); A[up] <- rbinom(sum(up), 1, Pr[up]); A <- A + t(A)
  dimnames(A) <- list(as.character(seq_len(n)), as.character(seq_len(n))); A
}
L <- lapply(1:5, function(t) build(mem, 0.3, 0.04))
fit <- fit_multilayer_identity_ties(L, algorithm = "leiden", omega = 1, seed = 123)
m <- extract_meta_membership(fit)
stopifnot(all(vapply(m, function(x) length(unique(x)), integer(1)) == K))   # K communities in every layer
stopifnot(length(unique(unlist(m))) == K)                                     # and they persist across layers
stopifnot(all(vapply(m, function(x) igraph::compare(x, mem, "nmi"), numeric(1)) > 0.95))
fit0 <- fit_multilayer_identity_ties(L, algorithm = "leiden", omega = 0, seed = 123)
stopifnot(length(unique(unlist(extract_meta_membership(fit0)))) == K * 5)   # omega = 0: layers decoupled
cat("multislice identity ties: ok\n")

# --- 3. weighted Jaccard fitter recovers the same static partition --------
fw <- fit_multilayer_weighted_jaccard(L, algorithm = "leiden", seed = 123)
mw <- extract_meta_membership(fw)
stopifnot(all(vapply(mw, function(x) length(unique(x)), integer(1)) == K))
stopifnot(mean(vapply(mw, function(x) igraph::compare(x, mem, "nmi"), numeric(1))) > 0.85)  # (before 1.2.1: ~1 community per layer)
cat("weighted Jaccard fitter: ok\n")

# --- 4. partition_stability() (1.3.0) ------------------------------------
boot <- bootstrap_multilayer(L, fit_type = "jaccard", algorithm = "leiden", n_boot = 10, seed = 123)
stopifnot(!is.null(boot$stability_samples), nrow(boot$stability_samples$nmi) == boot$n_boot)
ps <- partition_stability(boot)
stopifnot(ps$stability > 0.9, ps$floor > 0.9, ps$bin == "[0.9, 1.0)",
          abs(sum(ps$pair_summary) - 1) < 1e-12, length(ps$pairs) == length(L))
stopifnot(inherits(tryCatch(co_assignment_ci(boot, method = "calibrated"), error = function(e) e), "error"))
cat("partition_stability: ok\nALL PACKAGE REGRESSION TESTS PASSED\n")
