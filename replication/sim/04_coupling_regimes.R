# =============================================================================
# replication/sim/04_coupling_regimes.R
#
# WHEN DOES JACCARD COUPLING BEAT OVERLAP COUPLING (AND VICE VERSA)?
#
# DynMux links community c (layer t) to community c' (layer t') with weight
#   Jaccard : |c n c'| / |c u c'|            penalises size mismatch
#   Overlap : |c n c'| / min(|c|, |c'|)      nesting scores 1 regardless of size
# and runs second-stage detection on that community graph. The two weights
# agree when communities keep their size and membership and diverge under
# merges, splits, size skew, shrinkage and per-layer over-splitting. This
# script generates six regimes that isolate those mechanisms and scores every
# DynMux coupling on the SAME simulated series (fully paired, 30 reps).
#
# Regimes (T = 10 layers, fresh SBM draw each layer, adjacent links):
#   balanced : K = 4 equal communities, node-level switching           (control)
#   split    : K0 = 3; after t* = 5 one (low) or all (high) communities split
#              into two equal, persistent children
#   merge    : K0 = 6; after t* = 5 one (low) or three (high) disjoint pairs
#              merge into a single community
#   sizeskew : K = 5 with one giant community holding 40% (low) / 60% (high)
#              of the nodes and four small ones; node-level switching with
#              size-proportional destinations
#   shrink   : open population, K = 4; communities 1 and 2 lose 10% (low) /
#              20% (high) of their members to inactivity every layer, mild
#              switching among the rest; metrics on active nodes
#   nested   : K = 3 coarse communities, each with two sub-blocks whose
#              internal density is 1.5x (low) / 2.5x (high) the coarse
#              within-density; truth is the COARSE partition, so per-layer
#              detection may over-split and the coupling has to reassemble
#
# Two ground-truth conventions for the TRACKED (cross-layer) partition, since
# "which coupling is right" depends on what counts as the same community:
#   break    : a split produces two NEW lineages, a merge produces one NEW
#              lineage (identity ends when membership changes wholesale)
#   continue : one child of a split keeps the parent's lineage id, a merged
#              community keeps the first parent's id (identity persists
#              through nesting)
# The two conventions coincide for balanced / sizeskew / shrink / nested.
#
# Metrics per (config, rep, method), on active nodes:
#   nmi_layer, k_mae, comembership_acc, mean_n_comm, total_n_comm, runtime_s
#                                (as in sim/01, computed against truth_break;
#                                 per-layer metrics do not depend on lineage)
#   nmi_joint_break / nmi_joint_continue   NMI of the node x layer tracked
#                                partition vs each truth convention
#   purity_{break,continue}      size-weighted share of each tracked
#                                meta-community that comes from ONE true
#                                lineage (low = over-linking / chaining)
#   completeness_{break,continue} size-weighted share of each true lineage
#                                that lands in ONE meta-community (low =
#                                under-linking / fractured lineages)
#   n_true_break, n_true_continue  number of true lineages
#
# Density: dens_params(r, K_fine) with rho = 0.10 as in sim/01, where K_fine
# is the finest block count the regime ever uses (so children / sub-blocks
# stay detectable). Separation r in {4, 8}.
#
# CONFIG GRID (48): regime x N{100,200} x r{4,8} x intensity{low,high},
# shuffled with set.seed(123) so a partial run samples every condition.
# Per-rep seed = 31000 + TASK*1000 + rep.
#
# Usage (local smoke): CMP_MINI=1 CMP_CFG=1 DM_ROOT=. Rscript replication/sim/04_coupling_regimes.R
# Array: one SLURM_ARRAY_TASK_ID per config (1..48); see slurm/04_coupling.sbatch.
# Output: $DM_ROOT/output/coupling/coup_cfg%02d.csv
# Sequential by design (plain rep loop) so it stays debuggable.
# =============================================================================

ROOT <- Sys.getenv("DM_ROOT", unset = getwd())
suppressMessages(pkgload::load_all(file.path(ROOT, "r_code"), quiet = TRUE))
suppressPackageStartupMessages(library(igraph))

TASK <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", Sys.getenv("CMP_CFG", "1")))
MINI <- identical(Sys.getenv("CMP_MINI", "0"), "1")
REPS <- if (MINI) 2L else 30L
RHO  <- 0.10
T_   <- 10L
TSTAR <- 5L

outdir <- file.path(ROOT, "output", "coupling")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# SHARED BUILDING BLOCKS (copied from sim/01 so this file is self-contained)
# =============================================================================
dens_params <- function(r, K, rho = RHO) {
  p_in  <- rho * r * K / (r + K - 1)
  p_out <- p_in / r
  list(p_in = min(max(p_in, 0), 1), p_out = min(max(p_out, 0), 1))
}

# One undirected 0/1 SBM layer from a block-membership vector; `sub` optionally
# marks a finer block whose within-density is p_sub (nested regime).
build_layer <- function(m, p_in, p_out, sub = NULL, p_sub = NULL) {
  n  <- length(m)
  Pr <- ifelse(outer(m, m, "=="), p_in, p_out)
  if (!is.null(sub)) Pr[outer(sub, sub, "==") & outer(m, m, "==")] <- p_sub
  diag(Pr) <- 0
  A <- matrix(0, n, n); up <- upper.tri(Pr)
  A[up] <- rbinom(sum(up), 1, Pr[up]); A + t(A)
}

alive_idx <- function(sim, t) {
  if (!is.null(sim$active)) which(as.logical(sim$active[[t]])) else seq_len(nrow(sim$layers[[t]]))
}
AL <- function(sim) lapply(seq_along(sim$layers), function(t) {
  idx <- alive_idx(sim, t)
  A <- sim$layers[[t]][idx, idx, drop = FALSE]
  dimnames(A) <- list(as.character(idx), as.character(idx)); A
})
EXPn <- function(m, sim, n) lapply(seq_along(m), function(t) {
  idx <- alive_idx(sim, t); v <- rep(NA_integer_, n); v[idx] <- as.integer(m[[t]]); v
})
dm <- function(fitfun) function(sim) {
  n <- nrow(sim$layers[[1]])
  f <- fitfun(AL(sim), algorithm = "leiden", layer_links = sim$links, allow_unequal_nodes = TRUE)
  EXPn(extract_meta_membership(f), sim, n)
}

METHODS <- list(
  "DynMux Jaccard"          = dm(fit_multilayer_jaccard),
  "DynMux Overlap"          = dm(fit_multilayer_overlap),
  "DynMux weighted Jaccard" = dm(fit_multilayer_weighted_jaccard),
  "DynMux weighted Overlap" = dm(fit_multilayer_weighted_overlap),
  "DynMux multislice"       = dm(fit_multilayer_identity_ties)
)

adjacent_links <- function(T_) data.frame(from = seq_len(T_ - 1L), to = 2:T_, weight = 1)

# =============================================================================
# GENERATORS  (each returns list(layers, truth_break, truth_continue, links[, active]))
# =============================================================================

# 1. balanced: K equal communities, node-level switching. Control regime.
sim_balanced <- function(n, r, intensity, seed) {
  set.seed(seed); K <- 4L
  dp <- dens_params(r, K)
  p_switch <- if (intensity == "high") 0.10 else 0.03
  mem <- sample(seq_len(K), n, replace = TRUE)
  truth <- vector("list", T_); layers <- vector("list", T_)
  for (t in seq_len(T_)) {
    if (t > 1L) { sw <- runif(n) < p_switch
      if (any(sw)) mem[sw] <- sample(seq_len(K), sum(sw), replace = TRUE) }
    layers[[t]] <- build_layer(mem, dp$p_in, dp$p_out); truth[[t]] <- mem
  }
  list(layers = layers, truth_break = truth, truth_continue = truth, links = adjacent_links(T_))
}

# 2. split: after TSTAR, n_split of the K0 communities split into two equal
#    children that persist. Layers are drawn from the FINE (post-split) blocks.
sim_split <- function(n, r, intensity, seed) {
  set.seed(seed); K0 <- 3L
  dp <- dens_params(r, 2L * K0)
  n_split <- if (intensity == "high") 3L else 1L
  mem  <- sample(seq_len(K0), n, replace = TRUE)
  half <- integer(n)
  for (k in seq_len(K0)) { idx <- which(mem == k); half[idx] <- sample(rep(1:2, length.out = length(idx))) }
  splitters <- sort(sample(seq_len(K0), n_split))
  tb <- vector("list", T_); tc <- vector("list", T_); layers <- vector("list", T_)
  for (t in seq_len(T_)) {
    b <- mem; cc <- mem
    if (t > TSTAR) for (k in splitters) {
      i1 <- which(mem == k & half == 1L); i2 <- which(mem == k & half == 2L)
      b[i1] <- K0 + 2L * k - 1L; b[i2] <- K0 + 2L * k      # break: two new lineages
      cc[i2] <- K0 + k                                     # continue: child 1 keeps k
    }
    layers[[t]] <- build_layer(b, dp$p_in, dp$p_out); tb[[t]] <- b; tc[[t]] <- cc
  }
  list(layers = layers, truth_break = tb, truth_continue = tc, links = adjacent_links(T_))
}

# 3. merge: K0 = 6 communities in disjoint pairs (1,2) (3,4) (5,6); after TSTAR
#    n_merge pairs merge. Layers are drawn from the COARSE (post-merge) blocks.
sim_merge <- function(n, r, intensity, seed) {
  set.seed(seed); K0 <- 6L
  dp <- dens_params(r, K0)
  n_merge <- if (intensity == "high") 3L else 1L
  mem <- sample(seq_len(K0), n, replace = TRUE)
  mergers <- sort(sample(1:3, n_merge))                    # pair j = communities (2j-1, 2j)
  tb <- vector("list", T_); tc <- vector("list", T_); layers <- vector("list", T_)
  for (t in seq_len(T_)) {
    b <- mem; cc <- mem
    if (t > TSTAR) for (j in mergers) {
      idx <- which(mem %in% c(2L * j - 1L, 2L * j))
      b[idx] <- K0 + j                                     # break: one new lineage
      cc[idx] <- 2L * j - 1L                               # continue: keeps first parent's id
    }
    layers[[t]] <- build_layer(cc, dp$p_in, dp$p_out); tb[[t]] <- b; tc[[t]] <- cc
  }
  list(layers = layers, truth_break = tb, truth_continue = tc, links = adjacent_links(T_))
}

# 4. sizeskew: one giant community plus four small ones; switching destinations
#    are size-proportional so the skew is stable over layers.
sim_sizeskew <- function(n, r, intensity, seed) {
  set.seed(seed); K <- 5L
  dp <- dens_params(r, K)
  giant <- if (intensity == "high") 0.60 else 0.40
  shares <- c(giant, rep((1 - giant) / 4, 4)); p_switch <- 0.05
  mem <- sample(seq_len(K), n, replace = TRUE, prob = shares)
  truth <- vector("list", T_); layers <- vector("list", T_)
  for (t in seq_len(T_)) {
    if (t > 1L) { sw <- runif(n) < p_switch
      if (any(sw)) mem[sw] <- sample(seq_len(K), sum(sw), replace = TRUE, prob = shares) }
    layers[[t]] <- build_layer(mem, dp$p_in, dp$p_out); truth[[t]] <- mem
  }
  list(layers = layers, truth_break = truth, truth_continue = truth, links = adjacent_links(T_))
}

# 5. shrink: open population. Communities 1 and 2 lose `frac` of their active
#    members to inactivity every layer (never return); everyone else switches
#    with p = 0.03 among the four communities. Inactive nodes get unique truth
#    ids (ignored by the active-node metrics), as in sim/01 birthdeath.
sim_shrink <- function(n, r, intensity, seed) {
  set.seed(seed); K <- 4L
  dp <- dens_params(r, K)
  frac <- if (intensity == "high") 0.20 else 0.10; p_switch <- 0.03
  mem <- sample(seq_len(K), n, replace = TRUE); state <- rep(TRUE, n)
  truth <- vector("list", T_); layers <- vector("list", T_); active <- vector("list", T_)
  for (t in seq_len(T_)) {
    if (t > 1L) {
      for (k in 1:2) { idx <- which(state & mem == k)
        n_out <- floor(frac * length(idx)); if (n_out >= 1L && length(idx) - n_out >= 3L)
          state[sample(idx, n_out)] <- FALSE }
      sw <- state & runif(n) < p_switch
      if (any(sw)) mem[sw] <- sample(seq_len(K), sum(sw), replace = TRUE)
    }
    active[[t]] <- state
    A <- matrix(0, n, n); ai <- which(state)
    sub <- build_layer(mem[ai], dp$p_in, dp$p_out); A[ai, ai] <- sub
    layers[[t]] <- A
    tr <- mem; inact <- which(!state); if (length(inact)) tr[inact] <- K + seq_along(inact)
    truth[[t]] <- tr
  }
  list(layers = layers, truth_break = truth, truth_continue = truth, links = adjacent_links(T_), active = active)
}

# 6. nested: K = 3 coarse communities x 2 sub-blocks. Truth is coarse. Coarse
#    switching p = 0.05 (random new sub-block); within-coarse sub-block flip
#    p = 0.05. Sub-block within-density = mult x coarse within-density.
sim_nested <- function(n, r, intensity, seed) {
  set.seed(seed); K <- 3L
  dp <- dens_params(r, K)
  mult <- if (intensity == "high") 2.5 else 1.5
  p_sub <- min(1, dp$p_in * mult); p_switch <- 0.05; p_flip <- 0.05
  mem <- sample(seq_len(K), n, replace = TRUE); sub <- sample(1:2, n, replace = TRUE)
  truth <- vector("list", T_); layers <- vector("list", T_)
  for (t in seq_len(T_)) {
    if (t > 1L) {
      sw <- runif(n) < p_switch
      if (any(sw)) { mem[sw] <- sample(seq_len(K), sum(sw), replace = TRUE); sub[sw] <- sample(1:2, sum(sw), replace = TRUE) }
      fl <- !sw & runif(n) < p_flip; sub[fl] <- 3L - sub[fl]
    }
    layers[[t]] <- build_layer(mem, dp$p_in, dp$p_out, sub = sub, p_sub = p_sub); truth[[t]] <- mem
  }
  list(layers = layers, truth_break = truth, truth_continue = truth, links = adjacent_links(T_))
}

simulate_regime <- function(cfg, seed) switch(cfg$regime,
  balanced = sim_balanced(cfg$n, cfg$r, cfg$intensity, seed),
  split    = sim_split   (cfg$n, cfg$r, cfg$intensity, seed),
  merge    = sim_merge   (cfg$n, cfg$r, cfg$intensity, seed),
  sizeskew = sim_sizeskew(cfg$n, cfg$r, cfg$intensity, seed),
  shrink   = sim_shrink  (cfg$n, cfg$r, cfg$intensity, seed),
  nested   = sim_nested  (cfg$n, cfg$r, cfg$intensity, seed),
  stop("unknown regime: ", cfg$regime))

# =============================================================================
# METRICS  (active nodes only)
# =============================================================================
stack_active <- function(det, truth, active) {
  d <- integer(0); tr <- integer(0)
  for (t in seq_along(det)) {
    idx <- if (is.null(active)) seq_along(det[[t]]) else which(active[[t]])
    d <- c(d, as.integer(det[[t]][idx])); tr <- c(tr, as.integer(truth[[t]][idx]))
  }
  list(d = d, tr = tr)
}
lineage_scores <- function(det, truth, active) {
  s <- stack_active(det, truth, active)
  ct <- table(s$d, s$tr)
  c(nmi_joint = igraph::compare(s$d, s$tr, method = "nmi"),
    purity = sum(apply(ct, 1, max)) / sum(ct),
    completeness = sum(apply(ct, 2, max)) / sum(ct),
    n_true = length(unique(s$tr)))
}
eval_method <- function(det, sim) {
  active <- sim$active
  per <- lapply(seq_along(det), function(t) {
    idx <- if (is.null(active)) seq_along(det[[t]]) else which(active[[t]])
    list(d = as.integer(det[[t]][idx]), tr = as.integer(sim$truth_break[[t]][idx]))
  })
  nmi_layer <- mean(vapply(per, function(x) igraph::compare(x$d, x$tr, method = "nmi"), numeric(1)))
  k_mae <- mean(vapply(per, function(x) abs(length(unique(x$d)) - length(unique(x$tr))), numeric(1)))
  comemb <- mean(vapply(per, function(x) {
    Dd <- outer(x$d, x$d, "=="); Dt <- outer(x$tr, x$tr, "=="); up <- upper.tri(Dd); mean(Dd[up] == Dt[up]) }, numeric(1)))
  mean_n_comm <- mean(vapply(per, function(x) length(unique(x$d)), numeric(1)))
  total_n_comm <- length(unique(unlist(lapply(per, `[[`, "d"))))
  lb <- lineage_scores(det, sim$truth_break, active)
  lc <- lineage_scores(det, sim$truth_continue, active)
  c(nmi_layer = nmi_layer, k_mae = k_mae, comembership_acc = comemb,
    mean_n_comm = mean_n_comm, total_n_comm = total_n_comm,
    nmi_joint_break = lb[["nmi_joint"]], purity_break = lb[["purity"]],
    completeness_break = lb[["completeness"]], n_true_break = lb[["n_true"]],
    nmi_joint_continue = lc[["nmi_joint"]], purity_continue = lc[["purity"]],
    completeness_continue = lc[["completeness"]], n_true_continue = lc[["n_true"]])
}

# =============================================================================
# CONFIG GRID (48) + fixed shuffle
# =============================================================================
cfgs <- expand.grid(regime = c("balanced", "split", "merge", "sizeskew", "shrink", "nested"),
                    n = c(100L, 200L), r = c(4, 8), intensity = c("low", "high"),
                    stringsAsFactors = FALSE)
stopifnot(nrow(cfgs) == 48L)
set.seed(123); ord <- sample(nrow(cfgs))
if (TASK < 1L || TASK > nrow(cfgs)) stop("task out of range")
cfg <- cfgs[ord[TASK], ]
outfile <- file.path(outdir, sprintf("coup_cfg%02d.csv", TASK))
if (nzchar(Sys.getenv("SLURM_ARRAY_TASK_ID")) && file.exists(outfile)) {
  cat("[skip] task", TASK, "already complete\n"); quit(save = "no")
}
cat(sprintf("[coupling] task=%d regime=%s n=%d r=%g intensity=%s reps=%d\n",
            TASK, cfg$regime, cfg$n, cfg$r, cfg$intensity, REPS))

# =============================================================================
# MAIN LOOP  (one sim per rep, every method on the same sim)
# =============================================================================
rows <- list(); t0 <- proc.time()[["elapsed"]]
for (rep in seq_len(REPS)) {
  seed <- 31000L + TASK * 1000L + rep
  sim <- simulate_regime(cfg, seed)
  stopifnot(length(sim$layers) == T_, length(sim$truth_break) == T_, length(sim$truth_continue) == T_)
  for (mname in names(METHODS)) {
    tm <- system.time(det <- tryCatch(METHODS[[mname]](sim), error = function(e) {
      message(sprintf("rep %d %s failed: %s", rep, mname, conditionMessage(e))); NULL }))[["elapsed"]]
    if (is.null(det)) next
    stopifnot(length(det) == T_)
    m <- eval_method(det, sim)
    rows[[length(rows) + 1L]] <- data.frame(task = TASK, regime = cfg$regime, n = cfg$n, r = cfg$r,
      intensity = cfg$intensity, rep = rep, seed = seed, method = mname, runtime_s = tm,
      as.list(m), stringsAsFactors = FALSE)
  }
  if (rep %% 5L == 0L || MINI) cat(sprintf("  rep %d/%d  %.1f min\n", rep, REPS, (proc.time()[["elapsed"]] - t0) / 60))
}
out <- do.call(rbind, rows)
stopifnot(nrow(out) >= length(METHODS) * (if (MINI) 1L else 20L))
write.csv(out, outfile, row.names = FALSE)
cat(sprintf("[coupling] task %d done: %d rows, %.1f min\n", TASK, nrow(out), (proc.time()[["elapsed"]] - t0) / 60))
print(aggregate(cbind(nmi_layer, nmi_joint_break, nmi_joint_continue, purity_break, completeness_break) ~ method,
                data = out, FUN = function(x) round(mean(x), 3)))
