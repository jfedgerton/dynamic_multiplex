# =============================================================================
# replication/extended/09_dynamic_regimes.R
#
# HEAD-TO-HEAD METHOD COMPARISON ACROSS FOUR DYNAMIC-COMMUNITY REGIMES,
# focused on RECURRING and CHANGE-POINT structure (where DynMux's interlayer
# coupling should pay off) with a transition-layer metric (nmi_change).
#
# Four synthetic generators exercise qualitatively different kinds of dynamic
# community structure, and every available method is scored on the SAME
# simulated networks, on the TRACKED (meta / cross-layer) partition:
#
#   seasonality : a small bank of latent partitions recur on a period, each
#                 lightly drifting when revisited (period-lagged coupling) --
#                 DynMux's designed win. Ported from 06 so it is scored with the
#                 same specs (incl. weighted) and metrics (incl. nmi_change).
#   churnswitch : edge churn PLUS Markov switching between latent partitions.
#   birthdeath  : open population where whole communities are born / die over
#                 the series (membership enters and exits existence).
#   regimeshift : an abrupt CHANGE-POINT at t* -- era-1 K1 communities give way
#                 to era-2 K2 DISJOINT communities, so a single pooled partition
#                 cannot fit both eras (the Cold-War-end / bipolar->multipolar
#                 case).
#
# Density parameterization (all regimes): given separation ratio r, target
# density rho = 0.10 and community count K,
#     p_in  = rho * r * K / (r + K - 1)
#     p_out = p_in / r                         (both clipped to [0, 1]).
#
# CONFIG GRID (72): regime x N{50,100,200} x r{1.5,3,6} x intensity{low,high}.
# The SLURM array index selects a config AFTER a fixed shuffle so that a
# partial run (early array indices) already samples every regime / N / r /
# intensity condition:
#     set.seed(123); ord <- sample(nrow(cfgs)); cfg <- cfgs[ord[TASK], ]
#
# Methods (named list METHODS, each takes the sim object; ALG = "leiden" for
# the main set, "louvain" for the appendix set):
#   DynMux multislice (adjacent)   identity coupling, default adjacent links
#   DynMux multislice (custom)     identity coupling, generator's layer_links
#   DynMux Jaccard                 Jaccard coupling + generator's layer_links
#   DynMux Overlap                 overlap coupling + generator's layer_links
#   DynMux weighted Jaccard        degree-weighted Jaccard coupling
#   DynMux weighted Overlap        degree-weighted overlap coupling
#   Pooled Leiden                  one partition on Reduce("+", L), replicated
#   Cross-sectional + Hungarian    per-layer Leiden + optimal (Hungarian) label
#                                  matching across consecutive layers
#   multinet GLouvain              multinet generalized Louvain (GUARDED)
#   Dynamic SBM                    dynsbm, ICL-selected over Qmin..Qmax (GUARDED)
#   ... plus Louvain-algorithm appendix variants ("(Louvain)").
#
# Metrics per (config, rep, method), computed on ACTIVE nodes where the
# generator supplies sim$active, else on all nodes:
#   nmi_layer         mean per-layer NMI vs truth
#   nmi_joint         NMI of the full node x layer partition vs truth
#   k_mae             mean |K_hat_t - K_true_t|
#   comembership_acc  mean per-layer fraction of node-pairs whose same/different
#                     -community status matches truth
#   runtime_s         wall-clock (system.time elapsed) of the method call
#   mean_n_comm       mean per-layer number of detected communities
#   total_n_comm      distinct communities across all layers (tracked/meta count)
#
# Reps: fast methods 100 (2 when CMP_MINI=1), Dynamic SBM 50 (2 when CMP_MINI).
# Per-rep seed = 9000 + TASK*1000 + rep. Each rep simulates once and runs every
# available method on that one sim.
#
# Usage (local smoke, one config, fast methods + dynsbm if installed):
#   CMP_MINI=1 CMP_CFG=1 Rscript replication/extended/09_dynamic_regimes.R
# Array: one SLURM_ARRAY_TASK_ID per config (1..72). See slurm/09_dynamic_regimes.sbatch.
# The code is intentionally sequential (a plain rep loop) to stay debuggable.
# =============================================================================

suppressMessages(pkgload::load_all("r_code", quiet = TRUE))
suppressPackageStartupMessages(library(igraph))

HAVE_DYNSBM  <- requireNamespace("dynsbm",  quietly = TRUE)
HAVE_MULTINET <- requireNamespace("multinet", quietly = TRUE)
HAVE_CLUE    <- requireNamespace("clue",    quietly = TRUE)

TASK <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID",
                              Sys.getenv("CMP_CFG", "1")))
MINI <- identical(Sys.getenv("CMP_MINI", "0"), "1")
REPS_FAST   <- if (MINI) 2L else 100L
REPS_DYNSBM <- if (MINI) 2L else 50L
QMIN <- 2L; QMAX <- 8L; NSTART <- 3L
RHO  <- 0.10

outdir <- Sys.getenv("DYN_OUTROOT", "replication/extended/output/dynamic")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# SHARED BUILDING BLOCKS
# =============================================================================

# Density from separation ratio r and community count K (clipped to [0, 1]).
dens_params <- function(r, K, rho = RHO) {
  p_in  <- rho * r * K / (r + K - 1)
  p_out <- p_in / r
  list(p_in  = min(max(p_in,  0), 1),
       p_out = min(max(p_out, 0), 1))
}

# A single undirected 0/1 symmetric SBM layer from a membership vector.
build_layer <- function(m, p_in, p_out) {
  n  <- length(m)
  P  <- outer(m, m, "==")
  Pr <- ifelse(P, p_in, p_out); diag(Pr) <- 0
  A  <- matrix(0, n, n); up <- upper.tri(Pr)
  A[up] <- rbinom(sum(up), 1, Pr[up]); A + t(A)
}

# =============================================================================
# GENERATORS  (each returns list(layers, truth, links[, active]))
# =============================================================================

# 1. churnswitch: edgechurn + Markov membership SWITCHING. Each layer is a
#    fresh SBM draw (edge churn) from the CURRENT membership, and a fraction of
#    nodes switch community every layer -> the partition genuinely changes over
#    time (unlike plain edgechurn, whose membership is static).
sim_churnswitch <- function(n, K, r, intensity, seed) {
  set.seed(seed)
  dp <- dens_params(r, K)
  T_       <- if (identical(intensity, "high")) 14L else 8L
  p_switch <- if (identical(intensity, "high")) 0.15 else 0.05
  mem <- sample(seq_len(K), n, replace = TRUE)
  truth <- vector("list", T_); layers <- vector("list", T_)
  for (t in seq_len(T_)) {
    if (t > 1L) { sw <- runif(n) < p_switch
      if (any(sw)) mem[sw] <- sample(seq_len(K), sum(sw), replace = TRUE) }
    layers[[t]] <- build_layer(mem, dp$p_in, dp$p_out)
    truth[[t]]  <- mem
  }
  links <- data.frame(from = seq_len(T_ - 1L), to = 2:T_, weight = 1)
  list(layers = layers, truth = truth, links = links)
}

# 2. birthdeath: openpop node turnover + community BIRTH/DEATH. Every 3 layers
#    the oldest community dies (members dispersed to living communities) and a
#    new community is born (recruits a fresh cohort). Communities therefore
#    exist only within time windows -- a single pooled partition cannot
#    represent them. Metrics use active nodes (sim$active).
sim_birthdeath <- function(n, K0, r, intensity, seed) {
  set.seed(seed)
  dp <- dens_params(r, K0)
  T_      <- 12L
  frac    <- if (identical(intensity, "high")) 0.30 else 0.10   # node turnover
  recruit <- if (identical(intensity, "high")) 0.25 else 0.15   # birth cohort
  mem <- sample(seq_len(K0), n, replace = TRUE)
  alive <- seq_len(K0); nxt <- K0 + 1L; state <- rep(TRUE, n)
  truth <- vector("list", T_); layers <- vector("list", T_); active <- vector("list", T_)
  for (t in seq_len(T_)) {
    if (t > 1L) {
      toggle <- sample(seq_len(n), floor(frac * n)); state[toggle] <- !state[toggle]
      if (t %% 3L == 0L) {
        dead <- alive[1]; alive <- alive[-1]
        mv <- which(mem == dead)
        if (length(mv) && length(alive)) mem[mv] <- sample(alive, length(mv), replace = TRUE)
        born <- nxt; nxt <- nxt + 1L; alive <- c(alive, born)
        rec <- sample(seq_len(n), floor(recruit * n)); mem[rec] <- born
      }
    }
    act <- state; active[[t]] <- act
    A <- matrix(0, n, n); ai <- which(act)
    if (length(ai) >= 2L) {
      ma <- mem[ai]; P <- outer(ma, ma, "=="); Pr <- ifelse(P, dp$p_in, dp$p_out); diag(Pr) <- 0
      sub <- matrix(0, length(ai), length(ai)); up <- upper.tri(Pr)
      sub[up] <- rbinom(sum(up), 1, Pr[up]); sub <- sub + t(sub); A[ai, ai] <- sub
    }
    layers[[t]] <- A
    tr <- mem; inact <- which(!act); if (length(inact)) tr[inact] <- max(mem) + seq_along(inact)
    truth[[t]] <- tr
  }
  links <- data.frame(from = seq_len(T_ - 1L), to = 2:T_, weight = 1)
  list(layers = layers, truth = truth, links = links, active = active)
}

# 3. regimeshift: an abrupt CHANGE-POINT where the partition fully
#    reconfigures at t*. Era 1 (layers 1..t*) uses K1 communities; era 2
#    (t*+1..T) uses K2 *disjoint* communities (new ids) -- the two eras share
#    no community, so a single pooled partition cannot fit both. Clean SBM
#    draws (no edge churn); small within-era drift. This is the Cold-War-end /
#    bipolar->multipolar case, and the setting where coupling that can bridge
#    or break at t* should beat pooling on aggregate metrics.
sim_regimeshift <- function(n, r, intensity, seed) {
  set.seed(seed)
  K1 <- if (identical(intensity, "high")) 2L else 3L
  K2 <- if (identical(intensity, "high")) 6L else 4L
  dp <- dens_params(r, max(K1, K2))
  T_ <- 12L; tstar <- 6L; drift <- 0.03
  m1 <- sample(seq_len(K1), n, replace = TRUE)
  m2 <- sample(K1 + seq_len(K2), n, replace = TRUE)   # disjoint era-2 community ids
  truth <- vector("list", T_); layers <- vector("list", T_)
  cur <- m1
  for (t in seq_len(T_)) {
    if (t == tstar + 1L) {
      cur <- m2                                        # abrupt full reconfiguration
    } else if (t > 1L) {
      dr <- runif(n) < drift
      if (any(dr)) {
        pool <- if (t <= tstar) seq_len(K1) else K1 + seq_len(K2)
        cur[dr] <- sample(pool, sum(dr), replace = TRUE)
      }
    }
    layers[[t]] <- build_layer(cur, dp$p_in, dp$p_out)
    truth[[t]]  <- cur
  }
  links <- data.frame(from = seq_len(T_ - 1L), to = 2:T_, weight = 1)
  list(layers = layers, truth = truth, links = links)
}

# 4. seasonality: a small bank of latent partitions RECUR on a period (each
#    era returns to a prior configuration), with light drift on each visit and
#    period-lagged (non-adjacent) interlayer links -- DynMux's designed use
#    case, where similarity coupling that bridges non-adjacent recurrences
#    beats pooling and adjacent-only methods. Ported from 06 so seasonality is
#    scored with the same specs (incl. weighted) and metrics (incl. nmi_change)
#    as the change-point regimes.
sim_seasonality <- function(n, K, period, r, seed) {
  set.seed(seed)
  dp <- dens_params(r, K)
  cur <- lapply(seq_len(period), function(s) sample(seq_len(K), n, replace = TRUE))
  T_  <- 4L * period
  truth <- vector("list", T_); layers <- vector("list", T_)
  for (t in seq_len(T_)) {
    s <- ((t - 1L) %% period) + 1L
    mask <- runif(n) < 0.02                       # light drift each visit
    if (any(mask)) cur[[s]][mask] <- sample(seq_len(K), sum(mask), replace = TRUE)
    m <- cur[[s]]
    truth[[t]]  <- m
    layers[[t]] <- build_layer(m, dp$p_in, dp$p_out)
  }
  tt <- (period + 1L):T_
  links <- data.frame(from = tt - period, to = tt, weight = 1)
  list(layers = layers, truth = truth, links = links)
}

# Dispatch a config row to its generator (intensity / N / r wired per regime).
simulate_regime <- function(cfg, seed) {
  switch(cfg$regime,
    seasonality = sim_seasonality(cfg$n, K = 4L,
                                  period = if (cfg$intensity == "high") 4L else 2L,
                                  r = cfg$r, seed = seed),
    churnswitch = sim_churnswitch(cfg$n, K = 4L, r = cfg$r, intensity = cfg$intensity, seed = seed),
    birthdeath  = sim_birthdeath (cfg$n, K0 = 4L, r = cfg$r, intensity = cfg$intensity, seed = seed),
    regimeshift = sim_regimeshift(cfg$n, r = cfg$r, intensity = cfg$intensity, seed = seed),
    stop("unknown regime: ", cfg$regime)
  )
}

# =============================================================================
# METHOD HELPERS
# =============================================================================

# Per-layer Leiden (modularity objective, weighted). Isolated / empty layers
# collapse to a single community.
leiden_layer <- function(mat) {
  g <- igraph::graph_from_adjacency_matrix(mat, mode = "undirected",
                                           weighted = TRUE, diag = FALSE)
  if (igraph::ecount(g) == 0L) return(rep(1L, nrow(mat)))
  as.integer(igraph::membership(
    igraph::cluster_leiden(g, objective_function = "modularity",
                           weights = igraph::E(g)$weight)))
}

# Pooled: one partition on the summed adjacency, replicated across layers.
method_pooled <- function(sim, alg) {
  L   <- sim$layers
  agg <- Reduce("+", L)
  g   <- igraph::graph_from_adjacency_matrix(agg, mode = "undirected",
                                             weighted = TRUE, diag = FALSE)
  mem <- if (identical(alg, "leiden")) {
    igraph::membership(igraph::cluster_leiden(g, objective_function = "modularity",
                                              weights = igraph::E(g)$weight))
  } else {
    igraph::membership(igraph::cluster_louvain(g, weights = igraph::E(g)$weight))
  }
  replicate(length(L), as.integer(mem), simplify = FALSE)
}

# Greedy fallback assignment (used only if clue is not installed): for each row
# in decreasing best-overlap order, take its best still-free column.
.greedy_assign <- function(M) {
  nr <- nrow(M); asg <- rep(NA_integer_, nr); used <- integer(0)
  for (i in order(-apply(M, 1, max))) {
    cand <- order(-M[i, ]); cand <- cand[!(cand %in% used)]
    asg[i] <- cand[1]; used <- c(used, cand[1])
  }
  asg
}

# Optimal (Hungarian) label matching across consecutive layers: each layer's
# labels are relabelled to maximise overlap with the previous (already tracked)
# layer, so that persistent communities keep a stable id over time.
match_labels_hungarian <- function(mems) {
  out <- vector("list", length(mems))
  out[[1]] <- as.integer(mems[[1]])
  next_free <- max(out[[1]]) + 1L
  for (t in 2:length(mems)) {
    prev <- out[[t - 1L]]
    cur  <- as.integer(mems[[t]])
    cur_labs  <- sort(unique(cur))
    prev_labs <- sort(unique(prev))
    M <- matrix(0, length(cur_labs), length(prev_labs))
    for (i in seq_along(cur_labs))
      for (j in seq_along(prev_labs))
        M[i, j] <- sum(cur == cur_labs[i] & prev == prev_labs[j])
    d   <- max(nrow(M), ncol(M))
    Msq <- matrix(0, d, d); Msq[seq_len(nrow(M)), seq_len(ncol(M))] <- M
    asg <- if (HAVE_CLUE) as.integer(clue::solve_LSAP(max(Msq) - Msq))
           else .greedy_assign(Msq)
    mapping <- rep(NA_integer_, length(cur_labs))
    for (i in seq_along(cur_labs)) {
      col <- asg[i]
      if (col <= length(prev_labs) && M[i, col] > 0) mapping[i] <- prev_labs[col]
    }
    for (i in seq_along(cur_labs))
      if (is.na(mapping[i])) { mapping[i] <- next_free; next_free <- next_free + 1L }
    next_free <- max(next_free, max(mapping) + 1L)
    names(mapping) <- as.character(cur_labs)
    out[[t]] <- unname(mapping[as.character(cur)])
  }
  lapply(out, as.integer)
}

# multinet generalized Louvain. GUARDED: returns NULL if multinet is absent
# (recorded as NA). Untested in this sandbox (multinet not installed here).
method_multinet <- function(sim) {
  if (!HAVE_MULTINET) return(NULL)
  L <- sim$layers; n <- nrow(L[[1]]); T_ <- length(L)
  net    <- multinet::ml_empty()
  actors <- as.character(seq_len(n))
  for (t in seq_len(T_)) {
    g <- igraph::graph_from_adjacency_matrix(L[[t]], mode = "undirected", diag = FALSE)
    igraph::V(g)$name <- actors
    multinet::add_igraph_layer_ml(net, g, paste0("layer", t))
  }
  comm <- multinet::glouvain_ml(net, gamma = 1, omega = 1)
  # comm: data.frame(actor, layer, cid) -> per-layer integer membership vectors.
  lapply(seq_len(T_), function(t) {
    sub <- comm[comm$layer == paste0("layer", t), , drop = FALSE]
    v <- rep(NA_integer_, n)
    v[as.integer(sub$actor)] <- as.integer(as.factor(sub$cid))
    # actors absent from this layer's rows fall back to a fresh singleton id.
    if (anyNA(v)) v[is.na(v)] <- max(0L, v, na.rm = TRUE) + seq_len(sum(is.na(v)))
    v
  })
}

# Dynamic SBM on binarized layers; ICL selection over Qmin..Qmax. GUARDED.
# For open-population regimes the per-layer active mask is passed to dynsbm as
# its `present` matrix (n x T), so isolate/inactive nodes are handled rather
# than tripping dynsbm's "node never present" check.
method_dynsbm <- function(sim) {
  if (!HAVE_DYNSBM) return(NULL)
  L <- sim$layers; T_ <- length(L); n <- nrow(L[[1]])
  Y <- array(0, c(T_, n, n))
  for (t in seq_len(T_)) Y[t, , ] <- (L[[t]] > 0) * 1
  present <- NULL
  if (!is.null(sim$active)) {
    present <- vapply(sim$active, function(a) as.integer(a), integer(n))  # n x T
  }
  models <- dynsbm::select.dynsbm(Y, present = present, Qmin = QMIN, Qmax = QMAX,
                                  edge.type = "binary", nstart = NSTART,
                                  nb.cores = 1, plot = FALSE)
  icl  <- vapply(models, function(md) dynsbm:::compute.icl(md), numeric(1))
  best <- models[[which.max(icl)]]
  mem  <- best$membership                     # n x T
  lapply(seq_len(T_), function(t) as.integer(mem[, t]))
}

# METHODS registry. Each entry is function(sim) -> list of per-layer integer
# membership vectors (or NULL if unavailable). Main set uses Leiden; the
# appendix set at the end uses Louvain.
METHODS <- list(
  "DynMux multislice (adjacent)" = function(sim)
    extract_meta_membership(fit_multilayer_identity_ties(sim$layers, algorithm = "leiden")),
  "DynMux multislice (custom)"   = function(sim)
    extract_meta_membership(fit_multilayer_identity_ties(sim$layers, algorithm = "leiden",
                                                         layer_links = sim$links)),
  "DynMux Jaccard"               = function(sim)
    extract_meta_membership(fit_multilayer_jaccard(sim$layers, algorithm = "leiden",
                                                   layer_links = sim$links)),
  "DynMux Overlap"               = function(sim)
    extract_meta_membership(fit_multilayer_overlap(sim$layers, algorithm = "leiden",
                                                   layer_links = sim$links)),
  "DynMux weighted Jaccard"      = function(sim)
    extract_meta_membership(fit_multilayer_weighted_jaccard(sim$layers, algorithm = "leiden",
                                                            layer_links = sim$links)),
  "DynMux weighted Overlap"      = function(sim)
    extract_meta_membership(fit_multilayer_weighted_overlap(sim$layers, algorithm = "leiden",
                                                            layer_links = sim$links)),
  "Pooled Leiden"                = function(sim) method_pooled(sim, "leiden"),
  "Cross-sectional + Hungarian"  = function(sim)
    match_labels_hungarian(lapply(sim$layers, leiden_layer)),
  "multinet GLouvain"            = function(sim) method_multinet(sim),
  "Dynamic SBM"                  = function(sim) method_dynsbm(sim),
  # ---- appendix: Louvain-algorithm variants -------------------------------
  "DynMux multislice (adjacent, Louvain)" = function(sim)
    extract_meta_membership(fit_multilayer_identity_ties(sim$layers, algorithm = "louvain")),
  "DynMux Jaccard (Louvain)"     = function(sim)
    extract_meta_membership(fit_multilayer_jaccard(sim$layers, algorithm = "louvain",
                                                   layer_links = sim$links)),
  "DynMux Overlap (Louvain)"     = function(sim)
    extract_meta_membership(fit_multilayer_overlap(sim$layers, algorithm = "louvain",
                                                   layer_links = sim$links)),
  "DynMux weighted Jaccard (Louvain)" = function(sim)
    extract_meta_membership(fit_multilayer_weighted_jaccard(sim$layers, algorithm = "louvain",
                                                            layer_links = sim$links)),
  "DynMux weighted Overlap (Louvain)" = function(sim)
    extract_meta_membership(fit_multilayer_weighted_overlap(sim$layers, algorithm = "louvain",
                                                            layer_links = sim$links)),
  "Pooled Louvain"               = function(sim) method_pooled(sim, "louvain")
)

# =============================================================================
# METRICS  (computed on ACTIVE nodes where sim$active is present, else all)
# =============================================================================
eval_method <- function(det, sim) {
  T_ <- length(det)
  has_active <- !is.null(sim$active)
  per <- lapply(seq_len(T_), function(t) {
    idx <- if (has_active) which(sim$active[[t]]) else seq_along(det[[t]])
    list(d = as.integer(det[[t]][idx]), tr = as.integer(sim$truth[[t]][idx]))
  })
  nmi_layer <- mean(vapply(per, function(x)
    if (length(x$d) < 2L) NA_real_ else igraph::compare(x$d, x$tr, method = "nmi"),
    numeric(1)), na.rm = TRUE)
  d_all  <- unlist(lapply(per, `[[`, "d"))
  tr_all <- unlist(lapply(per, `[[`, "tr"))
  nmi_joint <- igraph::compare(d_all, tr_all, method = "nmi")
  k_mae <- mean(vapply(per, function(x)
    abs(length(unique(x$d)) - length(unique(x$tr))), numeric(1)))
  comemb <- mean(vapply(per, function(x) {
    if (length(x$d) < 2L) return(NA_real_)
    Dd <- outer(x$d, x$d, "=="); Dt <- outer(x$tr, x$tr, "==")
    up <- upper.tri(Dd); mean(Dd[up] == Dt[up])
  }, numeric(1)), na.rm = TRUE)
  mean_n_comm  <- mean(vapply(per, function(x) length(unique(x$d)), numeric(1)))
  total_n_comm <- length(unique(d_all))
  # transition-layer NMI: per-layer NMI averaged over layers where >5% of
  # (common active) nodes changed community vs the previous layer -- i.e. where
  # the structure actually reconfigures. NA if the truth never changes.
  chg <- vapply(seq_len(T_), function(t) {
    if (t == 1L) return(0)
    it <- if (has_active) which(sim$active[[t]])     else seq_along(sim$truth[[t]])
    ip <- if (has_active) which(sim$active[[t - 1L]]) else seq_along(sim$truth[[t - 1L]])
    cm <- intersect(it, ip); if (length(cm) < 2L) return(0)
    mean(as.integer(sim$truth[[t]][cm]) != as.integer(sim$truth[[t - 1L]][cm]))
  }, numeric(1))
  trans <- which(chg > 0.05)
  nmi_change <- if (!length(trans)) NA_real_ else
    mean(vapply(trans, function(t) {
      idx <- if (has_active) which(sim$active[[t]]) else seq_along(det[[t]])
      if (length(idx) < 2L) NA_real_ else
        igraph::compare(as.integer(det[[t]][idx]), as.integer(sim$truth[[t]][idx]), method = "nmi")
    }, numeric(1)), na.rm = TRUE)
  c(nmi_layer = nmi_layer, nmi_joint = nmi_joint, nmi_change = nmi_change, k_mae = k_mae,
    comembership_acc = comemb, mean_n_comm = mean_n_comm, total_n_comm = total_n_comm)
}

# =============================================================================
# CONFIG GRID  (72) + fixed shuffle so early array indices sample everything
# =============================================================================
cfgs <- expand.grid(
  regime    = c("seasonality", "churnswitch", "birthdeath", "regimeshift"),
  n         = c(50L, 100L, 200L),
  r         = c(1.5, 3, 6),
  intensity = c("low", "high"),
  stringsAsFactors = FALSE
)
stopifnot(nrow(cfgs) == 72L)

set.seed(123)
ord <- sample(nrow(cfgs))
stopifnot(TASK >= 1L, TASK <= nrow(cfgs))
cfg <- cfgs[ord[TASK], ]

outfile <- file.path(outdir, sprintf("dyn_cfg%02d.csv", TASK))
if (nzchar(Sys.getenv("SLURM_ARRAY_TASK_ID")) && file.exists(outfile)) {
  cat("[dyn skip] task", TASK, "already complete:", outfile, "\n")
  quit(save = "no")
}

if (!HAVE_MULTINET)
  message("[note] multinet not installed; 'multinet GLouvain' recorded as NA. ",
          "install.packages('multinet') on ROAR to enable it.")
if (!HAVE_DYNSBM)
  message("[note] dynsbm not installed; 'Dynamic SBM' skipped. ",
          "install.packages('dynsbm') on ROAR to enable it.")

cat(sprintf(paste0("[dyn] task=%d/%d -> regime=%s n=%d r=%.1f intensity=%s | ",
                   "REPS fast=%d dynsbm=%d | clue=%s dynsbm=%s multinet=%s\n"),
            TASK, nrow(cfgs), cfg$regime, cfg$n, cfg$r, cfg$intensity,
            REPS_FAST, REPS_DYNSBM, HAVE_CLUE, HAVE_DYNSBM, HAVE_MULTINET))

# =============================================================================
# RUN  (sequential: one rep at a time, all methods on the same sim)
# =============================================================================
metric_cols <- c("nmi_layer", "nmi_joint", "nmi_change", "k_mae", "comembership_acc",
                 "mean_n_comm", "total_n_comm")

run_rep <- function(rep) {
  seed <- 9000L + TASK * 1000L + rep
  sim  <- simulate_regime(cfg, seed)
  # Dynamic SBM is far slower and runs only on the first REPS_DYNSBM reps.
  names_this <- names(METHODS)
  if (rep > REPS_DYNSBM) names_this <- setdiff(names_this, "Dynamic SBM")
  if (!HAVE_DYNSBM)      names_this <- setdiff(names_this, "Dynamic SBM")
  rows <- vector("list", length(names_this)); ri <- 0L
  for (mname in names_this) {
    t0  <- proc.time()["elapsed"]
    det <- tryCatch(METHODS[[mname]](sim), error = function(e) {
      message(sprintf("  [rep %d] method '%s' errored: %s",
                      rep, mname, conditionMessage(e))); NULL
    })
    el  <- as.numeric(proc.time()["elapsed"] - t0)
    m   <- if (is.null(det)) setNames(rep(NA_real_, length(metric_cols)), metric_cols)
           else eval_method(det, sim)
    ri  <- ri + 1L
    rows[[ri]] <- data.frame(
      regime = cfg$regime, n = cfg$n, r = cfg$r, intensity = cfg$intensity,
      rep = rep, method = mname,
      nmi_layer        = round(m[["nmi_layer"]], 4),
      nmi_joint        = round(m[["nmi_joint"]], 4),
      nmi_change       = round(m[["nmi_change"]], 4),
      k_mae            = round(m[["k_mae"]], 4),
      comembership_acc = round(m[["comembership_acc"]], 4),
      runtime_s        = round(el, 3),
      mean_n_comm      = round(m[["mean_n_comm"]], 4),
      total_n_comm     = m[["total_n_comm"]],
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

t0  <- proc.time()["elapsed"]
all_rows <- vector("list", REPS_FAST)
for (rep in seq_len(REPS_FAST)) {
  all_rows[[rep]] <- tryCatch(run_rep(rep), error = function(e) {
    message(sprintf("rep %d failed: %s", rep, conditionMessage(e))); NULL
  })
  if (rep %% 10L == 0L || MINI)
    cat(sprintf("  ... rep %d/%d done (%.1f min elapsed)\n",
                rep, REPS_FAST, (proc.time()["elapsed"] - t0) / 60))
}
res <- do.call(rbind, Filter(Negate(is.null), all_rows))
write.csv(res, outfile, row.names = FALSE)

agg <- aggregate(cbind(nmi_layer, nmi_joint, comembership_acc) ~ method,
                 data = res, FUN = function(x) mean(x, na.rm = TRUE))
cat(sprintf("[dyn] task %d done: %d rows in %.1f min -> %s\n",
            TASK, nrow(res), (proc.time()["elapsed"] - t0) / 60, outfile))
print(agg[order(-agg$nmi_joint), ], row.names = FALSE)
