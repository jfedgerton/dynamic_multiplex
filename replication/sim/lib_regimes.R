# =============================================================================
# replication/sim/lib_regimes.R  --  shared building blocks for the regime grid
#
# Sourced by sim/01 (method comparison), sim/05 (multislice omega sweep) and
# sim/06 (stability as a selection rule). Defines:
#   dens_params(), build_layer()                density parameterisation
#   sim_churnswitch/birthdeath/regimeshift/seasonality(), simulate_regime()
#   alive_idx(), AL(), EXPn(), dm()             open-population helpers, DynMux wrapper
#   method_pooled(), method_hungarian(), method_multinet(), method_dynsbm()
#   METHODS                                     named list: every method on one sim
#   eval_method()                               joint/layer NMI, K MAE, co-membership
#   regime_grid()                               the 72-configuration grid, shuffled
# Expects the package to be loaded (pkgload::load_all(file.path(ROOT, "r_code")))
# and library(igraph). Globals RHO, QMIN/QMAX/NSTART, HAVE_* are set here if absent.
# =============================================================================
if (!exists("RHO"))    RHO <- 0.10
if (!exists("QMIN"))   { QMIN <- 2L; QMAX <- 8L; NSTART <- 3L }
if (!exists("HAVE_DYNSBM"))   HAVE_DYNSBM   <- requireNamespace("dynsbm",   quietly = TRUE)
if (!exists("HAVE_MULTINET")) HAVE_MULTINET <- requireNamespace("multinet", quietly = TRUE)
if (!exists("HAVE_CLUE"))     HAVE_CLUE     <- requireNamespace("clue",     quietly = TRUE)

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

# Alive-node handling for open-population regimes (birth/death). The generator
# supplies sim$active (per-layer logical); coupling and tracking methods must
# see only the alive nodes in each layer, with absent nodes dropped and
# alive-isolates kept. For closed-population regimes sim$active is NULL and
# these helpers are identity operations. This replaces the former post-hoc
# 06_bd_fix.R + splice_bdfix.R patch so a single run produces correct output.
alive_idx <- function(sim, t) {
  if (!is.null(sim$active)) which(as.logical(sim$active[[t]]))
  else seq_len(nrow(sim$layers[[t]]))
}
AL <- function(sim) lapply(seq_along(sim$layers), function(t) {
  idx <- alive_idx(sim, t)
  A   <- sim$layers[[t]][idx, idx, drop = FALSE]
  dimnames(A) <- list(as.character(idx), as.character(idx))
  A
})
EXPn <- function(m, sim, n) lapply(seq_along(m), function(t) {
  idx <- alive_idx(sim, t)
  v <- rep(NA_integer_, n); v[idx] <- as.integer(m[[t]]); v
})
# DynMux wrapper: fit on alive layers, expand memberships back to n.
dm <- function(fitfun, alg, links = TRUE) {
  function(sim) {
    n <- nrow(sim$layers[[1]])
    # AL() names every layer's nodes by their index in the full population, so
    # birth-death layers with different alive sets are matched by name.
    f <- if (links) fitfun(AL(sim), algorithm = alg, layer_links = sim$links, allow_unequal_nodes = TRUE)
         else        fitfun(AL(sim), algorithm = alg, allow_unequal_nodes = TRUE)
    EXPn(extract_meta_membership(f), sim, n)
  }
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

# Optimal (Hungarian) label matching across consecutive layers, restricted to
# nodes alive in both layers when sim$active is present.
method_hungarian <- function(sim) {
  n <- nrow(sim$layers[[1]]); T_ <- length(sim$layers)
  mems <- lapply(sim$layers, leiden_layer)
  out  <- vector("list", T_)
  a1 <- alive_idx(sim, 1); v1 <- rep(NA_integer_, n); v1[a1] <- as.integer(mems[[1]][a1])
  out[[1]] <- v1; next_free <- max(v1, na.rm = TRUE) + 1L
  for (t in 2:T_) {
    prev <- out[[t - 1L]]; cur <- as.integer(mems[[t]])
    ac <- alive_idx(sim, t); ap <- alive_idx(sim, t - 1L)
    avc <- logical(n); avc[ac] <- TRUE
    avp <- logical(n); avp[ap] <- TRUE
    both <- avc & avp
    cur_labs  <- sort(unique(cur[avc]))
    prev_labs <- sort(unique(prev[avp & !is.na(prev)]))
    M <- matrix(0, length(cur_labs), length(prev_labs))
    for (i in seq_along(cur_labs)) for (j in seq_along(prev_labs))
      M[i, j] <- sum(cur == cur_labs[i] & prev == prev_labs[j] & both, na.rm = TRUE)
    d <- max(nrow(M), ncol(M)); Ms <- matrix(0, d, d)
    Ms[seq_len(nrow(M)), seq_len(ncol(M))] <- M
    asg <- if (HAVE_CLUE) as.integer(clue::solve_LSAP(max(Ms) - Ms)) else .greedy_assign(Ms)
    mp <- rep(NA_integer_, length(cur_labs))
    for (i in seq_along(cur_labs)) {
      co <- asg[i]
      if (co <= length(prev_labs) && M[i, co] > 0) mp[i] <- prev_labs[co]
    }
    for (i in seq_along(cur_labs))
      if (is.na(mp[i])) { mp[i] <- next_free; next_free <- next_free + 1L }
    next_free <- max(next_free, max(mp, na.rm = TRUE) + 1L)
    names(mp) <- as.character(cur_labs)
    v <- rep(NA_integer_, n); v[avc] <- mp[as.character(cur[avc])]
    out[[t]] <- v
  }
  out
}

# multinet generalized Louvain on the alive node set per layer. GUARDED:
# returns NULL if multinet is absent (recorded as NA).
method_multinet <- function(sim) {
  if (!HAVE_MULTINET) return(NULL)
  n <- nrow(sim$layers[[1]]); T_ <- length(sim$layers)
  nn <- multinet::ml_empty()
  for (t in seq_len(T_)) {
    idx <- alive_idx(sim, t); if (!length(idx)) next
    A <- sim$layers[[t]][idx, idx, drop = FALSE]
    g <- igraph::graph_from_adjacency_matrix(A, mode = "undirected", diag = FALSE)
    igraph::V(g)$name <- as.character(idx)
    multinet::add_igraph_layer_ml(nn, g, paste0("layer", t))
  }
  cm <- multinet::glouvain_ml(nn, gamma = 1, omega = 1)
  lapply(seq_len(T_), function(t) {
    s <- cm[cm$layer == paste0("layer", t), , drop = FALSE]
    v <- rep(NA_integer_, n)
    if (nrow(s)) v[as.integer(s$actor)] <- as.integer(as.factor(s$cid))
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
  "DynMux multislice (adjacent)"          = dm(fit_multilayer_identity_ties,    "leiden", links = FALSE),
  "DynMux multislice (custom)"            = dm(fit_multilayer_identity_ties,    "leiden", links = TRUE),
  "DynMux Jaccard"                        = dm(fit_multilayer_jaccard,          "leiden"),
  "DynMux Overlap"                        = dm(fit_multilayer_overlap,          "leiden"),
  "DynMux weighted Jaccard"               = dm(fit_multilayer_weighted_jaccard, "leiden"),
  "DynMux weighted Overlap"               = dm(fit_multilayer_weighted_overlap, "leiden"),
  "Pooled Leiden"                         = function(sim) method_pooled(sim, "leiden"),
  "Cross-sectional + Hungarian"           = method_hungarian,
  "multinet GLouvain"                     = method_multinet,
  "Dynamic SBM"                           = method_dynsbm,
  # ---- appendix: Louvain-algorithm variants -------------------------------
  "DynMux multislice (adjacent, Louvain)" = dm(fit_multilayer_identity_ties,    "louvain", links = FALSE),
  "DynMux Jaccard (Louvain)"              = dm(fit_multilayer_jaccard,          "louvain"),
  "DynMux Overlap (Louvain)"              = dm(fit_multilayer_overlap,          "louvain"),
  "DynMux weighted Jaccard (Louvain)"     = dm(fit_multilayer_weighted_jaccard, "louvain"),
  "DynMux weighted Overlap (Louvain)"     = dm(fit_multilayer_weighted_overlap, "louvain"),
  "Pooled Louvain"                        = function(sim) method_pooled(sim, "louvain")
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
regime_grid <- function() {
  cfgs <- expand.grid(
    regime    = c("seasonality", "churnswitch", "birthdeath", "regimeshift"),
    n         = c(50L, 100L, 200L),
    r         = c(1.5, 3, 6),
    intensity = c("low", "high"),
    stringsAsFactors = FALSE)
  stopifnot(nrow(cfgs) == 72L)
  set.seed(123); ord <- sample(nrow(cfgs))
  cfgs[ord, ]
}
