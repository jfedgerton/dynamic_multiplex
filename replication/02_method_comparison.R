# =============================================================================
# replication/02_method_comparison.R
#
# UNIFIED HEAD-TO-HEAD METHOD COMPARISON (all methods, incl. Dynamic SBM).
#
# This is the exact code and seeds behind the Stage I method comparison in the
# paper. It merges what were two development scripts (the 10-method comparison
# and the separately-run Dynamic SBM) into one program with a single method
# registry, so every method is evaluated on the SAME simulated networks.
#
# Methods compared (11):
#   DynMux (Jaccard, Leiden)      - proposed, Jaccard coupling + Leiden
#   DynMux (Identity, Leiden)     - proposed, identity coupling + Leiden
#   DynMux (Jaccard)              - proposed, Jaccard coupling + Louvain
#   DynMux (Identity)             - proposed, identity coupling + Louvain
#   DynMux (filtered)             - proposed, two-layer rolling window (no
#                                   future information / leakage-free)
#   Independent Leiden            - per-layer Leiden, labels not linked
#   Independent Leiden (matched)  - per-layer Leiden + greedy adjacent-layer
#                                   label matching (the fair temporal baseline)
#   Aggregated Leiden             - single pooled partition over all layers
#   Multislice (Adjacent)         - Mucha et al. multislice modularity
#   multinet GLouvain             - multinet generalized Louvain
#   Dynamic SBM                   - dynsbm, fixed K across the series, ICL-
#                                   selected over Qmin..Qmax (generous range)
#
# Two arms (each a config grid; both run in every array task's config):
#   "borrow"     : layers individually near/below the per-layer detectability
#                  threshold but structure persistent across time. Tests
#                  whether joint detection borrows strength across layers.
#   "mergesplit" : the true number of communities CHANGES mid-series (6->3 or
#                  3->6). One fixed K cannot represent this; joint detection
#                  and per-layer K estimation should track it.
#
# Metrics per (config, rep, method):
#   nmi_layer : mean per-layer NMI vs truth (conventional per-snapshot metric)
#   nmi_joint : NMI of the full node x layer partition vs truth (credits
#               cross-layer label consistency)
#   k_mae     : mean |K_hat_t - K_true_t| (tracks time-varying K)
#   runtime_s
#
# Seeds are identical to the as-run study: seed = 9000 + cfg_i*1000 + rep.
# Dynamic SBM is ~1000x slower than the other methods; this is why the study
# runs one array task per config (all configs in parallel), with each task
# doing all reps x all methods. Fairness for Dynamic SBM: Qmax = 10 gives
# generous headroom above the true K (<= 6) in every cell; ICL selection via
# dynsbm:::compute.icl.
#
# Usage (smoke, one config, fast methods only unless dynsbm installed):
#   CMP_MINI=1 CMP_CFG=19 Rscript replication/02_method_comparison.R
# Array: one SLURM_ARRAY_TASK_ID per config (1..22). See slurm/02_*.sbatch.
# =============================================================================

set.seed(123)
options(dynamicmultiplex.skip_main = TRUE)
suppressPackageStartupMessages({ library(dynamicmultiplex); library(parallel) })

# Base method registry (the synthetic-study script defines base_methods and the
# per-layer / DynMux helper fitters used below). Sourced with the main block
# skipped so only the definitions load.
source("manuscript/01_synthetic_experiments.R", local = FALSE, echo = FALSE,
       max.deparse.length = Inf)

HAVE_DYNSBM <- requireNamespace("dynsbm", quietly = TRUE)

CFG_I <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", Sys.getenv("CMP_CFG", "1")))
CORES <- as.integer(Sys.getenv("CMP_CORES", "8"))
MINI  <- identical(Sys.getenv("CMP_MINI", "0"), "1")
REPS  <- if (MINI) 1L else 100L
QMIN  <- 2L; QMAX <- 10L; NSTART <- 10L
outdir <- Sys.getenv("CMP_OUTROOT", "replication/output/comparison")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# NEW METHODS  (greedy label matching; leakage-free filtered DynMux)
# =============================================================================
match_labels <- function(mems) {
  out <- mems
  for (t in 2:length(mems)) {
    prev <- out[[t - 1]]
    cur <- mems[[t]]
    tab <- table(cur, prev)
    cur_labs <- rownames(tab)
    ord <- order(-apply(tab, 1, max))
    mapping <- setNames(rep(NA_integer_, length(cur_labs)), cur_labs)
    used <- integer(0)
    for (i in ord) {
      cand <- as.integer(colnames(tab)[order(-tab[i, ])])
      cand <- cand[!(cand %in% used)]
      hit <- if (length(cand) > 0 && max(tab[i, !(as.integer(colnames(tab)) %in% used), drop = FALSE]) > 0)
        cand[1] else NA_integer_
      if (!is.na(hit)) { mapping[cur_labs[i]] <- hit; used <- c(used, hit) }
    }
    nxt <- max(c(prev, used, 0)) + 1L
    for (i in seq_along(mapping)) {
      if (is.na(mapping[i])) { mapping[i] <- nxt; nxt <- nxt + 1L }
    }
    out[[t]] <- unname(mapping[as.character(cur)])
  }
  out
}

method_independent_leiden_matched <- function(L) {
  match_labels(method_independent_leiden(L))
}

method_dynmux_filtered <- function(L) {
  # layer t fitted from layers (t-1, t) only: no future information
  T_ <- length(L)
  mems <- vector("list", T_)
  first <- extract_mem(fit_multilayer_jaccard(L[1:2], algorithm = "leiden"))
  mems[[1]] <- first[[1]]; mems[[2]] <- first[[2]]
  if (T_ >= 3) {
    for (t in 3:T_) {
      fit <- extract_mem(fit_multilayer_jaccard(L[(t - 1):t], algorithm = "leiden"))
      mems[[t]] <- fit[[2]]
    }
  }
  match_labels(mems)
}

# Dynamic SBM: T x N x N array -> ICL selection over Qmin..Qmax -> per-layer
# memberships. dynsbm holds one K fixed across the series (by design).
method_dynsbm <- function(L) {
  T_ <- length(L); n <- nrow(L[[1]])
  Y <- array(0, c(T_, n, n))
  for (t in seq_len(T_)) Y[t, , ] <- (L[[t]] > 0) * 1
  models <- dynsbm::select.dynsbm(Y, Qmin = QMIN, Qmax = QMAX,
                                  edge.type = "binary", nstart = NSTART,
                                  nb.cores = 1, plot = FALSE)
  icl <- vapply(models, function(md) dynsbm:::compute.icl(md), numeric(1))
  best <- models[[which.max(icl)]]
  mem <- best$membership                                   # n x T
  lapply(seq_len(T_), function(t) mem[, t])
}

METHODS <- list(
  "DynMux (Jaccard, Leiden)"     = base_methods[["DynMux (Jaccard, Leiden)"]],
  "DynMux (Identity, Leiden)"    = base_methods[["DynMux (Identity, Leiden)"]],
  "DynMux (Jaccard)"             = base_methods[["DynMux (Jaccard)"]],
  "DynMux (Identity)"            = base_methods[["DynMux (Identity)"]],
  "DynMux (filtered)"            = method_dynmux_filtered,
  "Independent Leiden"           = base_methods[["Independent Leiden"]],
  "Independent Leiden (matched)" = method_independent_leiden_matched,
  "Aggregated Leiden"            = base_methods[["Aggregated Leiden"]],
  "Multislice (Adjacent)"        = base_methods[["Multislice (Adjacent)"]],
  "multinet GLouvain"            = base_methods[["multinet GLouvain"]]
)
if (HAVE_DYNSBM) METHODS[["Dynamic SBM"]] <- method_dynsbm else
  message("dynsbm not installed; skipping Dynamic SBM. install.packages('dynsbm').")

# =============================================================================
# SIMULATORS  (identical to the as-run study)
# =============================================================================
sim_borrow <- function(n, K, p_in, p_out, p_switch, T_, seed) {
  set.seed(seed)
  memberships <- vector("list", T_)
  memberships[[1]] <- sample(seq_len(K), n, replace = TRUE)
  for (t in 2:T_) {
    prev <- memberships[[t - 1]]
    mask <- runif(n) < p_switch
    prev[mask] <- sample(seq_len(K), sum(mask), replace = TRUE)
    memberships[[t]] <- prev
  }
  layers <- lapply(memberships, function(mem) {
    P <- outer(mem, mem, "=="); Pr <- ifelse(P, p_in, p_out); diag(Pr) <- 0
    A <- matrix(0, n, n); up <- upper.tri(Pr)
    A[up] <- rbinom(sum(up), 1, Pr[up]); A + t(A)
  })
  list(layers = layers, truth = memberships)
}

sim_mergesplit <- function(n, T_, event, p_in, p_out, p_switch, seed) {
  set.seed(seed)
  t_star <- floor(T_ / 2) + 1L
  K_pre <- if (event == "merge") 6L else 3L
  memberships <- vector("list", T_)
  memberships[[1]] <- sample(seq_len(K_pre), n, replace = TRUE)
  for (t in 2:T_) {
    prev <- memberships[[t - 1]]
    if (t == t_star) {
      if (event == "merge") prev <- ceiling(prev / 2)
      else prev <- prev * 2L - rbinom(length(prev), 1, 0.5)
    }
    mask <- runif(n) < p_switch
    prev[mask] <- sample(sort(unique(prev)), sum(mask), replace = TRUE)
    memberships[[t]] <- prev
  }
  layers <- lapply(memberships, function(mem) {
    P <- outer(mem, mem, "=="); Pr <- ifelse(P, p_in, p_out); diag(Pr) <- 0
    A <- matrix(0, n, n); up <- upper.tri(Pr)
    A[up] <- rbinom(sum(up), 1, Pr[up]); A + t(A)
  })
  list(layers = layers, truth = memberships)
}

# =============================================================================
# METRICS  (identical to the as-run study)
# =============================================================================
eval_method <- function(det, truth) {
  nmi_layer <- mean(vapply(seq_along(det), function(t)
    igraph::compare(as.integer(det[[t]]), truth[[t]], method = "nmi"),
    numeric(1)))
  nmi_joint <- igraph::compare(as.integer(unlist(det)),
                               as.integer(unlist(truth)), method = "nmi")
  k_mae <- mean(vapply(seq_along(det), function(t)
    abs(length(unique(det[[t]])) - length(unique(truth[[t]]))), numeric(1)))
  c(nmi_layer = nmi_layer, nmi_joint = nmi_joint, k_mae = k_mae)
}

# =============================================================================
# CONFIG GRID  (18 borrow + 4 mergesplit = 22; identical to the as-run study)
# =============================================================================
cfgA <- expand.grid(arm = "borrow", n = 100, K = 5,
                    p_in = c(0.12, 0.15, 0.18), p_switch = c(0.02, 0.05),
                    T_ = c(10, 20, 40), event = "none",
                    stringsAsFactors = FALSE)
cfgA$p_out <- cfgA$p_in / 3
cfgB <- expand.grid(arm = "mergesplit", n = 100, K = NA,
                    p_in = 0.30, p_switch = 0.03, T_ = c(10, 20),
                    event = c("merge", "split"), stringsAsFactors = FALSE)
cfgB$p_out <- 0.05
cfgs <- rbind(cfgA, cfgB)
stopifnot(CFG_I >= 1, CFG_I <= nrow(cfgs))
cfg <- cfgs[CFG_I, ]

outfile <- file.path(outdir, sprintf("comparison_cfg%02d.csv", CFG_I))
if (nzchar(Sys.getenv("SLURM_ARRAY_TASK_ID")) && file.exists(outfile)) {
  cat("[compare skip] cfg", CFG_I, "already complete\n"); quit(save = "no")
}

cat(sprintf("[compare] cfg=%d/%d arm=%s T=%d event=%s p_in=%.2f sw=%.2f REPS=%d methods=%d dynsbm=%s\n",
            CFG_I, nrow(cfgs), cfg$arm, cfg$T_, cfg$event, cfg$p_in,
            cfg$p_switch, REPS, length(METHODS), HAVE_DYNSBM))

# =============================================================================
# RUN
# =============================================================================
run_rep <- function(rep) {
  seed <- 9000L + CFG_I * 1000L + rep                      # identical to as-run
  sim <- if (cfg$arm == "borrow") {
    sim_borrow(cfg$n, cfg$K, cfg$p_in, cfg$p_out, cfg$p_switch, cfg$T_, seed)
  } else {
    sim_mergesplit(cfg$n, cfg$T_, cfg$event, cfg$p_in, cfg$p_out,
                   cfg$p_switch, seed)
  }
  rows <- list(); ri <- 0L
  for (mname in names(METHODS)) {
    t0 <- proc.time()["elapsed"]
    det <- tryCatch(METHODS[[mname]](sim$layers), error = function(e) NULL)
    el <- proc.time()["elapsed"] - t0
    if (is.null(det)) next
    m <- eval_method(det, sim$truth)
    ri <- ri + 1L
    rows[[ri]] <- data.frame(
      arm = cfg$arm, n = cfg$n, K = cfg$K, p_in = cfg$p_in,
      p_out = cfg$p_out, p_switch = cfg$p_switch, T_layers = cfg$T_,
      event = cfg$event, rep = rep, method = mname,
      nmi_layer = round(m["nmi_layer"], 4),
      nmi_joint = round(m["nmi_joint"], 4),
      k_mae = round(m["k_mae"], 4),
      runtime_s = round(el, 3), stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

t0 <- proc.time()["elapsed"]
res <- mclapply(seq_len(REPS), function(r) {
  tryCatch(run_rep(r), error = function(e) {
    message(sprintf("rep %d failed: %s", r, conditionMessage(e))); NULL
  })
}, mc.cores = CORES, mc.preschedule = FALSE)
res <- do.call(rbind, Filter(Negate(is.null), res))
write.csv(res, outfile, row.names = FALSE)

agg <- aggregate(cbind(nmi_layer, nmi_joint, k_mae) ~ method, data = res,
                 FUN = mean)
cat(sprintf("[compare] cfg %d done: %d rows, %.1f min. means by method:\n",
            CFG_I, nrow(res), (proc.time()["elapsed"] - t0) / 60))
print(agg[order(-agg$nmi_joint), ], row.names = FALSE)
