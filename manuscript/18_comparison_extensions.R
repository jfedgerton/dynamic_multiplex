# =============================================================================
# manuscript/18_comparison_extensions.R
#
# Stage I extensions (Jared, 2026-07-21/22):
#   ARM A "borrowing strength": layers individually near/below the per-layer
#         detectability threshold but structure persistent across time. The
#         documented regime where JOINT multilayer detection beats per-layer
#         detection + post-hoc label matching (threshold drops ~ sqrt(T)).
#   ARM B "merge/split": true number of communities CHANGES mid-series.
#         One-to-one label matching has no coherent answer; joint detection
#         and per-layer K estimation should track it. (Also the regime a
#         fixed-K dynamic SBM cannot represent.)
#
# New methods added to the registry from 01_synthetic_experiments.R:
#   - Independent Leiden (matched): per-layer Leiden + greedy max-overlap
#     label matching between adjacent layers (the fair temporal baseline).
#   - DynMux (filtered, Jaccard Leiden): two-layer rolling-window fits so
#     layer t uses only information from layers <= t (no temporal leakage).
#
# Metrics per (arm, config, rep, method):
#   nmi_layer : mean per-layer NMI vs truth (the conventional metric)
#   nmi_joint : NMI of the FULL node x layer partition vs truth - credits
#               cross-layer label consistency; per-layer-arbitrary labels
#               score poorly here by construction
#   k_mae     : mean |K_hat_t - K_true_t| (tracks time-varying K in ARM B)
#   runtime_s
#
# Usage: COMP_MINI=1 Rscript manuscript/18_comparison_extensions.R  (smoke)
#        Rscript manuscript/18_comparison_extensions.R              (full)
# =============================================================================

set.seed(123)
options(dynamicmultiplex.skip_main = TRUE)
source("manuscript/_setup_parallel.R")
source("manuscript/01_synthetic_experiments.R", local = FALSE, echo = FALSE,
       max.deparse.length = Inf)
suppressPackageStartupMessages(library(parallel))

MINI  <- identical(Sys.getenv("COMP_MINI", "0"), "1")
REPS  <- if (MINI) 1L else 10L
CORES <- as.integer(Sys.getenv("COMP_CORES", "8"))
outdir <- "manuscript/output"

# --------------------- new method: greedy label matching ---------------------
match_labels <- function(mems) {
  out <- mems
  for (t in 2:length(mems)) {
    prev <- out[[t - 1]]
    cur <- mems[[t]]
    tab <- table(cur, prev)
    cur_labs <- rownames(tab)
    ord <- order(-apply(tab, 1, max))          # biggest communities first
    mapping <- setNames(rep(NA_integer_, length(cur_labs)), cur_labs)
    used <- integer(0)
    for (i in ord) {
      cand <- as.integer(colnames(tab)[order(-tab[i, ])])
      cand <- cand[!(cand %in% used)]
      hit <- if (length(cand) > 0 && max(tab[i, !(as.integer(colnames(tab)) %in% used), drop = FALSE]) > 0)
        cand[1] else NA_integer_
      if (!is.na(hit)) {
        mapping[cur_labs[i]] <- hit
        used <- c(used, hit)
      }
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
  # layer t fitted from layers (t-1, t) only: no future information.
  T_ <- length(L)
  mems <- vector("list", T_)
  first <- extract_mem(fit_multilayer_jaccard(L[1:2], algorithm = "leiden"))
  mems[[1]] <- first[[1]]
  mems[[2]] <- first[[2]]
  if (T_ >= 3) {
    for (t in 3:T_) {
      fit <- extract_mem(fit_multilayer_jaccard(L[(t - 1):t], algorithm = "leiden"))
      mems[[t]] <- fit[[2]]
    }
  }
  match_labels(mems)
}

METHODS2 <- list(
  "DynMux (Jaccard, Leiden)"      = base_methods[["DynMux (Jaccard, Leiden)"]],
  "DynMux (Identity, Leiden)"     = base_methods[["DynMux (Identity, Leiden)"]],
  "DynMux (Jaccard)"              = base_methods[["DynMux (Jaccard)"]],
  "DynMux (Identity)"             = base_methods[["DynMux (Identity)"]],
  "DynMux (filtered)"             = method_dynmux_filtered,
  "Independent Leiden"            = base_methods[["Independent Leiden"]],
  "Independent Leiden (matched)"  = method_independent_leiden_matched,
  "Aggregated Leiden"             = base_methods[["Aggregated Leiden"]],
  "Multislice (Adjacent)"         = base_methods[["Multislice (Adjacent)"]],
  "multinet GLouvain"             = base_methods[["multinet GLouvain"]]
)

# --------------------- simulators ------------------------------------------
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
      if (event == "merge") {
        prev <- ceiling(prev / 2)                       # 6 -> 3
      } else {
        prev <- prev * 2L - rbinom(length(prev), 1, 0.5) # 3 -> 6
      }
    }
    K_now <- length(unique(prev))
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

# --------------------- metrics ----------------------------------------------
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

# --------------------- run grid ---------------------------------------------
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
if (MINI) cfgs <- cfgs[c(1, nrow(cfgA) + 1), ]

jobs <- expand.grid(cfg_i = seq_len(nrow(cfgs)), rep = seq_len(REPS))
cat(sprintf("[comparison2] %d configs x %d reps = %d sims x %d methods\n",
            nrow(cfgs), REPS, nrow(jobs), length(METHODS2)))

run_job <- function(j) {
  cfg <- cfgs[jobs$cfg_i[j], ]
  seed <- 9000L + jobs$cfg_i[j] * 1000L + jobs$rep[j]
  sim <- if (cfg$arm == "borrow") {
    sim_borrow(cfg$n, cfg$K, cfg$p_in, cfg$p_out, cfg$p_switch, cfg$T_, seed)
  } else {
    sim_mergesplit(cfg$n, cfg$T_, cfg$event, cfg$p_in, cfg$p_out,
                   cfg$p_switch, seed)
  }
  rows <- list(); ri <- 0L
  for (mname in names(METHODS2)) {
    t0 <- proc.time()["elapsed"]
    det <- tryCatch(METHODS2[[mname]](sim$layers), error = function(e) NULL)
    el <- proc.time()["elapsed"] - t0
    if (is.null(det)) next
    m <- eval_method(det, sim$truth)
    ri <- ri + 1L
    rows[[ri]] <- data.frame(
      arm = cfg$arm, n = cfg$n, K = cfg$K, p_in = cfg$p_in,
      p_out = cfg$p_out, p_switch = cfg$p_switch, T_layers = cfg$T_,
      event = cfg$event, rep = jobs$rep[j], method = mname,
      nmi_layer = round(m["nmi_layer"], 4),
      nmi_joint = round(m["nmi_joint"], 4),
      k_mae = round(m["k_mae"], 4),
      runtime_s = round(el, 3), stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

res <- mclapply(seq_len(nrow(jobs)), function(j) {
  tryCatch(run_job(j), error = function(e) {
    message(sprintf("job %d failed: %s", j, conditionMessage(e))); NULL
  })
}, mc.cores = CORES, mc.preschedule = FALSE)
res <- do.call(rbind, Filter(Negate(is.null), res))
write.csv(res, file.path(outdir, "comparison2_results.csv"), row.names = FALSE)

agg <- aggregate(cbind(nmi_layer, nmi_joint, k_mae) ~ arm + method,
                 data = res, FUN = mean)
cat("\n[comparison2] means by arm x method:\n")
print(agg[order(agg$arm, -agg$nmi_joint), ], row.names = FALSE)
cat(sprintf("\n[comparison2] done: %d rows -> comparison2_results.csv\n",
            nrow(res)))
