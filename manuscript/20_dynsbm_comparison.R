# =============================================================================
# manuscript/20_dynsbm_comparison.R
#
# Adds Dynamic SBM (dynsbm) to the Stage I head-to-head (Jared, 2026-07-22).
# dynsbm is ~1000x slower than the other 10 methods (~160s+/fit), so it runs
# as its OWN array (one task per config, 100 reps each, reps parallelized over
# cores) rather than inside 18_. It uses the IDENTICAL config grid, seeds, and
# simulators as 18_comparison_extensions.R, so its rows rbind directly with
# comparison2_results.csv for a matched comparison on the same networks.
#
# Fairness: Qmax = 10 gives generous headroom above the true K (<= 6) in every
# cell, so dynsbm can never be blamed for too small a candidate range; ICL
# model selection via dynsbm:::compute.icl (verified to recover true K).
# dynsbm holds ONE K fixed across the whole series by design -> the merge/split
# arm is where that structural limit shows, with every advantage granted.
#
# Usage (smoke): DYN_MINI=1 DYN_CFG=23 Rscript manuscript/20_dynsbm_comparison.R
#   (cfg 23 = mergesplit merge T=10, the fast/interesting one)
# Array: one SLURM_ARRAY_TASK_ID per config (1..22).
# =============================================================================

suppressPackageStartupMessages({ library(dynsbm); library(parallel) })

CFG_I <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", Sys.getenv("DYN_CFG", "1")))
CORES <- as.integer(Sys.getenv("DYN_CORES", "8"))
MINI  <- identical(Sys.getenv("DYN_MINI", "0"), "1")
REPS  <- if (MINI) 1L else 100L
QMIN  <- 2L
QMAX  <- 10L
NSTART <- 10L
outdir <- file.path("manuscript", "output", "dynsbm")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ------------------- identical cfg grid as 18_ ------------------------------
cfgA <- expand.grid(arm = "borrow", n = 100, K = 5,
                    p_in = c(0.12, 0.15, 0.18), p_switch = c(0.02, 0.05),
                    T_ = c(10, 20, 40), event = "none",
                    stringsAsFactors = FALSE)
cfgA$p_out <- cfgA$p_in / 3
cfgB <- expand.grid(arm = "mergesplit", n = 100, K = NA,
                    p_in = 0.30, p_switch = 0.03, T_ = c(10, 20),
                    event = c("merge", "split"), stringsAsFactors = FALSE)
cfgB$p_out <- 0.05
cfgs <- rbind(cfgA, cfgB)                       # 18 borrow + 4 mergesplit = 22
stopifnot(CFG_I >= 1, CFG_I <= nrow(cfgs))

# skip-if-done guard (resubmission never recomputes a finished config) ----
outfile <- file.path(outdir, sprintf("dynsbm_cfg%02d.csv", CFG_I))
if (nzchar(Sys.getenv("SLURM_ARRAY_TASK_ID")) && file.exists(outfile)) {
  cat("[skip] cfg", CFG_I, "already complete\n"); quit(save = "no")
}
cfg <- cfgs[CFG_I, ]

cat(sprintf("[dynsbm] cfg=%d/%d arm=%s T=%d event=%s p_in=%.2f sw=%.2f REPS=%d Qmax=%d nstart=%d cores=%d\n",
            CFG_I, nrow(cfgs), cfg$arm, cfg$T_, cfg$event, cfg$p_in,
            cfg$p_switch, REPS, QMAX, NSTART, CORES))

# ------------------- identical simulators as 18_ ----------------------------
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
        prev <- ceiling(prev / 2)
      } else {
        prev <- prev * 2L - rbinom(length(prev), 1, 0.5)
      }
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

# ------------------- identical metrics as 18_ -------------------------------
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

# ------------------- dynsbm fitter: ICL over Qmin:Qmax -----------------------
fit_dynsbm <- function(L) {
  T_ <- length(L); n <- nrow(L[[1]])
  Y <- array(0, c(T_, n, n))
  for (t in seq_len(T_)) Y[t, , ] <- (L[[t]] > 0) * 1   # binarize
  models <- select.dynsbm(Y, Qmin = QMIN, Qmax = QMAX, edge.type = "binary",
                          nstart = NSTART, nb.cores = 1, plot = FALSE)
  icl <- vapply(models, function(md) dynsbm:::compute.icl(md), numeric(1))
  best <- models[[which.max(icl)]]
  mem <- best$membership                                # n x T
  lapply(seq_len(T_), function(t) mem[, t])
}

# ------------------- run 100 reps for this config ---------------------------
run_rep <- function(rep) {
  seed <- 9000L + CFG_I * 1000L + rep                   # identical to 18_
  sim <- if (cfg$arm == "borrow") {
    sim_borrow(cfg$n, cfg$K, cfg$p_in, cfg$p_out, cfg$p_switch, cfg$T_, seed)
  } else {
    sim_mergesplit(cfg$n, cfg$T_, cfg$event, cfg$p_in, cfg$p_out,
                   cfg$p_switch, seed)
  }
  t0 <- proc.time()["elapsed"]
  det <- tryCatch(fit_dynsbm(sim$layers), error = function(e) NULL)
  el <- proc.time()["elapsed"] - t0
  if (is.null(det)) return(NULL)
  m <- eval_method(det, sim$truth)
  data.frame(arm = cfg$arm, n = cfg$n, K = cfg$K, p_in = cfg$p_in,
             p_out = cfg$p_out, p_switch = cfg$p_switch, T_layers = cfg$T_,
             event = cfg$event, rep = rep, method = "Dynamic SBM",
             nmi_layer = round(m["nmi_layer"], 4),
             nmi_joint = round(m["nmi_joint"], 4),
             k_mae = round(m["k_mae"], 4),
             runtime_s = round(el, 2), stringsAsFactors = FALSE)
}

t0 <- proc.time()["elapsed"]
res <- mclapply(seq_len(REPS), function(r) {
  tryCatch(run_rep(r), error = function(e) {
    message(sprintf("rep %d failed: %s", r, conditionMessage(e))); NULL
  })
}, mc.cores = CORES, mc.preschedule = FALSE)
res <- do.call(rbind, Filter(Negate(is.null), res))
write.csv(res, outfile, row.names = FALSE)

cat(sprintf("[dynsbm] cfg %d done: %d/%d reps ok, %.1f min. nmi_layer=%.3f nmi_joint=%.3f k_mae=%.3f\n",
            CFG_I, nrow(res), REPS, (proc.time()["elapsed"] - t0) / 60,
            mean(res$nmi_layer), mean(res$nmi_joint), mean(res$k_mae)))
