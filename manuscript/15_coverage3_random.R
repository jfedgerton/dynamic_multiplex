# =============================================================================
# manuscript/15_coverage3_random.R
#
# Coverage validation v3: RANDOMIZED DESIGN, edge-resampling bootstrap.
#
# Each simulation draws its own configuration (network size, communities,
# switching, densities, layers, weights regime, fitter, algorithm), so every
# completed task samples the full design space and interim results are
# globally informative immediately (Jared's request, 2026-07-21).
#
# The bootstrap here is the PARAMETRIC EDGE-RESAMPLING scheme now implemented
# as bootstrap_multilayer(resample = "edges") on the dev-1.1.0 branch. Because
# the installed package on ROAR is 1.0.0 (and must not be touched while the
# v2 sweep runs), the resampling loop is implemented inline in this script
# with the identical algorithm: estimate within/between-community edge
# probabilities and weight pools from the observed network using the
# point-estimate partition, then redraw the full edge set each replicate.
#
# Estimands per simulation
#   K coverage : planted community count vs percentile CI from bootstrap
#                count samples (same construction as community_count_ci)
#   P coverage : Wilson 95% intervals for pairwise co-assignment propensity
#                vs ground truth from R_TRUTH fresh edge draws (same planted
#                membership path, same fitter)
#   calibration: 20-bin calibration of co-assignment probabilities
#
# Usage (mini test): COV_TASK=1 COV_MINI=1 Rscript manuscript/15_coverage3_random.R
# =============================================================================

suppressPackageStartupMessages({
  library(dynamicmultiplex)
  library(parallel)
})

TASKID <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID",
                                Sys.getenv("COV_TASK", "1")))
TASK   <- TASKID
CORES  <- as.integer(Sys.getenv("COV_CORES", "8"))
MINI   <- identical(Sys.getenv("COV_MINI", "0"), "1")

M_SIMS  <- if (MINI) 3L else 250L
B_BOOT  <- if (MINI) 8L else 100L
R_TRUTH <- if (MINI) 5L else 100L
ALPHA   <- 0.05
N_BINS  <- 20L
Z_ALPHA <- qnorm(1 - ALPHA / 2)

outdir <- file.path("manuscript", "output", "coverage3_random")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("[coverage3] task=%d M=%d B=%d R=%d (randomized design, edge-resampling bootstrap)\n",
            TASK, M_SIMS, B_BOOT, R_TRUTH))

FITTERS <- list(
  jaccard          = fit_multilayer_jaccard,
  overlap          = fit_multilayer_overlap,
  identity         = fit_multilayer_identity_ties,
  weighted_jaccard = fit_multilayer_weighted_jaccard,
  weighted_overlap = fit_multilayer_weighted_overlap
)

# --------------------------- random configuration ----------------------------
draw_config <- function() {
  n <- round(exp(runif(1, log(50), log(400))))          # log-uniform 50..400
  K <- sample(3:10, 1)
  p_switch <- runif(1, 0.02, 0.80)
  p_in  <- runif(1, 0.15, 0.50)
  p_out <- min(runif(1, 0.02, 0.15), p_in / 2)
  T_ <- sample(5:15, 1)
  weights <- sample(c("binary", "aligned", "orthogonal"), 1,
                    prob = c(0.6, 0.2, 0.2))
  fit_type <- if (weights == "binary") {
    sample(c("jaccard", "overlap", "identity"), 1)
  } else {
    sample(c("jaccard", "overlap", "identity",
             "weighted_jaccard", "weighted_overlap"), 1)
  }
  algorithm <- sample(c("louvain", "leiden"), 1)
  list(n = n, K = K, p_switch = p_switch, p_in = p_in, p_out = p_out,
       T_ = T_, weights = weights, fit_type = fit_type, algorithm = algorithm)
}

# --------------------------- simulator ---------------------------------------
gen_memberships <- function(cfg) {
  memberships <- vector("list", cfg$T_)
  memberships[[1]] <- sample(seq_len(cfg$K), cfg$n, replace = TRUE)
  for (t in 2:cfg$T_) {
    prev <- memberships[[t - 1]]
    mask <- runif(cfg$n) < cfg$p_switch
    prev[mask] <- sample(seq_len(cfg$K), sum(mask), replace = TRUE)
    memberships[[t]] <- prev
  }
  memberships
}

gen_layers <- function(memberships, cfg) {
  n <- cfg$n
  lapply(memberships, function(mem) {
    P  <- outer(mem, mem, "==")
    Pr <- ifelse(P, cfg$p_in, cfg$p_out); diag(Pr) <- 0
    A  <- matrix(0, n, n)
    up <- upper.tri(Pr)
    A[up] <- rbinom(sum(up), 1, Pr[up])
    A  <- A + t(A)
    if (cfg$weights == "aligned") {
      W <- matrix(0, n, n)
      mu <- ifelse(P, 1.0, 0.0)
      W[up] <- A[up] * rlnorm(sum(up), meanlog = mu[up], sdlog = 0.75)
      A <- W + t(W)
    } else if (cfg$weights == "orthogonal") {
      W <- matrix(0, n, n)
      W[up] <- A[up] * rlnorm(sum(up), meanlog = 0.5, sdlog = 0.75)
      A <- W + t(W)
    }
    A
  })
}

# ------------------- inline edge-resampling bootstrap ------------------------
# identical algorithm to bootstrap_multilayer(resample = "edges") on the
# dev-1.1.0 branch; implemented here so the installed 1.0.0 stays untouched
fit_edge_models <- function(layers, memberships_hat) {
  n <- nrow(layers[[1]])
  lapply(seq_along(layers), function(t) {
    A <- layers[[t]]
    mem <- memberships_hat[[t]]
    same <- outer(mem, mem, "==")
    sel <- upper.tri(A)
    E <- A > 0
    in_d <- sel & same; out_d <- sel & !same
    p_all <- mean(E[sel])
    p_in <- if (any(in_d)) mean(E[in_d]) else p_all
    p_out <- if (any(out_d)) mean(E[out_d]) else p_all
    w_all <- A[sel & E]; if (length(w_all) == 0) w_all <- 1
    w_in <- A[in_d & E]; if (length(w_in) == 0) w_in <- w_all
    w_out <- A[out_d & E]; if (length(w_out) == 0) w_out <- w_all
    list(same = same, sel = sel, p_in = p_in, p_out = p_out,
         w_in = w_in, w_out = w_out)
  })
}

redraw_layers <- function(edge_models, n) {
  lapply(edge_models, function(em) {
    idx <- which(em$sel)
    probs <- ifelse(em$same[idx], em$p_in, em$p_out)
    on <- idx[rbinom(length(idx), 1, probs) == 1]
    M <- matrix(0, n, n)
    if (length(on) > 0) {
      same_on <- em$same[on]
      w <- numeric(length(on))
      if (any(same_on)) w[same_on] <- sample(em$w_in, sum(same_on), replace = TRUE)
      if (any(!same_on)) w[!same_on] <- sample(em$w_out, sum(!same_on), replace = TRUE)
      M[on] <- w
    }
    M <- M + t(M); diag(M) <- 0
    M
  })
}

wilson_bounds <- function(phat, b) {
  z2 <- Z_ALPHA^2
  denom  <- 1 + z2 / b
  center <- (phat + z2 / (2 * b)) / denom
  half   <- Z_ALPHA * sqrt(phat * (1 - phat) / b + z2 / (4 * b^2)) / denom
  list(lo = pmax(center - half, 0), hi = pmin(center + half, 1))
}

get_memberships <- function(fit, T_) {
  lapply(seq_len(T_), function(t) fit$layer_communities[[t]]$membership)
}

# --------------------------- one simulation ----------------------------------
run_one <- function(sim_id) {
  seed <- (10000L + TASK) * 100000L + sim_id
  set.seed(seed)
  cfg <- draw_config()
  fitfun <- FITTERS[[cfg$fit_type]]
  n <- cfg$n; T_ <- cfg$T_

  truth  <- gen_memberships(cfg)
  layers <- gen_layers(truth, cfg)

  # point estimate ----
  fit0 <- fitfun(layers, algorithm = cfg$algorithm)
  mem_hat <- get_memberships(fit0, T_)

  # edge-resampling bootstrap ----
  edge_models <- fit_edge_models(layers, mem_hat)
  co_acc <- lapply(seq_len(T_), function(t) matrix(0, n, n))
  count_samples <- matrix(NA_integer_, nrow = B_BOOT, ncol = T_)
  b_ok <- 0L
  for (b in seq_len(B_BOOT)) {
    boot_layers <- redraw_layers(edge_models, n)
    bfit <- tryCatch(fitfun(boot_layers, algorithm = cfg$algorithm),
                     error = function(e) NULL)
    if (is.null(bfit)) next
    b_ok <- b_ok + 1L
    bmem <- get_memberships(bfit, T_)
    for (t in seq_len(T_)) {
      co_acc[[t]] <- co_acc[[t]] + outer(bmem[[t]], bmem[[t]], "==")
      count_samples[b, t] <- length(unique(bmem[[t]]))
    }
  }
  if (b_ok < 10L && !MINI) stop("too few completed bootstrap replicates")
  co_prob <- lapply(co_acc, function(m) m / b_ok)

  # K coverage: percentile CI on bootstrap count samples ----
  true_K <- vapply(truth, function(m) length(unique(m)), integer(1))
  cov_K <- vapply(seq_len(T_), function(t) {
    s <- count_samples[, t]; s <- s[!is.na(s)]
    qs <- quantile(s, probs = c(ALPHA / 2, 1 - ALPHA / 2), names = FALSE)
    (qs[1] <= true_K[t]) && (true_K[t] <= qs[2])
  }, logical(1))
  width_K <- vapply(seq_len(T_), function(t) {
    s <- count_samples[, t]; s <- s[!is.na(s)]
    diff(quantile(s, probs = c(ALPHA / 2, 1 - ALPHA / 2), names = FALSE))
  }, numeric(1))

  # ground-truth co-clustering propensity ----
  pstar_acc <- lapply(seq_len(T_), function(t) matrix(0, n, n))
  for (r in seq_len(R_TRUTH)) {
    fresh <- gen_layers(truth, cfg)
    tfit <- fitfun(fresh, algorithm = cfg$algorithm)
    tmem <- get_memberships(tfit, T_)
    for (t in seq_len(T_)) {
      pstar_acc[[t]] <- pstar_acc[[t]] + outer(tmem[[t]], tmem[[t]], "==")
    }
  }

  # Wilson coverage of p_star ----
  bin_total <- integer(N_BINS); bin_true <- integer(N_BINS)
  covP_num <- 0; covP_den <- 0
  covP_same_num <- 0; covP_same_den <- 0
  covP_diff_num <- 0; covP_diff_den <- 0
  widP_sum <- 0
  for (t in seq_len(T_)) {
    up    <- upper.tri(co_prob[[t]])
    phat  <- co_prob[[t]][up]
    pstar <- (pstar_acc[[t]] / R_TRUTH)[up]
    same  <- outer(truth[[t]], truth[[t]], "==")[up]
    wb    <- wilson_bounds(phat, b_ok)
    covp  <- (wb$lo <= pstar) & (pstar <= wb$hi)
    covP_num <- covP_num + sum(covp); covP_den <- covP_den + length(covp)
    covP_same_num <- covP_same_num + sum(covp[same])
    covP_same_den <- covP_same_den + sum(same)
    covP_diff_num <- covP_diff_num + sum(covp[!same])
    covP_diff_den <- covP_diff_den + sum(!same)
    widP_sum <- widP_sum + sum(wb$hi - wb$lo)
    bin <- pmin(N_BINS, 1L + floor(phat * N_BINS))
    bin_total <- bin_total + tabulate(bin, nbins = N_BINS)
    bin_true  <- bin_true + tabulate(bin[same], nbins = N_BINS)
  }

  row <- data.frame(
    task = TASK, sim = sim_id, resample = "edges",
    n = cfg$n, K = cfg$K, p_switch = cfg$p_switch,
    p_in = cfg$p_in, p_out = cfg$p_out, T_layers = cfg$T_,
    weights = cfg$weights, fit_type = cfg$fit_type,
    algorithm = cfg$algorithm, B = b_ok, R_truth = R_TRUTH,
    cov_K_mean = mean(cov_K),
    cov_P_mean = covP_num / covP_den,
    cov_P_same = ifelse(covP_same_den > 0, covP_same_num / covP_same_den, NA),
    cov_P_diff = ifelse(covP_diff_den > 0, covP_diff_num / covP_diff_den, NA),
    width_K_mean = mean(width_K),
    width_P_mean = widP_sum / covP_den,
    stringsAsFactors = FALSE
  )
  out <- list(row = row, bin_total = bin_total, bin_true = bin_true)
  rm(truth, layers, fit0, mem_hat, edge_models, co_acc, co_prob, pstar_acc)
  gc(verbose = FALSE)
  out
}

# --------------------------- run + write -------------------------------------
t0 <- proc.time()["elapsed"]
res <- mclapply(seq_len(M_SIMS), function(i) {
  out <- tryCatch(run_one(i), error = function(e) {
    message(sprintf("sim %d failed: %s", i, conditionMessage(e))); NULL
  })
  gc(verbose = FALSE)
  out
}, mc.cores = CORES, mc.preschedule = FALSE)
res <- Filter(Negate(is.null), res)

rows  <- do.call(rbind, lapply(res, `[[`, "row"))
btot  <- Reduce(`+`, lapply(res, `[[`, "bin_total"))
btrue <- Reduce(`+`, lapply(res, `[[`, "bin_true"))
calib <- data.frame(task = TASK, bin = seq_len(N_BINS),
                    bin_lo = (seq_len(N_BINS) - 1) / N_BINS,
                    bin_hi = seq_len(N_BINS) / N_BINS,
                    n_pairs = btot, n_true = btrue,
                    stringsAsFactors = FALSE)

write.csv(rows, file.path(outdir, sprintf("cov_task%05d.csv", TASK)),
          row.names = FALSE)
write.csv(calib, file.path(outdir, sprintf("calib_task%05d.csv", TASK)),
          row.names = FALSE)

cat(sprintf("[coverage3] task %d done: %d/%d sims ok, %.1f min. cov_K=%.3f cov_P=%.3f (same=%.3f diff=%.3f)\n",
            TASK, nrow(rows), M_SIMS, (proc.time()["elapsed"] - t0) / 60,
            mean(rows$cov_K_mean), mean(rows$cov_P_mean),
            mean(rows$cov_P_same, na.rm = TRUE),
            mean(rows$cov_P_diff, na.rm = TRUE)))
