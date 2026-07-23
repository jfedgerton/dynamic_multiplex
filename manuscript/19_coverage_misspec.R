# =============================================================================
# manuscript/19_coverage_misspec.R
#
# MISSPECIFICATION coverage arm (Jared, 2026-07-22): ship-or-die test for the
# CI reliability gate (width_K <= 1 AND width_P < 0.05).
#
# Question: the gate was validated on CLEAN planted-partition (SBM) networks,
# where narrow intervals were also correct. Does it still hold when the true
# structure is MISSPECIFIED relative to that assumption? The dangerous mode is
# "confidently wrong": tight intervals centered on a wrong partition, which a
# width gate cannot see. Degree heterogeneity (hubs) is where that could hide.
#
# DGP = degree-corrected SBM (DC-SBM): planted partition + per-node degree
# propensity theta_i (E[theta]=1) so Pr(edge_ij) = p_block(i,j) * theta_i *
# theta_j. Every node keeps a true block label, so cov_K / cov_P stay
# computable exactly as in 16_coverage3_grid.R.
#
# Two misspecification dials:
#   hetero  : none (theta==1, MUST reproduce sweep 16) / moderate (lognormal)
#             / severe (power-law hubs, Pareto)
#   balance : balanced communities / skewed (largest group ~40%)
#
# Everything else (edge-resampling bootstrap, K/P estimands, Wilson intervals,
# calibration bins, shuffled order, M/B/R_TRUTH) is identical to sweep 16.
#
# Usage (mini test): COV_TASK=1 COV_MINI=1 Rscript manuscript/19_coverage_misspec.R
# =============================================================================

suppressPackageStartupMessages({
  library(dynamicmultiplex)
  library(parallel)
})

OFFSET <- as.integer(Sys.getenv("COV_OFFSET", "0"))
TASKID <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID",
                                Sys.getenv("COV_TASK", "1")))
TASK   <- OFFSET + TASKID
CORES  <- as.integer(Sys.getenv("COV_CORES", "8"))
MINI   <- identical(Sys.getenv("COV_MINI", "0"), "1")

M_SIMS  <- if (MINI) 3L else 250L
B_BOOT  <- if (MINI) 8L else 100L
R_TRUTH <- if (MINI) 5L else 100L
ALPHA   <- 0.05
N_BINS  <- 20L
Z_ALPHA <- qnorm(1 - ALPHA / 2)

DENSITIES <- list(weak    = c(p_in = 0.20, p_out = 0.10),
                  default = c(p_in = 0.30, p_out = 0.05),
                  strong  = c(p_in = 0.50, p_out = 0.02))

outdir <- file.path("manuscript", "output", "coverage3_misspec")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ------------------- reduced grid, deterministically shuffled ----------------
# n dropped 50 (mostly refuse); T fixed at 10; switching two levels (not the
# driver); densities all 3 to span the accept/refuse boundary; the two misspec
# dials fully crossed. 2 specs (both leiden, the two coupling philosophies).
cells <- expand.grid(
  n = c(100, 200, 400), K = c(3, 5),
  p_switch = c(0.05, 0.20),
  density = names(DENSITIES), T_ = 10L,
  hetero = c("none", "moderate", "severe"),
  balance = c("balanced", "skewed"),
  stringsAsFactors = FALSE)
specs <- expand.grid(fit_type = c("jaccard", "identity"),
                     algorithm = "leiden",
                     stringsAsFactors = FALSE)
task_tab <- merge(cells, specs, by = NULL)
task_tab <- task_tab[order(task_tab$n, task_tab$K, task_tab$p_switch,
                           task_tab$density, task_tab$hetero,
                           task_tab$balance, task_tab$fit_type), ]
rownames(task_tab) <- NULL

# shuffle iteration order (deterministic, seed 123) so all misspec levels
# produce results from the first hours
set.seed(123)
perm <- sample.int(nrow(task_tab))

if (TASK < 1 || TASK > nrow(task_tab)) {
  stop(sprintf("Task %d out of range 1..%d", TASK, nrow(task_tab)))
}

# skip-if-done guard (resubmission never recomputes a finished task) ----
if (nzchar(Sys.getenv("SLURM_ARRAY_TASK_ID"))) {
  if (file.exists(file.path(outdir, sprintf("cov_task%05d.csv", TASK))) &&
      file.exists(file.path(outdir, sprintf("calib_task%05d.csv", TASK)))) {
    cat("[skip] task", TASK, "already complete\n"); quit(save = "no")
  }
}

cfg <- as.list(task_tab[perm[TASK], ])
cfg$p_in <- DENSITIES[[cfg$density]][["p_in"]]
cfg$p_out <- DENSITIES[[cfg$density]][["p_out"]]

cat(sprintf("[misspec] task=%d/%d cell=%d : n=%d K=%d sw=%.2f dens=%s hetero=%s bal=%s %s/%s M=%d B=%d R=%d\n",
            TASK, nrow(task_tab), perm[TASK], cfg$n, cfg$K, cfg$p_switch,
            cfg$density, cfg$hetero, cfg$balance, cfg$fit_type, cfg$algorithm,
            M_SIMS, B_BOOT, R_TRUTH))

FITTERS <- list(
  jaccard  = fit_multilayer_jaccard,
  overlap  = fit_multilayer_overlap,
  identity = fit_multilayer_identity_ties
)

# --------------------------- misspec generators ------------------------------
# community assignment probabilities: uniform, or geometric skew (ratio 0.6)
comm_probs <- function(K, balance) {
  if (balance == "balanced") return(rep(1 / K, K))
  w <- 0.6 ^ (seq_len(K) - 1)
  w / sum(w)
}

# per-node degree propensity, E[theta] = 1 in every regime
gen_theta <- function(n, hetero) {
  if (hetero == "none") return(rep(1, n))
  if (hetero == "moderate") {
    sdlog <- 0.6
    return(rlnorm(n, meanlog = -sdlog^2 / 2, sdlog = sdlog))   # mean 1
  }
  # severe: Pareto(alpha, xmin=1) normalized to mean 1 -> power-law hubs
  alpha <- 2.5
  raw <- runif(n) ^ (-1 / alpha)          # inverse-CDF Pareto, xmin = 1
  raw / (alpha / (alpha - 1))             # divide by theoretical mean
}

gen_memberships <- function(cfg, probs) {
  memberships <- vector("list", cfg$T_)
  memberships[[1]] <- sample(seq_len(cfg$K), cfg$n, replace = TRUE, prob = probs)
  for (t in 2:cfg$T_) {
    prev <- memberships[[t - 1]]
    mask <- runif(cfg$n) < cfg$p_switch
    prev[mask] <- sample(seq_len(cfg$K), sum(mask), replace = TRUE, prob = probs)
    memberships[[t]] <- prev
  }
  memberships
}

# DC-SBM layers: Pr(edge) = p_block * theta_i * theta_j, clipped to [0,1].
# theta is FIXED per simulation (a property of the network); fresh draws reuse
# the same theta with new Bernoulli edges. hetero="none" -> plain SBM = sweep16.
gen_layers <- function(memberships, cfg, theta) {
  n <- cfg$n
  th <- outer(theta, theta)
  lapply(memberships, function(mem) {
    P  <- outer(mem, mem, "==")
    Pr <- ifelse(P, cfg$p_in, cfg$p_out) * th
    Pr[Pr > 1] <- 1
    diag(Pr) <- 0
    A  <- matrix(0, n, n)
    up <- upper.tri(Pr)
    A[up] <- rbinom(sum(up), 1, Pr[up])
    A + t(A)
  })
}

# ------------------- inline edge-resampling bootstrap ------------------------
# identical to bootstrap_multilayer(resample = "edges") on dev-1.1.0
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
  seed <- (12000L + TASK) * 100000L + sim_id   # distinct band, no overflow
  set.seed(seed)
  fitfun <- FITTERS[[cfg$fit_type]]
  n <- cfg$n; T_ <- cfg$T_
  probs <- comm_probs(cfg$K, cfg$balance)
  theta <- gen_theta(n, cfg$hetero)             # fixed per simulation

  truth  <- gen_memberships(cfg, probs)
  layers <- gen_layers(truth, cfg, theta)
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

  # ground-truth co-clustering propensity (same theta, fresh edges) ----
  pstar_acc <- lapply(seq_len(T_), function(t) matrix(0, n, n))
  for (r in seq_len(R_TRUTH)) {
    fresh <- gen_layers(truth, cfg, theta)
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
    p_in = cfg$p_in, p_out = cfg$p_out, density = cfg$density,
    T_layers = cfg$T_, hetero = cfg$hetero, balance = cfg$balance,
    weights = "binary", fit_type = cfg$fit_type,
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

cat(sprintf("[misspec] task %d done: %d/%d sims ok, %.1f min. hetero=%s bal=%s cov_K=%.3f cov_P=%.3f gate_pass=%.0f%%\n",
            TASK, nrow(rows), M_SIMS, (proc.time()["elapsed"] - t0) / 60,
            cfg$hetero, cfg$balance,
            mean(rows$cov_K_mean), mean(rows$cov_P_mean),
            100 * mean(rows$width_K_mean <= 1 & rows$width_P_mean < 0.05)))
