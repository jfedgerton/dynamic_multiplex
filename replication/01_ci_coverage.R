# =============================================================================
# replication/01_ci_coverage.R
#
# UNIFIED CI-COVERAGE VALIDATION (all three arms in one script).
#
# This is the exact code and seeds used to validate the bootstrap confidence
# intervals reported in the paper. It combines what were three development
# scripts into one arm-dispatched program that shares a single bootstrap /
# estimand / calibration core. The three arms are:
#
#   ARM = "binary"  : planted-partition (SBM) binary networks, full factorial
#                     grid. The primary coverage validation.  (seed band 15000)
#   ARM = "valued"  : weighted networks (aligned / orthogonal edge weights),
#                     including the weighted fitters.           (seed band 20000)
#   ARM = "misspec" : degree-corrected SBM (DC-SBM) with degree heterogeneity
#                     and community-size skew - the structural-misspecification
#                     stress test of the reliability gate.      (seed band 12000)
#
# Shared across arms (defined once): the parametric edge-resampling bootstrap,
# the K coverage (percentile CI on bootstrap community counts), the pairwise
# co-assignment coverage (Wilson intervals vs a fresh-draw ground-truth
# propensity), and the 20-bin calibration. Each arm differs only in (i) its
# factorial grid, (ii) its data-generating process, and (iii) its seed band.
#
# THE GATE (adopted): accept an interval set iff
#     width_K_mean <= 1  AND  width_P_mean < 0.05
# Coverage is nominal inside this observable region and is refused outside it.
#
# Seeds are deterministic. Iteration order within each arm is a fixed shuffle
# (set.seed(123)) so that completed tasks span the whole design early. Per-sim
# seeds are (band + task) * 1e5 + sim_id, identical to the as-run study.
#
# ---------------------------------------------------------------------------
# HOW THE ARRAY MAPS TO ARMS
#   Global SLURM_ARRAY_TASK_ID 1..N_BINARY               -> binary  task i
#                              N_BINARY+1 .. +N_VALUED    -> valued  task i
#                              ...       +1 .. +N_MISSPEC -> misspec task i
# so a single sbatch array covers all three arms. Per-arm local task indices
# drive the shuffle and the seed exactly as in the original per-arm scripts.
#
# Usage (single arm, local index):
#   COV_ARM=binary  COV_TASK=1 COV_MINI=1 Rscript replication/01_ci_coverage.R
#   COV_ARM=valued  COV_TASK=1            Rscript replication/01_ci_coverage.R
#   COV_ARM=misspec COV_TASK=1            Rscript replication/01_ci_coverage.R
# Usage (array, global index -> arm auto-detected): see slurm/01_ci_coverage.sbatch
# =============================================================================

suppressPackageStartupMessages({
  library(dynamicmultiplex)
  library(parallel)
})

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

OUTROOT <- Sys.getenv("COV_OUTROOT", "replication/output")

# =============================================================================
# ARM GRIDS  (each block is verbatim the grid of the original per-arm script)
# =============================================================================

build_binary <- function() {
  cells <- expand.grid(
    n = c(50, 100, 200, 400), K = c(3, 5, 10),
    p_switch = c(0.02, 0.05, 0.10, 0.20, 0.50, 0.80),
    density = names(DENSITIES), T_ = c(5, 10, 15),
    stringsAsFactors = FALSE)
  cells <- cells[!(cells$K == 10 & cells$n == 50), ]
  specs <- expand.grid(fit_type = c("jaccard", "overlap", "identity"),
                       algorithm = c("louvain", "leiden"),
                       stringsAsFactors = FALSE)
  tt <- merge(cells, specs, by = NULL)
  tt <- tt[order(tt$n, tt$T_, tt$K, tt$p_switch, tt$density,
                 tt$fit_type, tt$algorithm), ]
  rownames(tt) <- NULL
  tt$weights <- "binary"
  tt$hetero <- "none"; tt$balance <- "balanced"      # inert for this arm
  list(task_tab = tt, band = 15000L, dgp = "sbm",
       outdir = file.path(OUTROOT, "coverage_binary"))
}

build_valued <- function() {
  cells <- expand.grid(
    n = c(100, 200, 400), K = c(3, 5),
    p_switch = c(0.05, 0.20, 0.80),
    density = c("default", "strong"), T_ = 10,
    weights = c("aligned", "orthogonal"),
    stringsAsFactors = FALSE)
  specs <- expand.grid(fit_type = c("jaccard", "overlap", "weighted_jaccard",
                                    "weighted_overlap", "identity"),
                       algorithm = c("louvain", "leiden"),
                       stringsAsFactors = FALSE)
  tt <- merge(cells, specs, by = NULL)
  tt <- tt[order(tt$n, tt$T_, tt$K, tt$p_switch, tt$density,
                 tt$weights, tt$fit_type, tt$algorithm), ]
  rownames(tt) <- NULL
  tt$hetero <- "none"; tt$balance <- "balanced"      # inert for this arm
  list(task_tab = tt, band = 20000L, dgp = "sbm",
       outdir = file.path(OUTROOT, "coverage_valued"))
}

build_misspec <- function() {
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
  tt <- merge(cells, specs, by = NULL)
  tt <- tt[order(tt$n, tt$K, tt$p_switch, tt$density, tt$hetero,
                 tt$balance, tt$fit_type), ]
  rownames(tt) <- NULL
  tt$weights <- "binary"
  list(task_tab = tt, band = 12000L, dgp = "dcsbm",
       outdir = file.path(OUTROOT, "coverage_misspec"))
}

ARMS <- list(binary = build_binary, valued = build_valued, misspec = build_misspec)

# --- resolve arm + local task index (single arm, or global array index) ------
resolve_task <- function() {
  arm_env <- Sys.getenv("COV_ARM", "")
  gid <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID",
                               Sys.getenv("COV_TASK", "1")))
  if (nzchar(arm_env)) {
    return(list(arm = arm_env, task = gid))
  }
  # global-index mode: walk arms in fixed order until gid lands in one
  order_arms <- c("binary", "valued", "misspec")
  offset <- 0L
  for (a in order_arms) {
    n_a <- nrow(ARMS[[a]]()$task_tab)
    if (gid <= offset + n_a) return(list(arm = a, task = gid - offset))
    offset <- offset + n_a
  }
  stop(sprintf("global task %d exceeds total %d", gid, offset))
}

sel <- resolve_task()
ARM <- sel$arm
TASK <- sel$task
spec <- ARMS[[ARM]]()
task_tab <- spec$task_tab
BAND <- spec$band
DGP <- spec$dgp
outdir <- spec$outdir
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# deterministic shuffle of iteration order (seed 123), identical to originals
set.seed(123)
perm <- sample.int(nrow(task_tab))

if (TASK < 1 || TASK > nrow(task_tab)) {
  stop(sprintf("[%s] task %d out of range 1..%d", ARM, TASK, nrow(task_tab)))
}

# skip-if-done guard (resubmission never recomputes a finished task)
if (nzchar(Sys.getenv("SLURM_ARRAY_TASK_ID"))) {
  if (file.exists(file.path(outdir, sprintf("cov_task%05d.csv", TASK))) &&
      file.exists(file.path(outdir, sprintf("calib_task%05d.csv", TASK)))) {
    cat(sprintf("[%s skip] task %d already complete\n", ARM, TASK))
    quit(save = "no")
  }
}

cfg <- as.list(task_tab[perm[TASK], ])
cfg$p_in <- DENSITIES[[cfg$density]][["p_in"]]
cfg$p_out <- DENSITIES[[cfg$density]][["p_out"]]

cat(sprintf("[ci:%s] task=%d/%d cell=%d : n=%d K=%d sw=%.2f dens=%s T=%d w=%s hetero=%s bal=%s %s/%s M=%d B=%d R=%d\n",
            ARM, TASK, nrow(task_tab), perm[TASK], cfg$n, cfg$K, cfg$p_switch,
            cfg$density, cfg$T_, cfg$weights, cfg$hetero, cfg$balance,
            cfg$fit_type, cfg$algorithm, M_SIMS, B_BOOT, R_TRUTH))

FITTERS <- list(
  jaccard          = fit_multilayer_jaccard,
  overlap          = fit_multilayer_overlap,
  identity         = fit_multilayer_identity_ties,
  weighted_jaccard = fit_multilayer_weighted_jaccard,
  weighted_overlap = fit_multilayer_weighted_overlap
)

# =============================================================================
# DATA-GENERATING PROCESSES
#   sbm   : planted partition, uniform community sizes (binary + valued arms).
#           weighted arms multiply Bernoulli edges by lognormal weights.
#   dcsbm : degree-corrected SBM - per-node degree propensity theta_i (mean 1)
#           and optional community-size skew (misspec arm).
# =============================================================================

# ---- shared: uniform-membership Markov switching ----------------------------
gen_memberships_sbm <- function(cfg) {
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

gen_layers_sbm <- function(memberships, cfg) {
  n <- cfg$n
  lapply(memberships, function(mem) {
    P  <- outer(mem, mem, "==")
    Pr <- ifelse(P, cfg$p_in, cfg$p_out); diag(Pr) <- 0
    A  <- matrix(0, n, n)
    up <- upper.tri(Pr)
    A[up] <- rbinom(sum(up), 1, Pr[up])
    A  <- A + t(A)
    if (identical(cfg$weights, "aligned")) {
      W <- matrix(0, n, n)
      mu <- ifelse(P, 1.0, 0.0)
      W[up] <- A[up] * rlnorm(sum(up), meanlog = mu[up], sdlog = 0.75)
      A <- W + t(W)
    } else if (identical(cfg$weights, "orthogonal")) {
      W <- matrix(0, n, n)
      W[up] <- A[up] * rlnorm(sum(up), meanlog = 0.5, sdlog = 0.75)
      A <- W + t(W)
    }
    A
  })
}

# ---- misspec: degree-corrected SBM ------------------------------------------
comm_probs <- function(K, balance) {
  if (balance == "balanced") return(rep(1 / K, K))
  w <- 0.6 ^ (seq_len(K) - 1)
  w / sum(w)
}

gen_theta <- function(n, hetero) {
  if (hetero == "none") return(rep(1, n))
  if (hetero == "moderate") {
    sdlog <- 0.6
    return(rlnorm(n, meanlog = -sdlog^2 / 2, sdlog = sdlog))   # mean 1
  }
  alpha <- 2.5                                # severe: Pareto power-law hubs
  raw <- runif(n) ^ (-1 / alpha)              # inverse-CDF Pareto, xmin = 1
  raw / (alpha / (alpha - 1))                 # normalize to mean 1
}

gen_memberships_dcsbm <- function(cfg) {
  probs <- comm_probs(cfg$K, cfg$balance)
  memberships <- vector("list", cfg$T_)
  memberships[[1]] <- sample(seq_len(cfg$K), cfg$n, replace = TRUE, prob = probs)
  for (t in 2:cfg$T_) {
    prev <- memberships[[t - 1]]
    mask <- runif(cfg$n) < cfg$p_switch
    prev[mask] <- sample(seq_len(cfg$K), sum(mask), replace = TRUE, prob = probs)
    memberships[[t]] <- prev
  }
  attr(memberships, "theta") <- gen_theta(cfg$n, cfg$hetero)   # fixed per sim
  memberships
}

gen_layers_dcsbm <- function(memberships, cfg) {
  n <- cfg$n
  theta <- attr(memberships, "theta")
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

gen_memberships <- if (DGP == "dcsbm") gen_memberships_dcsbm else gen_memberships_sbm
gen_layers      <- if (DGP == "dcsbm") gen_layers_dcsbm      else gen_layers_sbm

# =============================================================================
# SHARED BOOTSTRAP + ESTIMANDS  (identical across arms; parametric edge
# resampling, mirroring bootstrap_multilayer(resample = "edges") in 1.1.0)
# =============================================================================
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

# =============================================================================
# ONE SIMULATION
# =============================================================================
run_one <- function(sim_id) {
  seed <- (BAND + TASK) * 100000L + sim_id
  set.seed(seed)
  fitfun <- FITTERS[[cfg$fit_type]]
  n <- cfg$n; T_ <- cfg$T_

  truth  <- gen_memberships(cfg)
  layers <- gen_layers(truth, cfg)
  fit0 <- fitfun(layers, algorithm = cfg$algorithm)
  mem_hat <- get_memberships(fit0, T_)

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

  # K coverage: percentile CI on bootstrap community counts
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

  # fresh-draw ground-truth co-clustering propensity
  pstar_acc <- lapply(seq_len(T_), function(t) matrix(0, n, n))
  for (r in seq_len(R_TRUTH)) {
    fresh <- gen_layers(truth, cfg)
    tfit <- fitfun(fresh, algorithm = cfg$algorithm)
    tmem <- get_memberships(tfit, T_)
    for (t in seq_len(T_)) {
      pstar_acc[[t]] <- pstar_acc[[t]] + outer(tmem[[t]], tmem[[t]], "==")
    }
  }

  # Wilson coverage of the co-assignment propensity + calibration bins
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
    arm = ARM, task = TASK, sim = sim_id, resample = "edges",
    n = cfg$n, K = cfg$K, p_switch = cfg$p_switch,
    p_in = cfg$p_in, p_out = cfg$p_out, density = cfg$density,
    T_layers = cfg$T_, weights = cfg$weights,
    hetero = cfg$hetero, balance = cfg$balance,
    fit_type = cfg$fit_type, algorithm = cfg$algorithm,
    B = b_ok, R_truth = R_TRUTH,
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

# =============================================================================
# RUN + WRITE
# =============================================================================
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
calib <- data.frame(arm = ARM, task = TASK, bin = seq_len(N_BINS),
                    bin_lo = (seq_len(N_BINS) - 1) / N_BINS,
                    bin_hi = seq_len(N_BINS) / N_BINS,
                    n_pairs = btot, n_true = btrue,
                    stringsAsFactors = FALSE)

write.csv(rows, file.path(outdir, sprintf("cov_task%05d.csv", TASK)),
          row.names = FALSE)
write.csv(calib, file.path(outdir, sprintf("calib_task%05d.csv", TASK)),
          row.names = FALSE)

cat(sprintf("[ci:%s] task %d done: %d/%d sims ok, %.1f min. cov_K=%.3f cov_P=%.3f gate_pass=%.0f%%\n",
            ARM, TASK, nrow(rows), M_SIMS, (proc.time()["elapsed"] - t0) / 60,
            mean(rows$cov_K_mean), mean(rows$cov_P_mean),
            100 * mean(rows$width_K_mean <= 1 & rows$width_P_mean < 0.05)))
