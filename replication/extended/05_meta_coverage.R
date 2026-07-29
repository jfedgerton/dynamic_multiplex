# =============================================================================
# replication/extended/05_meta_coverage.R
#
# META CO-MEMBERSHIP INTERVAL COVERAGE -- re-validation on the NEW cross-layer
# (meta) communities, across ALL TEN DynMux deployments.
#
# replication/01_ci_coverage.R validated the co-assignment intervals when the
# package still exposed only per-layer communities. The package now computes
# co_assignment_ci() on the tracked META partition automatically:
# bootstrap_multilayer() accumulates co-assignment on fit$meta_communities
# (the cross-layer second-stage partition), and co_assignment_ci() turns that
# into Wilson intervals. So this harness measures META co-assignment coverage
# simply by using the package path -- it is the same estimand and the same gate
# logic as 01_ci_coverage.R, now applied to the tracked partition.
#
# WHAT ONE SIM DOES (mirrors 01_ci_coverage.R, but via the package):
#   1. simulate a dynamic network (planted-partition SBM or DC-SBM);
#   2. boot <- bootstrap_multilayer(layers, fit_type=<spec>, algorithm=<algo>,
#              n_boot=B) ; ci <- co_assignment_ci(boot)   [both META];
#   3. ground truth: R_TRUTH fresh draws, refit, take the META membership
#      (extract_meta_membership) and accumulate the co-clustering propensity
#      p* for each node pair;
#   4. coverage of the Wilson [lower, upper] interval for p*, exactly as 01.
#
# DEPLOYMENTS (10) = spec x algorithm:
#   specs {jaccard, overlap, weighted_jaccard, weighted_overlap, identity}
#   algos {leiden, louvain}
#
# COVERAGE CONFIG BLOCKS (15 per deployment):
#   BALANCED (9): N{50,100,200} x separation{weak,moderate,strong}, balanced
#                 planted-partition SBM (K=3, T=6, p_switch=0.05).
#   SKEWED   (3): N{50,100,200}, moderate separation, strongly skewed community
#                 sizes (geometric ratio 0.6), plain SBM (K=4).
#   MISSPEC  (3): N{50,100,200}, moderate separation, DC-SBM with severe degree
#                 heterogeneity (Pareto power-law hubs, theta mean 1), K=3.
#   The SKEWED and MISSPEC blocks are the ship-or-die check on the meta
#   communities: nominal cov_P must hold there, not just in the balanced grid.
#
# Total tasks = 10 deployments x 15 configs = 150. A fixed shuffle
# (set.seed(123); ord <- sample(150)) maps the SLURM array index to a task so
# that even a partial run samples every deployment / block / N / separation
# early. COV_TASK selects a task for local runs.
#
# THE GATE: width_P_mean < 0.05. Coverage is trusted (nominal ~95%) inside the
# gate and refused outside it; cov_P inside the gate is the ship criterion.
#
# Seeds are deterministic: per-sim seed = (BAND + TASK) * 1e5 + sim_id.
# The code is sequential (a plain sim loop) to stay debuggable.
#
# Usage (local smoke, one task, MINI):
#   COV_MINI=1 COV_TASK=1 Rscript replication/extended/05_meta_coverage.R
# Array: one SLURM_ARRAY_TASK_ID per task (1..150). See slurm/07_*.sbatch.
# =============================================================================

suppressMessages(pkgload::load_all("r_code", quiet = TRUE))

MINI    <- identical(Sys.getenv("COV_MINI", "0"), "1")
M_SIMS  <- if (MINI)  20L else 200L
B_BOOT  <- if (MINI)  30L else 100L
R_TRUTH <- if (MINI)  30L else 100L
ALPHA   <- 0.05
BAND    <- 30000L                      # seed band, distinct from 01's arms

OUTROOT <- Sys.getenv("COV_OUTROOT", "replication/extended/output/coverage")
dir.create(OUTROOT, showWarnings = FALSE, recursive = TRUE)

# Separation levels: within/between edge probabilities (moderate == 01's
# "default"). Larger p_in/p_out gap = easier detection.
SEP <- list(
  weak     = c(p_in = 0.20, p_out = 0.10),
  moderate = c(p_in = 0.30, p_out = 0.05),
  strong   = c(p_in = 0.50, p_out = 0.02)
)

# Deployment fitters (the package fit_type used by bootstrap_multilayer and
# for the fresh-draw ground-truth fits).
FITTERS <- list(
  jaccard          = fit_multilayer_jaccard,
  overlap          = fit_multilayer_overlap,
  weighted_jaccard = fit_multilayer_weighted_jaccard,
  weighted_overlap = fit_multilayer_weighted_overlap,
  identity         = fit_multilayer_identity_ties
)

# =============================================================================
# DATA-GENERATING PROCESSES  (verbatim structure from 01_ci_coverage.R)
#   sbm   : planted partition; community sizes uniform (balanced) or geometric
#           skew ratio 0.6 (skewed).
#   dcsbm : degree-corrected SBM; per-node degree propensity theta_i (mean 1),
#           fixed per simulation, reused across fresh draws.
# =============================================================================

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
    return(rlnorm(n, meanlog = -sdlog^2 / 2, sdlog = sdlog))    # mean 1
  }
  alpha <- 2.5                                # severe: Pareto power-law hubs
  raw <- runif(n) ^ (-1 / alpha)              # inverse-CDF Pareto, xmin = 1
  raw / (alpha / (alpha - 1))                 # normalize to mean 1
}

# membership sequence with Markov switching; skewed sizes via comm_probs.
# For DC-SBM the per-sim theta is attached as an attribute.
gen_memberships <- function(cfg) {
  probs <- comm_probs(cfg$K, cfg$balance)
  memberships <- vector("list", cfg$T_)
  memberships[[1]] <- sample(seq_len(cfg$K), cfg$n, replace = TRUE, prob = probs)
  for (t in 2:cfg$T_) {
    prev <- memberships[[t - 1]]
    mask <- runif(cfg$n) < cfg$p_switch
    prev[mask] <- sample(seq_len(cfg$K), sum(mask), replace = TRUE, prob = probs)
    memberships[[t]] <- prev
  }
  if (cfg$dgp == "dcsbm") attr(memberships, "theta") <- gen_theta(cfg$n, cfg$hetero)
  memberships
}

gen_layers_sbm <- function(memberships, cfg) {
  n <- cfg$n
  lapply(memberships, function(mem) {
    P  <- outer(mem, mem, "==")
    Pr <- ifelse(P, cfg$p_in, cfg$p_out); diag(Pr) <- 0
    A  <- matrix(0, n, n); up <- upper.tri(Pr)
    A[up] <- rbinom(sum(up), 1, Pr[up]); A + t(A)
  })
}

gen_layers_dcsbm <- function(memberships, cfg) {
  n <- cfg$n
  theta <- attr(memberships, "theta")
  th <- outer(theta, theta)
  lapply(memberships, function(mem) {
    P  <- outer(mem, mem, "==")
    Pr <- ifelse(P, cfg$p_in, cfg$p_out) * th
    Pr[Pr > 1] <- 1; diag(Pr) <- 0
    A  <- matrix(0, n, n); up <- upper.tri(Pr)
    A[up] <- rbinom(sum(up), 1, Pr[up]); A + t(A)
  })
}

gen_layers <- function(memberships, cfg) {
  if (cfg$dgp == "dcsbm") gen_layers_dcsbm(memberships, cfg)
  else                    gen_layers_sbm(memberships, cfg)
}

# =============================================================================
# TASK GRID  (10 deployments x 15 configs = 150) + fixed shuffle
# =============================================================================
deployments <- expand.grid(
  spec = c("jaccard", "overlap", "weighted_jaccard", "weighted_overlap", "identity"),
  algo = c("leiden", "louvain"),
  stringsAsFactors = FALSE)

cfg_balanced <- expand.grid(n = c(50L, 100L, 200L),
                            separation = c("weak", "moderate", "strong"),
                            stringsAsFactors = FALSE)
cfg_balanced$block <- "balanced"; cfg_balanced$dgp <- "sbm"
cfg_balanced$balance <- "balanced"; cfg_balanced$hetero <- "none"; cfg_balanced$K <- 3L

cfg_skewed <- data.frame(n = c(50L, 100L, 200L), separation = "moderate",
                         block = "skewed", dgp = "sbm",
                         balance = "skewed", hetero = "none", K = 4L,
                         stringsAsFactors = FALSE)

cfg_misspec <- data.frame(n = c(50L, 100L, 200L), separation = "moderate",
                          block = "misspec", dgp = "dcsbm",
                          balance = "balanced", hetero = "severe", K = 3L,
                          stringsAsFactors = FALSE)

cols <- c("n", "separation", "block", "dgp", "balance", "hetero", "K")
configs <- rbind(cfg_balanced[cols], cfg_skewed[cols], cfg_misspec[cols])
configs$T_ <- 6L
configs$p_switch <- 0.05
stopifnot(nrow(configs) == 15L)

tasks <- merge(deployments, configs, by = NULL)
tasks <- tasks[order(tasks$spec, tasks$algo, tasks$block, tasks$separation, tasks$n), ]
rownames(tasks) <- NULL
tasks$deployment <- paste(tasks$spec, tasks$algo, sep = "_")
stopifnot(nrow(tasks) == 150L)

TASK <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", Sys.getenv("COV_TASK", "1")))
stopifnot(TASK >= 1L, TASK <= nrow(tasks))

set.seed(123)
ord <- sample(nrow(tasks))
cfg <- as.list(tasks[ord[TASK], ])
cfg$p_in  <- SEP[[cfg$separation]][["p_in"]]
cfg$p_out <- SEP[[cfg$separation]][["p_out"]]

outfile <- file.path(OUTROOT, sprintf("cov_task%03d.csv", TASK))
if (nzchar(Sys.getenv("SLURM_ARRAY_TASK_ID")) && file.exists(outfile)) {
  cat("[meta-cov skip] task", TASK, "already complete:", outfile, "\n")
  quit(save = "no")
}

cat(sprintf(paste0("[meta-cov] task=%d/%d -> %s | block=%s n=%d sep=%s (p_in=%.2f p_out=%.2f) ",
                   "K=%d T=%d dgp=%s bal=%s het=%s | M=%d B=%d R=%d\n"),
            TASK, nrow(tasks), cfg$deployment, cfg$block, cfg$n, cfg$separation,
            cfg$p_in, cfg$p_out, cfg$K, cfg$T_, cfg$dgp, cfg$balance, cfg$hetero,
            M_SIMS, B_BOOT, R_TRUTH))

# =============================================================================
# ONE SIMULATION  (package bootstrap for meta co-assignment; fresh-draw p*)
# =============================================================================
run_one <- function(sim_id) {
  # compact, collision-free, within 32-bit int range (sim_id < 1000 per task)
  seed <- BAND + TASK * 1000L + sim_id
  set.seed(seed)
  fitfun <- FITTERS[[cfg$spec]]
  n <- cfg$n; T_ <- cfg$T_

  truth  <- gen_memberships(cfg)
  layers <- gen_layers(truth, cfg)

  # --- package bootstrap + Wilson intervals on the META co-assignment --------
  boot <- bootstrap_multilayer(layers, fit_type = cfg$spec, algorithm = cfg$algo,
                               n_boot = B_BOOT)
  if (boot$n_boot < 10L && !MINI) stop("too few completed bootstrap replicates")
  ci <- co_assignment_ci(boot, alpha = ALPHA)         # per-layer est/lower/upper

  # meta point-estimate community counts (per layer) vs the true counts
  meta_pt <- boot$point_estimate$meta_communities
  meta_K  <- vapply(meta_pt, function(m) length(unique(m)), integer(1))
  true_K  <- vapply(truth,   function(m) length(unique(m)), integer(1))

  # --- fresh-draw ground-truth META co-clustering propensity p* --------------
  pstar_acc <- lapply(seq_len(T_), function(t) matrix(0, n, n))
  for (r in seq_len(R_TRUTH)) {
    fresh <- gen_layers(truth, cfg)
    ffit  <- fitfun(fresh, algorithm = cfg$algo)
    fmem  <- extract_meta_membership(ffit)            # META membership
    for (t in seq_len(T_))
      pstar_acc[[t]] <- pstar_acc[[t]] + outer(fmem[[t]], fmem[[t]], "==")
  }

  # --- coverage of the Wilson interval for p* (upper-tri pairs) ---------------
  covP_num <- 0; covP_den <- 0
  covP_same_num <- 0; covP_same_den <- 0
  covP_diff_num <- 0; covP_diff_den <- 0
  widP_sum <- 0
  for (t in seq_len(T_)) {
    up    <- upper.tri(ci[[t]]$lower)
    lo    <- ci[[t]]$lower[up]
    hi    <- ci[[t]]$upper[up]
    pstar <- (pstar_acc[[t]] / R_TRUTH)[up]
    same  <- outer(truth[[t]], truth[[t]], "==")[up]
    covp  <- (lo <= pstar) & (pstar <= hi)
    covP_num <- covP_num + sum(covp); covP_den <- covP_den + length(covp)
    covP_same_num <- covP_same_num + sum(covp[same]); covP_same_den <- covP_same_den + sum(same)
    covP_diff_num <- covP_diff_num + sum(covp[!same]); covP_diff_den <- covP_diff_den + sum(!same)
    widP_sum <- widP_sum + sum(hi - lo)
  }

  width_P_mean <- widP_sum / covP_den
  data.frame(
    task = TASK, deployment = cfg$deployment, spec = cfg$spec, algorithm = cfg$algo,
    block = cfg$block, n = cfg$n, K = cfg$K, separation = cfg$separation,
    p_in = cfg$p_in, p_out = cfg$p_out, p_switch = cfg$p_switch, T_layers = cfg$T_,
    dgp = cfg$dgp, balance = cfg$balance, hetero = cfg$hetero,
    sim = sim_id, B = boot$n_boot, R_truth = R_TRUTH,
    cov_P_mean  = covP_num / covP_den,
    cov_P_same  = if (covP_same_den > 0) covP_same_num / covP_same_den else NA_real_,
    cov_P_diff  = if (covP_diff_den > 0) covP_diff_num / covP_diff_den else NA_real_,
    width_P_mean = width_P_mean,
    gate = isTRUE(width_P_mean < 0.05),
    meta_K_mean = mean(meta_K), true_K_mean = mean(true_K),
    k_mae = mean(abs(meta_K - true_K)),
    kcount_repro_mean = mean(boot$community_count_reproducibility, na.rm = TRUE),
    stringsAsFactors = FALSE)
}

# =============================================================================
# RUN + WRITE  (sequential)
# =============================================================================
t0 <- proc.time()["elapsed"]
rows <- vector("list", M_SIMS)
for (i in seq_len(M_SIMS)) {
  rows[[i]] <- tryCatch(run_one(i), error = function(e) {
    message(sprintf("sim %d failed: %s", i, conditionMessage(e))); NULL
  })
  if (i %% 25L == 0L || MINI)
    cat(sprintf("  ... sim %d/%d done (%.1f min elapsed)\n",
                i, M_SIMS, (proc.time()["elapsed"] - t0) / 60))
}
res <- do.call(rbind, Filter(Negate(is.null), rows))
write.csv(res, outfile, row.names = FALSE)

cat(sprintf(paste0("[meta-cov] task %d done: %d/%d sims ok in %.1f min. ",
                   "cov_P=%.3f width_P=%.4f gate_pass=%.0f%% -> %s\n"),
            TASK, nrow(res), M_SIMS, (proc.time()["elapsed"] - t0) / 60,
            mean(res$cov_P_mean), mean(res$width_P_mean),
            100 * mean(res$gate), outfile))
