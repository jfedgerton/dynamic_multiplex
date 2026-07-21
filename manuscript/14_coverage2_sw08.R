# =============================================================================
# manuscript/14_coverage2_sw08.R
#
# Coverage study v2 ADD-ON: p_switch = 0.80 stress arm (near-independent
# memberships across layers). Same machinery as 13_, single switching value,
# estimand: coverage of Wilson intervals for pairwise co-assignment
# propensities (the co_assignment_ci() feature planned for 1.1.0).
#
# Estimands per simulated dataset
#   K  : planted community count per layer vs community_count_ci   (kept)
#   Q  : planted-partition modularity vs modularity_ci             (kept as a
#        measured quantity; function slated for removal in 1.1.0 pending
#        full-grid evidence)
#   P  : NEW. For each node pair (i,j) and layer t, the bootstrap
#        co-assignment probability p_hat gives a Wilson 95% interval
#        (B replicates). Ground truth p_star = co-clustering propensity:
#        the frequency with which the SAME algorithm pairs (i,j) when the
#        edges are redrawn fresh (R_TRUTH draws) conditional on the same
#        planted membership path. Coverage = share of pairs whose Wilson
#        interval contains p_star; also split by true co-membership.
#   calibration : 20-bin calibration of p_hat vs binary truth (kept)
#
# Costs: M=500 sims (was 1000), each sim runs B=100 bootstrap fits plus
# R_TRUTH=100 truth fits, so per-task wall time is comparable to v1.
#
# Usage (mini test): COV_MODE=main COV_TASK=1 COV_MINI=1 Rscript manuscript/14_coverage2_sw08.R
# =============================================================================

suppressPackageStartupMessages({
  library(dynamicmultiplex)
  library(parallel)
})

MODE   <- Sys.getenv("COV_MODE", "main")
OFFSET <- as.integer(Sys.getenv("COV_OFFSET", "0"))
TASKID <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID",
                                Sys.getenv("COV_TASK", "1")))
TASK   <- OFFSET + TASKID
CORES  <- as.integer(Sys.getenv("COV_CORES", "8"))
MINI   <- identical(Sys.getenv("COV_MINI", "0"), "1")

M_SIMS  <- if (MINI) 4L else 500L
B_BOOT  <- if (MINI) 10L else 100L
R_TRUTH <- if (MINI) 5L  else 100L
ALPHA   <- 0.05
N_BINS  <- 20L
Z_ALPHA <- qnorm(1 - ALPHA / 2)

DENSITIES <- list(weak    = c(p_in = 0.20, p_out = 0.10),
                  default = c(p_in = 0.30, p_out = 0.05),
                  strong  = c(p_in = 0.50, p_out = 0.02))

outdir <- file.path("manuscript", "output",
                    sprintf("coverage2sw08_chunks_%s", MODE))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# --------------------------- task table (identical to v1) --------------------
if (MODE == "main") {
  cells <- expand.grid(
    n = c(50, 100, 200, 400), K = c(3, 5, 10),
    p_switch = c(0.80),
    density = names(DENSITIES), T_ = c(5, 10, 15),
    stringsAsFactors = FALSE)
  cells <- cells[!(cells$K == 10 & cells$n == 50), ]
  specs <- expand.grid(fit_type = c("jaccard", "overlap", "identity"),
                       algorithm = c("louvain", "leiden"),
                       stringsAsFactors = FALSE)
  specs$B <- B_BOOT
  cells$weights <- "binary"
} else if (MODE == "valued") {
  cells <- expand.grid(
    n = c(100, 200), K = c(3, 5), p_switch = c(0.05, 0.20),
    density = "default", T_ = 10,
    weights = c("aligned", "orthogonal"), stringsAsFactors = FALSE)
  specs <- expand.grid(
    fit_type = c("jaccard", "overlap", "weighted_jaccard",
                 "weighted_overlap", "identity"),
    algorithm = c("louvain", "leiden"), stringsAsFactors = FALSE)
  specs$B <- B_BOOT
} else if (MODE == "bsens") {
  cells <- data.frame(n = 100, K = 3, p_switch = 0.05, density = "default",
                      T_ = 10, weights = "binary", stringsAsFactors = FALSE)
  specs <- expand.grid(fit_type = "jaccard",
                       algorithm = c("louvain", "leiden"),
                       B = c(100L, 500L, 1000L), stringsAsFactors = FALSE)
} else stop("Unknown COV_MODE: ", MODE)

task_tab <- merge(cells, specs, by = NULL)
task_tab <- task_tab[order(task_tab$n, task_tab$T_, task_tab$K,
                           task_tab$p_switch, task_tab$density,
                           task_tab$weights, task_tab$fit_type,
                           task_tab$algorithm, task_tab$B), ]
rownames(task_tab) <- NULL
if (TASK < 1 || TASK > nrow(task_tab)) {
  stop(sprintf("Task %d out of range 1..%d for mode %s",
               TASK, nrow(task_tab), MODE))
}
cfg <- task_tab[TASK, ]
cat(sprintf("[coverage2sw08] mode=%s task=%d/%d : n=%d K=%d sw=%.2f dens=%s T=%d w=%s %s/%s B=%d M=%d R=%d\n",
            MODE, TASK, nrow(task_tab), cfg$n, cfg$K, cfg$p_switch,
            cfg$density, cfg$T_, cfg$weights, cfg$fit_type, cfg$algorithm,
            cfg$B, M_SIMS, R_TRUTH))

# --------------------------- simulator ---------------------------------------
gen_memberships <- function(cfg) {
  n <- cfg$n; K <- cfg$K; T_ <- cfg$T_
  memberships <- vector("list", T_)
  memberships[[1]] <- sample(seq_len(K), n, replace = TRUE)
  for (t in 2:T_) {
    prev <- memberships[[t - 1]]
    mask <- runif(n) < cfg$p_switch
    prev[mask] <- sample(seq_len(K), sum(mask), replace = TRUE)
    memberships[[t]] <- prev
  }
  memberships
}

gen_layers <- function(memberships, cfg) {
  dens <- DENSITIES[[cfg$density]]
  n <- cfg$n
  lapply(memberships, function(mem) {
    P  <- outer(mem, mem, "==")
    Pr <- ifelse(P, dens["p_in"], dens["p_out"]); diag(Pr) <- 0
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

true_modularity <- function(A, mem) {
  g <- igraph::graph_from_adjacency_matrix(A, mode = "undirected",
                                           weighted = TRUE, diag = FALSE)
  igraph::modularity(g, membership = mem, weights = igraph::E(g)$weight)
}

FITTERS <- list(
  jaccard          = fit_multilayer_jaccard,
  overlap          = fit_multilayer_overlap,
  identity         = fit_multilayer_identity_ties,
  weighted_jaccard = fit_multilayer_weighted_jaccard,
  weighted_overlap = fit_multilayer_weighted_overlap
)

wilson_bounds <- function(phat, b) {
  z2 <- Z_ALPHA^2
  denom  <- 1 + z2 / b
  center <- (phat + z2 / (2 * b)) / denom
  half   <- Z_ALPHA * sqrt(phat * (1 - phat) / b + z2 / (4 * b^2)) / denom
  list(lo = pmax(0, center - half), hi = pmin(1, center + half))
}

# --------------------------- one simulation ----------------------------------
run_one <- function(sim_id) {
  seed <- (3000L + TASK) * 550000L + sim_id
  set.seed(seed)
  truth  <- gen_memberships(cfg)
  layers <- gen_layers(truth, cfg)
  T_ <- cfg$T_

  boot <- bootstrap_multilayer(layers, fit_type = cfg$fit_type,
                               algorithm = cfg$algorithm,
                               n_boot = cfg$B, seed = seed)
  ci <- community_ci(boot, alpha = ALPHA)

  true_K <- vapply(truth, function(m) length(unique(m)), integer(1))
  true_Q <- vapply(seq_len(T_), function(t)
    true_modularity(layers[[t]], truth[[t]]), numeric(1))
  kci <- ci$community_count_ci
  qci <- ci$modularity_ci
  cov_K <- (kci$lower <= true_K) & (true_K <= kci$upper)
  cov_Q <- (qci$lower <= true_Q) & (true_Q <= qci$upper)

  # ---- ground-truth co-clustering propensity: R_TRUTH fresh edge draws ----
  fitfun <- FITTERS[[cfg$fit_type]]
  n <- cfg$n
  pstar_acc <- lapply(seq_len(T_), function(t) matrix(0, n, n))
  for (r in seq_len(R_TRUTH)) {
    set.seed((3000L + TASK) * 550000L + 1000L + sim_id * 1000L + r)
    fresh <- gen_layers(truth, cfg)
    fit <- fitfun(fresh, algorithm = cfg$algorithm)
    for (t in seq_len(T_)) {
      mem <- fit$layer_communities[[t]]$membership
      pstar_acc[[t]] <- pstar_acc[[t]] + outer(mem, mem, "==")
    }
    rm(fresh, fit)
  }

  # ---- Wilson interval coverage of p_star, pooled over layers -------------
  bin_total <- integer(N_BINS); bin_true <- integer(N_BINS)
  covP_num <- 0; covP_den <- 0
  covP_same_num <- 0; covP_same_den <- 0
  covP_diff_num <- 0; covP_diff_den <- 0
  widP_sum <- 0
  for (t in seq_len(T_)) {
    P     <- ci$co_assignment[[t]]
    up    <- upper.tri(P)
    phat  <- P[up]
    pstar <- (pstar_acc[[t]] / R_TRUTH)[up]
    same  <- outer(truth[[t]], truth[[t]], "==")[up]
    wb    <- wilson_bounds(phat, cfg$B)
    covp  <- (wb$lo <= pstar) & (pstar <= wb$hi)
    covP_num <- covP_num + sum(covp);        covP_den <- covP_den + length(covp)
    covP_same_num <- covP_same_num + sum(covp[same])
    covP_same_den <- covP_same_den + sum(same)
    covP_diff_num <- covP_diff_num + sum(covp[!same])
    covP_diff_den <- covP_diff_den + sum(!same)
    widP_sum <- widP_sum + sum(wb$hi - wb$lo)
    bin  <- pmin(N_BINS, 1L + floor(phat * N_BINS))
    bin_total <- bin_total + tabulate(bin, nbins = N_BINS)
    bin_true  <- bin_true  + tabulate(bin[same], nbins = N_BINS)
  }

  row <- data.frame(
    mode = MODE, task = TASK, sim = sim_id,
    n = cfg$n, K = cfg$K, p_switch = cfg$p_switch, density = cfg$density,
    T_layers = cfg$T_, weights = cfg$weights,
    fit_type = cfg$fit_type, algorithm = cfg$algorithm, B = cfg$B,
    R_truth = R_TRUTH,
    cov_K_mean = mean(cov_K), cov_Q_mean = mean(cov_Q),
    cov_P_mean = covP_num / covP_den,
    cov_P_same = ifelse(covP_same_den > 0, covP_same_num / covP_same_den, NA),
    cov_P_diff = ifelse(covP_diff_den > 0, covP_diff_num / covP_diff_den, NA),
    width_K_mean = mean(kci$upper - kci$lower),
    width_Q_mean = mean(qci$upper - qci$lower),
    width_P_mean = widP_sum / covP_den,
    stringsAsFactors = FALSE
  )
  out <- list(row = row, bin_total = bin_total, bin_true = bin_true)
  rm(truth, layers, boot, ci, kci, qci, cov_K, cov_Q, pstar_acc)
  gc(verbose = FALSE)
  out
}

# --------------------------- run + write chunks ------------------------------
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
calib <- data.frame(
  mode = MODE, task = TASK,
  n = cfg$n, K = cfg$K, p_switch = cfg$p_switch, density = cfg$density,
  T_layers = cfg$T_, weights = cfg$weights,
  fit_type = cfg$fit_type, algorithm = cfg$algorithm, B = cfg$B,
  bin = seq_len(N_BINS),
  bin_lo = (seq_len(N_BINS) - 1) / N_BINS, bin_hi = seq_len(N_BINS) / N_BINS,
  n_pairs = btot, n_true = btrue, stringsAsFactors = FALSE)

write.csv(rows, file.path(outdir, sprintf("cov_task%05d.csv", TASK)),
          row.names = FALSE)
write.csv(calib, file.path(outdir, sprintf("calib_task%05d.csv", TASK)),
          row.names = FALSE)

cat(sprintf("[coverage2sw08] task %d done: %d/%d sims ok, %.1f min. cov_K=%.3f cov_Q=%.3f cov_P=%.3f\n",
            TASK, nrow(rows), M_SIMS, (proc.time()["elapsed"] - t0) / 60,
            mean(rows$cov_K_mean), mean(rows$cov_Q_mean),
            mean(rows$cov_P_mean)))
