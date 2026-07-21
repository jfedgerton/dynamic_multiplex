# =============================================================================
# manuscript/12_bootstrap_coverage.R
#
# Coverage study for bootstrap_multilayer() / community_ci().
# For each (design cell x model spec), simulate M datasets from the switching
# DGP, run the B-replicate bootstrap on each, and record whether the nominal
# 95% CIs cover the true value. Also accumulates calibration bins for the
# pairwise co-assignment probabilities.
#
# Estimands
#   K  : true number of distinct planted communities per layer
#        vs community_count_ci
#   Q  : modularity of the PLANTED partition on the REALIZED graph
#        vs modularity_ci
#   co-assignment: binned calibration of Pr(same community | bootstrap)
#        against true co-membership (binary truth), 20 bins
#
# Modes (env var COV_MODE)
#   main   : binary networks, 6 specs ({jaccard, overlap, identity} x
#            {louvain, leiden}), full factorial grid minus K=10 & n=50
#   valued : log-normal edge weights (aligned + orthogonal regimes),
#            10 specs (5 similarities x 2 algorithms), reduced grid
#   bsens  : central cell, jaccard x {louvain, leiden}, B in {100, 500, 1000}
#
# Task selection (slurm array): task id = COV_OFFSET + SLURM_ARRAY_TASK_ID.
# Each task = one (cell x spec) combination; M sims parallelized over
# COV_CORES forked workers; results written as one CSV chunk per task.
#
# Memory: per-sim objects are removed and gc() runs every iteration; only
# summary rows and binned counts are accumulated.
#
# Usage (local mini-test):
#   COV_MODE=main COV_TASK=1 COV_MINI=1 Rscript manuscript/12_bootstrap_coverage.R
# =============================================================================

suppressPackageStartupMessages({
  library(dynamicmultiplex)
  library(parallel)
})

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
MODE   <- Sys.getenv("COV_MODE", "main")
OFFSET <- as.integer(Sys.getenv("COV_OFFSET", "0"))
TASKID <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID",
                                Sys.getenv("COV_TASK", "1")))
TASK   <- OFFSET + TASKID
CORES  <- as.integer(Sys.getenv("COV_CORES", "8"))
MINI   <- identical(Sys.getenv("COV_MINI", "0"), "1")

M_SIMS <- if (MINI) 4L else 1000L
B_BOOT <- if (MINI) 10L else 100L
ALPHA  <- 0.05
N_BINS <- 20L

DENSITIES <- list(weak    = c(p_in = 0.20, p_out = 0.10),
                  default = c(p_in = 0.30, p_out = 0.05),
                  strong  = c(p_in = 0.50, p_out = 0.02))

outdir <- file.path("manuscript", "output",
                    sprintf("coverage_chunks_%s", MODE))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Build the deterministic task table for each mode
# -----------------------------------------------------------------------------
if (MODE == "main") {
  cells <- expand.grid(
    n        = c(50, 100, 200, 400),
    K        = c(3, 5, 10),
    p_switch = c(0.02, 0.05, 0.10, 0.20, 0.50),
    density  = names(DENSITIES),
    T_       = c(5, 10, 15),
    stringsAsFactors = FALSE
  )
  cells <- cells[!(cells$K == 10 & cells$n == 50), ]
  specs <- expand.grid(
    fit_type  = c("jaccard", "overlap", "identity"),
    algorithm = c("louvain", "leiden"),
    stringsAsFactors = FALSE
  )
  specs$B <- B_BOOT
  cells$weights <- "binary"
} else if (MODE == "valued") {
  cells <- expand.grid(
    n        = c(100, 200),
    K        = c(3, 5),
    p_switch = c(0.05, 0.20),
    density  = "default",
    T_       = 10,
    weights  = c("aligned", "orthogonal"),
    stringsAsFactors = FALSE
  )
  specs <- expand.grid(
    fit_type  = c("jaccard", "overlap", "weighted_jaccard",
                  "weighted_overlap", "identity"),
    algorithm = c("louvain", "leiden"),
    stringsAsFactors = FALSE
  )
  specs$B <- B_BOOT
} else if (MODE == "bsens") {
  cells <- data.frame(n = 100, K = 3, p_switch = 0.05,
                      density = "default", T_ = 10, weights = "binary",
                      stringsAsFactors = FALSE)
  specs <- expand.grid(
    fit_type  = "jaccard",
    algorithm = c("louvain", "leiden"),
    B         = c(100L, 500L, 1000L),
    stringsAsFactors = FALSE
  )
} else stop("Unknown COV_MODE: ", MODE)

task_tab <- merge(cells, specs, by = NULL)  # full crossing, deterministic order
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
cat(sprintf("[coverage] mode=%s task=%d/%d : n=%d K=%d sw=%.2f dens=%s T=%d w=%s %s/%s B=%d M=%d\n",
            MODE, TASK, nrow(task_tab), cfg$n, cfg$K, cfg$p_switch,
            cfg$density, cfg$T_, cfg$weights, cfg$fit_type, cfg$algorithm,
            cfg$B, M_SIMS))

# -----------------------------------------------------------------------------
# Simulator: switching DGP, optionally with log-normal edge weights
# -----------------------------------------------------------------------------
simulate_cell <- function(cfg, seed) {
  set.seed(seed)
  dens <- DENSITIES[[cfg$density]]
  n <- cfg$n; K <- cfg$K; T_ <- cfg$T_
  memberships <- vector("list", T_)
  memberships[[1]] <- sample(seq_len(K), n, replace = TRUE)
  for (t in 2:T_) {
    prev <- memberships[[t - 1]]
    mask <- runif(n) < cfg$p_switch
    prev[mask] <- sample(seq_len(K), sum(mask), replace = TRUE)
    memberships[[t]] <- prev
  }
  layers <- lapply(memberships, function(mem) {
    P  <- outer(mem, mem, "==")
    Pr <- ifelse(P, dens["p_in"], dens["p_out"]); diag(Pr) <- 0
    A  <- matrix(0, n, n)
    up <- upper.tri(Pr)
    A[up] <- rbinom(sum(up), 1, Pr[up])
    A  <- A + t(A)
    if (cfg$weights == "aligned") {
      # within-community edges heavier (log-normal, trade-flow style)
      W <- matrix(0, n, n)
      mu <- ifelse(P, 1.0, 0.0)     # mu_in = 1, mu_out = 0 on log scale
      W[up] <- A[up] * rlnorm(sum(up), meanlog = mu[up], sdlog = 0.75)
      A <- W + t(W)
    } else if (cfg$weights == "orthogonal") {
      # weights carry no community signal (adversarial regime)
      W <- matrix(0, n, n)
      W[up] <- A[up] * rlnorm(sum(up), meanlog = 0.5, sdlog = 0.75)
      A <- W + t(W)
    }
    A
  })
  list(layers = layers, truth = memberships)
}

true_modularity <- function(A, mem) {
  g <- igraph::graph_from_adjacency_matrix(A, mode = "undirected",
                                           weighted = TRUE, diag = FALSE)
  w <- igraph::E(g)$weight
  igraph::modularity(g, membership = mem, weights = w)
}

# -----------------------------------------------------------------------------
# One simulation draw -> coverage row + calibration bin counts
# -----------------------------------------------------------------------------
run_one <- function(sim_id) {
  seed <- TASK * 1000000L + sim_id
  sim <- simulate_cell(cfg, seed)

  boot <- bootstrap_multilayer(
    sim$layers,
    fit_type  = cfg$fit_type,
    algorithm = cfg$algorithm,
    n_boot    = cfg$B,
    seed      = seed
  )
  ci <- community_ci(boot, alpha = ALPHA)

  T_ <- cfg$T_
  true_K <- vapply(sim$truth, function(m) length(unique(m)), integer(1))
  true_Q <- vapply(seq_len(T_), function(t)
    true_modularity(sim$layers[[t]], sim$truth[[t]]), numeric(1))

  kci <- ci$community_count_ci
  qci <- ci$modularity_ci
  cov_K <- (kci$lower <= true_K) & (true_K <= kci$upper)
  cov_Q <- (qci$lower <= true_Q) & (true_Q <= qci$upper)

  # Calibration bins for co-assignment probabilities (upper triangle, all
  # layers pooled): counts of pairs per bin and how many truly co-assigned
  bin_total <- integer(N_BINS)
  bin_true  <- integer(N_BINS)
  for (t in seq_len(T_)) {
    P    <- ci$co_assignment[[t]]
    up   <- upper.tri(P)
    phat <- P[up]
    same <- outer(sim$truth[[t]], sim$truth[[t]], "==")[up]
    bin  <- pmin(N_BINS, 1L + floor(phat * N_BINS))
    bin_total <- bin_total + tabulate(bin, nbins = N_BINS)
    bin_true  <- bin_true  + tabulate(bin[same], nbins = N_BINS)
  }

  row <- data.frame(
    mode = MODE, task = TASK, sim = sim_id,
    n = cfg$n, K = cfg$K, p_switch = cfg$p_switch, density = cfg$density,
    T_layers = cfg$T_, weights = cfg$weights,
    fit_type = cfg$fit_type, algorithm = cfg$algorithm, B = cfg$B,
    cov_K_mean = mean(cov_K), cov_Q_mean = mean(cov_Q),
    n_layers_cov_K = sum(cov_K), n_layers_cov_Q = sum(cov_Q),
    width_K_mean = mean(kci$upper - kci$lower),
    width_Q_mean = mean(qci$upper - qci$lower),
    stringsAsFactors = FALSE
  )

  out <- list(row = row, bin_total = bin_total, bin_true = bin_true)
  rm(sim, boot, ci, kci, qci, cov_K, cov_Q); gc(verbose = FALSE)
  out
}

# -----------------------------------------------------------------------------
# Run all sims for this task (forked parallelism), then write chunk files
# -----------------------------------------------------------------------------
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
  n_pairs = btot, n_true = btrue,
  stringsAsFactors = FALSE
)

write.csv(rows, file.path(outdir, sprintf("cov_task%05d.csv", TASK)),
          row.names = FALSE)
write.csv(calib, file.path(outdir, sprintf("calib_task%05d.csv", TASK)),
          row.names = FALSE)

cat(sprintf("[coverage] task %d done: %d/%d sims ok, %.1f min. cov_K=%.3f cov_Q=%.3f\n",
            TASK, nrow(rows) , M_SIMS,
            (proc.time()["elapsed"] - t0) / 60,
            mean(rows$cov_K_mean), mean(rows$cov_Q_mean)))
