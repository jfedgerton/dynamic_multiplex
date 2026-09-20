# =============================================================================
# replication/sim/02_stability.R  --  bootstrap stability vs accuracy (Section 4)
#
# Does the bootstrap STABILITY of a DynMux partition predict its ACCURACY?
# For every simulated series: fit DynMux (Jaccard, Leiden), redraw B bootstrap
# networks from the block densities of the fitted partition (the scheme in
# bootstrap_multilayer()), refit, and record
#   acc_nmi / acc_ari     mean over layers of NMI / ARI(point partition, truth)
#   stab_nmi / stab_ari   mean over replicates and layers of NMI / ARI(replicate, point)
#   node sample (first NODE_SIMS sims, NODE_PER_SIM nodes x layers):
#     stab = mean_b Jaccard(C_b(i), C_0(i)),  acc = Jaccard(C_0(i), C_true(i))
#   pair sample (same sims): co-assignment share, whether the pair is truly
#     together, so decided / undetermined shares can be validated
# post/12_stability.R turns the binary arm into the calibration table shipped
# with the packages and validates it on the held-out half and on the two
# robustness arms.
#
# Arms (STAB_ARM):
#   binary   planted partition, 594 cells: n{50,100,200,400} x K{3,5,10} x
#            p_switch{0.02..0.80} x density{weak,default,strong} x T{5,10,15}
#            (K=10 with n=50 excluded)                       -- calibration arm
#   dcsbm    degree-corrected generator (theta_i theta_j scaling; hetero none /
#            moderate lognormal / severe Pareto; balanced or skewed sizes),
#            216 cells; the bootstrap still redraws a plain SBM  -- misspecified
#   weighted lognormal edge weights aligned with or orthogonal to the blocks,
#            72 cells; weighted Jaccard coupling, weights resampled by block
# M = 50 simulations per cell, B = 100 replicates. Seeds: (base + TASK) * 100000
# + sim with base 12000 (binary), 13000 (dcsbm), 14000 (weighted).
#
# Usage (mini): STAB_ARM=binary COV_TASK=1 COV_MINI=1 DM_ROOT=. Rscript replication/sim/02_stability.R
# Arrays: binary 1..594, dcsbm 1..216, weighted 1..72 (slurm/02a-c).
# Output: $DM_ROOT/output/stability/<arm>_stab_task%05d.csv, <arm>_node_task%05d.csv,
#         <arm>_pair_task%05d.csv
# =============================================================================
suppressPackageStartupMessages({ library(dynamicmultiplex); library(parallel); library(igraph) })

ARM   <- match.arg(Sys.getenv("STAB_ARM", "binary"), c("binary", "dcsbm", "weighted"))
TASK  <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", Sys.getenv("COV_TASK", "1")))
CORES <- as.integer(Sys.getenv("COV_CORES", "8"))
MINI  <- identical(Sys.getenv("COV_MINI", "0"), "1")
M_SIMS       <- if (MINI) 3L else 50L
B_BOOT       <- if (MINI) 8L else 100L
NODE_SIMS    <- if (MINI) 2L else 10L
NODE_PER_SIM <- if (MINI) 10L else 40L
PAIR_PER_SIM <- if (MINI) 20L else 200L
SEED_BASE <- c(binary = 12000L, dcsbm = 13000L, weighted = 14000L)[[ARM]]   # keeps (base + task) * 1e5 + sim inside int32

DENSITIES <- list(weak = c(p_in = 0.20, p_out = 0.10), default = c(p_in = 0.30, p_out = 0.05),
                  strong = c(p_in = 0.50, p_out = 0.02))
ROOT   <- Sys.getenv("DM_ROOT", unset = getwd())
outdir <- file.path(ROOT, "output", "stability"); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- cell grids (one row per array task, shuffled with seed 123) ------------
if (ARM == "binary") {
  cells <- expand.grid(n = c(50, 100, 200, 400), K = c(3, 5, 10), p_switch = c(0.02, 0.05, 0.10, 0.20, 0.50, 0.80),
                       density = names(DENSITIES), T_ = c(5, 10, 15), stringsAsFactors = FALSE)
  cells <- cells[!(cells$K == 10 & cells$n == 50), ]
  cells$hetero <- "none"; cells$balance <- "balanced"; cells$weights <- "none"
  cells <- cells[order(cells$n, cells$T_, cells$K, cells$p_switch, cells$density), ]
  stopifnot(nrow(cells) == 594)
} else if (ARM == "dcsbm") {
  cells <- expand.grid(n = c(100, 200, 400), K = c(3, 5), p_switch = c(0.05, 0.20), density = names(DENSITIES), T_ = 10,
                       hetero = c("none", "moderate", "severe"), balance = c("balanced", "skewed"), stringsAsFactors = FALSE)
  cells$weights <- "none"
  cells <- cells[order(cells$n, cells$K, cells$p_switch, cells$density, cells$hetero, cells$balance), ]
  stopifnot(nrow(cells) == 216)
} else {
  cells <- expand.grid(n = c(100, 200, 400), K = c(3, 5), p_switch = c(0.05, 0.20, 0.80), density = c("default", "strong"), T_ = 10,
                       weights = c("aligned", "orthogonal"), stringsAsFactors = FALSE)
  cells$hetero <- "none"; cells$balance <- "balanced"
  cells <- cells[order(cells$n, cells$K, cells$p_switch, cells$density, cells$weights), ]
  stopifnot(nrow(cells) == 72)
}
rownames(cells) <- NULL
set.seed(123); perm <- sample.int(nrow(cells))
if (TASK < 1 || TASK > nrow(cells)) stop("task out of range for arm ", ARM)
outfiles <- file.path(outdir, sprintf(c("%s_stab_task%05d.csv", "%s_node_task%05d.csv", "%s_pair_task%05d.csv"), ARM, TASK))
if (nzchar(Sys.getenv("SLURM_ARRAY_TASK_ID")) && all(file.exists(outfiles))) {
  cat("[skip] task", TASK, "already complete\n"); quit(save = "no")
}
cfg <- as.list(cells[perm[TASK], ])
cfg$p_in <- DENSITIES[[cfg$density]][["p_in"]]; cfg$p_out <- DENSITIES[[cfg$density]][["p_out"]]
cat(sprintf("[stab:%s] task=%d/%d cell=%d : n=%d K=%d sw=%.2f dens=%s T=%d hetero=%s bal=%s weights=%s M=%d B=%d\n",
            ARM, TASK, nrow(cells), perm[TASK], cfg$n, cfg$K, cfg$p_switch, cfg$density, cfg$T_, cfg$hetero, cfg$balance,
            cfg$weights, M_SIMS, B_BOOT))
fitfun <- if (ARM == "weighted") fit_multilayer_weighted_jaccard else fit_multilayer_jaccard

# ---- generators ---------------------------------------------------------------
size_probs <- function(K, balance) { if (balance == "balanced") return(rep(1 / K, K)); w <- 0.6 ^ (seq_len(K) - 1); w / sum(w) }
gen_theta <- function(n, hetero) {                       # per-node degree propensity, E[theta] = 1
  if (hetero == "none") return(rep(1, n))
  if (hetero == "moderate") return(rlnorm(n, meanlog = -0.6^2 / 2, sdlog = 0.6))
  raw <- runif(n) ^ (-1 / 2.5); raw / (2.5 / 1.5)          # Pareto(2.5, xmin 1) / mean
}
gen_memberships <- function(cfg, probs) {
  m <- vector("list", cfg$T_); m[[1]] <- sample(seq_len(cfg$K), cfg$n, replace = TRUE, prob = probs)
  for (t in 2:cfg$T_) { prev <- m[[t - 1]]; mask <- runif(cfg$n) < cfg$p_switch
    prev[mask] <- sample(seq_len(cfg$K), sum(mask), replace = TRUE, prob = probs); m[[t]] <- prev }
  m
}
gen_layers <- function(memberships, cfg, theta) lapply(memberships, function(mem) {
  n <- cfg$n; P <- outer(mem, mem, "=="); Pr <- ifelse(P, cfg$p_in, cfg$p_out) * outer(theta, theta)
  Pr[Pr > 1] <- 1; diag(Pr) <- 0
  A <- matrix(0, n, n); up <- upper.tri(Pr); A[up] <- rbinom(sum(up), 1, Pr[up])
  if (cfg$weights == "aligned") {            # lognormal weights, mean log 1 within / 0 between blocks
    W <- matrix(0, n, n); mu <- ifelse(P, 1.0, 0.0); W[up] <- A[up] * rlnorm(sum(up), meanlog = mu[up], sdlog = 0.75); A <- W
  } else if (cfg$weights == "orthogonal") {  # weights carry no block signal
    W <- matrix(0, n, n); W[up] <- A[up] * rlnorm(sum(up), meanlog = 0.5, sdlog = 0.75); A <- W
  }
  A + t(A) })
get_memberships <- function(fit, T_) lapply(seq_len(T_), function(t) as.integer(fit$meta_communities[[t]]))

# ---- bootstrap: block-density redraw from the detected partition (as in the
# package); weights resampled from the observed within / between pools --------
fit_models <- function(layers, mem_hat) lapply(seq_along(layers), function(t) {
  A <- layers[[t]]; mem <- mem_hat[[t]]; same <- outer(mem, mem, "=="); sel <- upper.tri(A); E <- A > 0
  p_all <- mean(E[sel]); in_d <- sel & same; out_d <- sel & !same
  w_all <- A[sel & E]; if (!length(w_all)) w_all <- 1
  w_in <- A[in_d & E]; if (!length(w_in)) w_in <- w_all
  w_out <- A[out_d & E]; if (!length(w_out)) w_out <- w_all
  list(same = same, p_in = if (any(in_d)) mean(E[in_d]) else p_all, p_out = if (any(out_d)) mean(E[out_d]) else p_all,
       w_in = w_in, w_out = w_out) })
redraw <- function(models, n, weighted) lapply(models, function(em) {
  probs <- ifelse(em$same, em$p_in, em$p_out); up <- upper.tri(probs)
  M <- matrix(0, n, n); on <- rbinom(sum(up), 1, probs[up]) == 1
  if (weighted) { w <- numeric(sum(up)); s <- em$same[up]
    w[on & s] <- sample(em$w_in, sum(on & s), replace = TRUE); w[on & !s] <- sample(em$w_out, sum(on & !s), replace = TRUE)
    M[up] <- w } else M[up] <- as.numeric(on)
  M + t(M) })

nmi <- function(a, b) igraph::compare(a, b, method = "nmi")
ari <- function(a, b) igraph::compare(a, b, method = "adjusted.rand")
node_jaccard <- function(mem_a, mem_b) {
  same_a <- outer(mem_a, mem_a, "=="); same_b <- outer(mem_b, mem_b, "=="); rowSums(same_a & same_b) / rowSums(same_a | same_b)
}

run_one <- function(sim_id) {
  seed <- (SEED_BASE + TASK) * 100000L + sim_id; set.seed(seed)
  n <- cfg$n; T_ <- cfg$T_
  theta <- gen_theta(n, cfg$hetero)
  truth <- gen_memberships(cfg, size_probs(cfg$K, cfg$balance)); layers <- gen_layers(truth, cfg, theta)
  fit0 <- fitfun(layers, algorithm = "leiden"); mem0 <- get_memberships(fit0, T_)
  acc_nmi <- mean(vapply(seq_len(T_), function(t) nmi(mem0[[t]], truth[[t]]), numeric(1)))
  acc_ari <- mean(vapply(seq_len(T_), function(t) ari(mem0[[t]], truth[[t]]), numeric(1)))
  node_acc <- lapply(seq_len(T_), function(t) node_jaccard(mem0[[t]], truth[[t]]))
  models <- fit_models(layers, mem0)
  s_nmi <- numeric(0); s_ari <- numeric(0); b_ok <- 0L
  node_stab_acc <- lapply(seq_len(T_), function(t) numeric(n))
  coassign <- lapply(seq_len(T_), function(t) matrix(0, n, n))
  for (b in seq_len(B_BOOT)) {
    bfit <- tryCatch(fitfun(redraw(models, n, ARM == "weighted"), algorithm = "leiden"), error = function(e) NULL)
    if (is.null(bfit)) next
    b_ok <- b_ok + 1L; bmem <- get_memberships(bfit, T_)
    s_nmi <- c(s_nmi, mean(vapply(seq_len(T_), function(t) nmi(bmem[[t]], mem0[[t]]), numeric(1))))
    s_ari <- c(s_ari, mean(vapply(seq_len(T_), function(t) ari(bmem[[t]], mem0[[t]]), numeric(1))))
    for (t in seq_len(T_)) { node_stab_acc[[t]] <- node_stab_acc[[t]] + node_jaccard(bmem[[t]], mem0[[t]])
      if (sim_id <= NODE_SIMS) coassign[[t]] <- coassign[[t]] + outer(bmem[[t]], bmem[[t]], "==") }
  }
  if (b_ok < 10L && !MINI) stop("too few completed bootstrap replicates")
  row <- data.frame(arm = ARM, task = TASK, sim = sim_id, n = n, K = cfg$K, p_switch = cfg$p_switch, density = cfg$density,
                    T_layers = T_, hetero = cfg$hetero, balance = cfg$balance, weights = cfg$weights, B = b_ok,
                    acc_nmi = acc_nmi, acc_ari = acc_ari, stab_nmi = mean(s_nmi), stab_ari = mean(s_ari),
                    stab_nmi_sd = sd(s_nmi), K_true_mean = mean(vapply(truth, function(m) length(unique(m)), integer(1))),
                    K_hat_mean = mean(vapply(mem0, function(m) length(unique(m)), integer(1))), stringsAsFactors = FALSE)
  nodes <- NULL; pairs <- NULL
  if (sim_id <= NODE_SIMS) {
    pick <- sample.int(n, min(NODE_PER_SIM, n))
    nodes <- do.call(rbind, lapply(seq_len(T_), function(t) data.frame(arm = ARM, task = TASK, sim = sim_id, layer = t, node = pick,
      n = n, K = cfg$K, p_switch = cfg$p_switch, density = cfg$density, T_layers = T_, hetero = cfg$hetero, balance = cfg$balance,
      weights = cfg$weights, stab = node_stab_acc[[t]][pick] / b_ok, acc = node_acc[[t]][pick],
      degree = rowSums(layers[[t]] > 0)[pick], stringsAsFactors = FALSE)))
    pairs <- do.call(rbind, lapply(seq_len(T_), function(t) {
      ij <- cbind(sample.int(n, PAIR_PER_SIM, replace = TRUE), sample.int(n, PAIR_PER_SIM, replace = TRUE)); ij <- ij[ij[, 1] != ij[, 2], , drop = FALSE]
      data.frame(arm = ARM, task = TASK, sim = sim_id, layer = t, n = n, K = cfg$K, p_switch = cfg$p_switch, density = cfg$density,
                 T_layers = T_, hetero = cfg$hetero, balance = cfg$balance, weights = cfg$weights,
                 share = coassign[[t]][ij] / b_ok, point_same = as.integer(mem0[[t]][ij[, 1]] == mem0[[t]][ij[, 2]]),
                 true_same = as.integer(truth[[t]][ij[, 1]] == truth[[t]][ij[, 2]]), stringsAsFactors = FALSE) }))
  }
  list(row = row, nodes = nodes, pairs = pairs)
}

t0 <- proc.time()[["elapsed"]]
res <- mclapply(seq_len(M_SIMS), function(i) tryCatch(run_one(i), error = function(e) {
  message(sprintf("sim %d failed: %s", i, conditionMessage(e))); NULL }), mc.cores = CORES, mc.preschedule = FALSE)
res <- Filter(Negate(is.null), res)
stopifnot(length(res) >= (if (MINI) 1L else 10L))
rows <- do.call(rbind, lapply(res, `[[`, "row")); nodes <- do.call(rbind, lapply(res, `[[`, "nodes")); pairs <- do.call(rbind, lapply(res, `[[`, "pairs"))
stopifnot(all(rows$B >= 1), all(is.finite(rows$stab_nmi)), all(is.finite(rows$acc_nmi)))
write.csv(rows, outfiles[1], row.names = FALSE); write.csv(nodes, outfiles[2], row.names = FALSE); write.csv(pairs, outfiles[3], row.names = FALSE)
cat(sprintf("[stab:%s] task %d done: %d sims, %.1f min; acc_nmi=%.3f stab_nmi=%.3f cor=%.3f\n", ARM, TASK, nrow(rows),
            (proc.time()[["elapsed"]] - t0) / 60, mean(rows$acc_nmi), mean(rows$stab_nmi), suppressWarnings(cor(rows$acc_nmi, rows$stab_nmi))))
