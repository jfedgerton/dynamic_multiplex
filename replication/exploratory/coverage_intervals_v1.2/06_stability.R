# =============================================================================
# replication/sim/06_stability.R   (EXPLORATORY ARM: does bootstrap STABILITY
# predict ACCURACY at the partition and node level?)
#
# Decision rule pre-registered 2026-09-20 (see post/16_stability.R):
#   partition level: Spearman(stability, accuracy) >= 0.7 on validation
#     configurations AND calibrated 5th percentile of accuracy given
#     stability is monotone and >= 0.8 when stability >= 0.9
#   node level: the same, for node stability vs node accuracy
# If neither holds, the bootstrap has no inferential product beyond the
# pair-level decided/undetermined split and the uncertainty section is cut.
#
# Design: the 594-cell grid of sim/02, jaccard/leiden, M = 50 sims per cell,
# B = 100 bootstrap replicates (block-density redraw from the detected
# partition, the scheme in the package). Per sim:
#   acc_nmi   mean over layers of NMI(point-estimate partition, planted truth)
#   acc_ari   same with ARI
#   stab_nmi  mean over replicates and layers of NMI(replicate, point estimate)
#   stab_ari  same with ARI
#   K_hat, K_true, K_cov (percentile interval on replicate counts covers K_true)
#   node sample (first NODE_SIMS sims, NODE_PER_SIM nodes x all layers):
#     s_i = mean_b Jaccard(C_b(i), C_0(i))   node stability
#     a_i = Jaccard(C_0(i), C_true(i))       node accuracy
#     plus the node's degree, detected community size, design columns
#
# Outputs per task under $DM_ROOT/output/stability/:
#   stab_task%05d.csv   one row per sim
#   node_task%05d.csv   node sample
# Usage (mini): COV_TASK=1 COV_MINI=1 DM_ROOT=. Rscript replication/sim/06_stability.R
# Array: 1..594; seeds (20000 + TASK) * 100000 + sim.
# =============================================================================
suppressPackageStartupMessages({ library(dynamicmultiplex); library(parallel); library(igraph) })

TASK  <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", Sys.getenv("COV_TASK", "1")))
CORES <- as.integer(Sys.getenv("COV_CORES", "8"))
MINI  <- identical(Sys.getenv("COV_MINI", "0"), "1")
M_SIMS       <- if (MINI) 3L else 50L
B_BOOT       <- if (MINI) 8L else 100L
NODE_SIMS    <- if (MINI) 2L else 10L
NODE_PER_SIM <- if (MINI) 10L else 40L
ALPHA <- 0.05

DENSITIES <- list(weak = c(p_in = 0.20, p_out = 0.10), default = c(p_in = 0.30, p_out = 0.05),
                  strong = c(p_in = 0.50, p_out = 0.02))
ROOT   <- Sys.getenv("DM_ROOT", unset = getwd())
outdir <- file.path(ROOT, "output", "stability")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
cells <- expand.grid(n = c(50, 100, 200, 400), K = c(3, 5, 10),
                     p_switch = c(0.02, 0.05, 0.10, 0.20, 0.50, 0.80),
                     density = names(DENSITIES), T_ = c(5, 10, 15), stringsAsFactors = FALSE)
cells <- cells[!(cells$K == 10 & cells$n == 50), ]
cells <- cells[order(cells$n, cells$T_, cells$K, cells$p_switch, cells$density), ]
rownames(cells) <- NULL; stopifnot(nrow(cells) == 594)
set.seed(123); perm <- sample.int(nrow(cells))
if (TASK < 1 || TASK > nrow(cells)) stop("task out of range")
outfiles <- file.path(outdir, sprintf(c("stab_task%05d.csv", "node_task%05d.csv"), TASK))
if (nzchar(Sys.getenv("SLURM_ARRAY_TASK_ID")) && all(file.exists(outfiles))) {
  cat("[skip] task", TASK, "already complete\n"); quit(save = "no")
}
cfg <- as.list(cells[perm[TASK], ])
cfg$p_in <- DENSITIES[[cfg$density]][["p_in"]]; cfg$p_out <- DENSITIES[[cfg$density]][["p_out"]]
cat(sprintf("[stab] task=%d/594 cell=%d : n=%d K=%d sw=%.2f dens=%s T=%d M=%d B=%d\n",
            TASK, perm[TASK], cfg$n, cfg$K, cfg$p_switch, cfg$density, cfg$T_, M_SIMS, B_BOOT))
fitfun <- fit_multilayer_jaccard

gen_memberships <- function(cfg) {
  m <- vector("list", cfg$T_); m[[1]] <- sample(seq_len(cfg$K), cfg$n, replace = TRUE)
  for (t in 2:cfg$T_) { prev <- m[[t - 1]]; mask <- runif(cfg$n) < cfg$p_switch
    prev[mask] <- sample(seq_len(cfg$K), sum(mask), replace = TRUE); m[[t]] <- prev }
  m
}
gen_layers <- function(memberships, cfg) lapply(memberships, function(mem) {
  n <- cfg$n; P <- outer(mem, mem, "=="); Pr <- ifelse(P, cfg$p_in, cfg$p_out); diag(Pr) <- 0
  A <- matrix(0, n, n); up <- upper.tri(Pr); A[up] <- rbinom(sum(up), 1, Pr[up]); A + t(A) })
get_memberships <- function(fit, T_) lapply(seq_len(T_), function(t) as.integer(fit$layer_communities[[t]]$membership))
fit_sbm <- function(layers, mem_hat) lapply(seq_along(layers), function(t) {
  A <- layers[[t]]; mem <- mem_hat[[t]]; same <- outer(mem, mem, "=="); sel <- upper.tri(A); E <- A > 0
  p_all <- mean(E[sel]); in_d <- sel & same; out_d <- sel & !same
  list(same = same, p_in = if (any(in_d)) mean(E[in_d]) else p_all, p_out = if (any(out_d)) mean(E[out_d]) else p_all) })
redraw_sbm <- function(models, n) lapply(models, function(em) {
  probs <- ifelse(em$same, em$p_in, em$p_out); up <- upper.tri(probs)
  M <- matrix(0, n, n); M[up] <- rbinom(sum(up), 1, probs[up]); M + t(M) })
nmi <- function(a, b) igraph::compare(a, b, method = "nmi")
ari <- function(a, b) igraph::compare(a, b, method = "adjusted.rand")
node_jaccard <- function(mem_a, mem_b) {           # per node: Jaccard of its community in a vs in b
  n <- length(mem_a); out <- numeric(n)
  same_a <- outer(mem_a, mem_a, "=="); same_b <- outer(mem_b, mem_b, "==")
  inter <- rowSums(same_a & same_b); uni <- rowSums(same_a | same_b)
  inter / uni
}

run_one <- function(sim_id) {
  seed <- (20000L + TASK) * 100000L + sim_id; set.seed(seed)
  n <- cfg$n; T_ <- cfg$T_
  truth <- gen_memberships(cfg); layers <- gen_layers(truth, cfg)
  fit0 <- fitfun(layers, algorithm = "leiden"); mem0 <- get_memberships(fit0, T_)
  acc_nmi <- mean(vapply(seq_len(T_), function(t) nmi(mem0[[t]], truth[[t]]), numeric(1)))
  acc_ari <- mean(vapply(seq_len(T_), function(t) ari(mem0[[t]], truth[[t]]), numeric(1)))
  node_acc <- lapply(seq_len(T_), function(t) node_jaccard(mem0[[t]], truth[[t]]))
  models <- fit_sbm(layers, mem0)
  s_nmi <- numeric(0); s_ari <- numeric(0); kcount <- matrix(NA_integer_, B_BOOT, T_)
  node_stab_acc <- lapply(seq_len(T_), function(t) numeric(n)); b_ok <- 0L
  for (b in seq_len(B_BOOT)) {
    bfit <- tryCatch(fitfun(redraw_sbm(models, n), algorithm = "leiden"), error = function(e) NULL)
    if (is.null(bfit)) next
    b_ok <- b_ok + 1L; bmem <- get_memberships(bfit, T_)
    s_nmi <- c(s_nmi, mean(vapply(seq_len(T_), function(t) nmi(bmem[[t]], mem0[[t]]), numeric(1))))
    s_ari <- c(s_ari, mean(vapply(seq_len(T_), function(t) ari(bmem[[t]], mem0[[t]]), numeric(1))))
    for (t in seq_len(T_)) { kcount[b, t] <- length(unique(bmem[[t]]))
      node_stab_acc[[t]] <- node_stab_acc[[t]] + node_jaccard(bmem[[t]], mem0[[t]]) }
  }
  if (b_ok < 10L && !MINI) stop("too few completed bootstrap replicates")
  K_true <- vapply(truth, function(m) length(unique(m)), integer(1))
  K_hat  <- vapply(mem0, function(m) length(unique(m)), integer(1))
  K_cov <- mean(vapply(seq_len(T_), function(t) { s <- kcount[, t]; s <- s[!is.na(s)]
    q <- quantile(s, c(ALPHA / 2, 1 - ALPHA / 2), names = FALSE); q[1] <= K_true[t] && K_true[t] <= q[2] }, logical(1)))
  row <- data.frame(task = TASK, sim = sim_id, n = n, K = cfg$K, p_switch = cfg$p_switch, density = cfg$density,
                    T_layers = T_, B = b_ok, acc_nmi = acc_nmi, acc_ari = acc_ari,
                    stab_nmi = mean(s_nmi), stab_ari = mean(s_ari), stab_nmi_sd = sd(s_nmi),
                    K_true_mean = mean(K_true), K_hat_mean = mean(K_hat), K_cov = K_cov, stringsAsFactors = FALSE)
  nodes <- NULL
  if (sim_id <= NODE_SIMS) {
    pick <- sample.int(n, min(NODE_PER_SIM, n))
    nodes <- do.call(rbind, lapply(seq_len(T_), function(t) {
      csize <- as.integer(table(mem0[[t]]))[match(mem0[[t]], sort(unique(mem0[[t]])))]
      data.frame(task = TASK, sim = sim_id, layer = t, node = pick, n = n, K = cfg$K, p_switch = cfg$p_switch,
                 density = cfg$density, T_layers = T_,
                 stab = node_stab_acc[[t]][pick] / b_ok, acc = node_acc[[t]][pick],
                 degree = rowSums(layers[[t]] > 0)[pick], csize = csize[pick], stringsAsFactors = FALSE) }))
  }
  list(row = row, nodes = nodes)
}

t0 <- proc.time()["elapsed"]
res <- mclapply(seq_len(M_SIMS), function(i) tryCatch(run_one(i), error = function(e) {
  message(sprintf("sim %d failed: %s", i, conditionMessage(e))); NULL }), mc.cores = CORES, mc.preschedule = FALSE)
res <- Filter(Negate(is.null), res)
stopifnot(length(res) >= (if (MINI) 1L else 10L))
rows <- do.call(rbind, lapply(res, `[[`, "row")); nodes <- do.call(rbind, lapply(res, `[[`, "nodes"))
write.csv(rows, outfiles[1], row.names = FALSE); write.csv(nodes, outfiles[2], row.names = FALSE)
cat(sprintf("[stab] task %d done: %d sims, %.1f min; acc_nmi=%.3f stab_nmi=%.3f cor=%.3f K_cov=%.3f\n",
            TASK, nrow(rows), (proc.time()["elapsed"] - t0) / 60, mean(rows$acc_nmi), mean(rows$stab_nmi),
            suppressWarnings(cor(rows$acc_nmi, rows$stab_nmi)), mean(rows$K_cov)))
