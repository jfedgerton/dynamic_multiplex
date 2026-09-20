# =============================================================================
# replication/sim/05_alt_bootstrap.R   (EXPLORATORY ARM, not used by the paper yet)
#
# Tests alternatives to the co-assignment interval on the same 594-cell grid as
# sim/02 (jaccard/leiden spec only), M = 50 sims per cell:
#
#   Bootstrap variants, each producing its own phat and (phat bin, p* bin) joint
#   counts + Wilson coverage by bin, on the SAME simulated networks:
#     sbm    the existing scheme: per-layer block densities from the detected
#            partition, fresh Bernoulli redraw (sim/02's "edge-resampling")
#     dcsbm  degree-corrected redraw: P_ij = theta_i theta_j p_block, theta
#            from observed degrees within the detected block
#     rewire nonparametric perturbation: 10% of observed edges moved to random
#            non-edges, no model fitted
#
#   Per-pair FEATURE sample (sbm variant only, first FEAT_SIMS sims), for
#   feature-conditioned calibration in post: phat, p*, layer, degrees,
#   detected-community sizes, same_detected, node co-assignment entropy,
#   the sim's mean Wilson width (polarisation), and the design columns.
#   Mid-range pairs (0.1 < phat < 0.9) are oversampled because they are the
#   ones the pooled interval cannot resolve.
#
# Outputs per task under $DM_ROOT/output/alt_bootstrap/:
#   alt_joint_task%05d.csv   variant, bin_phat, bin_pstar, count       (3 x 2500 rows)
#   alt_calib_task%05d.csv   variant, bin, n_pairs, n_wilson_covered   (3 x 50 rows)
#   alt_cov_task%05d.csv     one row per sim x variant: cov_P_mean, width_P_mean
#   alt_pairs_task%05d.csv   the feature sample
#
# Usage (mini): COV_TASK=1 COV_MINI=1 DM_ROOT=. Rscript replication/sim/05_alt_bootstrap.R
# Array: 1..594 (one cell per task); seeds (19000 + TASK) * 100000 + sim.
# =============================================================================
suppressPackageStartupMessages({ library(dynamicmultiplex); library(parallel) })

TASK  <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", Sys.getenv("COV_TASK", "1")))
CORES <- as.integer(Sys.getenv("COV_CORES", "8"))
MINI  <- identical(Sys.getenv("COV_MINI", "0"), "1")

M_SIMS    <- if (MINI) 3L else 50L
FEAT_SIMS <- if (MINI) 2L else 10L          # sims that also write the pair sample
B_BOOT    <- if (MINI) 8L else 100L
R_TRUTH   <- if (MINI) 5L else 100L
ALPHA     <- 0.05
N_BINS    <- 50L
Z_ALPHA   <- qnorm(1 - ALPHA / 2)
REWIRE_FRAC <- 0.10
PAIRS_MID <- if (MINI) 20L else 100L        # per layer: mid-range pairs sampled
PAIRS_EXT <- if (MINI) 10L else 30L         # per layer: extreme pairs sampled
VARIANTS  <- c("sbm", "dcsbm", "rewire")

DENSITIES <- list(weak    = c(p_in = 0.20, p_out = 0.10),
                  default = c(p_in = 0.30, p_out = 0.05),
                  strong  = c(p_in = 0.50, p_out = 0.02))
ROOT   <- Sys.getenv("DM_ROOT", unset = getwd())
outdir <- file.path(ROOT, "output", "alt_bootstrap")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cells <- expand.grid(n = c(50, 100, 200, 400), K = c(3, 5, 10),
                     p_switch = c(0.02, 0.05, 0.10, 0.20, 0.50, 0.80),
                     density = names(DENSITIES), T_ = c(5, 10, 15), stringsAsFactors = FALSE)
cells <- cells[!(cells$K == 10 & cells$n == 50), ]
cells <- cells[order(cells$n, cells$T_, cells$K, cells$p_switch, cells$density), ]
rownames(cells) <- NULL
stopifnot(nrow(cells) == 594)
set.seed(123); perm <- sample.int(nrow(cells))
if (TASK < 1 || TASK > nrow(cells)) stop("task out of range")
outfiles <- file.path(outdir, sprintf(c("alt_joint_task%05d.csv", "alt_calib_task%05d.csv",
                                        "alt_cov_task%05d.csv", "alt_pairs_task%05d.csv"), TASK))
if (nzchar(Sys.getenv("SLURM_ARRAY_TASK_ID")) && all(file.exists(outfiles))) {
  cat("[skip] task", TASK, "already complete\n"); quit(save = "no")
}
cfg <- as.list(cells[perm[TASK], ])
cfg$p_in <- DENSITIES[[cfg$density]][["p_in"]]; cfg$p_out <- DENSITIES[[cfg$density]][["p_out"]]
cfg$fit_type <- "jaccard"; cfg$algorithm <- "leiden"
cat(sprintf("[alt] task=%d/594 cell=%d : n=%d K=%d sw=%.2f dens=%s T=%d M=%d B=%d R=%d\n",
            TASK, perm[TASK], cfg$n, cfg$K, cfg$p_switch, cfg$density, cfg$T_, M_SIMS, B_BOOT, R_TRUTH))
fitfun <- fit_multilayer_jaccard

# --------------------------- simulator (as sim/02) ---------------------------
gen_memberships <- function(cfg) {
  m <- vector("list", cfg$T_); m[[1]] <- sample(seq_len(cfg$K), cfg$n, replace = TRUE)
  for (t in 2:cfg$T_) { prev <- m[[t - 1]]; mask <- runif(cfg$n) < cfg$p_switch
    prev[mask] <- sample(seq_len(cfg$K), sum(mask), replace = TRUE); m[[t]] <- prev }
  m
}
gen_layers <- function(memberships, cfg) {
  n <- cfg$n
  lapply(memberships, function(mem) {
    P <- outer(mem, mem, "=="); Pr <- ifelse(P, cfg$p_in, cfg$p_out); diag(Pr) <- 0
    A <- matrix(0, n, n); up <- upper.tri(Pr); A[up] <- rbinom(sum(up), 1, Pr[up]); A + t(A) })
}
get_memberships <- function(fit, T_) lapply(seq_len(T_), function(t) as.integer(fit$layer_communities[[t]]$membership))
wilson_bounds <- function(phat, b) {
  z2 <- Z_ALPHA^2; denom <- 1 + z2 / b
  center <- (phat + z2 / (2 * b)) / denom
  half <- Z_ALPHA * sqrt(phat * (1 - phat) / b + z2 / (4 * b^2)) / denom
  list(lo = pmax(center - half, 0), hi = pmin(center + half, 1))
}

# --------------------------- bootstrap variants ------------------------------
# sbm: block densities from detected partition, Bernoulli redraw
fit_sbm <- function(layers, mem_hat) lapply(seq_along(layers), function(t) {
  A <- layers[[t]]; mem <- mem_hat[[t]]; same <- outer(mem, mem, "=="); sel <- upper.tri(A); E <- A > 0
  p_all <- mean(E[sel]); in_d <- sel & same; out_d <- sel & !same
  list(same = same, p_in = if (any(in_d)) mean(E[in_d]) else p_all,
       p_out = if (any(out_d)) mean(E[out_d]) else p_all)
})
redraw_sbm <- function(models, n) lapply(models, function(em) {
  probs <- ifelse(em$same, em$p_in, em$p_out); up <- upper.tri(probs)
  M <- matrix(0, n, n); M[up] <- rbinom(sum(up), 1, probs[up]); M + t(M) })

# dcsbm: theta_i = deg_i / mean degree of i's detected block; P_ij = min(1, theta_i theta_j p_block)
fit_dcsbm <- function(layers, mem_hat) lapply(seq_along(layers), function(t) {
  A <- layers[[t]]; mem <- mem_hat[[t]]; deg <- rowSums(A > 0)
  blockmean <- ave(deg, mem, FUN = mean); theta <- ifelse(blockmean > 0, deg / blockmean, 1)
  same <- outer(mem, mem, "=="); sel <- upper.tri(A); E <- A > 0
  p_all <- mean(E[sel]); in_d <- sel & same; out_d <- sel & !same
  p_in <- if (any(in_d)) mean(E[in_d]) else p_all; p_out <- if (any(out_d)) mean(E[out_d]) else p_all
  probs <- outer(theta, theta) * ifelse(same, p_in, p_out)
  probs[probs > 1] <- 1                       # (pmin(1, M) would drop the matrix dims)
  list(probs = probs)
})
redraw_dcsbm <- function(models, n) lapply(models, function(em) {
  up <- upper.tri(em$probs); M <- matrix(0, n, n); M[up] <- rbinom(sum(up), 1, em$probs[up]); M + t(M) })

# rewire: move REWIRE_FRAC of observed edges to random non-edges (no model)
redraw_rewire <- function(layers, n) lapply(layers, function(A) {
  up <- which(upper.tri(A)); on <- up[A[up] > 0]; off <- up[A[up] == 0]
  k <- round(REWIRE_FRAC * length(on)); if (k == 0 || length(off) == 0) return(A)
  drop <- sample(on, k); add <- sample(off, min(k, length(off)))
  M <- A; M[drop] <- 0; M[add] <- 1; M[lower.tri(M)] <- 0; M + t(M) })

boot_coassign <- function(draw_fun, n, T_) {
  co_acc <- lapply(seq_len(T_), function(t) matrix(0, n, n)); b_ok <- 0L
  for (b in seq_len(B_BOOT)) {
    bl <- draw_fun()
    bfit <- tryCatch(fitfun(bl, algorithm = cfg$algorithm), error = function(e) NULL)
    if (is.null(bfit)) next
    b_ok <- b_ok + 1L; bmem <- get_memberships(bfit, T_)
    for (t in seq_len(T_)) co_acc[[t]] <- co_acc[[t]] + outer(bmem[[t]], bmem[[t]], "==")
  }
  if (b_ok < 10L && !MINI) stop("too few completed bootstrap replicates")
  list(co_prob = lapply(co_acc, function(m) m / b_ok), b_ok = b_ok)
}

# --------------------------- one simulation ----------------------------------
run_one <- function(sim_id) {
  seed <- (19000L + TASK) * 100000L + sim_id   # < 2^31; disjoint from sim/02 (15000 + task)
  set.seed(seed)
  n <- cfg$n; T_ <- cfg$T_
  truth  <- gen_memberships(cfg)
  layers <- gen_layers(truth, cfg)
  fit0 <- fitfun(layers, algorithm = cfg$algorithm)
  mem_hat <- get_memberships(fit0, T_)

  # ground truth propensity p* (fresh draws from the TRUE memberships)
  pstar_acc <- lapply(seq_len(T_), function(t) matrix(0, n, n))
  for (r in seq_len(R_TRUTH)) {
    tfit <- fitfun(gen_layers(truth, cfg), algorithm = cfg$algorithm); tmem <- get_memberships(tfit, T_)
    for (t in seq_len(T_)) pstar_acc[[t]] <- pstar_acc[[t]] + outer(tmem[[t]], tmem[[t]], "==")
  }
  pstar_mat <- lapply(pstar_acc, function(m) m / R_TRUTH)

  # three bootstraps on the same network
  sbm_models <- fit_sbm(layers, mem_hat); dc_models <- fit_dcsbm(layers, mem_hat)
  boots <- list(
    sbm    = boot_coassign(function() redraw_sbm(sbm_models, n), n, T_),
    dcsbm  = boot_coassign(function() redraw_dcsbm(dc_models, n), n, T_),
    rewire = boot_coassign(function() redraw_rewire(layers, n), n, T_))

  joint <- list(); calib <- list(); covrows <- list()
  for (v in VARIANTS) {
    cp <- boots[[v]]$co_prob; b_ok <- boots[[v]]$b_ok
    jnt <- integer(N_BINS * N_BINS); btot <- integer(N_BINS); bwc <- integer(N_BINS)
    num <- 0; den <- 0; wsum <- 0
    for (t in seq_len(T_)) {
      up <- upper.tri(cp[[t]]); phat <- cp[[t]][up]; pstar <- pstar_mat[[t]][up]
      wb <- wilson_bounds(phat, b_ok); covp <- (wb$lo <= pstar) & (pstar <= wb$hi)
      num <- num + sum(covp); den <- den + length(covp); wsum <- wsum + sum(wb$hi - wb$lo)
      bin <- pmin(N_BINS, 1L + floor(phat * N_BINS)); binS <- pmin(N_BINS, 1L + floor(pstar * N_BINS))
      btot <- btot + tabulate(bin, N_BINS); bwc <- bwc + tabulate(bin[covp], N_BINS)
      jnt <- jnt + tabulate(bin + N_BINS * (binS - 1L), N_BINS * N_BINS)
    }
    joint[[v]] <- jnt; calib[[v]] <- list(btot = btot, bwc = bwc)
    covrows[[v]] <- data.frame(task = TASK, sim = sim_id, variant = v, n = n, K = cfg$K,
      p_switch = cfg$p_switch, density = cfg$density, T_layers = T_, B = b_ok,
      cov_P_mean = num / den, width_P_mean = wsum / den, stringsAsFactors = FALSE)
  }

  # per-pair feature sample (sbm variant), mid-range oversampled
  pairs <- NULL
  if (sim_id <= FEAT_SIMS) {
    cp <- boots$sbm$co_prob
    mean_width <- covrows$sbm$width_P_mean
    rows <- list()
    for (t in seq_len(T_)) {
      A <- layers[[t]]; deg <- rowSums(A > 0); mem <- mem_hat[[t]]
      csize <- as.integer(table(mem))[match(mem, sort(unique(mem)))]
      # node entropy of its co-assignment row (how undecided the node is)
      P <- cp[[t]]; Pc <- pmin(pmax(P, 1e-6), 1 - 1e-6)
      ent <- rowMeans(-(Pc * log(Pc) + (1 - Pc) * log(1 - Pc)))
      idx <- which(upper.tri(P), arr.ind = TRUE)
      ph <- P[idx]; mid <- which(ph > 0.1 & ph < 0.9); ext <- which(!(ph > 0.1 & ph < 0.9))
      pick <- c(if (length(mid)) mid[sample.int(length(mid), min(PAIRS_MID, length(mid)))],
                if (length(ext)) ext[sample.int(length(ext), min(PAIRS_EXT, length(ext)))])
      if (!length(pick)) next
      i <- idx[pick, 1]; j <- idx[pick, 2]
      rows[[t]] <- data.frame(task = TASK, sim = sim_id, layer = t, n = n, K = cfg$K,
        p_switch = cfg$p_switch, density = cfg$density, T_layers = T_,
        phat = ph[pick], pstar = pstar_mat[[t]][cbind(i, j)],
        same_true = as.integer(truth[[t]][i] == truth[[t]][j]),
        same_detected = as.integer(mem[i] == mem[j]),
        deg_i = deg[i], deg_j = deg[j], csize_i = csize[i], csize_j = csize[j],
        ent_i = ent[i], ent_j = ent[j], mean_width = mean_width,
        phat_dcsbm = boots$dcsbm$co_prob[[t]][cbind(i, j)],
        phat_rewire = boots$rewire$co_prob[[t]][cbind(i, j)],
        stringsAsFactors = FALSE)
    }
    pairs <- do.call(rbind, rows)
  }
  list(joint = joint, calib = calib, cov = do.call(rbind, covrows), pairs = pairs)
}

t0 <- proc.time()["elapsed"]
res <- mclapply(seq_len(M_SIMS), function(i) tryCatch(run_one(i), error = function(e) {
  message(sprintf("sim %d failed: %s", i, conditionMessage(e))); NULL }),
  mc.cores = CORES, mc.preschedule = FALSE)
res <- Filter(Negate(is.null), res)
stopifnot(length(res) >= (if (MINI) 1L else 10L))

design <- data.frame(task = TASK, n = cfg$n, K = cfg$K, p_switch = cfg$p_switch,
                     density = cfg$density, T_ = cfg$T_, stringsAsFactors = FALSE)
joint_df <- do.call(rbind, lapply(VARIANTS, function(v) {
  jnt <- Reduce(`+`, lapply(res, function(r) r$joint[[v]]))
  cbind(design, variant = v, bin_phat = rep(seq_len(N_BINS), times = N_BINS),
        bin_pstar = rep(seq_len(N_BINS), each = N_BINS), count = jnt) }))
calib_df <- do.call(rbind, lapply(VARIANTS, function(v) {
  cbind(design, variant = v, bin = seq_len(N_BINS),
        n_pairs = Reduce(`+`, lapply(res, function(r) r$calib[[v]]$btot)),
        n_wilson_covered = Reduce(`+`, lapply(res, function(r) r$calib[[v]]$bwc))) }))
cov_df   <- do.call(rbind, lapply(res, `[[`, "cov"))
pairs_df <- do.call(rbind, lapply(res, `[[`, "pairs"))
stopifnot(sum(joint_df$count[joint_df$variant == "sbm"]) == sum(calib_df$n_pairs[calib_df$variant == "sbm"]))

write.csv(joint_df, outfiles[1], row.names = FALSE)
write.csv(calib_df, outfiles[2], row.names = FALSE)
write.csv(cov_df,   outfiles[3], row.names = FALSE)
write.csv(pairs_df, outfiles[4], row.names = FALSE)
cat(sprintf("[alt] task %d done: %d/%d sims, %.1f min; Wilson cov sbm=%.3f dcsbm=%.3f rewire=%.3f; pair rows=%d\n",
            TASK, length(res), M_SIMS, (proc.time()["elapsed"] - t0) / 60,
            mean(cov_df$cov_P_mean[cov_df$variant == "sbm"]), mean(cov_df$cov_P_mean[cov_df$variant == "dcsbm"]),
            mean(cov_df$cov_P_mean[cov_df$variant == "rewire"]), nrow(pairs_df)))
