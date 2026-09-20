# =============================================================================
# replication/sim/03_mechanism_tests.R  --  where should DynMux beat multislice?
#
# Pre-specified (2026-09-20, before running) stress tests of the ONE assumption
# that separates the two couplings: multislice ties node i at t to node i at
# t+1 (identity, weight omega); DynMux ties COMMUNITIES by membership overlap
# and never looks at node identity. Three regimes violate the identity
# assumption in different ways. Each is a HYPOTHESIS about who wins; the sim
# decides. Every method runs on the same simulated series.
#
#   turnover       birth/death with heavy node turnover (frac 0.5 / 0.7 of nodes
#                  toggle presence every layer; base regime uses 0.1 / 0.3).
#                  Absent nodes have no identity tie to attach to.
#   coreperiphery  K communities that persist as SETS: a stable core never
#                  moves, the periphery re-draws its community every layer with
#                  prob p_rot. Node-level stickiness is wrong for the periphery,
#                  set-level overlap (the core) is right.
#   longT          long series (T = 30 / 60), many small communities (K = 6 /
#                  12), slow drift (p_switch = 0.02). Identity ties accumulate
#                  over T; with small per-layer degree, omega = 1 is a strong
#                  smoothing prior relative to the within-layer signal.
#
# Methods: DynMux Jaccard (Leiden), multislice identity ties (adjacent links,
# omega = 1, refined generalized Louvain), cross-sectional + Hungarian.
# Metrics: as sim/01 (eval_method). Per-rep seed = 15000 + TASK*1000 + rep.
#
# Usage (smoke): CMP_MINI=1 CMP_CFG=1 DM_ROOT=. Rscript replication/sim/03_mechanism_tests.R
# Array: 1..50 (slurm/03_mechanism.sbatch). Output: output/mechanism/mech_cfg%02d.csv
# =============================================================================
ROOT <- Sys.getenv("DM_ROOT", unset = getwd())
suppressMessages(pkgload::load_all(file.path(ROOT, "r_code"), quiet = TRUE))
suppressPackageStartupMessages(library(igraph))
TASK <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", Sys.getenv("CMP_CFG", "1")))
MINI <- identical(Sys.getenv("CMP_MINI", "0"), "1")
REPS <- if (MINI) 2L else 10L
QMIN <- 2L; QMAX <- 8L; NSTART <- 3L; RHO <- 0.10
source(file.path(ROOT, "replication", "sim", "lib_regimes.R"))   # dens_params, build_layer, dm(), method_hungarian, eval_method

outdir <- file.path(ROOT, "output", "mechanism"); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- generators (hard-coded copies, not shared with sim/01) ------------------
# turnover: sim_birthdeath with the presence-toggle fraction as a parameter.
sim_turnover <- function(n, K0, r, frac, seed) {
  set.seed(seed)
  dp <- dens_params(r, K0); T_ <- 12L; recruit <- 0.15
  mem <- sample(seq_len(K0), n, replace = TRUE)
  alive <- seq_len(K0); nxt <- K0 + 1L; state <- rep(TRUE, n)
  truth <- vector("list", T_); layers <- vector("list", T_); active <- vector("list", T_)
  for (t in seq_len(T_)) {
    if (t > 1L) {
      toggle <- sample(seq_len(n), floor(frac * n)); state[toggle] <- !state[toggle]
      if (t %% 3L == 0L) {
        dead <- alive[1]; alive <- alive[-1]
        mv <- which(mem == dead)
        if (length(mv) && length(alive)) mem[mv] <- sample(alive, length(mv), replace = TRUE)
        born <- nxt; nxt <- nxt + 1L; alive <- c(alive, born)
        rec <- sample(seq_len(n), floor(recruit * n)); mem[rec] <- born
      }
    }
    act <- state; active[[t]] <- act
    A <- matrix(0, n, n); ai <- which(act)
    if (length(ai) >= 2L) {
      ma <- mem[ai]; P <- outer(ma, ma, "=="); Pr <- ifelse(P, dp$p_in, dp$p_out); diag(Pr) <- 0
      sub <- matrix(0, length(ai), length(ai)); up <- upper.tri(Pr)
      sub[up] <- rbinom(sum(up), 1, Pr[up]); sub <- sub + t(sub); A[ai, ai] <- sub
    }
    layers[[t]] <- A
    tr <- mem; inact <- which(!act); if (length(inact)) tr[inact] <- max(mem) + seq_along(inact)
    truth[[t]] <- tr
  }
  stopifnot(all(vapply(active, sum, numeric(1)) >= 2))
  list(layers = layers, truth = truth, links = data.frame(from = seq_len(T_ - 1L), to = 2:T_, weight = 1), active = active)
}

# coreperiphery: fixed K = 4 sets; core share `core` of each community never
# moves; every periphery node re-draws a uniform community with prob p_rot.
sim_coreperiphery <- function(n, K, r, core, p_rot, seed) {
  set.seed(seed)
  dp <- dens_params(r, K); T_ <- 12L
  mem <- sample(seq_len(K), n, replace = TRUE)
  is_core <- runif(n) < core
  truth <- vector("list", T_); layers <- vector("list", T_)
  for (t in seq_len(T_)) {
    if (t > 1L) { rot <- !is_core & runif(n) < p_rot
      if (any(rot)) mem[rot] <- sample(seq_len(K), sum(rot), replace = TRUE) }
    layers[[t]] <- build_layer(mem, dp$p_in, dp$p_out); truth[[t]] <- mem
  }
  list(layers = layers, truth = truth, links = data.frame(from = seq_len(T_ - 1L), to = 2:T_, weight = 1))
}

# longT: churnswitch with T, K and p_switch as parameters.
sim_longT <- function(n, K, r, T_, p_switch, seed) {
  set.seed(seed)
  dp <- dens_params(r, K)
  mem <- sample(seq_len(K), n, replace = TRUE)
  truth <- vector("list", T_); layers <- vector("list", T_)
  for (t in seq_len(T_)) {
    if (t > 1L) { sw <- runif(n) < p_switch
      if (any(sw)) mem[sw] <- sample(seq_len(K), sum(sw), replace = TRUE) }
    layers[[t]] <- build_layer(mem, dp$p_in, dp$p_out); truth[[t]] <- mem
  }
  list(layers = layers, truth = truth, links = data.frame(from = seq_len(T_ - 1L), to = 2:T_, weight = 1))
}

# ---- grid (50 configs) --------------------------------------------------------
mech_grid <- function() {
  g1 <- expand.grid(regime = "turnover", n = c(50L, 100L, 200L), r = c(1.5, 3, 6), frac = c(0.5, 0.7),
                    core = NA_real_, p_rot = NA_real_, K = 4L, T_ = 12L, p_switch = NA_real_, stringsAsFactors = FALSE)
  g2 <- expand.grid(regime = "coreperiphery", n = c(100L, 200L), r = c(3, 6), frac = NA_real_,
                    core = c(0.5, 0.7), p_rot = c(0.3, 0.5), K = 4L, T_ = 12L, p_switch = NA_real_, stringsAsFactors = FALSE)
  g3 <- expand.grid(regime = "longT", n = c(100L, 200L), r = c(3, 6), frac = NA_real_,
                    core = NA_real_, p_rot = NA_real_, K = c(6L, 12L), T_ = c(30L, 60L), p_switch = 0.02, stringsAsFactors = FALSE)
  g <- rbind(g1, g2, g3); stopifnot(nrow(g) == 50L)
  set.seed(123); g[sample(nrow(g)), ]
}
simulate_mech <- function(cfg, seed) switch(cfg$regime,
  turnover      = sim_turnover(cfg$n, K0 = cfg$K, r = cfg$r, frac = cfg$frac, seed = seed),
  coreperiphery = sim_coreperiphery(cfg$n, K = cfg$K, r = cfg$r, core = cfg$core, p_rot = cfg$p_rot, seed = seed),
  longT         = sim_longT(cfg$n, K = cfg$K, r = cfg$r, T_ = cfg$T_, p_switch = cfg$p_switch, seed = seed),
  stop("unknown regime: ", cfg$regime))

MECH_METHODS <- list(
  "DynMux Jaccard"              = dm(fit_multilayer_jaccard,       "leiden"),
  "Multislice adjacent"         = dm(fit_multilayer_identity_ties, "leiden", links = FALSE),
  "Cross-sectional + Hungarian" = method_hungarian)

cfgs <- mech_grid(); stopifnot(TASK >= 1L, TASK <= nrow(cfgs)); cfg <- cfgs[TASK, ]
outfile <- file.path(outdir, sprintf("mech_cfg%02d.csv", TASK))
if (nzchar(Sys.getenv("SLURM_ARRAY_TASK_ID")) && file.exists(outfile)) { cat("[mech skip]", outfile, "\n"); quit(save = "no") }
cat(sprintf("[mech] task=%d/%d regime=%s n=%d r=%.1f frac=%s core=%s p_rot=%s K=%d T=%d p_switch=%s reps=%d\n",
            TASK, nrow(cfgs), cfg$regime, cfg$n, cfg$r, cfg$frac, cfg$core, cfg$p_rot, cfg$K, cfg$T_, cfg$p_switch, REPS))

metric_cols <- c("nmi_layer", "nmi_joint", "nmi_change", "k_mae", "comembership_acc", "mean_n_comm", "total_n_comm")
rows <- list(); t0 <- proc.time()[["elapsed"]]
for (rep in seq_len(REPS)) {
  seed <- 15000L + TASK * 1000L + rep
  sim <- simulate_mech(cfg, seed)
  for (mname in names(MECH_METHODS)) {
    t1  <- proc.time()[["elapsed"]]
    det <- tryCatch(MECH_METHODS[[mname]](sim), error = function(e) { message(sprintf("  [rep %d] %s errored: %s", rep, mname, conditionMessage(e))); NULL })
    el  <- proc.time()[["elapsed"]] - t1
    m   <- if (is.null(det)) setNames(rep(NA_real_, length(metric_cols)), metric_cols) else eval_method(det, sim)
    rows[[length(rows) + 1L]] <- data.frame(regime = cfg$regime, n = cfg$n, r = cfg$r, frac = cfg$frac, core = cfg$core, p_rot = cfg$p_rot,
      K = cfg$K, T = cfg$T_, p_switch = cfg$p_switch, rep = rep, method = mname,
      nmi_layer = round(m[["nmi_layer"]], 4), nmi_joint = round(m[["nmi_joint"]], 4), nmi_change = round(m[["nmi_change"]], 4),
      k_mae = round(m[["k_mae"]], 4), comembership_acc = round(m[["comembership_acc"]], 4), runtime_s = round(el, 3),
      mean_n_comm = round(m[["mean_n_comm"]], 4), total_n_comm = m[["total_n_comm"]], stringsAsFactors = FALSE)
  }
  if (rep %% 5L == 0L || MINI) cat(sprintf("  ... rep %d/%d (%.1f min)\n", rep, REPS, (proc.time()[["elapsed"]] - t0) / 60))
}
res <- do.call(rbind, rows)
stopifnot(nrow(res) == REPS * length(MECH_METHODS), all(metric_cols %in% names(res)))
write.csv(res, outfile, row.names = FALSE)
agg <- aggregate(cbind(nmi_layer, nmi_joint, k_mae, runtime_s) ~ method, data = res, FUN = function(x) mean(x, na.rm = TRUE))
cat(sprintf("[mech] task %d done in %.1f min -> %s\n", TASK, (proc.time()[["elapsed"]] - t0) / 60, outfile))
print(agg[order(-agg$nmi_joint), ], row.names = FALSE)
