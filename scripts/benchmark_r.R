#!/usr/bin/env Rscript
#
# Benchmark dynamic_multiplex against competing temporal community detection
# approaches in R.
#
# Competitors:
#   1. dynamic_multiplex (adjacent-layer coupling) -- our method
#   2. Independent Louvain (no temporal coupling)
#   3. Full-coupling Louvain (all layer pairs connected)
#   4. Aggregated Louvain (collapse layers into one network)
#   5. Multislice Louvain -- all-to-all (Mucha et al. 2010)
#   6. Multislice Louvain -- adjacent-only (supra-adjacency, adjacent coupling)
#   7. Spectral SBM (knows true k, informational advantage)
#
# Extended simulation dimensions:
#   - p_switch: community switching rate (temporal persistence)
#   - n_layers: tests long-range pooling degradation
#
# Usage:
#   Rscript scripts/benchmark_r.R

library(dynamic_multiplex)

# ---------------------------------------------------------------------------
# Simulation with evolving communities
# ---------------------------------------------------------------------------

simulate_evolving_networks <- function(n_nodes, n_layers, n_communities,
                                       p_in = 0.3, p_out = 0.05,
                                       p_switch = 0.05, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  memberships <- list()
  memberships[[1]] <- sample(seq_len(n_communities), size = n_nodes, replace = TRUE)

  for (t in 2:n_layers) {
    prev <- memberships[[t - 1]]
    switch_mask <- runif(n_nodes) < p_switch
    prev[switch_mask] <- sample(seq_len(n_communities), sum(switch_mask), replace = TRUE)
    memberships[[t]] <- prev
  }

  layers <- lapply(seq_len(n_layers), function(t) {
    mat <- matrix(0, nrow = n_nodes, ncol = n_nodes)
    for (i in seq_len(n_nodes - 1)) {
      for (j in (i + 1):n_nodes) {
        prob <- if (memberships[[t]][i] == memberships[[t]][j]) p_in else p_out
        tie <- rbinom(1, 1, prob)
        mat[i, j] <- tie
        mat[j, i] <- tie
      }
    }
    mat
  })

  list(layers = layers, true_memberships = memberships)
}


# ---------------------------------------------------------------------------
# NMI computation (igraph::compare)
# ---------------------------------------------------------------------------

compute_nmi <- function(detected, true_mem) {
  mean(vapply(seq_along(detected), function(i) {
    igraph::compare(detected[[i]], true_mem[[i]], method = "nmi")
  }, numeric(1)))
}

compute_ari <- function(detected, true_mem) {
  mean(vapply(seq_along(detected), function(i) {
    igraph::compare(detected[[i]], true_mem[[i]], method = "adjusted.rand")
  }, numeric(1)))
}


# ---------------------------------------------------------------------------
# Competing methods
# ---------------------------------------------------------------------------

extract_memberships <- function(fit) {
  lapply(fit$layer_communities, function(x) x$membership)
}

run_dynmux <- function(layers, fit_type = "jaccard") {
  fit_fn <- switch(
    fit_type,
    jaccard          = fit_multilayer_jaccard,
    overlap          = fit_multilayer_overlap,
    weighted_jaccard = fit_multilayer_weighted_jaccard,
    weighted_overlap = fit_multilayer_weighted_overlap,
    identity         = fit_multilayer_identity_ties
  )
  fit <- fit_fn(layers, algorithm = "louvain")
  extract_memberships(fit)
}

run_independent <- function(layers) {
  lapply(layers, function(mat) {
    g <- igraph::graph_from_adjacency_matrix(mat, mode = "undirected",
                                              weighted = TRUE, diag = FALSE)
    cl <- igraph::cluster_louvain(g, weights = igraph::E(g)$weight)
    igraph::membership(cl)
  })
}

run_full_coupling <- function(layers) {
  n <- length(layers)
  links <- do.call(rbind, lapply(seq_len(n - 1), function(i) {
    do.call(rbind, lapply((i + 1):n, function(j) {
      data.frame(from = i, to = j, weight = 1, stringsAsFactors = FALSE)
    }))
  }))
  fit <- fit_multilayer_jaccard(layers, algorithm = "louvain", layer_links = links)
  extract_memberships(fit)
}

run_aggregated <- function(layers) {
  agg <- Reduce("+", layers)
  g <- igraph::graph_from_adjacency_matrix(agg, mode = "undirected",
                                            weighted = TRUE, diag = FALSE)
  cl <- igraph::cluster_louvain(g, weights = igraph::E(g)$weight)
  mem <- igraph::membership(cl)
  replicate(length(layers), mem, simplify = FALSE)
}

run_multislice_all_to_all <- function(layers, omega = 1.0) {
  # Multislice modularity (Mucha et al. 2010).
  # Builds supra-adjacency with identity coupling between ALL layer pairs.
  # The temporal problem: community structure at t=1 influences t=T.
  n_nodes <- nrow(layers[[1]])
  n_layers <- length(layers)
  N <- n_nodes * n_layers

  supra <- matrix(0, N, N)

  # Diagonal blocks: within-layer adjacency
  for (t in seq_len(n_layers)) {
    r0 <- (t - 1) * n_nodes
    supra[(r0 + 1):(r0 + n_nodes), (r0 + 1):(r0 + n_nodes)] <- layers[[t]]
  }

  # Off-diagonal: identity coupling between ALL layer pairs
  for (t1 in seq_len(n_layers - 1)) {
    for (t2 in (t1 + 1):n_layers) {
      for (i in seq_len(n_nodes)) {
        r1 <- (t1 - 1) * n_nodes + i
        r2 <- (t2 - 1) * n_nodes + i
        supra[r1, r2] <- omega
        supra[r2, r1] <- omega
      }
    }
  }

  g <- igraph::graph_from_adjacency_matrix(supra, mode = "undirected",
                                            weighted = TRUE, diag = FALSE)
  cl <- igraph::cluster_louvain(g, weights = igraph::E(g)$weight)
  supra_mem <- igraph::membership(cl)

  # Extract per-layer memberships
  lapply(seq_len(n_layers), function(t) {
    r0 <- (t - 1) * n_nodes
    supra_mem[(r0 + 1):(r0 + n_nodes)]
  })
}

run_multislice_adjacent <- function(layers, omega = 1.0) {
  # Multislice with ADJACENT-ONLY coupling.
  # Same supra-adjacency but interlayer ties only connect t to t+1.
  n_nodes <- nrow(layers[[1]])
  n_layers <- length(layers)
  N <- n_nodes * n_layers

  supra <- matrix(0, N, N)

  for (t in seq_len(n_layers)) {
    r0 <- (t - 1) * n_nodes
    supra[(r0 + 1):(r0 + n_nodes), (r0 + 1):(r0 + n_nodes)] <- layers[[t]]
  }

  # Adjacent-only coupling
  for (t in seq_len(n_layers - 1)) {
    for (i in seq_len(n_nodes)) {
      r1 <- (t - 1) * n_nodes + i
      r2 <- t * n_nodes + i
      supra[r1, r2] <- omega
      supra[r2, r1] <- omega
    }
  }

  g <- igraph::graph_from_adjacency_matrix(supra, mode = "undirected",
                                            weighted = TRUE, diag = FALSE)
  cl <- igraph::cluster_louvain(g, weights = igraph::E(g)$weight)
  supra_mem <- igraph::membership(cl)

  lapply(seq_len(n_layers), function(t) {
    r0 <- (t - 1) * n_nodes
    supra_mem[(r0 + 1):(r0 + n_nodes)]
  })
}

run_spectral_sbm <- function(layers, n_communities) {
  # Spectral SBM: adjacency spectral embedding + k-means per layer.
  # Uses the leading eigenvectors of the adjacency matrix (Rohe et al. 2011).
  # Knows the true number of communities (informational advantage).
  lapply(layers, function(mat) {
    k <- min(n_communities, nrow(mat))
    eig <- eigen(mat, symmetric = TRUE)
    ord <- order(abs(eig$values), decreasing = TRUE)[seq_len(k)]
    embedding <- eig$vectors[, ord, drop = FALSE]
    km <- kmeans(embedding, centers = k, nstart = 10)
    km$cluster
  })
}


# ---------------------------------------------------------------------------
# Main benchmark
# ---------------------------------------------------------------------------

# Extended configurations: more n_layers and p_switch values
# to test the temporal pooling problem
configs <- expand.grid(
  n_nodes       = c(50, 100),
  n_layers      = c(5, 10, 20),
  n_communities = c(3, 5),
  p_switch      = c(0.0, 0.02, 0.05, 0.10, 0.20),
  stringsAsFactors = FALSE
)

base_methods <- list(
  DynMux_Jaccard         = function(L) run_dynmux(L, "jaccard"),
  DynMux_Overlap         = function(L) run_dynmux(L, "overlap"),
  DynMux_WtJaccard       = function(L) run_dynmux(L, "weighted_jaccard"),
  DynMux_WtOverlap       = function(L) run_dynmux(L, "weighted_overlap"),
  DynMux_Identity        = function(L) run_dynmux(L, "identity"),
  Independent            = run_independent,
  FullCoupling           = run_full_coupling,
  Aggregated             = run_aggregated,
  Multislice_AllToAll     = run_multislice_all_to_all,
  Multislice_Adjacent    = run_multislice_adjacent
)

n_reps <- 5
results <- list()
idx <- 0

cat(sprintf("Running %d configs x %d reps x %d methods ...\n",
            nrow(configs), n_reps, length(base_methods) + 1L))

for (cfg_i in seq_len(nrow(configs))) {
  cfg <- configs[cfg_i, ]
  # SBM needs the true k; capture via local scope
  k_true <- cfg$n_communities
  methods <- c(base_methods, list(
    SBM = function(L) run_spectral_sbm(L, k_true)
  ))

  for (rep in seq_len(n_reps)) {
    seed <- cfg_i * 1000 + rep
    sim <- simulate_evolving_networks(
      n_nodes = cfg$n_nodes, n_layers = cfg$n_layers,
      n_communities = cfg$n_communities,
      p_in = 0.3, p_out = 0.05, p_switch = cfg$p_switch, seed = seed
    )

    for (method_name in names(methods)) {
      t0 <- proc.time()["elapsed"]
      detected <- tryCatch(
        methods[[method_name]](sim$layers),
        error = function(e) { message("  SKIP ", method_name, ": ", e$message); NULL }
      )
      elapsed <- proc.time()["elapsed"] - t0

      if (is.null(detected)) next

      nmi_val <- compute_nmi(detected, sim$true_memberships)
      ari_val <- compute_ari(detected, sim$true_memberships)

      idx <- idx + 1
      results[[idx]] <- data.frame(
        n_nodes      = cfg$n_nodes,
        n_layers     = cfg$n_layers,
        n_communities = cfg$n_communities,
        p_switch     = cfg$p_switch,
        rep          = rep,
        method       = method_name,
        nmi          = round(nmi_val, 4),
        ari          = round(ari_val, 4),
        runtime_s    = round(elapsed, 4),
        stringsAsFactors = FALSE
      )
    }
  }

  if (cfg_i %% 4 == 0 || cfg_i == nrow(configs)) {
    cat(sprintf("  [%d/%d configs done]\n", cfg_i, nrow(configs)))
  }
}

df <- do.call(rbind, results)
outpath <- "scripts/benchmark_r_results.csv"
write.csv(df, outpath, row.names = FALSE)
cat(sprintf("\nResults saved to %s\n", outpath))

# Summary
summary_df <- aggregate(
  cbind(nmi, ari, runtime_s) ~ method + p_switch,
  data = df, FUN = mean
)
summary_df <- summary_df[order(summary_df$p_switch, -summary_df$nmi), ]

cat("\n=== Summary (mean across configs and reps) ===\n")
print(summary_df, row.names = FALSE)


# ---------------------------------------------------------------------------
# Bootstrap CI demo
# ---------------------------------------------------------------------------

cat("\n=== Bootstrap CI Demo ===\n")
sim_boot <- simulate_evolving_networks(
  n_nodes = 50, n_layers = 5, n_communities = 3,
  p_in = 0.3, p_out = 0.05, p_switch = 0.05, seed = 42
)

cat("Running bootstrap (n_boot=50) ...\n")
boot <- bootstrap_multilayer(
  sim_boot$layers, fit_type = "jaccard", algorithm = "louvain",
  n_boot = 50, seed = 42
)
ci <- community_ci(boot, alpha = 0.05)

cat("\nModularity 95% CIs:\n")
print(ci$modularity_ci, row.names = FALSE)

cat("\nCommunity count 95% CIs:\n")
print(ci$community_count_ci, row.names = FALSE)

cat("\nMean node stability per layer:\n")
print(ci$mean_node_stability, row.names = FALSE)

write.csv(ci$modularity_ci, "scripts/bootstrap_modularity_ci_r.csv", row.names = FALSE)
write.csv(ci$community_count_ci, "scripts/bootstrap_count_ci_r.csv", row.names = FALSE)
write.csv(ci$mean_node_stability, "scripts/bootstrap_stability_r.csv", row.names = FALSE)
cat("\nBootstrap results saved to scripts/bootstrap_*_r.csv\n")
