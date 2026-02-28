#!/usr/bin/env Rscript
#
# Manuscript Section: Synthetic Experiments
#
# Monte Carlo comparison of dynamic_multiplex against competing temporal
# community detection strategies on planted-partition networks with evolving
# community structure.
#
# Outputs (manuscript/output/):
#   - synthetic_results.csv        : raw simulation results
#   - fig_nmi_by_switch.pdf        : NMI vs community switching rate
#   - fig_ari_by_switch.pdf        : ARI vs community switching rate
#   - fig_nmi_by_layers.pdf        : NMI vs number of temporal layers
#   - fig_runtime_comparison.pdf   : runtime comparison across methods
#   - tab_main_results.csv         : summary table for LaTeX
#
# Usage:
#   Rscript manuscript/01_synthetic_experiments.R

library(dynamic_multiplex)
library(ggplot2)

outdir <- "manuscript/output"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)


# ===========================================================================
# 1. Simulation engine (evolving planted partition)
# ===========================================================================

simulate_evolving <- function(n_nodes, n_layers, n_communities,
                               p_in = 0.3, p_out = 0.05,
                               p_switch = 0.05, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  memberships <- list()
  memberships[[1]] <- sample(seq_len(n_communities), n_nodes, replace = TRUE)
  for (t in 2:n_layers) {
    prev <- memberships[[t - 1]]
    mask <- runif(n_nodes) < p_switch
    prev[mask] <- sample(seq_len(n_communities), sum(mask), replace = TRUE)
    memberships[[t]] <- prev
  }

  layers <- lapply(seq_len(n_layers), function(t) {
    mat <- matrix(0, n_nodes, n_nodes)
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


# ===========================================================================
# 2. Methods
# ===========================================================================

extract_mem <- function(fit) lapply(fit$layer_communities, function(x) x$membership)

method_dynmux_jaccard <- function(L) {
  extract_mem(fit_multilayer_jaccard(L, algorithm = "louvain"))
}
method_dynmux_overlap <- function(L) {
  extract_mem(fit_multilayer_overlap(L, algorithm = "louvain"))
}
method_dynmux_wt_jaccard <- function(L) {
  extract_mem(fit_multilayer_weighted_jaccard(L, algorithm = "louvain"))
}
method_dynmux_wt_overlap <- function(L) {
  extract_mem(fit_multilayer_weighted_overlap(L, algorithm = "louvain"))
}
method_dynmux_identity <- function(L) {
  extract_mem(fit_multilayer_identity_ties(L, algorithm = "louvain"))
}

method_independent <- function(L) {
  lapply(L, function(mat) {
    g <- igraph::graph_from_adjacency_matrix(mat, mode = "undirected",
                                              weighted = TRUE, diag = FALSE)
    igraph::membership(igraph::cluster_louvain(g, weights = igraph::E(g)$weight))
  })
}

method_full_coupling <- function(L) {
  n <- length(L)
  links <- do.call(rbind, lapply(seq_len(n - 1), function(i) {
    do.call(rbind, lapply((i + 1):n, function(j) {
      data.frame(from = i, to = j, weight = 1, stringsAsFactors = FALSE)
    }))
  }))
  extract_mem(fit_multilayer_jaccard(L, algorithm = "louvain", layer_links = links))
}

method_aggregated <- function(L) {
  agg <- Reduce("+", L)
  g <- igraph::graph_from_adjacency_matrix(agg, mode = "undirected",
                                            weighted = TRUE, diag = FALSE)
  mem <- igraph::membership(igraph::cluster_louvain(g, weights = igraph::E(g)$weight))
  replicate(length(L), mem, simplify = FALSE)
}

method_spectral_sbm <- function(L, n_communities) {
  # Spectral SBM: adjacency spectral embedding + k-means per layer.
  # Uses the leading eigenvectors of the adjacency matrix (Rohe et al. 2011).
  # Knows the true number of communities (informational advantage).
  lapply(L, function(mat) {
    k <- min(n_communities, nrow(mat))
    eig <- eigen(mat, symmetric = TRUE)
    ord <- order(abs(eig$values), decreasing = TRUE)[seq_len(k)]
    embedding <- eig$vectors[, ord, drop = FALSE]
    km <- kmeans(embedding, centers = k, nstart = 10)
    km$cluster
  })
}

base_methods <- list(
  "DynMux (Jaccard)"          = method_dynmux_jaccard,
  "DynMux (Overlap)"          = method_dynmux_overlap,
  "DynMux (Wt. Jaccard)"     = method_dynmux_wt_jaccard,
  "DynMux (Wt. Overlap)"     = method_dynmux_wt_overlap,
  "DynMux (Identity)"         = method_dynmux_identity,
  "Independent Louvain"       = method_independent,
  "Full Coupling"             = method_full_coupling,
  "Aggregated Louvain"        = method_aggregated
)


# ===========================================================================
# 3. Evaluation
# ===========================================================================

compute_nmi <- function(det, tru) {
  mean(vapply(seq_along(det), function(i) {
    igraph::compare(det[[i]], tru[[i]], method = "nmi")
  }, numeric(1)))
}

compute_ari <- function(det, tru) {
  mean(vapply(seq_along(det), function(i) {
    igraph::compare(det[[i]], tru[[i]], method = "adjusted.rand")
  }, numeric(1)))
}


# ===========================================================================
# 4. Monte Carlo sweep
# ===========================================================================

configs <- expand.grid(
  n_nodes       = c(50, 100, 200),
  n_layers      = c(3, 5, 10, 20),
  n_communities = c(3, 5),
  p_switch      = c(0.0, 0.02, 0.05, 0.10, 0.20),
  stringsAsFactors = FALSE
)

n_reps <- 10
results <- list()
idx <- 0

cat(sprintf("Monte Carlo: %d configs x %d reps x %d methods\n",
            nrow(configs), n_reps, length(base_methods) + 1L))

for (cfg_i in seq_len(nrow(configs))) {
  cfg <- configs[cfg_i, ]
  # SBM needs the true k; capture via local scope
  k_true <- cfg$n_communities
  all_methods <- c(base_methods, list(
    "Spectral SBM" = function(L) method_spectral_sbm(L, k_true)
  ))
  for (rep in seq_len(n_reps)) {
    seed <- cfg_i * 1000 + rep
    sim <- simulate_evolving(
      n_nodes = cfg$n_nodes, n_layers = cfg$n_layers,
      n_communities = cfg$n_communities,
      p_switch = cfg$p_switch, seed = seed
    )

    for (mname in names(all_methods)) {
      t0 <- proc.time()["elapsed"]
      det <- tryCatch(all_methods[[mname]](sim$layers),
                      error = function(e) NULL)
      elapsed <- proc.time()["elapsed"] - t0
      if (is.null(det)) next

      idx <- idx + 1
      results[[idx]] <- data.frame(
        n_nodes       = cfg$n_nodes,
        n_layers      = cfg$n_layers,
        n_communities = cfg$n_communities,
        p_switch      = cfg$p_switch,
        rep           = rep,
        method        = mname,
        nmi           = round(compute_nmi(det, sim$true_memberships), 4),
        ari           = round(compute_ari(det, sim$true_memberships), 4),
        runtime_s     = round(elapsed, 4),
        stringsAsFactors = FALSE
      )
    }
  }
  if (cfg_i %% 10 == 0 || cfg_i == nrow(configs)) {
    cat(sprintf("  [%d/%d]\n", cfg_i, nrow(configs)))
  }
}

df <- do.call(rbind, results)
write.csv(df, file.path(outdir, "synthetic_results.csv"), row.names = FALSE)
cat("Raw results saved.\n")


# ===========================================================================
# 5. Summary table
# ===========================================================================

summary_tab <- aggregate(
  cbind(nmi, ari, runtime_s) ~ method + p_switch,
  data = df, FUN = mean
)
summary_tab <- summary_tab[order(summary_tab$p_switch, -summary_tab$nmi), ]
write.csv(summary_tab, file.path(outdir, "tab_main_results.csv"), row.names = FALSE)
cat("Summary table saved.\n")


# ===========================================================================
# 6. Figures
# ===========================================================================

theme_paper <- theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

# Classify methods into groups for cleaner plotting
df$method_group <- ifelse(
  grepl("DynMux", df$method), "Dynamic Multiplex",
  ifelse(df$method == "Independent Louvain", "Independent",
         ifelse(df$method == "Full Coupling", "Full Coupling", "Aggregated"))
)

# Fig 1: NMI vs switching rate
agg1 <- aggregate(nmi ~ method + p_switch, data = df, FUN = mean)
p1 <- ggplot(agg1, aes(x = p_switch, y = nmi, color = method, linetype = method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  labs(x = "Community switching rate (p_switch)",
       y = "Mean NMI",
       color = "Method", linetype = "Method") +
  theme_paper +
  guides(color = guide_legend(ncol = 2))
ggsave(file.path(outdir, "fig_nmi_by_switch.pdf"), p1, width = 8, height = 5)

# Fig 2: ARI vs switching rate
agg2 <- aggregate(ari ~ method + p_switch, data = df, FUN = mean)
p2 <- ggplot(agg2, aes(x = p_switch, y = ari, color = method, linetype = method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  labs(x = "Community switching rate (p_switch)",
       y = "Mean ARI",
       color = "Method", linetype = "Method") +
  theme_paper +
  guides(color = guide_legend(ncol = 2))
ggsave(file.path(outdir, "fig_ari_by_switch.pdf"), p2, width = 8, height = 5)

# Fig 3: NMI vs number of layers (at moderate switching)
df_mod <- df[df$p_switch == 0.05, ]
agg3 <- aggregate(nmi ~ method + n_layers, data = df_mod, FUN = mean)
p3 <- ggplot(agg3, aes(x = n_layers, y = nmi, color = method, linetype = method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  labs(x = "Number of temporal layers",
       y = "Mean NMI (p_switch = 0.05)",
       color = "Method", linetype = "Method") +
  theme_paper +
  guides(color = guide_legend(ncol = 2))
ggsave(file.path(outdir, "fig_nmi_by_layers.pdf"), p3, width = 8, height = 5)

# Fig 4: Runtime comparison
agg4 <- aggregate(runtime_s ~ method + n_nodes, data = df, FUN = mean)
p4 <- ggplot(agg4, aes(x = factor(n_nodes), y = runtime_s, fill = method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  labs(x = "Network size (nodes)", y = "Mean runtime (seconds)", fill = "Method") +
  theme_paper +
  guides(fill = guide_legend(ncol = 2))
ggsave(file.path(outdir, "fig_runtime_comparison.pdf"), p4, width = 9, height = 5)

cat("Figures saved to manuscript/output/\n")
cat("Done.\n")
