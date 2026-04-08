#!/usr/bin/env Rscript
#
# Manuscript Section: Empirical Application
#
# Applies dynamic_multiplex to temporal international-relations networks
# constructed from the peacesciencer R package:
#   - Trade networks     (COW bilateral trade)
#   - IGO networks       (shared IGO memberships)
#   - Alliance networks  (ATOP defense pacts)
#
# Each network type is treated as a sequence of yearly layers, and
# dynamic_multiplex identifies evolving communities over time.
#
# Extended analysis:
#   - Fits all 5 interlayer weighting strategies for robustness
#   - Compares against multislice Louvain (Mucha et al. 2010)
#   - Computes cross-model agreement (NMI between specifications)
#   - Measures community persistence across time
#
# Outputs (manuscript/output/):
#   - empirical_fits_all.rds                 : all fit objects
#   - empirical_community_summary.csv        : community sizes and modularity
#   - empirical_robustness_agreement.csv     : cross-spec NMI agreement
#   - empirical_persistence.csv              : community persistence over time
#   - empirical_runtime.csv                  : runtime comparison
#   - fig_community_count_over_time.pdf      : community count by year
#   - fig_modularity_over_time.pdf           : modularity by year
#   - fig_robustness_agreement.pdf           : cross-spec agreement heatmap
#   - fig_persistence_over_time.pdf          : community persistence
#   - fig_*_communities.pdf                  : alluvial plots
#
# Usage:
#   install.packages("peacesciencer")
#   Rscript manuscript/02_empirical_application.R

suppressPackageStartupMessages({
  library(dynamicmultiplex)
  library(ggplot2)
})

outdir <- "manuscript/output"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# ===========================================================================
# Configuration
# ===========================================================================

YEAR_RANGE <- 1965:2000
ALGORITHM  <- "louvain"

# All interlayer tie strategies for robustness comparison
FIT_TYPES <- c("jaccard", "overlap", "weighted_jaccard",
               "weighted_overlap", "identity")


# ===========================================================================
# 1. Load and prepare data from peacesciencer
# ===========================================================================

cat("Loading peacesciencer data ...\n")

if (!requireNamespace("peacesciencer", quietly = TRUE)) {
  stop("Install peacesciencer:  install.packages('peacesciencer')", call. = FALSE)
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop("Install dplyr:  install.packages('dplyr')", call. = FALSE)
}
if (!requireNamespace("tidyr", quietly = TRUE)) {
  stop("Install tidyr:  install.packages('tidyr')", call. = FALSE)
}

library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)

# Build undirected dyad-year data with trade, IGO, and alliance variables.
dyads <- peacesciencer::create_dyadyears(
  subset_years = YEAR_RANGE,
  directed = FALSE
) %>%
  peacesciencer::add_cow_trade() %>%
  peacesciencer::add_cow_alliance()

dyads <- dyads %>%
  peacesciencer::add_igo_dyadic()

cat(sprintf("  Dyad-years loaded: %d rows, years %d-%d\n",
            nrow(dyads), min(YEAR_RANGE), max(YEAR_RANGE)))


# ===========================================================================
# 2. Build yearly adjacency matrices
# ===========================================================================

build_yearly_matrices <- function(dyads, year_range, value_col,
                                   state_col1 = "ccode1",
                                   state_col2 = "ccode2",
                                   year_col = "year") {
  all_states <- sort(unique(c(dyads[[state_col1]], dyads[[state_col2]])))
  n <- length(all_states)
  state_idx <- setNames(seq_along(all_states), as.character(all_states))

  layers <- lapply(year_range, function(yr) {
    sub <- dyads[dyads[[year_col]] == yr, ]
    mat <- matrix(0, nrow = n, ncol = n)

    if (nrow(sub) > 0 && value_col %in% names(sub)) {
      vals <- sub[[value_col]]
      vals[is.na(vals)] <- 0
      i_idx <- state_idx[as.character(sub[[state_col1]])]
      j_idx <- state_idx[as.character(sub[[state_col2]])]
      valid <- !is.na(i_idx) & !is.na(j_idx)
      for (k in which(valid)) {
        mat[i_idx[k], j_idx[k]] <- vals[k]
        mat[j_idx[k], i_idx[k]] <- vals[k]
      }
    }
    mat
  })

  names(layers) <- as.character(year_range)
  list(layers = layers, states = all_states)
}

cat("Building yearly adjacency matrices ...\n")

# Trade: sum of imports + exports
if (all(c("flow1", "flow2") %in% names(dyads))) {
  dyads$trade_total <- ifelse(is.na(dyads$flow1), 0, dyads$flow1) +
                       ifelse(is.na(dyads$flow2), 0, dyads$flow2)
} else {
  warning("Trade flow columns not found; creating empty trade matrices.")
  dyads$trade_total <- 0
}

trade_data <- build_yearly_matrices(dyads, YEAR_RANGE, "trade_total")

# Alliance: binary defense pact indicator
alliance_col <- if ("cow_defense" %in% names(dyads)) "cow_defense" else {
  alliance_cols <- grep("defense|alliance", names(dyads), value = TRUE, ignore.case = TRUE)
  if (length(alliance_cols) > 0) alliance_cols[1] else {
    warning("No alliance column found; creating empty alliance matrices.")
    dyads$alliance_dummy <- 0
    "alliance_dummy"
  }
}
alliance_data <- build_yearly_matrices(dyads, YEAR_RANGE, alliance_col)

# IGO: shared IGO memberships
igo_col <- intersect(c("dyadigos", "IGOcount", "igo_count", "SharedIGO"), names(dyads))
if (length(igo_col) == 0) {
  igo_candidates <- grep("igo", names(dyads), value = TRUE, ignore.case = TRUE)
  igo_col <- if (length(igo_candidates) > 0) igo_candidates[1] else {
    warning("No IGO column found; creating empty IGO matrices.")
    dyads$igo_dummy <- 0
    "igo_dummy"
  }
} else {
  igo_col <- igo_col[1]
}
igo_data <- build_yearly_matrices(dyads, YEAR_RANGE, igo_col)

cat(sprintf("  Trade matrices: %d layers, %d states\n",
            length(trade_data$layers), length(trade_data$states)))
cat(sprintf("  Alliance matrices: %d layers, %d states\n",
            length(alliance_data$layers), length(alliance_data$states)))
cat(sprintf("  IGO matrices: %d layers, %d states\n",
            length(igo_data$layers), length(igo_data$states)))


# ===========================================================================
# 3. Fit all interlayer weighting specifications (robustness)
# ===========================================================================

fit_network <- function(data, network_name, fit_type, algorithm) {
  cat(sprintf("  Fitting %s / %s ...\n", network_name, fit_type))

  fit_fn <- switch(
    fit_type,
    jaccard          = fit_multilayer_jaccard,
    overlap          = fit_multilayer_overlap,
    weighted_jaccard = fit_multilayer_weighted_jaccard,
    weighted_overlap = fit_multilayer_weighted_overlap,
    identity         = fit_multilayer_identity_ties
  )

  t0 <- proc.time()["elapsed"]
  fit <- fit_fn(data$layers, algorithm = algorithm)
  elapsed <- proc.time()["elapsed"] - t0

  list(fit = fit, runtime = elapsed)
}

# Multislice all-to-all (Mucha et al. 2010) for comparison
fit_multislice <- function(data, network_name, omega = 1.0) {
  cat(sprintf("  Fitting %s / multislice_all_to_all ...\n", network_name))
  n_nodes <- nrow(data$layers[[1]])
  n_layers <- length(data$layers)
  N <- n_nodes * n_layers

  t0 <- proc.time()["elapsed"]

  supra <- matrix(0, N, N)
  for (t in seq_len(n_layers)) {
    r0 <- (t - 1) * n_nodes
    supra[(r0 + 1):(r0 + n_nodes), (r0 + 1):(r0 + n_nodes)] <- data$layers[[t]]
  }
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
  elapsed <- proc.time()["elapsed"] - t0

  # Convert supra partition to per-layer fit structure
  layer_communities <- lapply(seq_len(n_layers), function(t) {
    r0 <- (t - 1) * n_nodes
    mem <- supra_mem[(r0 + 1):(r0 + n_nodes)]
    list(
      membership = mem,
      modularity = NA_real_,
      communities = split(seq_along(mem), mem)
    )
  })

  list(
    fit = list(
      algorithm = "louvain_multislice",
      layer_communities = layer_communities,
      directed = FALSE
    ),
    runtime = elapsed
  )
}

cat("\nFitting all specifications across all networks ...\n")

# Store all fits in a nested list: all_fits[[network]][[fit_type]]
networks <- list(Trade = trade_data, Alliance = alliance_data, IGO = igo_data)
all_fits <- list()
runtime_records <- list()
ri <- 0

for (net_name in names(networks)) {
  all_fits[[net_name]] <- list()

  for (ft in FIT_TYPES) {
    res <- fit_network(networks[[net_name]], net_name, ft, ALGORITHM)
    all_fits[[net_name]][[ft]] <- res$fit
    ri <- ri + 1
    runtime_records[[ri]] <- data.frame(
      network = net_name, method = ft, runtime_s = round(res$runtime, 2),
      stringsAsFactors = FALSE
    )
  }

  # Multislice all-to-all comparison
  res_ms <- fit_multislice(networks[[net_name]], net_name)
  all_fits[[net_name]][["multislice_all_to_all"]] <- res_ms$fit
  ri <- ri + 1
  runtime_records[[ri]] <- data.frame(
    network = net_name, method = "multislice_all_to_all",
    runtime_s = round(res_ms$runtime, 2), stringsAsFactors = FALSE
  )
}

saveRDS(all_fits, file.path(outdir, "empirical_fits_all.rds"))
runtime_df <- do.call(rbind, runtime_records)
write.csv(runtime_df, file.path(outdir, "empirical_runtime.csv"), row.names = FALSE)
cat("All fits and runtimes saved.\n")


# ===========================================================================
# 4. Community summary table
# ===========================================================================

summarize_fit <- function(fit, network_name, fit_type, years) {
  do.call(rbind, lapply(seq_along(fit$layer_communities), function(i) {
    lc <- fit$layer_communities[[i]]
    n_comms <- length(unique(lc$membership))
    mod_val <- lc$modularity
    mod <- if (is.null(mod_val) || is.na(mod_val)) NA_real_ else round(mod_val, 4)
    sizes <- table(lc$membership)
    data.frame(
      network    = network_name,
      fit_type   = fit_type,
      year       = years[i],
      n_communities = n_comms,
      modularity = mod,
      largest_community = max(sizes),
      smallest_community = min(sizes),
      stringsAsFactors = FALSE
    )
  }))
}

comm_summary_list <- list()
si <- 0
for (net_name in names(all_fits)) {
  for (ft in names(all_fits[[net_name]])) {
    si <- si + 1
    comm_summary_list[[si]] <- summarize_fit(
      all_fits[[net_name]][[ft]], net_name, ft, YEAR_RANGE
    )
  }
}
comm_summary <- do.call(rbind, comm_summary_list)

write.csv(comm_summary, file.path(outdir, "empirical_community_summary.csv"),
          row.names = FALSE)
cat("Community summary table saved.\n")


# ===========================================================================
# 5. Cross-specification agreement (NMI between fit types)
# ===========================================================================

# For each network and year, compute NMI between all pairs of specifications.
# This measures whether substantive conclusions are robust to weighting choice.

all_methods <- c(FIT_TYPES, "multislice_all_to_all")

agreement_records <- list()
ai <- 0

for (net_name in names(all_fits)) {
  for (yr_i in seq_along(YEAR_RANGE)) {
    for (m1_i in seq_len(length(all_methods) - 1)) {
      for (m2_i in (m1_i + 1):length(all_methods)) {
        m1 <- all_methods[m1_i]
        m2 <- all_methods[m2_i]

        mem1 <- all_fits[[net_name]][[m1]]$layer_communities[[yr_i]]$membership
        mem2 <- all_fits[[net_name]][[m2]]$layer_communities[[yr_i]]$membership

        nmi_val <- igraph::compare(mem1, mem2, method = "nmi")

        ai <- ai + 1
        agreement_records[[ai]] <- data.frame(
          network = net_name, year = YEAR_RANGE[yr_i],
          method1 = m1, method2 = m2, nmi = round(nmi_val, 4),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

agreement_df <- do.call(rbind, agreement_records)
write.csv(agreement_df, file.path(outdir, "empirical_robustness_agreement.csv"),
          row.names = FALSE)
cat("Cross-specification agreement saved.\n")


# ===========================================================================
# 6. Community persistence over time
# ===========================================================================

# For each network and fit_type, compute NMI between adjacent years.
# High NMI = communities are stable over time.

persistence_records <- list()
pi <- 0

for (net_name in names(all_fits)) {
  for (ft in names(all_fits[[net_name]])) {
    fit <- all_fits[[net_name]][[ft]]
    for (yr_i in seq_len(length(YEAR_RANGE) - 1)) {
      mem_t  <- fit$layer_communities[[yr_i]]$membership
      mem_t1 <- fit$layer_communities[[yr_i + 1]]$membership

      nmi_val <- igraph::compare(mem_t, mem_t1, method = "nmi")

      pi <- pi + 1
      persistence_records[[pi]] <- data.frame(
        network = net_name, fit_type = ft,
        year_from = YEAR_RANGE[yr_i], year_to = YEAR_RANGE[yr_i + 1],
        persistence_nmi = round(nmi_val, 4),
        stringsAsFactors = FALSE
      )
    }
  }
}

persistence_df <- do.call(rbind, persistence_records)
write.csv(persistence_df, file.path(outdir, "empirical_persistence.csv"),
          row.names = FALSE)
cat("Community persistence saved.\n")


# ===========================================================================
# 7. Figures
# ===========================================================================

theme_paper <- theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

# --- Fig: Number of communities over time (Jaccard baseline) ---
jaccard_summary <- comm_summary[comm_summary$fit_type == "jaccard", ]
p_count <- ggplot(jaccard_summary, aes(x = year, y = n_communities,
                                       color = network, linetype = network)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.5) +
  labs(x = "Year", y = "Number of communities",
       color = "Network", linetype = "Network",
       title = "Community count over time (Jaccard specification)") +
  theme_paper
ggsave(file.path(outdir, "fig_community_count_over_time.pdf"),
       p_count, width = 8, height = 4.5)

# --- Fig: Modularity over time ---
mod_data <- jaccard_summary[!is.na(jaccard_summary$modularity), ]
if (nrow(mod_data) > 0) {
  p_mod <- ggplot(mod_data, aes(x = year, y = modularity,
                                 color = network, linetype = network)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.5) +
    labs(x = "Year", y = "Modularity",
         color = "Network", linetype = "Network") +
    theme_paper
  ggsave(file.path(outdir, "fig_modularity_over_time.pdf"),
         p_mod, width = 8, height = 4.5)
}

# --- Fig: Robustness — mean cross-specification agreement ---
agg_agree <- aggregate(nmi ~ network + method1 + method2,
                        data = agreement_df, FUN = mean)
# Create symmetric matrix labels for display
p_robust <- ggplot(agg_agree, aes(x = method1, y = method2, fill = nmi)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(nmi, 2)), size = 2.5) +
  scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 1)) +
  facet_wrap(~network) +
  labs(x = "", y = "", fill = "Mean NMI",
       title = "Cross-specification agreement (mean NMI across years)") +
  theme_paper +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7))
ggsave(file.path(outdir, "fig_robustness_agreement.pdf"),
       p_robust, width = 12, height = 5)

# --- Fig: Community persistence over time ---
persist_jaccard <- persistence_df[persistence_df$fit_type == "jaccard", ]
p_persist <- ggplot(persist_jaccard, aes(x = year_from, y = persistence_nmi,
                                          color = network, linetype = network)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.5) +
  labs(x = "Year", y = "Persistence (NMI with next year)",
       color = "Network", linetype = "Network",
       title = "Community persistence: year-to-year stability") +
  ylim(0, 1) +
  theme_paper
ggsave(file.path(outdir, "fig_persistence_over_time.pdf"),
       p_persist, width = 8, height = 4.5)

# --- Fig: Community count comparison across specifications ---
p_count_robust <- ggplot(comm_summary, aes(x = year, y = n_communities,
                                            color = fit_type, linetype = fit_type)) +
  geom_line(linewidth = 0.6, alpha = 0.8) +
  facet_wrap(~network) +
  labs(x = "Year", y = "Number of communities",
       color = "Specification", linetype = "Specification",
       title = "Robustness: community count across interlayer weighting specs") +
  theme_paper +
  guides(color = guide_legend(ncol = 3))
ggsave(file.path(outdir, "fig_community_count_robustness.pdf"),
       p_count_robust, width = 10, height = 5)

# --- Alluvial plots (if ggalluvial available) ---
has_alluvial <- requireNamespace("ggalluvial", quietly = TRUE) &&
                requireNamespace("RColorBrewer", quietly = TRUE)

if (has_alluvial) {
  # Save alluvial for Jaccard baseline
  save_alluvial <- function(fit, filename, max_nodes = 80) {
    p <- plot_multilayer_alluvial(fit, max_nodes = max_nodes, palette = "Dark2")
    ggsave(file.path(outdir, filename), p, width = 10, height = 6)
  }

  save_alluvial(all_fits$Trade$jaccard,    "fig_trade_communities.pdf")
  save_alluvial(all_fits$Alliance$jaccard, "fig_alliance_communities.pdf")
  save_alluvial(all_fits$IGO$jaccard,      "fig_igo_communities.pdf")
  cat("Alluvial plots saved.\n")
} else {
  cat("Skipping alluvial plots (install ggalluvial and RColorBrewer).\n")
}

# ===========================================================================
# 8. Bootstrap confidence intervals for empirical fits
# ===========================================================================

cat("\n=== Bootstrap CIs for empirical data ===\n")

boot_empirical_records <- list()
bi <- 0

for (net_name in names(networks)) {
  cat(sprintf("  Bootstrap %s (n_boot=50) ...\n", net_name))
  boot <- bootstrap_multilayer(
    networks[[net_name]]$layers,
    fit_type = "jaccard",
    algorithm = ALGORITHM,
    n_boot = 50,
    seed = 123
  )
  ci <- community_ci(boot, alpha = 0.05)

  ci$modularity_ci$network <- net_name
  ci$community_count_ci$network <- net_name
  ci$mean_node_stability$network <- net_name
  ci$modularity_ci$year <- YEAR_RANGE
  ci$community_count_ci$year <- YEAR_RANGE
  ci$mean_node_stability$year <- YEAR_RANGE

  bi <- bi + 1
  boot_empirical_records[[bi]] <- list(
    modularity = ci$modularity_ci,
    count = ci$community_count_ci,
    stability = ci$mean_node_stability
  )
}

boot_mod_df <- do.call(rbind, lapply(boot_empirical_records, `[[`, "modularity"))
boot_count_df <- do.call(rbind, lapply(boot_empirical_records, `[[`, "count"))
boot_stab_df <- do.call(rbind, lapply(boot_empirical_records, `[[`, "stability"))

write.csv(boot_mod_df, file.path(outdir, "empirical_bootstrap_modularity_ci.csv"),
          row.names = FALSE)
write.csv(boot_count_df, file.path(outdir, "empirical_bootstrap_count_ci.csv"),
          row.names = FALSE)
write.csv(boot_stab_df, file.path(outdir, "empirical_bootstrap_stability.csv"),
          row.names = FALSE)

# Fig: Modularity CIs over time by network
p_boot_mod <- ggplot(boot_mod_df, aes(x = year, y = estimate,
                                        color = network, fill = network)) +
  geom_line(linewidth = 0.8) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  labs(x = "Year", y = "Modularity",
       color = "Network", fill = "Network",
       title = "Bootstrap 95% CI for modularity over time") +
  theme_paper
ggsave(file.path(outdir, "fig_empirical_bootstrap_modularity.pdf"),
       p_boot_mod, width = 8, height = 4.5)

# Fig: Community count CIs over time
p_boot_count <- ggplot(boot_count_df, aes(x = year, y = estimate,
                                            color = network, fill = network)) +
  geom_line(linewidth = 0.8) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  labs(x = "Year", y = "Number of communities",
       color = "Network", fill = "Network",
       title = "Bootstrap 95% CI for community count over time") +
  theme_paper
ggsave(file.path(outdir, "fig_empirical_bootstrap_count.pdf"),
       p_boot_count, width = 8, height = 4.5)

# Fig: Mean stability over time
p_boot_stab <- ggplot(boot_stab_df, aes(x = year, y = mean_stability,
                                          color = network, linetype = network)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  labs(x = "Year", y = "Mean node stability",
       color = "Network", linetype = "Network",
       title = "Bootstrap node stability over time") +
  ylim(0, 1) +
  theme_paper
ggsave(file.path(outdir, "fig_empirical_bootstrap_stability.pdf"),
       p_boot_stab, width = 8, height = 4.5)

cat("Empirical bootstrap results and figures saved.\n")

cat("\nEmpirical application complete.\n")
cat(sprintf("Outputs in: %s/\n", outdir))
