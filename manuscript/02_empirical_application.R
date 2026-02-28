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
# Outputs (manuscript/output/):
#   - empirical_trade_fit.rds          : trade community fit object
#   - empirical_igo_fit.rds            : IGO community fit object
#   - empirical_alliance_fit.rds       : alliance community fit object
#   - empirical_community_summary.csv  : community sizes and modularity over time
#   - fig_trade_communities.pdf        : trade community alluvial plot
#   - fig_igo_communities.pdf          : IGO community alluvial plot
#   - fig_alliance_communities.pdf     : alliance community alluvial plot
#   - fig_community_count_over_time.pdf: number of communities over time
#
# Usage:
#   install.packages("peacesciencer")
#   Rscript manuscript/02_empirical_application.R

suppressPackageStartupMessages({
  library(dynamic_multiplex)
  library(ggplot2)
})

outdir <- "manuscript/output"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# ===========================================================================
# Configuration
# ===========================================================================

YEAR_RANGE <- 1965:2000
FIT_TYPE   <- "jaccard"       # interlayer tie strategy
ALGORITHM  <- "louvain"       # community detection algorithm


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
# peacesciencer provides helper functions that add columns to a dyad-year frame.
dyads <- peacesciencer::create_dyadyears(
  subset_years = YEAR_RANGE,
  directed = FALSE
) %>%
  peacesciencer::add_cow_trade() %>%
  peacesciencer::add_cow_alliance()

# For IGOs, peacesciencer provides add_igo_dyadic() which adds the number of
# shared IGO memberships between two states.
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
  # Get the universe of states present across all years
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

# Trade: use the sum of imports + exports (flow1 + flow2) as edge weight.
# peacesciencer add_cow_trade() creates flow1 and flow2 columns.
if (all(c("flow1", "flow2") %in% names(dyads))) {
  dyads$trade_total <- ifelse(is.na(dyads$flow1), 0, dyads$flow1) +
                       ifelse(is.na(dyads$flow2), 0, dyads$flow2)
} else {
  warning("Trade flow columns not found; creating empty trade matrices.")
  dyads$trade_total <- 0
}

trade_data  <- build_yearly_matrices(dyads, YEAR_RANGE, "trade_total")

# Alliance: binary defense pact indicator.
# add_cow_alliance() creates cow_defense (among others).
alliance_col <- if ("cow_defense" %in% names(dyads)) "cow_defense" else {
  # Fallback: look for any alliance column
  alliance_cols <- grep("defense|alliance", names(dyads), value = TRUE, ignore.case = TRUE)
  if (length(alliance_cols) > 0) alliance_cols[1] else {
    warning("No alliance column found; creating empty alliance matrices.")
    dyads$alliance_dummy <- 0
    "alliance_dummy"
  }
}
alliance_data <- build_yearly_matrices(dyads, YEAR_RANGE, alliance_col)

# IGO: number of shared IGO memberships.
# add_igo_dyadic() should create a column for shared IGO count.
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
# 3. Fit dynamic_multiplex models
# ===========================================================================

fit_network <- function(data, network_name, fit_type, algorithm) {
  cat(sprintf("Fitting %s network (fit_type=%s, algorithm=%s) ...\n",
              network_name, fit_type, algorithm))

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
  cat(sprintf("  Done in %.1f seconds.\n", elapsed))

  fit
}

trade_fit    <- fit_network(trade_data,    "trade",    FIT_TYPE, ALGORITHM)
alliance_fit <- fit_network(alliance_data, "alliance", FIT_TYPE, ALGORITHM)
igo_fit      <- fit_network(igo_data,      "IGO",      FIT_TYPE, ALGORITHM)

# Save fit objects
saveRDS(trade_fit,    file.path(outdir, "empirical_trade_fit.rds"))
saveRDS(alliance_fit, file.path(outdir, "empirical_alliance_fit.rds"))
saveRDS(igo_fit,      file.path(outdir, "empirical_igo_fit.rds"))
cat("Fit objects saved.\n")


# ===========================================================================
# 4. Community summary table
# ===========================================================================

summarize_fit <- function(fit, network_name, years) {
  do.call(rbind, lapply(seq_along(fit$layer_communities), function(i) {
    lc <- fit$layer_communities[[i]]
    n_comms <- length(unique(lc$membership))
    mod <- if (is.na(lc$modularity)) NA_real_ else round(lc$modularity, 4)
    sizes <- table(lc$membership)
    data.frame(
      network    = network_name,
      year       = years[i],
      n_communities = n_comms,
      modularity = mod,
      largest_community = max(sizes),
      smallest_community = min(sizes),
      stringsAsFactors = FALSE
    )
  }))
}

comm_summary <- rbind(
  summarize_fit(trade_fit,    "Trade",    YEAR_RANGE),
  summarize_fit(alliance_fit, "Alliance", YEAR_RANGE),
  summarize_fit(igo_fit,      "IGO",      YEAR_RANGE)
)

write.csv(comm_summary, file.path(outdir, "empirical_community_summary.csv"),
          row.names = FALSE)
cat("Community summary table saved.\n")


# ===========================================================================
# 5. Figures
# ===========================================================================

theme_paper <- theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# --- Fig: Number of communities over time ---
p_count <- ggplot(comm_summary, aes(x = year, y = n_communities,
                                     color = network, linetype = network)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.5) +
  labs(x = "Year", y = "Number of communities",
       color = "Network", linetype = "Network") +
  theme_paper
ggsave(file.path(outdir, "fig_community_count_over_time.pdf"),
       p_count, width = 8, height = 4.5)

# --- Fig: Modularity over time ---
mod_data <- comm_summary[!is.na(comm_summary$modularity), ]
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

# --- Alluvial plots (if ggalluvial is available) ---
has_alluvial <- requireNamespace("ggalluvial", quietly = TRUE) &&
                requireNamespace("RColorBrewer", quietly = TRUE)

if (has_alluvial) {
  save_alluvial <- function(fit, filename, max_nodes = 80) {
    p <- plot_multilayer_alluvial(fit, max_nodes = max_nodes, palette = "Dark2")
    ggsave(file.path(outdir, filename), p, width = 10, height = 6)
  }

  save_alluvial(trade_fit,    "fig_trade_communities.pdf")
  save_alluvial(alliance_fit, "fig_alliance_communities.pdf")
  save_alluvial(igo_fit,      "fig_igo_communities.pdf")
  cat("Alluvial plots saved.\n")
} else {
  cat("Skipping alluvial plots (install ggalluvial and RColorBrewer).\n")
}

cat("\nEmpirical application complete.\n")
cat(sprintf("Outputs in: %s/\n", outdir))
