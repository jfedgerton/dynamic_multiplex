# =============================================================================
# manuscript/replot_from_cache.R  (v2)
#
# Regenerates manuscript figures FROM CACHED RESULTS in manuscript/output/.
# Design system:
#   - DynMux variants highlighted in color (Okabe-Ito, one color per
#     similarity measure); Louvain = solid, Leiden = dashed.
#   - All competitor methods drawn as thin grey lines ("Competitors").
#   - Legends at bottom, 2 columns.
#   - Facet strip text enlarged.
#   - No underscores in titles / axis / strip / legend text.
#   - Log axes in 10^x notation; large-count axes comma-formatted.
#   - Widths ~7.5 in.
# Sections needing ROAR-only files are guarded with file.exists().
#
# Usage (from repo root):  Rscript manuscript/replot_from_cache.R
# =============================================================================

set.seed(123)
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

outdir <- "manuscript/output"
if (!dir.exists(outdir)) stop("Run from the repo root: manuscript/output not found.")

# --- Style system ------------------------------------------------------------
SIM_COLS <- c(
  "Jaccard"      = "#E69F00",
  "Overlap"      = "#56B4E9",
  "Wt. Jaccard"  = "#009E73",
  "Wt. Overlap"  = "#D55E00",
  "Identity"     = "#CC79A7",
  "Competitors"  = "grey55"
)
LT_VALS <- c("Louvain" = "solid", "Leiden" = "dashed", "Competitors" = "solid")

# Classify a method name into the highlight scheme.
classify_methods <- function(df) {
  m <- df$method
  is_dyn <- grepl("^DynMux", m)
  sim <- rep("Competitors", length(m))
  sim[is_dyn & grepl("Wt\\. Jaccard", m)] <- "Wt. Jaccard"
  sim[is_dyn & grepl("Wt\\. Overlap", m)] <- "Wt. Overlap"
  sim[is_dyn & grepl("\\(Jaccard", m)]    <- "Jaccard"
  sim[is_dyn & grepl("\\(Overlap", m)]    <- "Overlap"
  sim[is_dyn & grepl("Identity", m)]      <- "Identity"
  alg <- ifelse(!is_dyn, "Competitors",
                ifelse(grepl("Leiden", m), "Leiden", "Louvain"))
  df$sim <- factor(sim, levels = names(SIM_COLS))
  df$alg <- factor(alg, levels = names(LT_VALS))
  # Draw competitors first (background), DynMux on top
  df[order(df$sim != "Competitors"), ]
}

de_underscore <- function(x) gsub("_", " ", x)
DGP_LABS <- c(drift = "Drift", regime_shift = "Regime shift",
              switching = "Switching")
FIT_LABS <- c(jaccard = "Jaccard", overlap = "Overlap",
              weighted_jaccard = "Weighted Jaccard",
              weighted_overlap = "Weighted Overlap",
              identity = "Identity",
              multislice_all_to_all = "Multislice (all-to-all)")

theme_paper <- theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 12)
  )

guides_2col <- guides(
  color = guide_legend(ncol = 3, order = 1),
  fill = guide_legend(ncol = 3, order = 1),
  linetype = guide_legend(ncol = 1, order = 2)
)

scale_dynmux <- list(
  scale_color_manual(values = SIM_COLS, name = NULL, drop = FALSE),
  scale_linetype_manual(values = LT_VALS, name = NULL,
                        breaks = c("Louvain", "Leiden"))
)

msg <- function(...) cat(sprintf(...), "\n")

# Shared builder for method-line figures.
# The DynMux variants often sit exactly on top of one another; point markers
# are staggered slightly along x so every variant stays visible, while the
# lines remain at their true positions.
dynmux_lineplot <- function(d, xvar, yvar, xlab, ylab, title_str,
                            subtitle_str = NULL, logx = FALSE) {
  d <- classify_methods(d)
  dyn <- subset(d, sim != "Competitors")
  cmp <- subset(d, sim == "Competitors")
  if (nrow(dyn) > 0) {
    k <- as.integer(dyn$sim) - 3L   # -2..2 across the 5 similarities
    if (logx) {
      dyn$x_pt <- dyn[[xvar]] * (1 + 0.012 * k)
    } else {
      rng <- diff(range(d[[xvar]], na.rm = TRUE))
      dyn$x_pt <- dyn[[xvar]] + 0.006 * rng * k
    }
  }
  ggplot(d, aes(x = .data[[xvar]], y = .data[[yvar]],
                color = sim, linetype = alg, group = method)) +
    geom_line(data = cmp, linewidth = 0.65, alpha = 0.9) +
    geom_line(data = dyn, linewidth = 0.9, alpha = 0.9) +
    geom_point(data = dyn, aes(x = x_pt), size = 1.8) +
    scale_dynmux +
    labs(x = xlab, y = ylab, title = title_str, subtitle = subtitle_str) +
    theme_paper + guides_2col
}

# =============================================================================
# 1. Synthetic experiments (source: synthetic_results.csv)
# =============================================================================
f <- file.path(outdir, "synthetic_results.csv")
if (file.exists(f)) {
  df <- read.csv(f, stringsAsFactors = FALSE)

  agg1 <- aggregate(nmi ~ method + p_switch, data = df, FUN = mean)
  p1 <- dynmux_lineplot(agg1, "p_switch", "nmi",
                        "Community switching rate", "Mean NMI", NULL)
  ggsave(file.path(outdir, "fig_nmi_by_switch.pdf"), p1, width = 7.5, height = 5.5)

  agg2 <- aggregate(ari ~ method + p_switch, data = df, FUN = mean)
  p2 <- dynmux_lineplot(agg2, "p_switch", "ari",
                        "Community switching rate", "Mean ARI", NULL)
  ggsave(file.path(outdir, "fig_ari_by_switch.pdf"), p2, width = 7.5, height = 5.5)

  df_mod <- df[df$p_switch == 0.05, ]
  agg3 <- aggregate(nmi ~ method + n_layers, data = df_mod, FUN = mean)
  p3 <- dynmux_lineplot(agg3, "n_layers", "nmi",
                        "Number of temporal layers",
                        "Mean NMI (switching rate 0.05)", NULL)
  ggsave(file.path(outdir, "fig_nmi_by_layers.pdf"), p3, width = 7.5, height = 5.5)

  # Runtime comparison: grouped bars, DynMux colored / competitors grey
  agg4 <- classify_methods(aggregate(runtime_s ~ method + n_nodes, data = df,
                                     FUN = mean))
  p4 <- ggplot(agg4, aes(x = factor(n_nodes), y = runtime_s, fill = sim,
                         group = method)) +
    geom_col(position = position_dodge(width = 0.85), width = 0.8,
             color = "white", linewidth = 0.2) +
    scale_fill_manual(values = SIM_COLS, name = NULL, drop = FALSE) +
    scale_y_log10(labels = label_log(digits = 2)) +
    labs(x = "Network size (nodes)", y = "Mean runtime (seconds, log scale)") +
    theme_paper + guides_2col
  ggsave(file.path(outdir, "fig_runtime_comparison.pdf"), p4,
         width = 7.5, height = 5)

  msg("[1] synthetic figures done")
} else msg("[1] SKIP: %s not found", f)

# =============================================================================
# 2. Scalability (source: scalability_results.csv)
# =============================================================================
f <- file.path(outdir, "scalability_results.csv")
if (file.exists(f)) {
  df <- read.csv(f, stringsAsFactors = FALSE)
  agg <- df %>%
    filter(status == "ok") %>%
    group_by(n, method) %>%
    summarise(time_mean = mean(time), time_sd = sd(time),
              nmi_mean = mean(nmi), nmi_sd = sd(nmi), .groups = "drop") %>%
    as.data.frame()

  p_time <- dynmux_lineplot(agg, "n", "time_mean",
                            "Nodes per layer (log scale)",
                            "Wall-clock fit time (seconds, log scale)",
                            "Runtime scaling across network size",
                            logx = TRUE) +
    scale_x_log10(breaks = c(50, 100, 200, 400, 800)) +
    scale_y_log10(labels = label_log(digits = 2))
  ggsave(file.path(outdir, "fig_scalability.pdf"), p_time,
         width = 7.5, height = 6)

  p_acc <- dynmux_lineplot(agg, "n", "nmi_mean",
                           "Nodes per layer (log scale)",
                           "Mean NMI vs ground truth",
                           "Recovery accuracy across network size",
                           logx = TRUE) +
    scale_x_log10(breaks = c(50, 100, 200, 400, 800))
  ggsave(file.path(outdir, "fig_scalability_acc.pdf"), p_acc,
         width = 7.5, height = 6)

  msg("[2] scalability figures done")
} else msg("[2] SKIP: %s not found", f)

# =============================================================================
# 3. DGP misspecification (sources: dgp_misspec_results.csv, tab_dgp_misspec.csv)
# =============================================================================
f  <- file.path(outdir, "dgp_misspec_results.csv")
ft <- file.path(outdir, "tab_dgp_misspec.csv")
if (file.exists(f) && file.exists(ft)) {
  df  <- read.csv(f,  stringsAsFactors = FALSE)
  agg <- read.csv(ft, stringsAsFactors = FALSE)
  df_ok <- df %>% filter(status == "ok")
  fix_group <- function(g, m) ifelse(m %in% c("Aggregated Leiden", "Independent Leiden"),
                                     "Temporally agnostic", g)
  df_ok$group <- fix_group(df_ok$group, df_ok$method)
  agg$group   <- fix_group(agg$group, agg$method)

  # Headline boxplot: facets stacked (one DGP per row) so the method axis has
  # room; x labels at 90 degrees; fill by temporal-awareness group.
  grp_cols <- c("Temporally aware" = "#E69F00",
                "Temporally agnostic" = "#56B4E9",
                "Generative (SBM)" = "#009E73",
                "Other" = "#999999")
  p_nmi <- ggplot(df_ok, aes(x = method, y = nmi, fill = group)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.9, linewidth = 0.3) +
    facet_wrap(~ dgp, ncol = 1, labeller = labeller(dgp = DGP_LABS)) +
    scale_fill_manual(values = grp_cols, name = NULL) +
    labs(x = NULL, y = "NMI vs ground truth",
         title = "Recovery accuracy across data-generating processes") +
    theme_paper +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                     size = 7)) +
    guides(fill = guide_legend(ncol = 2))
  ggsave(file.path(outdir, "fig_dgp_misspec_nmi.pdf"), p_nmi,
         width = 7.5, height = 9)

  # Temporal-coherence diagnostic with repelled labels
  agg$dgp_lab <- DGP_LABS[agg$dgp]
  # Break labels onto two lines at a natural point (before the parenthetical
  # or between words) so they occupy less horizontal space.
  wrap_method <- function(s) {
    s <- sub(" \\(", "\n(", s)
    vapply(strsplit(s, "\n"), function(parts) {
      paste(unlist(lapply(parts, function(p)
        paste(strwrap(p, width = 16), collapse = "\n"))), collapse = "\n")
    }, character(1))
  }
  agg$method_lab <- wrap_method(agg$method)
  # Label each method ONCE (at its right-most point) instead of once per DGP;
  # leader lines connect the label to that point.
  lab_df <- agg %>% group_by(method) %>% slice_max(tc_truth, n = 1) %>% ungroup()
  label_layer <- if (requireNamespace("ggrepel", quietly = TRUE)) {
    ggrepel::geom_text_repel(data = lab_df, aes(label = method_lab),
                             size = 2.4, lineheight = 0.85,
                             max.overlaps = 40, seed = 123,
                             box.padding = 0.75, point.padding = 0.55,
                             force = 6, force_pull = 0.15,
                             min.segment.length = 0.1,
                             segment.color = "grey70", segment.size = 0.25)
  } else {
    geom_text(data = lab_df, aes(label = method_lab), vjust = -0.9,
              size = 2.4, lineheight = 0.85, check_overlap = TRUE)
  }
  p_tc <- ggplot(agg, aes(x = tc_truth, y = tc_method, color = group,
                          shape = dgp_lab)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "grey50") +
    geom_point(size = 2.6, alpha = 0.9) +
    label_layer +
    scale_color_manual(values = grp_cols, name = NULL) +
    scale_shape_discrete(name = NULL) +
    coord_cartesian(clip = "off") +
    scale_x_continuous(expand = expansion(mult = c(0.14, 0.16))) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.14))) +
    labs(x = "Ground-truth temporal coherence",
         y = "Method's temporal coherence",
         title = "Does the method preserve the right amount of temporal structure?",
         subtitle = "On the dashed line = matches truth. Above = over-pools. Below = under-couples.") +
    theme_paper +
    guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2))
  ggsave(file.path(outdir, "fig_dgp_misspec_temporal.pdf"), p_tc,
         width = 7.5, height = 7.5)

  msg("[3] dgp misspec figures done")
} else msg("[3] SKIP: dgp files not found")

# =============================================================================
# 4. Parameter sensitivity + shuffle test
# =============================================================================
f <- file.path(outdir, "tab_param_sensitivity.csv")
if (file.exists(f)) {
  agg <- read.csv(f, stringsAsFactors = FALSE)

  make_sens_plot <- function(d, x_lab, file_base, title_str) {
    p <- dynmux_lineplot(as.data.frame(d), "param_x", "nmi_mean",
                         x_lab, "NMI vs ground truth", title_str) +
      geom_errorbar(aes(ymin = nmi_mean - nmi_se, ymax = nmi_mean + nmi_se),
                    width = 0.02, alpha = 0.6, linetype = "solid")
    ggsave(file.path(outdir, paste0(file_base, ".pdf")), p,
           width = 7.5, height = 6)
  }
  make_sens_plot(agg %>% filter(sweep == "density"),
                 "Density ratio (within / between)", "fig_sens_density",
                 "Sensitivity to network density ratio (switching DGP)")
  make_sens_plot(agg %>% filter(sweep == "drift_rate"),
                 "Drift rate per layer", "fig_sens_drift",
                 "Sensitivity to drift rate (drift DGP)")
  make_sens_plot(agg %>% filter(sweep == "break_position"),
                 "Break occurs after layer k", "fig_sens_break",
                 "Sensitivity to break position (regime-shift DGP)")

  agg_break <- agg %>% filter(sweep == "break_position") %>% as.data.frame()
  p_bc <- dynmux_lineplot(agg_break, "param_x", "break_contrast_mean",
                          "Break occurs after layer k",
                          "Within-regime minus across-regime NMI",
                          "Break detection across break positions") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")
  ggsave(file.path(outdir, "fig_sens_break_contrast.pdf"), p_bc,
         width = 7.5, height = 6)

  msg("[4a] sensitivity figures done")
} else msg("[4a] SKIP: %s not found", f)

f <- file.path(outdir, "shuffle_test_results.csv")
if (file.exists(f)) {
  shuffle_df <- classify_methods(read.csv(f, stringsAsFactors = FALSE))
  p_shuffle <- ggplot(shuffle_df,
                      aes(x = tc_shuffled, y = tc_original, color = sim)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "grey50") +
    geom_point(data = subset(shuffle_df, sim == "Competitors"),
               alpha = 0.4, size = 1.4) +
    geom_point(data = subset(shuffle_df, sim != "Competitors"),
               aes(shape = alg), alpha = 0.75, size = 1.8) +
    facet_wrap(~ dgp, ncol = 3, labeller = labeller(dgp = DGP_LABS)) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    scale_color_manual(values = SIM_COLS, name = NULL, drop = FALSE) +
    scale_shape_manual(values = c(Louvain = 16, Leiden = 17), name = NULL) +
    labs(x = "Temporal coherence after shuffling layer order",
         y = "Temporal coherence on original layer order",
         title = "Convergent validity: does the method use layer order?",
         subtitle = "On the diagonal = ignores layer order. Above = exploits it.") +
    theme_paper +
    guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2))
  ggsave(file.path(outdir, "fig_shuffle_test.pdf"), p_shuffle,
         width = 8, height = 4.8)
  msg("[4b] shuffle test figure done")
} else msg("[4b] SKIP: %s not found", f)

# =============================================================================
# 5. Multislice omega sweep (source: tab_multislice_omega.csv)
# =============================================================================
f <- file.path(outdir, "tab_multislice_omega.csv")
if (file.exists(f)) {
  agg <- read.csv(f, stringsAsFactors = FALSE)
  p <- ggplot(agg, aes(x = factor(omega), y = mean_pct_largest,
                       color = method, group = method)) +
    geom_line(linewidth = 0.8) + geom_point(size = 2.2) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
    facet_wrap(~dgp, ncol = 3, labeller = labeller(dgp = DGP_LABS)) +
    scale_color_manual(values = c("#E69F00", "#0072B2"), name = NULL) +
    labs(x = expression(omega),
         y = "Mean share of nodes in largest community",
         title = "Multislice partition coarseness across omega",
         subtitle = "High values indicate the giant-component pathology") +
    theme_paper + guides(color = guide_legend(ncol = 2))
  ggsave(file.path(outdir, "fig_multislice_omega_sweep.pdf"), p,
         width = 7.5, height = 4.2)
  msg("[5] omega sweep figure done")
} else msg("[5] SKIP: %s not found", f)

# =============================================================================
# 6. Ideal-point external validity (source: idealpoint_alignment.csv)
# =============================================================================
f <- file.path(outdir, "idealpoint_alignment.csv")
if (file.exists(f)) {
  df <- read.csv(f, stringsAsFactors = FALSE)
  p <- ggplot(df, aes(x = year, y = ari, color = method)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("DynMux (Jaccard)" = "#E69F00",
                                  "multinet GLouvain" = "grey55"),
                       name = NULL) +
    labs(x = NULL, y = "Adjusted Rand Index vs UN-voting clusters",
         title = "External validity: alignment with UN-voting ideal-point clusters",
         subtitle = "Higher = alliance communities more closely match voting blocs") +
    theme_paper + guides(color = guide_legend(ncol = 2))
  ggsave(file.path(outdir, "fig_idealpoint_alignment.pdf"), p,
         width = 7.5, height = 5)
  msg("[6] idealpoint figure done")
} else msg("[6] SKIP: %s not found", f)

# =============================================================================
# 7. Pareto frontier (NEW figure; source: tab_dgp_misspec.csv)
# =============================================================================
f <- file.path(outdir, "tab_dgp_misspec.csv")
if (file.exists(f)) {
  agg <- read.csv(f, stringsAsFactors = FALSE)
  pf <- agg %>%
    group_by(method, group) %>%
    summarise(nmi = mean(nmi_mean, na.rm = TRUE),
              tc_abs_gap = mean(abs(tc_gap), na.rm = TRUE), .groups = "drop")
  pf$pareto <- vapply(seq_len(nrow(pf)), function(i) {
    !any(pf$tc_abs_gap < pf$tc_abs_gap[i] & pf$nmi > pf$nmi[i])
  }, logical(1))
  frontier <- pf %>% filter(pareto) %>% arrange(tc_abs_gap)
  grp_cols <- c("Temporally aware" = "#E69F00",
                "Temporally agnostic" = "#56B4E9",
                "Generative (SBM)" = "#009E73",
                "Other" = "#999999")
  wrap_method <- function(s) {
    s <- sub(" \\(", "\n(", s)
    vapply(strsplit(s, "\n"), function(parts) {
      paste(unlist(lapply(parts, function(p)
        paste(strwrap(p, width = 16), collapse = "\n"))), collapse = "\n")
    }, character(1))
  }
  pf$method_lab <- wrap_method(pf$method)
  label_layer <- if (requireNamespace("ggrepel", quietly = TRUE)) {
    ggrepel::geom_text_repel(aes(label = method_lab), size = 2.3,
                             lineheight = 0.85,
                             max.overlaps = 40, seed = 123,
                             box.padding = 0.45, point.padding = 0.35,
                             force = 2.5, force_pull = 0.4,
                             min.segment.length = 0.15,
                             segment.color = "grey70", segment.size = 0.25)
  } else {
    geom_text(aes(label = method_lab), vjust = -0.9, size = 2.3,
              lineheight = 0.85, check_overlap = TRUE)
  }
  p_pareto <- ggplot(pf, aes(x = tc_abs_gap, y = nmi)) +
    geom_step(data = frontier, direction = "hv", linetype = "dashed",
              color = "grey40") +
    geom_point(aes(color = group, shape = pareto), size = 2.8) +
    label_layer +
    scale_color_manual(values = grp_cols, name = NULL) +
    scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 16),
                       labels = c(`TRUE` = "On frontier",
                                  `FALSE` = "Dominated"), name = NULL) +
    labs(x = "Mean absolute temporal-coherence gap (0 = matches truth)",
         y = "Mean NMI vs ground truth",
         title = "Accuracy vs temporal fidelity across methods",
         subtitle = "Upper-left is better. Dashed line traces the Pareto frontier.") +
    theme_paper +
    guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2))
  ggsave(file.path(outdir, "fig_pareto_nmi_tc.pdf"), p_pareto,
         width = 7.5, height = 6.5)
  msg("[7] pareto figure done (NEW - review before use)")
} else msg("[7] SKIP: %s not found", f)

# =============================================================================
# 8. Link prediction (source: linkpred_results.csv)
# =============================================================================
f <- file.path(outdir, "linkpred_results.csv")
if (file.exists(f)) {
  df <- read.csv(f, stringsAsFactors = FALSE)
  N_FOLDS <- length(unique(df$fold))

  overall <- df %>% filter(!is.na(auc_overall)) %>%
    group_by(network, method) %>%
    summarise(auc_mean = mean(auc_overall),
              auc_se = sd(auc_overall) / sqrt(n()), .groups = "drop") %>%
    as.data.frame()
  overall <- classify_methods(overall)
  p_auc <- ggplot(overall, aes(x = method, y = auc_mean, fill = sim)) +
    geom_col(color = "white", linewidth = 0.2) +
    geom_errorbar(aes(ymin = auc_mean - auc_se, ymax = auc_mean + auc_se),
                  width = 0.25, linewidth = 0.3) +
    geom_hline(yintercept = 0.5, linetype = "dashed") +
    facet_wrap(~network, ncol = 3) +
    scale_fill_manual(values = SIM_COLS, name = NULL, drop = FALSE) +
    labs(x = NULL, y = "Held-out AUC",
         title = "Out-of-sample link prediction (10% edge hold-out)",
         subtitle = sprintf("Mean +/- SE over %d folds", N_FOLDS)) +
    theme_paper +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                     size = 7)) +
    guides(fill = guide_legend(ncol = 2))
  ggsave(file.path(outdir, "fig_linkpred_auc.pdf"), p_auc,
         width = 7.5, height = 6)

  per_layer <- df %>% filter(!is.na(auc_layer)) %>%
    group_by(network, method, layer) %>%
    summarise(auc_mean = mean(auc_layer), .groups = "drop") %>%
    as.data.frame()
  p_layers <- dynmux_lineplot(per_layer, "layer", "auc_mean",
                              "Layer (year index)", "Held-out AUC",
                              "Link-prediction AUC over time") +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
    facet_wrap(~network, ncol = 3, scales = "free_x")
  ggsave(file.path(outdir, "fig_linkpred_layers.pdf"), p_layers,
         width = 8, height = 5)

  msg("[8] link prediction figures done")
} else msg("[8] SKIP: %s not found", f)

# =============================================================================
# 9. Empirical application figures (ROAR only)
# =============================================================================
f_sum <- file.path(outdir, "empirical_community_summary.csv")
if (file.exists(f_sum)) {
  comm_summary <- read.csv(f_sum, stringsAsFactors = FALSE)
  comm_summary$fit_lab <- ifelse(comm_summary$fit_type %in% names(FIT_LABS),
                                 FIT_LABS[comm_summary$fit_type],
                                 de_underscore(comm_summary$fit_type))
  net_cols <- c("#E69F00", "#0072B2", "#009E73")
  jaccard_summary <- comm_summary[comm_summary$fit_type == "jaccard", ]

  p_count <- ggplot(jaccard_summary, aes(x = year, y = n_communities,
                                         color = network, linetype = network)) +
    geom_line(linewidth = 0.9) + geom_point(size = 1.5) +
    scale_color_manual(values = net_cols, name = NULL) +
    scale_linetype_discrete(name = NULL) +
    labs(x = "Year", y = "Number of communities",
         title = "Community count over time (Jaccard specification)") +
    theme_paper + guides_2col
  ggsave(file.path(outdir, "fig_community_count_over_time.pdf"),
         p_count, width = 7.5, height = 4.5)

  mod_data <- jaccard_summary[!is.na(jaccard_summary$modularity), ]
  if (nrow(mod_data) > 0) {
    p_mod <- ggplot(mod_data, aes(x = year, y = modularity,
                                  color = network, linetype = network)) +
      geom_line(linewidth = 0.9) + geom_point(size = 1.5) +
      scale_color_manual(values = net_cols, name = NULL) +
      scale_linetype_discrete(name = NULL) +
      labs(x = "Year", y = "Modularity") +
      theme_paper + guides_2col
    ggsave(file.path(outdir, "fig_modularity_over_time.pdf"),
           p_mod, width = 7.5, height = 4.5)
  }

  p_count_robust <- ggplot(comm_summary,
                           aes(x = year, y = n_communities,
                               color = fit_lab, linetype = fit_lab)) +
    geom_line(linewidth = 0.6, alpha = 0.85) +
    facet_wrap(~network) +
    scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73",
                                  "#D55E00", "#CC79A7", "grey55"),
                       name = NULL) +
    scale_linetype_discrete(name = NULL) +
    labs(x = "Year", y = "Number of communities",
         title = "Robustness: community count across interlayer weighting specifications") +
    theme_paper + guides_2col
  ggsave(file.path(outdir, "fig_community_count_robustness.pdf"),
         p_count_robust, width = 7.5, height = 4.5)
  msg("[9a] empirical summary figures done")
} else msg("[9a] SKIP (ROAR only): %s not found", f_sum)

f_agree <- file.path(outdir, "empirical_robustness_agreement.csv")
if (file.exists(f_agree)) {
  agreement_df <- read.csv(f_agree, stringsAsFactors = FALSE)
  agg_agree <- aggregate(nmi ~ network + method1 + method2,
                         data = agreement_df, FUN = mean)
  agg_agree$method1 <- ifelse(agg_agree$method1 %in% names(FIT_LABS),
                              FIT_LABS[agg_agree$method1],
                              de_underscore(agg_agree$method1))
  agg_agree$method2 <- ifelse(agg_agree$method2 %in% names(FIT_LABS),
                              FIT_LABS[agg_agree$method2],
                              de_underscore(agg_agree$method2))
  p_robust <- ggplot(agg_agree, aes(x = method1, y = method2, fill = nmi)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(nmi, 2)), size = 2.5) +
    scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 1)) +
    facet_wrap(~network) +
    labs(x = NULL, y = NULL, fill = "Mean NMI",
         title = "Cross-specification agreement (mean NMI across years)") +
    theme_paper +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                     size = 7),
          axis.text.y = element_text(size = 7))
  ggsave(file.path(outdir, "fig_robustness_agreement.pdf"),
         p_robust, width = 7.5, height = 3.6)
  msg("[9b] robustness agreement heatmap done")
} else msg("[9b] SKIP (ROAR only): %s not found", f_agree)

f_persist <- file.path(outdir, "empirical_persistence.csv")
if (file.exists(f_persist)) {
  persistence_df <- read.csv(f_persist, stringsAsFactors = FALSE)
  persist_jaccard <- persistence_df[persistence_df$fit_type == "jaccard", ]
  p_persist <- ggplot(persist_jaccard,
                      aes(x = year_from, y = persistence_nmi,
                          color = network, linetype = network)) +
    geom_line(linewidth = 0.9) + geom_point(size = 1.5) +
    scale_color_manual(values = c("#E69F00", "#0072B2", "#009E73"),
                       name = NULL) +
    scale_linetype_discrete(name = NULL) +
    labs(x = "Year", y = "Persistence (NMI with next year)",
         title = "Community persistence: year-to-year stability") +
    ylim(0, 1) + theme_paper + guides_2col
  ggsave(file.path(outdir, "fig_persistence_over_time.pdf"),
         p_persist, width = 7.5, height = 4.5)
  msg("[9c] persistence figure done")
} else msg("[9c] SKIP (ROAR only): %s not found", f_persist)

for (spec in list(
  list(csv = "empirical_bootstrap_modularity_ci.csv", y = "estimate",
       ylab = "Modularity", out = "fig_empirical_bootstrap_modularity.pdf",
       title = "Bootstrap 95% CI for modularity over time", ribbon = TRUE),
  list(csv = "empirical_bootstrap_count_ci.csv", y = "estimate",
       ylab = "Number of communities", out = "fig_empirical_bootstrap_count.pdf",
       title = "Bootstrap 95% CI for community count over time", ribbon = TRUE),
  list(csv = "empirical_bootstrap_stability.csv", y = "mean_stability",
       ylab = "Mean node stability", out = "fig_empirical_bootstrap_stability.pdf",
       title = "Bootstrap node stability over time", ribbon = FALSE)
)) {
  f_b <- file.path(outdir, spec$csv)
  if (!file.exists(f_b)) { msg("[9d] SKIP (ROAR only): %s not found", f_b); next }
  b_df <- read.csv(f_b, stringsAsFactors = FALSE)
  net_cols <- c("#E69F00", "#0072B2", "#009E73")
  p_b <- ggplot(b_df, aes(x = year, y = .data[[spec$y]],
                          color = network, fill = network)) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = net_cols, name = NULL) +
    scale_fill_manual(values = net_cols, name = NULL) +
    labs(x = "Year", y = spec$ylab, title = spec$title) +
    theme_paper + guides(color = guide_legend(ncol = 2),
                         fill = guide_legend(ncol = 2))
  if (spec$ribbon) {
    p_b <- p_b + geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2,
                             color = NA)
  } else {
    p_b <- p_b + geom_point(size = 1.5) + ylim(0, 1)
  }
  ggsave(file.path(outdir, spec$out), p_b, width = 7.5, height = 4.5)
  msg("[9d] %s done", spec$out)
}

# =============================================================================
# 10. Co-assignment figures from SAVED bootstrap array (ROAR only).
#     Requires coassign_prob.rds written by the updated
#     08_coassignment_bootstrap.R (saveRDS of prob + metadata).
# =============================================================================
f_prob <- file.path(outdir, "coassign_prob.rds")
if (file.exists(f_prob)) {
  pr <- readRDS(f_prob)   # list(prob = array [layers x n x n], years, states, B, network)
  HIGHLIGHT_YEARS <- c(1985, 1989, 1991, 1995)
  NAMED_STATES <- c("2" = "USA", "200" = "UK", "220" = "France",
                    "255" = "Germany (West / unified)", "265" = "East Germany",
                    "290" = "Poland", "345" = "Yugoslavia",
                    "365" = "USSR / Russia", "710" = "China", "750" = "India",
                    "40" = "Cuba")

  hist_df <- list()
  for (t in seq_along(pr$years)) {
    P <- pr$prob[t, , ]
    hist_df[[t]] <- data.frame(year = pr$years[t], p_same = P[upper.tri(P)])
  }
  hist_df <- do.call(rbind, hist_df)
  p_hist <- ggplot(hist_df %>% filter(year %in% HIGHLIGHT_YEARS),
                   aes(x = p_same)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "white") +
    facet_wrap(~ year, ncol = 2) +
    scale_y_continuous(labels = label_comma()) +
    labs(x = "Pr(same community | bootstrap)", y = "Number of state pairs",
         title = "Distribution of pairwise co-membership probabilities",
         subtitle = sprintf("%s, B = %d bootstrap replicates",
                            pr$network, pr$B)) +
    theme_paper
  ggsave(file.path(outdir, "fig_coassign_uncertainty.pdf"), p_hist,
         width = 7.5, height = 6)

  all_states <- as.integer(pr$states)
  ix_named <- match(as.integer(names(NAMED_STATES)), all_states)
  keep_named <- !is.na(ix_named)
  ix_named <- ix_named[keep_named]
  labs_named <- unname(NAMED_STATES[keep_named])
  named_long <- list(); rx <- 0L
  for (y in HIGHLIGHT_YEARS) {
    t <- match(y, pr$years)
    if (is.na(t)) next
    P <- pr$prob[t, ix_named, ix_named]
    rx <- rx + 1L
    named_long[[rx]] <- data.frame(
      year = y,
      a = rep(labs_named, each = length(labs_named)),
      b = rep(labs_named, times = length(labs_named)),
      p = as.vector(P), stringsAsFactors = FALSE)
  }
  named_df <- do.call(rbind, named_long)
  wrap10 <- function(s) vapply(s, function(x)
    paste(strwrap(x, width = 10), collapse = "\n"), character(1))
  named_df$a <- factor(wrap10(named_df$a), levels = wrap10(labs_named))
  named_df$b <- factor(wrap10(named_df$b), levels = wrap10(labs_named))
  p_named <- ggplot(named_df, aes(x = a, y = b, fill = p)) +
    geom_tile(color = "white") +
    facet_wrap(~ year, ncol = 2) +
    scale_fill_viridis_c(limits = c(0, 1), option = "magma", direction = -1) +
    labs(x = NULL, y = NULL, fill = "Pr(same\ncommunity)",
         title = "Bootstrap co-membership probability among named states",
         subtitle = sprintf("%s, B = %d", pr$network, pr$B)) +
    theme_paper +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 6.5))
  ggsave(file.path(outdir, "fig_coassign_named.pdf"), p_named,
         width = 7.5, height = 8)
  msg("[10] coassign figures done (2x2 facets)")
} else msg("[10] SKIP: %s not found (run updated 08_coassignment_bootstrap.R)", f_prob)

msg("\nreplot_from_cache v2 complete. Figures written to %s", outdir)
