#!/usr/bin/env Rscript
# =============================================================================
# replication/post/11_appendix_regimes.R
# Appendix material from the regime comparison (sim/01):
#
#   fig_regime_paired_nmi.pdf   paired Delta joint NMI, 4 regimes x 2 intensities
#   fig_regime_paired_kmae.pdf  paired Delta K MAE,     4 regimes x 2 intensities
#   tab_app_paired_ci_nmi.tex   paired 95% CIs, joint NMI (booktabs tabular)
#   tab_app_paired_ci_kmae.tex  paired 95% CIs, K MAE    (booktabs tabular)
#   tab_paired_wilcoxon.tex     paired Wilcoxon signed-rank tests, by regime
#
#   Delta Joint NMI = DynMux - baseline   (positive favors DynMux)
#   Delta K MAE     = baseline - DynMux   (positive favors DynMux; lower is better)
#
# Pairing unit = one simulated network (config file x rep); every method saw
# the same network, so the difference is within-network.
#
# Formerly paper_scripts/15b_sim_paired_intensity.R (figures + CSV),
# 33_paired_ci_table.R (longtable, now two tabulars) and
# 21_paired_wilcoxon.R (tests). The intermediate CSV that 33 read from 15b is
# gone: everything is computed once here from the raw dyn_cfg files.
#
# Tables are bare booktabs tabulars, no caption or label; the manuscript wraps
# each in its own float. Requires \usepackage{booktabs}.
#
# Usage:  DM_ROOT=/path/to/dynamic_multiplex Rscript replication/post/11_appendix_regimes.R
# =============================================================================
set.seed(123)
suppressMessages({ library(ggplot2) })

ROOT <- Sys.getenv("DM_ROOT", unset = getwd())
REG  <- file.path(ROOT, "output", "regime")
TAB  <- file.path(ROOT, "manuscript", "tables")
FIG  <- file.path(ROOT, "manuscript", "figures")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

write_tex <- function(header, body, align, path) {
  stopifnot(length(body) > 0, nchar(align) > 0)
  writeLines(c(sprintf("\\begin{tabular}{%s}", align),
               "\\toprule", header, "\\midrule",
               body,
               "\\bottomrule", "\\end{tabular}"), path)
  cat("wrote", basename(path), "(", length(body), "rows )\n")
}

# --- read ----------------------------------------------------------------
fs <- list.files(REG, "^dyn_cfg.*csv$", full.names = TRUE)
if (!length(fs)) stop("No regime output in ", REG, " -- run sim/01 first.", call. = FALSE)
d <- do.call(rbind, lapply(fs, function(f) { x <- read.csv(f); x$cfg <- basename(f); x }))
stopifnot(all(c("regime", "intensity", "method", "rep", "nmi_joint", "k_mae") %in% names(d)))
cat("files:", length(fs), " rows:", nrow(d), "\n")

relab <- c("DynMux multislice (adjacent)" = "Multislice (adjacent)",
           "DynMux multislice (custom)"   = "Multislice (same links)",
           "Cross-sectional + Hungarian"  = "Hungarian matching")
d$method <- ifelse(d$method %in% names(relab), relab[d$method], d$method)
d$unit   <- paste(d$cfg, d$rep, sep = "_")

rn <- c(birthdeath  = "Births & Deaths",
        churnswitch = "Gradual Switching",
        regimeshift = "Abrupt Rewiring",
        seasonality = "Recurring Structure")
rn_tex <- c(birthdeath  = "Births \\& deaths",
            churnswitch = "Gradual switching",
            regimeshift = "Abrupt rewiring",
            seasonality = "Recurring structure")
refs      <- c("DynMux Jaccard", "DynMux Overlap")
baselines <- c("Hungarian matching", "Multislice (adjacent)", "Multislice (same links)",
               "Dynamic SBM", "Pooled Leiden", "multinet GLouvain")
blab <- c("Hungarian matching"    = "Hungarian\nmatching",
          "Multislice (adjacent)" = "Multislice\nadjacent",
          "Multislice (same links)" = "Multislice\nsame links",
          "Dynamic SBM"           = "Dynamic\nSBM",
          "Pooled Leiden"         = "Pooled\nLeiden",
          "multinet GLouvain"     = "Multislice full\n(multinet)")
bl_tex <- c("Hungarian matching"    = "Hungarian matching",
            "Multislice (adjacent)" = "Multislice (adjacent links)",
            "Multislice (same links)" = "Multislice (same links as DynMux)",
            "Dynamic SBM"           = "Dynamic SBM",
            "Pooled Leiden"         = "Pooled Leiden",
            "multinet GLouvain"     = "Multislice full (\\texttt{multinet})")
stopifnot(all(c(refs, baselines) %in% unique(d$method)))

# --- paired differences within regime x intensity (figures + CI tables) ---
rows <- list()
for (rg in names(rn)) for (it in c("low", "high")) for (mt in c("nmi_joint", "k_mae")) {
  sub <- d[d$regime == rg & d$intensity == it, c("unit", "method", mt)]
  stopifnot(nrow(sub) > 0)
  w <- reshape(sub, idvar = "unit", timevar = "method", direction = "wide")
  names(w) <- sub(paste0(mt, "."), "", names(w), fixed = TRUE)
  for (rf in refs) for (bl in baselines) {
    stopifnot(all(c(rf, bl) %in% names(w)))
    dd <- if (mt == "nmi_joint") w[[rf]] - w[[bl]] else w[[bl]] - w[[rf]]
    dd <- dd[is.finite(dd)]
    stopifnot(length(dd) > 1)
    n <- length(dd); se <- sd(dd) / sqrt(n)
    rows[[length(rows) + 1]] <- data.frame(
      regime = rg, intensity = it, metric = mt, ref = rf, baseline = bl,
      n = n, diff = mean(dd), lo = mean(dd) - 1.96 * se, hi = mean(dd) + 1.96 * se,
      stringsAsFactors = FALSE)
  }
}
r <- do.call(rbind, rows)
stopifnot(nrow(r) == 160)                        # 4 x 2 x 2 x 2 x 5
cat("paired comparisons:", nrow(r), " n per comparison:",
    paste(unique(r$n), collapse = ","), "\n")
cat("CIs excluding zero:", sum(r$lo > 0 | r$hi < 0), "of", nrow(r), "\n")

# --- figures: one per metric; fixed x-scale across all panels ------------
r$regime_lab    <- factor(unname(rn[r$regime]), levels = unname(rn))
r$intensity_lab <- factor(ifelse(r$intensity == "low", "Low intensity", "High intensity"),
                          levels = c("Low intensity", "High intensity"))
r$baseline_lab  <- factor(unname(blab[r$baseline]), levels = rev(unname(blab[baselines])))
r$ref_short     <- sub("DynMux ", "", r$ref)

mk <- function(mt, xlab) {
  s <- r[r$metric == mt, ]
  ggplot(s, aes(diff, baseline_lab, colour = ref_short, shape = ref_short)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
    geom_errorbar(aes(xmin = lo, xmax = hi),
                  position = position_dodge(width = 0.55), width = 0.25) +
    geom_point(position = position_dodge(width = 0.55), size = 1.8) +
    facet_grid(regime_lab ~ intensity_lab) +
    scale_colour_manual(values = c(Jaccard = "#1b9e77", Overlap = "#d95f02")) +
    scale_shape_manual(values  = c(Jaccard = 16, Overlap = 17)) +
    labs(x = xlab, y = NULL, colour = "DynMux coupling", shape = "DynMux coupling") +
    theme_bw(base_size = 9) +
    theme(legend.position  = "bottom",
          panel.grid.minor = element_blank(),
          strip.text       = element_text(face = "bold", size = 8.5),
          axis.text.y      = element_text(size = 8))
}
p1 <- mk("nmi_joint", expression(paste(Delta, " Joint NMI (positive favors DynMux), 95% CI")))
p2 <- mk("k_mae",     expression(paste(Delta, " ", italic(K), " MAE (positive favors DynMux), 95% CI")))
ggsave(file.path(FIG, "fig_regime_paired_nmi.pdf"),  p1, width = 6.5, height = 8.0)
ggsave(file.path(FIG, "fig_regime_paired_nmi.png"),  p1, width = 6.5, height = 8.0, dpi = 300)
ggsave(file.path(FIG, "fig_regime_paired_kmae.pdf"), p2, width = 6.5, height = 8.0)
ggsave(file.path(FIG, "fig_regime_paired_kmae.png"), p2, width = 6.5, height = 8.0, dpi = 300)
cat("wrote fig_regime_paired_nmi / _kmae (.pdf/.png)\n")

# --- CI tables: one tabular per metric, Jaccard and Overlap side by side --
fmt3 <- function(x) sub("-", "$-$", sprintf("%.3f", round(x, 3) + 0), fixed = TRUE)  # +0 kills "-0.000"
r$cell <- sprintf("%s [%s, %s]", fmt3(r$diff), fmt3(r$lo), fmt3(r$hi))
for (mt in c("nmi_joint", "k_mae")) {
  s <- r[r$metric == mt, ]
  w <- reshape(s[, c("regime", "intensity", "baseline", "ref_short", "cell")],
               idvar = c("regime", "intensity", "baseline"),
               timevar = "ref_short", direction = "wide")
  names(w) <- sub("cell.", "", names(w), fixed = TRUE)
  stopifnot(all(c("Jaccard", "Overlap") %in% names(w)), nrow(w) == 8 * length(baselines))
  w <- w[order(factor(w$regime,    levels = names(rn)),
               factor(w$intensity, levels = c("low", "high")),
               factor(w$baseline,  levels = baselines)), ]
  grp  <- sprintf("%s (%s)", unname(rn_tex[w$regime]), ifelse(w$intensity == "low", "Low", "High"))
  show <- c(TRUE, grp[-1] != grp[-length(grp)])          # blank repeated group labels
  body <- sprintf("%s & %s & %s & %s \\\\",
                  ifelse(show, grp, ""), unname(bl_tex[w$baseline]), w$Jaccard, w$Overlap)
  write_tex(header = "Regime (intensity) & Baseline & Jaccard & Overlap \\\\",
            body   = body, align = "llcc",
            path   = file.path(TAB, sprintf("tab_app_paired_ci_%s.tex",
                                            if (mt == "nmi_joint") "nmi" else "kmae")))
}

# --- paired Wilcoxon signed-rank tests, pooled over intensity ------------
# One-sample signed-rank test on the within-network differences per regime
# (n = 540 networks per regime: 18 configs x 30 reps).
res <- list()
for (rg in names(rn)) for (mt in c("nmi_joint", "k_mae")) {
  sub <- d[d$regime == rg, c("unit", "method", mt)]
  w <- reshape(sub, idvar = "unit", timevar = "method", direction = "wide")
  names(w) <- sub(paste0(mt, "."), "", names(w), fixed = TRUE)
  for (rf in refs) for (bl in baselines) {
    dd <- if (mt == "nmi_joint") w[[rf]] - w[[bl]] else w[[bl]] - w[[rf]]
    dd <- dd[is.finite(dd)]
    wt <- suppressWarnings(wilcox.test(dd, mu = 0, exact = FALSE))
    res[[length(res) + 1]] <- data.frame(regime = rg, metric = mt, ref = rf,
      baseline = bl, n = length(dd), mean_diff = mean(dd),
      V = unname(wt$statistic), p = wt$p.value, stringsAsFactors = FALSE)
  }
}
wx <- do.call(rbind, res)
stopifnot(nrow(wx) == 80)                        # 4 x 2 x 2 x 5
wx$p_lab      <- ifelse(wx$p < 1e-16, "$<10^{-16}$", sprintf("%.2g", wx$p))
wx$metric_lab <- ifelse(wx$metric == "nmi_joint", "Joint NMI", "$K$ MAE")
write_tex(
  header = "Regime & Metric & DynMux & Baseline & $n$ & Mean diff. & $p$ \\\\",
  body   = sprintf("%s & %s & %s & %s & %d & %s & %s \\\\",
                   unname(rn_tex[wx$regime]), wx$metric_lab, sub("DynMux ", "", wx$ref),
                   unname(bl_tex[wx$baseline]), wx$n,
                   fmt3(wx$mean_diff), wx$p_lab),
  align  = "llllrrl",
  path   = file.path(TAB, "tab_paired_wilcoxon.tex"))
cat("max Wilcoxon p:", format(max(wx$p), digits = 3), "\n")

cat("done 11_appendix_regimes\n")
