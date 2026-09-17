#!/usr/bin/env Rscript
# 15b_sim_paired_intensity.R -- Paired-difference figures split by regime AND
# change intensity, one figure per metric so that every panel within a figure
# shares a fixed x-scale (Political Analysis requires identical panel scales).
#
#   Delta Joint NMI = DynMux - baseline   (positive favors DynMux)
#   Delta K MAE     = baseline - DynMux   (positive favors DynMux; lower is better)
#
# facet_wrap: 4 regimes x 2 intensities = 8 panels, 2 columns.
# Produces: manuscript/figures/fig_regime_paired_nmi.pdf/.png
#           manuscript/figures/fig_regime_paired_kmae.pdf/.png
#           manuscript/tables/tab_regime_paired_intensity.csv
set.seed(123)
suppressMessages({ library(ggplot2) })

# --- paths -------------------------------------------------------------
# DM_ROOT is the project root; defaults to the working directory.
ROOT <- Sys.getenv("DM_ROOT", unset = getwd())
DYN  <- file.path(ROOT, "replication", "extended", "output", "dynamic")
TAB  <- file.path(ROOT, "manuscript", "tables")
FIG  <- file.path(ROOT, "manuscript", "figures")
if (!dir.exists(DYN)) {
  stop("Simulation output not found at: ", DYN,
       "\n  Set DM_ROOT to the project root, or run from that directory.",
       call. = FALSE)
}
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

fs <- list.files(DYN, "dyn_cfg.*csv$", full.names = TRUE)
stopifnot(length(fs) > 0)
d <- do.call(rbind, lapply(fs, function(f) { x <- read.csv(f); x$cfg <- basename(f); x }))
stopifnot(all(c("regime", "intensity", "method", "rep",
                "nmi_joint", "k_mae") %in% names(d)))

relab <- c("DynMux multislice (adjacent)" = "Multislice (adjacent)",
           "Cross-sectional + Hungarian"  = "Hungarian matching")
d$method <- ifelse(d$method %in% names(relab), relab[d$method], d$method)
d$unit <- paste(d$cfg, d$rep, sep = "_")

# Regime labels match the main text.
rn <- c(birthdeath  = "Births & Deaths",
        churnswitch = "Gradual Switching",
        regimeshift = "Abrupt Rewiring",
        seasonality = "Recurring Structure")
refs      <- c("DynMux Jaccard", "DynMux Overlap")
baselines <- c("Hungarian matching", "Multislice (adjacent)",
               "Dynamic SBM", "Pooled Leiden", "multinet GLouvain")
blab <- c("Hungarian matching"    = "Hungarian\nmatching",
          "Multislice (adjacent)" = "Multislice\nadjacent",
          "Dynamic SBM"           = "Dynamic\nSBM",
          "Pooled Leiden"         = "Pooled\nLeiden",
          "multinet GLouvain"     = "Multislice full\n(multinet)")

# --- paired differences within regime x intensity ----------------------
rows <- list()
for (rg in unique(d$regime)) for (it in c("low", "high")) {
  for (mt in c("nmi_joint", "k_mae")) {
    sub <- d[d$regime == rg & d$intensity == it, c("unit", "method", mt)]
    if (!nrow(sub)) next
    w <- reshape(sub, idvar = "unit", timevar = "method", direction = "wide")
    names(w) <- sub(paste0(mt, "."), "", names(w), fixed = TRUE)
    for (rf in refs) for (bl in baselines) {
      if (!all(c(rf, bl) %in% names(w))) next
      dd <- if (mt == "nmi_joint") w[[rf]] - w[[bl]] else w[[bl]] - w[[rf]]
      dd <- dd[is.finite(dd)]
      if (!length(dd)) next
      n <- length(dd); se <- sd(dd) / sqrt(n)
      rows[[length(rows) + 1]] <- data.frame(
        regime = rg, intensity = it, metric = mt, ref = rf, baseline = bl,
        n = n, diff = mean(dd), lo = mean(dd) - 1.96 * se,
        hi = mean(dd) + 1.96 * se)
    }
  }
}
r <- do.call(rbind, rows)
stopifnot(nrow(r) > 0)
write.csv(r, file.path(TAB, "tab_regime_paired_intensity.csv"), row.names = FALSE)

r$regime_lab    <- factor(unname(rn[r$regime]), levels = unname(rn))
r$intensity_lab <- factor(ifelse(r$intensity == "low", "Low intensity",
                                 "High intensity"),
                          levels = c("Low intensity", "High intensity"))
r$baseline_lab  <- factor(unname(blab[r$baseline]),
                          levels = rev(unname(blab[baselines])))
r$ref <- sub("DynMux ", "", r$ref)

# --- one figure per metric; fixed x-scale across all panels ------------
mk <- function(mt, xlab) {
  s <- r[r$metric == mt, ]
  ggplot(s, aes(diff, baseline_lab, colour = ref, shape = ref)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
    geom_errorbarh(aes(xmin = lo, xmax = hi),
                   position = position_dodge(width = 0.55), height = 0.25) +
    geom_point(position = position_dodge(width = 0.55), size = 1.8) +
    facet_grid(regime_lab ~ intensity_lab) +               # fixed scales
    scale_colour_manual(values = c(Jaccard = "#1b9e77", Overlap = "#d95f02")) +
    scale_shape_manual(values  = c(Jaccard = 16, Overlap = 17)) +
    labs(x = xlab, y = NULL, colour = "DynMux coupling",
         shape = "DynMux coupling") +
    theme_bw(base_size = 9) +
    theme(legend.position  = "bottom",
          panel.grid.minor = element_blank(),
          strip.text       = element_text(face = "bold", size = 8.5),
          axis.text.y      = element_text(size = 8))
}

p1 <- mk("nmi_joint",
         expression(paste(Delta, " Joint NMI (positive favors DynMux), 95% CI")))
p2 <- mk("k_mae",
         expression(paste(Delta, " ", italic(K), " MAE (positive favors DynMux), 95% CI")))

ggsave(file.path(FIG, "fig_regime_paired_nmi.pdf"),  p1, width = 6.5, height = 8.0)
ggsave(file.path(FIG, "fig_regime_paired_nmi.png"),  p1, width = 6.5, height = 8.0, dpi = 300)
ggsave(file.path(FIG, "fig_regime_paired_kmae.pdf"), p2, width = 6.5, height = 8.0)
ggsave(file.path(FIG, "fig_regime_paired_kmae.png"), p2, width = 6.5, height = 8.0, dpi = 300)

cat("done 15b: 8 panels per figure (4 regimes x 2 intensities), 2 figures\n")
cat("rows:", nrow(r), " comparisons per panel:",
    nrow(r) / length(unique(paste(r$regime, r$intensity, r$metric))), "\n")
