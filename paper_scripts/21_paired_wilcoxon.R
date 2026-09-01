#!/usr/bin/env Rscript
# 21_paired_wilcoxon.R -- Paired Wilcoxon signed-rank tests for every
# DynMux-vs-baseline comparison shown in fig_regime_paired (both metrics).
# Delta Joint NMI = DynMux - baseline; Delta K MAE = baseline - DynMux,
# so positive always favors DynMux. One-sample signed-rank test on the
# within-network differences (n = 540 networks per regime).
# Outputs: manuscript/tables/tab_paired_wilcoxon.csv/.tex (appendix)
set.seed(123)
D <- "replication/extended/output/dynamic"
d <- do.call(rbind, lapply(list.files(D, "dyn_cfg.*csv$", full.names = TRUE),
  function(f) { x <- read.csv(f); x$cfg <- basename(f); x }))
relab <- c("DynMux multislice (adjacent)" = "Multislice adjacent",
           "Cross-sectional + Hungarian"  = "Hungarian matching",
           "multinet GLouvain"            = "Multislice full (multinet)")
d$method <- ifelse(d$method %in% names(relab), relab[d$method], d$method)
d$unit <- paste(d$cfg, d$rep, sep = "_")
rn <- c(birthdeath = "Births \\& Deaths", churnswitch = "Churn \\& Switch",
        regimeshift = "Regime Shift", seasonality = "Seasonality")
refs <- c("DynMux Jaccard", "DynMux Overlap")
baselines <- c("Hungarian matching", "Multislice adjacent", "Dynamic SBM",
               "Pooled Leiden", "Multislice full (multinet)")
res <- list()
for (rg in unique(d$regime)) for (mt in c("nmi_joint", "k_mae")) {
  sub <- d[d$regime == rg, c("unit", "method", mt)]
  w <- reshape(sub, idvar = "unit", timevar = "method", direction = "wide")
  names(w) <- sub(paste0(mt, "."), "", names(w), fixed = TRUE)
  for (rf in refs) for (bl in baselines) {
    if (mt == "nmi_joint") dd <- w[[rf]] - w[[bl]] else dd <- w[[bl]] - w[[rf]]
    dd <- dd[is.finite(dd)]
    wt <- suppressWarnings(wilcox.test(dd, mu = 0, exact = FALSE))
    res[[length(res) + 1]] <- data.frame(regime = rg, metric = mt, ref = rf,
      baseline = bl, n = length(dd), mean_diff = mean(dd),
      V = unname(wt$statistic), p = wt$p.value)
  }
}
r <- do.call(rbind, res)
write.csv(r, "manuscript/tables/tab_paired_wilcoxon.csv", row.names = FALSE)
r$p_lab <- ifelse(r$p < 1e-16, "$<10^{-16}$", sprintf("%.2g", r$p))
r$metric_lab <- ifelse(r$metric == "nmi_joint", "Joint NMI", "$K$ MAE")
tex <- c("\\begin{tabular}{llllrrl}", "\\toprule",
  "Regime & Metric & DynMux & Baseline & $n$ & Mean diff. & $p$ \\\\",
  "\\midrule",
  sprintf("%s & %s & %s & %s & %d & %.3f & %s \\\\",
          unname(rn[r$regime]), r$metric_lab, sub("DynMux ", "", r$ref),
          r$baseline, r$n, r$mean_diff, r$p_lab),
  "\\bottomrule", "\\end{tabular}")
writeLines(tex, "manuscript/tables/tab_paired_wilcoxon.tex")
cat("done 21; comparisons:", nrow(r), "; max p =", format(max(r$p), digits = 3), "\n")
