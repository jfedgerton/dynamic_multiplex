#!/usr/bin/env Rscript
# =============================================================================
# replication/post/14_coupling.R  --  Jaccard vs overlap coupling by regime
#
# Reads output/coupling/coup_cfg*.csv (sim/04) and writes
#   manuscript/tables/tab_app_coupling_regimes.tex   regime x metric: Jaccard,
#                     Overlap, paired difference (Jaccard - Overlap) with 95% CI
#   manuscript/tables/tab_app_coupling_methods.tex   regime x method means for
#                     all five DynMux couplings (NMI break / continue, purity,
#                     completeness, K MAE)
#   manuscript/figures/fig_coupling_purity.pdf       purity vs completeness by
#                     regime, one point per method x intensity
#   output/coupling_summary/*.csv                    the numbers behind them
# Paired CIs: per (config, rep) difference between the two methods on the same
# simulated series; mean and t-interval across all reps in the regime.
# Usage: DM_ROOT=... Rscript replication/post/14_coupling.R
# =============================================================================
ROOT <- Sys.getenv("DM_ROOT", unset = getwd())
suppressPackageStartupMessages(library(ggplot2))
IN  <- file.path(ROOT, "output", "coupling")
SUM <- file.path(ROOT, "output", "coupling_summary"); dir.create(SUM, showWarnings = FALSE, recursive = TRUE)
TAB <- file.path(ROOT, "manuscript", "tables");  dir.create(TAB, showWarnings = FALSE, recursive = TRUE)
FIG <- file.path(ROOT, "manuscript", "figures"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

d <- do.call(rbind, lapply(list.files(IN, "^coup_cfg.*csv$", full.names = TRUE), read.csv))
cat("rows:", nrow(d), " configs:", length(unique(d$task)), " methods:", paste(unique(d$method), collapse = "; "), "\n")
stopifnot(nrow(d) > 0, all(c("nmi_joint_break", "nmi_joint_continue", "purity_break", "completeness_break") %in% names(d)))
REGIMES <- c(balanced = "Balanced (control)", split = "Split", merge = "Merge", sizeskew = "Size skew",
             shrink = "Shrink", nested = "Nested")
d$regime_lab <- factor(REGIMES[d$regime], levels = REGIMES)
d$method <- sub("^DynMux ", "", d$method)
METRICS <- c(nmi_joint_break = "Joint NMI (break)", nmi_joint_continue = "Joint NMI (continue)",
             purity_break = "Purity", completeness_break = "Completeness", k_mae = "$K$ MAE")

# ---- table 1: regime x method means ------------------------------------------
means <- aggregate(d[, names(METRICS)], by = list(regime = d$regime_lab, method = d$method), FUN = mean)
write.csv(means, file.path(SUM, "coupling_means_by_regime_method.csv"), row.names = FALSE)
fmt <- function(x) sprintf("%.3f", x)
lines <- c("\\begin{tabular}{llccccc}", "\\toprule",
           paste("Regime & Coupling &", paste(METRICS, collapse = " & "), "\\\\"), "\\midrule")
for (rg in levels(means$regime)) {
  sub <- means[means$regime == rg, ]
  if (!nrow(sub)) next
  sub <- sub[order(match(sub$method, c("Jaccard", "Overlap", "weighted Jaccard", "weighted Overlap", "multislice"))), ]
  for (i in seq_len(nrow(sub)))
    lines <- c(lines, paste(if (i == 1) rg else "", "&", sub$method[i], "&",
                            paste(fmt(unlist(sub[i, names(METRICS)])), collapse = " & "), "\\\\"))
  lines <- c(lines, "\\addlinespace")
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}")
writeLines(lines, file.path(TAB, "tab_app_coupling_methods.tex"))

# ---- table 2: paired Jaccard - Overlap differences with 95% CI ---------------
key <- paste(d$task, d$rep)
J <- d[d$method == "Jaccard", ]; O <- d[d$method == "Overlap", ]
J <- J[order(J$task, J$rep), ]; O <- O[order(O$task, O$rep), ]
stopifnot(identical(paste(J$task, J$rep), paste(O$task, O$rep)))
rows <- list()
for (rg in levels(d$regime_lab)) for (m in names(METRICS)) {
  sel <- J$regime_lab == rg
  if (!any(sel)) next
  diff <- J[[m]][sel] - O[[m]][sel]
  ci <- tryCatch(t.test(diff)$conf.int, error = function(e) c(NA_real_, NA_real_))   # constant diffs (mini runs)
  rows[[length(rows) + 1]] <- data.frame(regime = rg, metric = METRICS[[m]], jaccard = mean(J[[m]][sel]),
    overlap = mean(O[[m]][sel]), diff = mean(diff), lo = ci[1], hi = ci[2], n = sum(sel),
    p_wilcoxon = if (sd(diff) > 0) suppressWarnings(wilcox.test(diff)$p.value) else NA_real_, stringsAsFactors = FALSE)
}
pd <- do.call(rbind, rows)
write.csv(pd, file.path(SUM, "coupling_paired_jaccard_minus_overlap.csv"), row.names = FALSE)
lines <- c("\\begin{tabular}{llcccc}", "\\toprule",
           "Regime & Metric & Jaccard & Overlap & Jaccard $-$ Overlap & 95\\% CI \\\\", "\\midrule")
for (rg in levels(d$regime_lab)) {
  sub <- pd[pd$regime == rg, ]
  if (!nrow(sub)) next
  for (i in seq_len(nrow(sub))) {
    star <- if (!is.na(sub$lo[i]) && (sub$lo[i] > 0 || sub$hi[i] < 0)) "$^{*}$" else ""
    lines <- c(lines, sprintf("%s & %s & %.3f & %.3f & %s%.3f%s & [%.3f, %.3f] \\\\", if (i == 1) rg else "",
                              sub$metric[i], sub$jaccard[i], sub$overlap[i], if (sub$diff[i] > 0) "$+$" else "",
                              sub$diff[i], star, sub$lo[i], sub$hi[i]))
  }
  lines <- c(lines, "\\addlinespace")
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}")
writeLines(lines, file.path(TAB, "tab_app_coupling_regimes.tex"))

# ---- figure: purity vs completeness -----------------------------------------
pf <- aggregate(cbind(purity_break, completeness_break) ~ regime_lab + method + intensity, data = d, FUN = mean)
pf$method <- factor(pf$method, levels = c("Jaccard", "Overlap", "weighted Jaccard", "weighted Overlap", "multislice"))
g <- ggplot(pf, aes(x = completeness_break, y = purity_break, colour = method, shape = intensity)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey80", linetype = 2) +
  geom_point(size = 2.6) +
  facet_wrap(~ regime_lab, ncol = 3) +
  scale_colour_manual(values = c("#0072B2", "#D55E00", "#56B4E9", "#E69F00", "#009E73")) +
  scale_shape_manual(values = c(low = 16, high = 17)) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Completeness (share of each true lineage in one tracked community)",
       y = "Purity (share of each tracked community\nfrom one true lineage)", colour = "Coupling", shape = "Intensity") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom", strip.background = element_rect(fill = "grey95"))
ggsave(file.path(FIG, "fig_coupling_purity.pdf"), g, width = 8.5, height = 6.8)

cat("\n--- Jaccard - Overlap by regime (joint NMI, break / continue) ---\n")
print(pd[pd$metric %in% METRICS[1:2], c("regime", "metric", "jaccard", "overlap", "diff", "lo", "hi")], row.names = FALSE, digits = 3)
cat("\nwrote", file.path(TAB, "tab_app_coupling_regimes.tex"), file.path(TAB, "tab_app_coupling_methods.tex"),
    file.path(FIG, "fig_coupling_purity.pdf"), "\n")
