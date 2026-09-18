#!/usr/bin/env Rscript
# 34_misspec_breakdown.R -- Coverage in the degree-corrected (DC-SBM) arm,
# broken out by the two misspecification dials rather than pooled.
#
# Why this exists: the pooled misspecification number in 20_coverage_gate.R
# averages over hetero = "none", which sets theta == 1 and is therefore NOT
# misspecified -- it is the plain SBM of the main sweep. Pooling it makes the
# method look more robust than it is. The "none" row here doubles as a
# regression check: it should reproduce the main sweep's gated coverage.
#
# Data: cov_task*.csv in manuscript/output/coverage3_misspec, written by
#   19_coverage_misspec.R. Columns used: cov_P_mean, width_P_mean, n,
#   hetero (none/moderate/severe), balance (balanced/skewed).
# Gate: width_P_mean < 0.05 AND n >= 100  (identical to 20_coverage_gate.R)
#
# Outputs:
#   manuscript/tables/tab_coverage_misspec_hetero.tex   (appendix, 3 rows)
#   manuscript/tables/tab_coverage_misspec_balance.tex  (appendix, 2 rows)
#   manuscript/figures/fig_misspec_curve.pdf/.png       (appendix)
# Tables are bare booktabs tabulars, no caption or label; the manuscript
# wraps each in its own float. Requires \usepackage{booktabs}.
set.seed(123)
suppressMessages({ library(ggplot2) })

# --- paths -------------------------------------------------------------
ROOT <- Sys.getenv("DM_ROOT", unset = getwd())
OUT  <- file.path(ROOT, "manuscript", "output")
TAB  <- file.path(ROOT, "manuscript", "tables")
FIG  <- file.path(ROOT, "manuscript", "figures")
SUB  <- file.path(OUT, "coverage3_misspec")
if (!dir.exists(SUB)) {
  stop("Misspecification output not found at: ", SUB,
       "\n  Set DM_ROOT to the project root, or run 19_coverage_misspec.R first.",
       call. = FALSE)
}
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
commafmt <- function(x) format(x, big.mark = ",", trim = TRUE)

# --- read --------------------------------------------------------------
fs <- list.files(SUB, "^cov_task.*csv$", full.names = TRUE)
stopifnot(length(fs) > 0)
x <- do.call(rbind, lapply(fs, read.csv))
cat("files:", length(fs), " rows:", nrow(x), "\n")
stopifnot(all(c("cov_P_mean", "width_P_mean", "n", "hetero", "balance") %in% names(x)))

x$gated <- x$width_P_mean < 0.05 & x$n >= 100

# --- assertions on the design -------------------------------------------
het_levels <- c("none", "moderate", "severe")
bal_levels <- c("balanced", "skewed")
stopifnot(setequal(unique(x$hetero),  het_levels))
stopifnot(setequal(unique(x$balance), bal_levels))
cat("\ncell counts (hetero x balance):\n"); print(table(x$hetero, x$balance))

# --- breakdown ----------------------------------------------------------
# One row per level: ungated coverage, gated coverage, share retained, sims.
breakdown <- function(col, levels_in_order) {
  do.call(rbind, lapply(levels_in_order, function(g) {
    y <- x[x[[col]] == g, ]
    stopifnot(nrow(y) > 0)
    data.frame(level    = g,
               ungated  = mean(y$cov_P_mean),
               gated    = mean(y$cov_P_mean[y$gated]),
               retained = mean(y$gated),
               n_sims   = nrow(y),
               stringsAsFactors = FALSE)
  }))
}

het <- breakdown("hetero",  het_levels)
bal <- breakdown("balance", bal_levels)

cat("\n--- by degree heterogeneity ---\n"); print(het)
cat("\n--- by community balance ---\n");    print(bal)
cat("\n--- pooled (cross-check against 20_coverage_gate.R) ---\n")
cat("ungated:", sprintf("%.4f", mean(x$cov_P_mean)),
    " gated:", sprintf("%.4f", mean(x$cov_P_mean[x$gated])),
    " retained:", sprintf("%.4f", mean(x$gated)), "\n")

het_lab <- c(none     = "None ($\\theta \\equiv 1$)",
             moderate = "Moderate (lognormal)",
             severe   = "Severe (Pareto hubs)")
bal_lab <- c(balanced = "Balanced", skewed = "Skewed (largest group $\\approx$ 40\\%)")

write_tex(
  header = "Degree heterogeneity & Ungated & Gated & Share retained & Simulations \\\\",
  body   = sprintf("%s & %.3f & %.3f & %.3f & %s \\\\",
                   het_lab[het$level], het$ungated, het$gated,
                   het$retained, commafmt(het$n_sims)),
  align  = "lcccc",
  path   = file.path(TAB, "tab_coverage_misspec_hetero.tex"))

write_tex(
  header = "Community sizes & Ungated & Gated & Share retained & Simulations \\\\",
  body   = sprintf("%s & %.3f & %.3f & %.3f & %s \\\\",
                   bal_lab[bal$level], bal$ungated, bal$gated,
                   bal$retained, commafmt(bal$n_sims)),
  align  = "lcccc",
  path   = file.path(TAB, "tab_coverage_misspec_balance.tex"))

# --- coverage vs interval width, by heterogeneity -----------------------
# Binning and styling deliberately match 20_coverage_gate.R so this figure
# can be read against fig_coverage_curve on identical axes (PA requires
# identical scales across comparable panels).
x$wbin <- cut(x$width_P_mean, breaks = c(seq(0, 0.15, 0.01), Inf), right = FALSE)
agg <- aggregate(cov_P_mean ~ wbin + hetero, data = x, FUN = mean)
cnt <- aggregate(cbind(nsims = cov_P_mean) ~ wbin + hetero, data = x, FUN = length)
agg <- merge(agg, cnt)
agg <- agg[agg$nsims >= 200, ]           # drop unstable bins, as in script 20
stopifnot(nrow(agg) > 0)
mids <- seq(0.005, 0.155, 0.01)
agg$wmid <- mids[as.integer(agg$wbin)]
agg$hetero <- factor(agg$hetero, levels = het_levels,
                     labels = c("None", "Moderate", "Severe"))

cat("\n--- curve: where each level crosses nominal 0.95 ---\n")
for (lv in levels(agg$hetero)) {
  a <- agg[agg$hetero == lv, ]
  a <- a[order(a$wmid), ]
  below <- a$wmid[a$cov_P_mean < 0.95]
  cat(sprintf("%-9s first bin below 0.95: %s\n", lv,
              if (length(below)) sprintf("%.3f", min(below)) else "never"))
}

p <- ggplot(agg, aes(wmid, cov_P_mean, colour = hetero,
                     linetype = hetero, shape = hetero)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = 0.05, linetype = "dotted", colour = "grey60") +
  geom_line(linewidth = 0.6) + geom_point(size = 1.8) +
  scale_colour_brewer(palette = "Dark2", name = "Degree heterogeneity") +
  scale_linetype_manual(values = c("solid", "longdash", "dotdash"),
                        name = "Degree heterogeneity") +
  scale_shape_manual(values = c(16, 17, 15), name = "Degree heterogeneity") +
  labs(x = "Mean interval width", y = "Empirical coverage (nominal 0.95)") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(FIG, "fig_misspec_curve.pdf"), p, width = 6.5, height = 3.6)
ggsave(file.path(FIG, "fig_misspec_curve.png"), p, width = 6.5, height = 3.6, dpi = 300)
cat("wrote fig_misspec_curve.pdf/.png\n")
cat("done 34\n")
