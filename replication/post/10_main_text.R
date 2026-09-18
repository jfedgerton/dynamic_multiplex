#!/usr/bin/env Rscript
# =============================================================================
# replication/post/10_main_text.R
# Every generated artifact that the MAIN TEXT of the manuscript includes.
#
#   Table 2  manuscript/tables/tab_metrics_wide.tex
#            rank-shaded regime comparison (joint NMI, K MAE) by regime.
#            Source: $DM_ROOT/output/regime/dyn_cfg*.csv   (sim/01)
#            Formerly paper_scripts/19_sim_metrics_heat.R.
#
#   Figure 2 manuscript/figures/fig_coverage_by_config.pdf
#            gated-free coverage by n and separation, faceted by K.
#            Source: $DM_ROOT/output/coverage_grid/cov_task*.csv (sim/02)
#            Formerly panel (1) of paper_scripts/23_coverage_alt_figs.R, which
#            had the jfe4_collab storage path hard-coded.
#
#   Figure 3 manuscript/figures/fig_order_recovery.pdf
#            precision vs recall of Braumoeller-coded orders, one panel per
#            empirical network, one point per order x method (mean over years).
#            Source: $DM_ROOT/output/empirical/order_recovery.csv (empirical/07)
#            No plotting script existed for this figure before; the manuscript
#            figure was produced interactively from 32_setlevel_allorders.R.
#
# Tables are full floats here because the caption of Table 2 documents the
# shading rule and belongs with the generator. Requires
# \usepackage[table]{xcolor} and \usepackage{booktabs}.
#
# Usage:  DM_ROOT=/path/to/dynamic_multiplex Rscript replication/post/10_main_text.R
# =============================================================================
set.seed(123)
suppressMessages({ library(ggplot2) })

ROOT <- Sys.getenv("DM_ROOT", unset = getwd())
REG  <- file.path(ROOT, "output", "regime")
COV  <- file.path(ROOT, "output", "coverage_grid")
EMP  <- file.path(ROOT, "output", "empirical")
TAB  <- file.path(ROOT, "manuscript", "tables")
FIG  <- file.path(ROOT, "manuscript", "figures")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# TABLE 2: tab_metrics_wide.tex
# Within-column ORDINAL shading by rank:
#   1 -> blue!30, 2 -> blue!20, 3 -> blue!10, 4 -> white,
#   5 -> red!10,  6 -> red!20,  7 -> red!30.
# Ranks on 4-dp column means; ties share a shade; displayed to hundredths.
# Rows sorted by overall joint NMI (mean of the four regime means).
# =============================================================================
fs <- list.files(REG, "^dyn_cfg.*csv$", full.names = TRUE)
if (!length(fs)) stop("No regime output in ", REG, " -- run sim/01 first.", call. = FALSE)
d <- do.call(rbind, lapply(fs, read.csv))
cat("regime files:", length(fs), " rows:", nrow(d), "\n")
if (length(fs) != 72L)                             # 4 regimes x 3 n x 3 r x 2 intensity
  warning(sprintf("expected 72 regime files, found %d: table built on a partial design", length(fs)))

mainset <- c("DynMux Jaccard", "DynMux Overlap", "Cross-sectional + Hungarian",
             "DynMux multislice (adjacent)", "Dynamic SBM", "Pooled Leiden",
             "multinet GLouvain")
stopifnot(all(mainset %in% unique(d$method)))
disp <- c("DynMux Jaccard"               = "DynMux Jaccard",
          "DynMux Overlap"               = "DynMux Overlap",
          "Cross-sectional + Hungarian"  = "Hungarian matching",
          "DynMux multislice (adjacent)" = "Multislice adjacent",
          "Dynamic SBM"                  = "Dynamic SBM",
          "Pooled Leiden"                = "Pooled Leiden",
          "multinet GLouvain"            = "Multislice full (\\texttt{multinet})")

agg <- aggregate(cbind(nmi_joint, k_mae) ~ method + regime, data = d, FUN = mean)
agg <- agg[agg$method %in% mainset, ]
overall <- aggregate(nmi_joint ~ method, data = agg, FUN = mean)
ord  <- overall$method[order(-overall$nmi_joint)]
regs <- c("birthdeath", "churnswitch", "regimeshift", "seasonality")
stopifnot(setequal(unique(agg$regime), regs))

shade <- c("\\cellcolor{blue!30}", "\\cellcolor{blue!20}", "\\cellcolor{blue!10}",
           "", "\\cellcolor{red!10}", "\\cellcolor{red!20}", "\\cellcolor{red!30}")
cellv <- function(vals, higher_better) {
  v  <- round(vals, 4)
  rk <- if (higher_better) rank(-v, ties.method = "min") else rank(v, ties.method = "min")
  paste0(shade[rk], sprintf("%.2f", vals))
}
M <- matrix("", nrow = length(ord), ncol = 8)
for (j in seq_along(regs)) {
  a <- agg[agg$regime == regs[j], ]
  a <- a[match(ord, a$method), ]
  stopifnot(!anyNA(a$nmi_joint), !anyNA(a$k_mae))
  M[, 2 * j - 1] <- cellv(a$nmi_joint, TRUE)
  M[, 2 * j]     <- cellv(a$k_mae, FALSE)
}
body <- paste0(unname(disp[ord]), " & ", apply(M, 1, paste, collapse = " & "), " \\\\")
tex <- c(
"\\begin{table}[t]",
"\\scriptsize",
"\\centering",
"\\caption{Method performance by dynamic regime (18 configurations $\\times$ 30 replicates per regime): mean joint NMI on the tracked partition and mean absolute error of the community count ($K$ MAE). Cells are colored by within-column rank, from blue (best) through white (median) to red (worst); shades are ordinal, ties share a shade, and colors are not comparable in magnitude across columns. DynMux specifications use the generator's layer links, including the period-lagged seasonal links in the seasonality regime; multislice adjacent uses default adjacent links. Multislice full (\\texttt{multinet}) is the generalized Louvain implementation in the \\texttt{multinet} R package; multislice adjacent is the \\texttt{dynamicmultiplex} multislice specification with adjacent-layer identity links. Methods sorted by overall joint NMI.}",
"\\label{tab:metrics_wide}",
"\\begin{tabular}{lcccccccc}",
"\\toprule",
"& \\multicolumn{2}{c}{Births \\& Deaths} & \\multicolumn{2}{c}{Churn \\& Switch} & \\multicolumn{2}{c}{Regime Shift} & \\multicolumn{2}{c}{Seasonality} \\\\",
"\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7} \\cmidrule(lr){8-9}",
"Method & Joint NMI & $K$ MAE & Joint NMI & $K$ MAE & Joint NMI & $K$ MAE & Joint NMI & $K$ MAE \\\\",
"\\midrule",
body,
"\\bottomrule",
"\\end{tabular}",
"\\end{table}")
writeLines(tex, file.path(TAB, "tab_metrics_wide.tex"))
cat("wrote tab_metrics_wide.tex; row order:", paste(unname(disp[ord]), collapse = " | "), "\n")

# =============================================================================
# FIGURE 2: fig_coverage_by_config.pdf
# Ungated coverage by n and community separation, faceted by K.
# =============================================================================
fs <- list.files(COV, "^cov_task.*csv$", full.names = TRUE)
if (!length(fs)) stop("No coverage output in ", COV, " -- run sim/02 first.", call. = FALSE)
cv <- do.call(rbind, lapply(fs, read.csv))
cat("coverage files:", length(fs), " rows:", nrow(cv), "\n")
stopifnot(all(c("n", "K", "density", "cov_P_mean") %in% names(cv)))

cv$sep <- factor(cv$density, levels = c("weak", "default", "strong"),
                 labels = c("Weak (0.20 / 0.10)", "Moderate (0.30 / 0.05)", "Strong (0.50 / 0.02)"))
stopifnot(!anyNA(cv$sep))
agg2 <- aggregate(cov_P_mean ~ n + K + sep, data = cv, FUN = mean)
agg2$K_lab <- factor(paste0("K = ", agg2$K), levels = paste0("K = ", sort(unique(agg2$K))))
LEG <- expression("Community separation (" * p['in'] * " / " * p[out] * ")")
p2 <- ggplot(agg2, aes(factor(n), cov_P_mean, group = sep, color = sep,
                       linetype = sep, shape = sep)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey40") +
  geom_line(linewidth = 0.6) + geom_point(size = 2.2) +
  facet_wrap(~K_lab) +
  scale_color_brewer(palette = "Dark2", name = LEG) +
  scale_linetype_manual(values = c("solid", "longdash", "dotdash"), name = LEG) +
  scale_shape_manual(values = c(16, 17, 15), name = LEG) +
  scale_y_continuous(limits = c(0.4, 1.0)) +
  labs(x = "Number of nodes", y = "Empirical coverage") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(FIG, "fig_coverage_by_config.pdf"), p2, width = 7.5, height = 3.6)
ggsave(file.path(FIG, "fig_coverage_by_config.png"), p2, width = 7.5, height = 3.6, dpi = 300)
cat("wrote fig_coverage_by_config.pdf/.png\n")

# =============================================================================
# FIGURE 3: fig_order_recovery.pdf
# Precision vs recall of coded-order recovery. One panel per network, one
# point per (order, method), averaged over the years the order is active.
# Dashed lines at 0.5 split the plane into quadrants; the upper-right
# quadrant is "the community both contains most of the order and is mostly
# the order".
# =============================================================================
orf <- file.path(EMP, "order_recovery.csv")
if (!file.exists(orf)) stop("Missing ", orf, " -- run empirical/07 first.", call. = FALSE)
o <- read.csv(orf, stringsAsFactors = FALSE)
stopifnot(all(c("net", "order", "kind", "year", "method", "J", "prec", "rec", "nB") %in% names(o)))
cat("order-recovery rows:", nrow(o), " nets:", paste(sort(unique(o$net)), collapse = ","), "\n")

net_lab <- c(atop = "Alliances (ATOP)", dca = "Defense cooperation (DCA)",
             igo = "IGO co-membership", trade = "Trade")
meth_lab <- c(Jaccard = "DynMux Jaccard", Overlap = "DynMux Overlap",
              multislice = "Multislice adjacent", multinet = "Multislice full (multinet)",
              Hungarian = "Hungarian matching", Pooled = "Pooled Leiden")
stopifnot(all(o$net %in% names(net_lab)), all(o$method %in% names(meth_lab)))

om <- aggregate(cbind(prec, rec, J) ~ net + order + kind + method, data = o, FUN = mean)
om$net_lab  <- factor(unname(net_lab[om$net]),   levels = unname(net_lab))
om$meth_lab <- factor(unname(meth_lab[om$method]), levels = unname(meth_lab))

p3 <- ggplot(om, aes(rec, prec, colour = meth_lab, shape = meth_lab)) +
  annotate("rect", xmin = 0.5, xmax = 1, ymin = 0.5, ymax = 1, fill = "grey92", colour = NA) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey55") +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey55") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey70") +
  geom_point(size = 2, alpha = 0.85) +
  facet_wrap(~net_lab, ncol = 2) +
  scale_colour_brewer(palette = "Dark2", name = "Method") +
  scale_shape_manual(values = c(16, 17, 15, 18, 3, 4), name = "Method") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  coord_equal() +
  labs(x = "Recall (share of the coded order inside its best-matching community)",
       y = "Precision (share of that community belonging to the order)") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
  guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))
ggsave(file.path(FIG, "fig_order_recovery.pdf"), p3, width = 6.5, height = 7.2)
ggsave(file.path(FIG, "fig_order_recovery.png"), p3, width = 6.5, height = 7.2, dpi = 300)
cat("wrote fig_order_recovery.pdf/.png (", nrow(om), "points )\n")

cat("done 10_main_text\n")
