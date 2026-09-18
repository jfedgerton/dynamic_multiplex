#!/usr/bin/env Rscript
# 19_sim_metrics_heat.R -- Rank-shaded main simulation table (published
# Table: tab_metrics_wide). Within-column ORDINAL shading by rank:
#   1 -> blue!30, 2 -> blue!20, 3 -> blue!10, 4 -> (white),
#   5 -> red!10, 6 -> red!20, 7 -> red!30.
# Ranks computed on 4-dp column means; ties share a shade. Displayed values
# rounded to hundredths. Rows sorted by overall joint NMI (mean of the four
# regime means), matching the caption. Requires \usepackage[table]{xcolor}.
# Output: manuscript/tables/tab_metrics_wide.tex
set.seed(123)
D <- "replication/extended/output/dynamic"
d <- do.call(rbind, lapply(list.files(D, "dyn_cfg.*csv$", full.names = TRUE), read.csv))
mainset <- c("DynMux Jaccard", "DynMux Overlap", "Cross-sectional + Hungarian",
             "DynMux multislice (adjacent)", "Dynamic SBM", "Pooled Leiden",
             "multinet GLouvain")
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
ord <- overall$method[order(-overall$nmi_joint)]
regs <- c("birthdeath", "churnswitch", "regimeshift", "seasonality")
shade <- c("\\cellcolor{blue!30}", "\\cellcolor{blue!20}", "\\cellcolor{blue!10}",
           "", "\\cellcolor{red!10}", "\\cellcolor{red!20}", "\\cellcolor{red!30}")
cellv <- function(vals, higher_better) {
  v <- round(vals, 4)
  rk <- if (higher_better) rank(-v, ties.method = "min") else rank(v, ties.method = "min")
  paste0(shade[rk], sprintf("%.2f", vals))
}
M <- matrix("", nrow = length(ord), ncol = 8)
for (j in seq_along(regs)) {
  a <- agg[agg$regime == regs[j], ]
  a <- a[match(ord, a$method), ]
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
writeLines(tex, "manuscript/tables/tab_metrics_wide.tex")
cat("done 19; row order:", paste(unname(disp[ord]), collapse = " | "), "\n")
