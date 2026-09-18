#!/usr/bin/env Rscript
# ============================================================================
# 16_sim_metrics_table.R -- Multi-metric summary table for the regime study.
# Pools replication/extended/output/dynamic/dyn_cfg*.csv across all four
# regimes (72 configs x 30 reps) and reports, per method: mean joint NMI,
# community-count MAE, co-membership accuracy, and runtime.
# Produces: manuscript/tables/tab_metrics.tex      (8 methods, main text)
#           manuscript/tables/tab_metrics_full.tex (all methods, appendix)
#           manuscript/tables/tab_metrics.csv
# ============================================================================
set.seed(123)
D <- "replication/extended/output/dynamic"
files <- list.files(D, "dyn_cfg.*csv$", full.names = TRUE)
d <- do.call(rbind, lapply(files, read.csv))
agg <- aggregate(cbind(nmi_joint, k_mae, comembership_acc, runtime_s) ~ method,
                 data = d, FUN = mean)
agg <- agg[order(-agg$nmi_joint), ]
write.csv(agg, "manuscript/tables/tab_metrics.csv", row.names = FALSE)
main_set <- c("DynMux weighted Overlap", "DynMux Jaccard", "DynMux Overlap",
              "DynMux multislice (adjacent)", "Cross-sectional + Hungarian",
              "Pooled Leiden", "multinet GLouvain", "Dynamic SBM")
esc <- function(s) gsub("([&_%])", "\\\\\\1", s)
mk <- function(sub, file, lab, cap) {
  con <- c("\\begin{table}[htbp]", "\\centering", sprintf("\\caption{%s}", cap),
           sprintf("\\label{%s}", lab), "\\begin{tabular}{lrrrr}", "\\toprule",
           "Method & Joint NMI & $K$ error (MAE) & Co-membership acc. & Runtime (s) \\\\",
           "\\midrule")
  for (i in seq_len(nrow(sub))) { r <- sub[i, ]
    con <- c(con, sprintf("%s & %.3f & %.2f & %.3f & %.3f \\\\",
      esc(r$method), r$nmi_joint, r$k_mae, r$comembership_acc, r$runtime_s)) }
  con <- c(con, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(con, file)
}
mk(agg[agg$method %in% main_set, ], "manuscript/tables/tab_metrics.tex",
   "tab:metrics",
   "Method performance pooled across all four regimes (72 configurations $\\times$ 30 replicates): mean joint NMI on the tracked partition, mean absolute error of the community count, co-membership accuracy, and mean runtime per network.")
mk(agg, "manuscript/tables/tab_metrics_full.tex", "tab:metrics_full",
   "Method performance pooled across all four regimes, all specifications.")
cat("done\n"); print(agg[agg$method %in% main_set,], row.names = FALSE, digits = 3)
