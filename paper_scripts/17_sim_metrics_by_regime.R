#!/usr/bin/env Rscript
# ============================================================================
# 17_sim_metrics_by_regime.R -- Multi-metric table broken out by regime.
# Within each regime block, methods are sorted by mean joint NMI (descending).
# Produces: manuscript/tables/tab_metrics_regime.tex (8 methods per regime)
#           manuscript/tables/tab_metrics_regime.csv
# ============================================================================
set.seed(123)
D <- "replication/extended/output/dynamic"
files <- list.files(D, "dyn_cfg.*csv$", full.names = TRUE)
d <- do.call(rbind, lapply(files, read.csv))
main_set <- c("DynMux weighted Overlap", "DynMux Jaccard", "DynMux Overlap",
              "DynMux multislice (adjacent)", "Cross-sectional + Hungarian",
              "Pooled Leiden", "multinet GLouvain", "Dynamic SBM")
d <- d[d$method %in% main_set, ]
agg <- aggregate(cbind(nmi_joint, k_mae, comembership_acc, runtime_s)
                 ~ regime + method, data = d, FUN = mean)
write.csv(agg, "manuscript/tables/tab_metrics_regime.csv", row.names = FALSE)
rn <- c(birthdeath = "Births and Deaths", churnswitch = "Churn and Switching",
        regimeshift = "Regime Shift", seasonality = "Seasonality")
esc <- function(s) gsub("([&_%])", "\\\\\\1", s)
con <- c("\\begin{table}[t]", "\\footnotesize", "\\centering",
  "\\caption{Method performance by dynamic regime (18 configurations $\\times$ 30 replicates per regime): mean joint NMI on the tracked partition, mean absolute error of the community count, co-membership accuracy, and mean runtime per network. Methods sorted by joint NMI within regime.}",
  "\\label{tab:metrics_regime}", "\\begin{tabular}{lcccc}", "\\toprule",
  "Method & NMI & $K$ MAE & Co-mem acc. & Time (s) \\\\")
for (rg in c("birthdeath", "churnswitch", "regimeshift", "seasonality")) {
  sub <- agg[agg$regime == rg, ]
  sub <- sub[order(-sub$nmi_joint), ]
  con <- c(con, "\\midrule",
    sprintf("\\multicolumn{5}{l}{\\textbf{%s}} \\\\", rn[rg]), "\\addlinespace[2pt]")
  for (i in seq_len(nrow(sub))) { r <- sub[i, ]
    con <- c(con, sprintf("%s & %.3f & %.2f & %.3f & %.3f \\\\",
      esc(r$method), r$nmi_joint, r$k_mae, r$comembership_acc, r$runtime_s)) }
}
con <- c(con, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(con, "manuscript/tables/tab_metrics_regime.tex")
cat("done\n")
for (rg in c("birthdeath","churnswitch","regimeshift","seasonality")) {
  s <- agg[agg$regime==rg,]; s <- s[order(-s$nmi_joint),]
  cat("\n==", rn[rg], "==\n")
  print(s[, c("method","nmi_joint","k_mae","comembership_acc","runtime_s")],
        row.names = FALSE, digits = 3) }
