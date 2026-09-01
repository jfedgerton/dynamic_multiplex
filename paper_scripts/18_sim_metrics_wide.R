#!/usr/bin/env Rscript
# 18_sim_metrics_wide.R -- Main + appendix simulation tables.
# MAIN (tab_metrics_wide.tex): 7 methods introduced in the main text
#   (DynMux Jaccard/Overlap with the generator's layer links, incl. the
#   period-lagged seasonal links; multislice relabeled as the benchmark;
#   Hungarian; DynSBM; Pooled; multinet). NMI + K MAE per regime + one
#   overall Time column. Best cell per column shaded (ties within 0.0005).
# APPENDIX (tab_metrics_appendix.tex): all specifications, NMI/MAE/Acc.
set.seed(123)
D <- "replication/extended/output/dynamic"
d <- do.call(rbind, lapply(list.files(D, "dyn_cfg.*csv$", full.names=TRUE), read.csv))
relab <- c("DynMux multislice (adjacent)" = "Multislice (adjacent)",
           "DynMux multislice (adjacent, Louvain)" = "Multislice (adjacent, Louvain)",
           "DynMux multislice (custom)" = "Multislice (custom links)",
           "Cross-sectional + Hungarian" = "Hungarian matching")
d$method <- ifelse(d$method %in% names(relab), relab[d$method], d$method)
regs <- c("birthdeath","churnswitch","regimeshift","seasonality")
rn <- c(birthdeath="Births/Deaths", churnswitch="Churn/Switch",
        regimeshift="Regime Shift", seasonality="Seasonality")
agg <- aggregate(cbind(nmi_joint,k_mae,comembership_acc,runtime_s)~regime+method, d, mean)
ov  <- aggregate(cbind(nmi_joint,runtime_s)~method, d, mean)
esc <- function(s) gsub("([&_%])","\\\\\\1",s)
mainset <- c("DynMux Jaccard","DynMux Overlap","Multislice (adjacent)",
             "Hungarian matching","Dynamic SBM","Pooled Leiden","multinet GLouvain")
mo <- ov[ov$method %in% mainset,]; mord <- mo$method[order(-mo$nmi_joint)]
val <- function(rg,m,col) agg[agg$regime==rg & agg$method==m, col]
tex <- c("\\begin{table}[t]","\\footnotesize","\\setlength{\\tabcolsep}{4pt}","\\centering",
 "\\caption{Method performance by dynamic regime (18 configurations $\\times$ 30 replicates per regime): mean joint NMI on the tracked partition and mean absolute error of the community count (MAE). Time is mean runtime per network in seconds, pooled across regimes. DynMux specifications use the generator's layer links, including the period-lagged seasonal links in the seasonality regime; multislice (adjacent) uses default adjacent links. Best value per column shaded (ties broken at three decimals). Methods sorted by overall joint NMI.}",
 "\\label{tab:metrics_wide}","\\begin{tabular}{lccccccccc}","\\toprule",
 paste0("& ", paste(sprintf("\\multicolumn{2}{c}{%s}", rn[regs]), collapse=" & "), " & \\\\"),
 "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7} \\cmidrule(lr){8-9}",
 paste0("Method & ", paste(rep("NMI & MAE",4), collapse=" & "), " & Time (s) \\\\"),
 "\\midrule")
fmt <- function(x,best,dec=2,lo=FALSE){ v <- sprintf(paste0("%.",dec,"f"),x)
  hit <- if (lo) x <= best+5e-4 else x >= best-5e-4
  ifelse(hit, paste0("\\best{",v,"}"), v) }
bn <- sapply(regs, function(rg) max(agg$nmi_joint[agg$regime==rg & agg$method %in% mainset]))
bm <- sapply(regs, function(rg) min(agg$k_mae[agg$regime==rg & agg$method %in% mainset]))
bt <- min(ov$runtime_s[ov$method %in% mainset])
for (m in mord) {
  cells <- unlist(lapply(regs, function(rg) c(
    fmt(val(rg,m,"nmi_joint"), bn[rg]),
    fmt(val(rg,m,"k_mae"), bm[rg], lo=TRUE))))
  tt <- fmt(ov$runtime_s[ov$method==m], bt, lo=TRUE)
  tex <- c(tex, sprintf("%s & %s & %s \\\\", esc(m), paste(cells,collapse=" & "), tt))
}
tex <- c(tex,"\\bottomrule","\\end{tabular}","\\end{table}")
writeLines(tex,"manuscript/tables/tab_metrics_wide.tex")
aord <- ov$method[order(-ov$nmi_joint)]
at <- c("\\begin{table}[t]","\\scriptsize","\\setlength{\\tabcolsep}{2.5pt}","\\centering",
 "\\caption{All specifications, by dynamic regime: mean joint NMI, community-count MAE, and raw pairwise co-membership accuracy (Acc.). Acc.\\ rewards predicting non-co-membership and therefore favors methods that fragment the community series; it is reported for completeness. Time is mean runtime per network in seconds.}",
 "\\label{tab:metrics_appendix}","\\begin{tabular}{lccccccccccccc}","\\toprule",
 paste0("& ", paste(sprintf("\\multicolumn{3}{c}{%s}", rn[regs]), collapse=" & "), " & \\\\"),
 "\\cmidrule(lr){2-4} \\cmidrule(lr){5-7} \\cmidrule(lr){8-10} \\cmidrule(lr){11-13}",
 paste0("Method & ", paste(rep("NMI & MAE & Acc.",4), collapse=" & "), " & Time (s) \\\\"),
 "\\midrule")
for (m in aord) {
  cells <- unlist(lapply(regs, function(rg) sprintf("%.2f",
    c(val(rg,m,"nmi_joint"), val(rg,m,"k_mae"), val(rg,m,"comembership_acc")))))
  at <- c(at, sprintf("%s & %s & %.2f \\\\", esc(m), paste(cells,collapse=" & "),
                      ov$runtime_s[ov$method==m]))
}
at <- c(at,"\\bottomrule","\\end{tabular}","\\end{table}")
writeLines(at,"manuscript/tables/tab_metrics_appendix.tex")
cat("written both tables\n")
