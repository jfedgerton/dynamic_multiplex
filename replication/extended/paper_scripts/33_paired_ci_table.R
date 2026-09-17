#!/usr/bin/env Rscript
# 33_paired_ci_table.R -- Appendix table of paired 95% confidence intervals
# for every DynMux-vs-baseline comparison, by regime, intensity, and metric.
# Reads the summary written by 15b_sim_paired_intensity.R; run that first.
#
#   Delta Joint NMI = DynMux - baseline   (positive favors DynMux)
#   Delta K MAE     = baseline - DynMux   (positive favors DynMux)
#
# Outputs: manuscript/tables/tab_app_paired_ci.tex  (longtable, appendix)
#          manuscript/tables/tab_app_paired_ci.csv
set.seed(123)

ROOT <- Sys.getenv("DM_ROOT", unset = getwd())
TAB  <- file.path(ROOT, "manuscript", "tables")
SRC  <- file.path(TAB, "tab_regime_paired_intensity.csv")
if (!file.exists(SRC)) {
  stop("Missing ", SRC, "\n  Run 15b_sim_paired_intensity.R first.", call. = FALSE)
}
r <- read.csv(SRC)
stopifnot(all(c("regime","intensity","metric","ref","baseline",
                "n","diff","lo","hi") %in% names(r)))
stopifnot(nrow(r) == 160)                       # 4 x 2 x 2 x 2 x 5

rn <- c(birthdeath  = "Births \\& deaths",
        churnswitch = "Gradual switching",
        regimeshift = "Abrupt rewiring",
        seasonality = "Recurring structure")
bl <- c("Hungarian matching"    = "Hungarian matching",
        "Multislice (adjacent)" = "Multislice (adjacent)",
        "Dynamic SBM"           = "Dynamic SBM",
        "Pooled Leiden"         = "Pooled Leiden",
        "multinet GLouvain"     = "Multislice full (\\texttt{multinet})")

r$regime_lab <- unname(rn[r$regime])
r$base_lab   <- unname(bl[r$baseline])
r$int_lab    <- ifelse(r$intensity == "low", "Low", "High")
r$ref_short  <- sub("DynMux ", "", r$ref)
r$cell <- sprintf("%.3f [%.3f, %.3f]", r$diff, r$lo, r$hi)

w <- reshape(r[, c("metric","regime_lab","int_lab","base_lab","ref_short","cell")],
             idvar = c("metric","regime_lab","int_lab","base_lab"),
             timevar = "ref_short", direction = "wide")
names(w) <- sub("cell.", "", names(w), fixed = TRUE)
w$metric_lab <- ifelse(w$metric == "nmi_joint", "Joint NMI", "$K$ MAE")
w <- w[order(factor(w$metric_lab, levels = c("Joint NMI","$K$ MAE")),
             factor(w$regime_lab, levels = unname(rn)),
             factor(w$int_lab,    levels = c("Low","High")),
             factor(w$base_lab,   levels = unname(bl))), ]
write.csv(w, file.path(TAB, "tab_app_paired_ci.csv"), row.names = FALSE)

# --- LaTeX longtable; blank repeated group labels for readability ----------
hdr <- c(
 "\\begin{longtable}{lllcc}",
 "\\caption{Paired differences and 95\\% confidence intervals for every",
 "DynMux-versus-baseline comparison, by metric, network regime, and change",
 "intensity. Differences are computed within each simulated network",
 sprintf("($n = %d$ per comparison). Positive values favor DynMux for both", unique(r$n)[1]),
 "metrics: joint NMI is DynMux minus baseline, and $K$ MAE is baseline minus",
 "DynMux, since a lower community-count error is better.}",
 "\\label{tab:app_paired_ci} \\\\",
 "\\toprule",
 "Metric & Regime (intensity) & Baseline & Jaccard & Overlap \\\\",
 "\\midrule", "\\endfirsthead",
 "\\multicolumn{5}{l}{\\textit{Table \\ref{tab:app_paired_ci} continued}} \\\\",
 "\\toprule",
 "Metric & Regime (intensity) & Baseline & Jaccard & Overlap \\\\",
 "\\midrule", "\\endhead",
 "\\bottomrule", "\\endlastfoot")

body <- character(0); pm <- ""; pg <- ""
for (i in seq_len(nrow(w))) {
  m <- if (w$metric_lab[i] == pm) "" else w$metric_lab[i]
  g <- sprintf("%s (%s)", w$regime_lab[i], w$int_lab[i])
  gg <- if (g == pg) "" else g
  if (m != "" && i > 1) body <- c(body, "\\addlinespace")
  body <- c(body, sprintf("%s & %s & %s & %s & %s \\\\",
                          m, gg, w$base_lab[i], w$Jaccard[i], w$Overlap[i]))
  pm <- w$metric_lab[i]; pg <- g
}
writeLines(c(hdr, body, "\\end{longtable}"),
           file.path(TAB, "tab_app_paired_ci.tex"))

cat("wrote tab_app_paired_ci.tex /.csv  (", nrow(w), "rows )\n")
cat("comparisons with CI excluding zero:",
    sum(r$lo > 0 | r$hi < 0), "of", nrow(r), "\n")
