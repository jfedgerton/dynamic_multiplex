#!/usr/bin/env Rscript
# ============================================================================
# 10_make_tables.R  --  Export LaTeX (booktabs) tables for the paper.
# Run from the repo root on ROAR. Reads:
#   replication/extended/output/dynamic/dyn_cfg*.csv              (regime comparison)
#   replication/extended/output/coverage/*.csv                   (meta-coverage)
#   replication/extended/output/empirical/empirical_summary.rds  (ATOP + DCA)
# Writes: manuscript/tables/tab_{regime,coverage,empirical}.tex
# No external package dependencies (hand-rolled booktabs writer).
# ============================================================================
set.seed(123)
OUTDIR <- "manuscript/tables"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# --- minimal booktabs writer ------------------------------------------------
latex_table <- function(df, file, caption, label, digits = 3, align = NULL) {
  cells <- as.data.frame(lapply(df, function(col)
    if (is.numeric(col)) formatC(round(col, digits), format = "f", digits = digits)
    else as.character(col)), stringsAsFactors = FALSE)
  rows <- apply(cells, 1, paste, collapse = " & ")
  al   <- if (is.null(align)) paste0("l", paste(rep("r", ncol(df) - 1), collapse = "")) else align
  head <- paste(gsub("_", "\\\\_", names(df)), collapse = " & ")
  con  <- c("\\begin{table}[htbp]", "\\centering",
            sprintf("\\caption{%s}", caption), sprintf("\\label{%s}", label),
            sprintf("\\begin{tabular}{%s}", al), "\\toprule",
            paste(head, "\\\\"), "\\midrule",
            paste0(rows, " \\\\"),
            "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(con, file.path(OUTDIR, file)); cat("wrote", file.path(OUTDIR, file), "\n")
}

# order methods with DynMux specs first, baselines last
method_order <- function(m) {
  base <- c("Pooled Leiden", "Pooled Louvain", "Cross-sectional + Hungarian",
            "multinet GLouvain", "Dynamic SBM")
  factor(m, levels = c(setdiff(unique(m), base), intersect(base, unique(m))))
}

# ============================================================================
# TABLE 1 -- regime comparison: mean nmi_joint by method x regime (wide)
# ============================================================================
rf <- list.files("replication/extended/output/dynamic", "dyn_cfg.*csv$", full.names = TRUE)
if (length(rf)) {
  d <- do.call(rbind, lapply(rf, read.csv))
  agg <- aggregate(nmi_joint ~ method + regime, d, mean)
  wide <- reshape(agg, idvar = "method", timevar = "regime", direction = "wide")
  names(wide) <- sub("nmi_joint\\.", "", names(wide))
  wide <- wide[order(method_order(wide$method)), ]
  latex_table(wide, "tab_regime.tex",
    caption = "Mean joint NMI by method and dynamic regime (72 configs $\\times$ 30 reps).",
    label = "tab:regime", digits = 3)
} else cat("skip regime: no dyn_cfg CSVs\n")

# ============================================================================
# TABLE 2 -- inside-gate coverage: mean cov_P_mean by spec x block (wide) + overall
# ============================================================================
cf <- list.files("replication/extended/output/coverage", "\\.csv$", full.names = TRUE)
if (length(cf)) {
  cc <- do.call(rbind, lapply(cf, read.csv))
  g  <- cc[cc$gate == TRUE, ]
  cov <- aggregate(cov_P_mean ~ spec + block, g, mean)
  wide <- reshape(cov, idvar = "spec", timevar = "block", direction = "wide")
  names(wide) <- sub("cov_P_mean\\.", "", names(wide))
  ov <- aggregate(cov_P_mean ~ spec, g, mean); names(ov)[2] <- "overall"
  wide <- merge(wide, ov, by = "spec")
  latex_table(wide, "tab_coverage.tex",
    caption = "Inside-gate empirical coverage of co-assignment intervals (nominal 0.95), by specification and block.",
    label = "tab:coverage", digits = 3)
} else cat("skip coverage: no coverage CSVs\n")

# ============================================================================
# TABLE 3 -- empirical: per-method modularity / persistence / communities
#            for ATOP and DCA side by side
# ============================================================================
sf <- "replication/extended/output/empirical/empirical_summary.rds"
if (file.exists(sf)) {
  s <- readRDS(sf)
  a <- s$atop$df; d <- s$dca$df
  names(a)[-1] <- paste0("atop_", c("mod", "persist", "ncomm", "totcomm"))
  names(d)[-1] <- paste0("dca_",  c("mod", "persist", "ncomm", "totcomm"))
  emp <- merge(a[, c("method","atop_mod","atop_persist","atop_ncomm")],
               d[, c("method","dca_mod","dca_persist","dca_ncomm")], by = "method", all = TRUE)
  emp <- emp[order(method_order(emp$method)), ]
  latex_table(emp, "tab_empirical.tex",
    caption = "Empirical partition quality on ATOP alliances (1816--2018) and Kinne DCAD (1980--2010): mean per-year modularity, year-to-year persistence (NMI), and communities per year.",
    label = "tab:empirical", digits = 3)
} else cat("skip empirical: no empirical_summary.rds\n")

# ============================================================================
# TABLE 4 -- omega robustness: multislice modularity vs interlayer omega
# ============================================================================
of <- "replication/extended/output/empirical/omega_sweep.csv"
if (file.exists(of)) {
  om <- read.csv(of)
  wide <- reshape(om[, c("omega", "network", "mod")], idvar = "omega",
                  timevar = "network", direction = "wide")
  names(wide) <- sub("mod\\.", "", names(wide))
  latex_table(wide, "tab_omega.tex",
    caption = "Multislice (identity) modularity across interlayer coupling $\\omega$ (Leiden resolution 1). Even the best $\\omega$ stays below the Overlap r1 reference (ATOP 0.421, DCA 0.330).",
    label = "tab:omega", digits = 3)
} else cat("skip omega: no omega_sweep.csv\n")

cat("DONE tables\n")
