set.seed(123)
OUTDIR <- "manuscript/tables"; dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
latex_table <- function(df, file, caption, label, digits = 3, align = NULL) {
  cells <- as.data.frame(lapply(df, function(col)
    if (is.numeric(col)) formatC(round(col, digits), format = "f", digits = digits)
    else as.character(col)), stringsAsFactors = FALSE)
  rows <- apply(cells, 1, paste, collapse = " & ")
  al <- if (is.null(align)) paste0("l", paste(rep("r", ncol(df) - 1), collapse = "")) else align
  head <- paste(gsub("_", "\\\\_", names(df)), collapse = " & ")
  con <- c("\\begin{table}[htbp]", "\\centering", sprintf("\\caption{%s}", caption),
    sprintf("\\label{%s}", label), sprintf("\\begin{tabular}{%s}", al), "\\toprule",
    paste(head, "\\\\"), "\\midrule", paste0(rows, " \\\\"), "\\bottomrule",
    "\\end{tabular}", "\\end{table}")
  writeLines(con, file.path(OUTDIR, file)); cat("wrote", file.path(OUTDIR, file), "\n") }
morder <- function(m) { base <- c("Pooled Leiden","Pooled Louvain","Cross-sectional + Hungarian","multinet GLouvain","Dynamic SBM")
  factor(m, levels = c(setdiff(unique(m), base), intersect(base, unique(m)))) }
rf <- list.files("replication/extended/output/dynamic", "dyn_cfg.*csv$", full.names = TRUE)
if (length(rf)) { d <- do.call(rbind, lapply(rf, read.csv))
  agg <- aggregate(nmi_joint ~ method + regime, d, mean)
  w <- reshape(agg, idvar = "method", timevar = "regime", direction = "wide")
  names(w) <- sub("nmi_joint\\.", "", names(w)); w <- w[order(morder(w$method)), ]
  latex_table(w, "tab_regime.tex", "Mean joint NMI by method and dynamic regime (72 configs $\\times$ 30 reps).", "tab:regime") } else cat("skip regime\n")
cf <- list.files("replication/extended/output/coverage", "\\.csv$", full.names = TRUE)
if (length(cf)) { cc <- do.call(rbind, lapply(cf, read.csv)); g <- cc[cc$gate == TRUE, ]
  cov <- aggregate(cov_P_mean ~ spec + block, g, mean)
  w <- reshape(cov, idvar = "spec", timevar = "block", direction = "wide"); names(w) <- sub("cov_P_mean\\.", "", names(w))
  ov <- aggregate(cov_P_mean ~ spec, g, mean); names(ov)[2] <- "overall"; w <- merge(w, ov, by = "spec")
  latex_table(w, "tab_coverage.tex", "Inside-gate coverage of co-assignment intervals (nominal 0.95), by specification and block.", "tab:coverage") } else cat("skip coverage\n")
sf <- "replication/extended/output/empirical/empirical_summary.rds"
if (file.exists(sf)) { s <- readRDS(sf); a <- s$atop$df; d <- s$dca$df
  names(a)[-1] <- paste0("atop_", c("mod","persist","ncomm","totcomm")); names(d)[-1] <- paste0("dca_", c("mod","persist","ncomm","totcomm"))
  emp <- merge(a[,c("method","atop_mod","atop_persist","atop_ncomm")], d[,c("method","dca_mod","dca_persist","dca_ncomm")], by = "method", all = TRUE)
  emp <- emp[order(morder(emp$method)), ]
  latex_table(emp, "tab_empirical.tex", "Empirical partition quality on ATOP (1816--2018) and DCAD (1980--2010).", "tab:empirical") } else cat("skip empirical\n")
of <- "replication/extended/output/empirical/omega_sweep.csv"
if (file.exists(of)) { om <- read.csv(of)
  w <- reshape(om[,c("omega","network","mod")], idvar = "omega", timevar = "network", direction = "wide"); names(w) <- sub("mod\\.", "", names(w))
  latex_table(w, "tab_omega.tex", "Multislice modularity across interlayer $\\omega$ (resolution 1); Overlap r1 reference ATOP 0.421, DCA 0.330.", "tab:omega") } else cat("skip omega\n")
cat("DONE tables\n")
