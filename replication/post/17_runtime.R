#!/usr/bin/env Rscript
# replication/post/17_runtime.R -- runtime table (appendix, tab_app_runtime.tex)
# Sources: output/regime/dyn_cfg*.csv (runtime_s per method per simulated series, sim/01)
#          replication/slurm/logs/08_emp_fit_*.out (per-method wall time on the empirical networks)
# Output:  manuscript/tables/tab_app_runtime.tex (bare tabular)
# Usage:   DM_ROOT=. Rscript replication/post/17_runtime.R

ROOT <- Sys.getenv("DM_ROOT", unset = getwd())
outdir <- file.path(ROOT, "manuscript", "tables"); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- 1. simulated runtimes -------------------------------------------------
files <- list.files(file.path(ROOT, "output", "regime"), pattern = "^dyn_cfg[0-9]+\\.csv$", full.names = TRUE)
stopifnot(length(files) == 72)
reg <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
stopifnot(all(c("regime", "n", "method", "runtime_s") %in% names(reg)))
stopifnot(nrow(reg) == 34560)  # 72 cfg x 30 reps x 16 method rows

main_methods <- c(
  "DynMux multislice (custom)"   = "Multislice, same links as DynMux",
  "DynMux multislice (adjacent)" = "Multislice, adjacent links",
  "DynMux Jaccard"               = "DynMux (Jaccard)",
  "Cross-sectional + Hungarian"  = "Hungarian matching",
  "Dynamic SBM"                  = "Dynamic SBM",
  "Pooled Leiden"                = "Pooled Leiden",
  "multinet GLouvain"            = "Multislice full (\\texttt{multinet})"
)
missing <- setdiff(names(main_methods), unique(reg$method))
if (length(missing)) stop("method names not found in output/regime: ", paste(missing, collapse = "; "),
                          "\nfound: ", paste(unique(reg$method), collapse = "; "))
reg <- reg[reg$method %in% names(main_methods), ]
ns <- sort(unique(reg$n)); stopifnot(length(ns) == 3)
sim <- aggregate(runtime_s ~ method + n, reg, mean)
cat("simulated series per method x n:", nrow(reg) / (length(main_methods) * length(ns)), "\n")

# ---- 2. empirical wall times from the fit logs ------------------------------
logs <- list.files(file.path(ROOT, "replication", "slurm", "logs"), pattern = "^08_emp_fit_[0-9]+_[0-9]+\\.out$", full.names = TRUE)
stopifnot(length(logs) >= 4)
logs <- logs[order(file.info(logs)$mtime)]  # later logs override earlier ones
emp <- list()
for (lg in logs) {
  L <- readLines(lg, warn = FALSE)
  net <- NA_character_
  for (z in L) {
    m <- regmatches(z, regexec("^\\[([a-z]+)\\] timing summary", z))[[1]]
    if (length(m) == 2) { net <- m[2]; next }
    m <- regmatches(z, regexec("^\\s+(.+?)\\s+(ok|FAILED)\\s+([0-9.]+) s\\s+ncomm=", z))[[1]]
    if (length(m) == 4 && !is.na(net) && m[3] == "ok") emp[[paste(net, m[2])]] <- as.numeric(m[4])
  }
}
stopifnot(length(emp) > 0)
emp_methods <- c(  # resolution 1 fits; DynSBM is not run on the empirical networks
  "DynMux multislice (custom)"   = NA,
  "DynMux multislice (adjacent)" = "DynMux multislice r1",
  "DynMux Jaccard"               = "DynMux Jaccard r1",
  "Cross-sectional + Hungarian"  = "Cross-sectional + Hungarian",
  "Dynamic SBM"                  = NA,
  "Pooled Leiden"                = "Pooled Leiden",
  "multinet GLouvain"            = "multinet GLouvain"
)
nets <- c(atop = "ATOP", dca = "DCA", igo = "IGO", trade = "Trade")
for (nn in names(nets)) for (mm in na.omit(emp_methods))
  if (is.null(emp[[paste(nn, mm)]])) stop("no wall time for ", nn, " / ", mm, " in the fit logs")

# ---- 3. table ----------------------------------------------------------------
fmt <- function(x) ifelse(is.na(x), "--", ifelse(x < 1, sprintf("%.2f", x), ifelse(x < 10, sprintf("%.1f", x), sprintf("%.0f", x))))
rows <- character(0)
for (mk in names(main_methods)) {
  s <- sapply(ns, function(k) sim$runtime_s[sim$method == mk & sim$n == k])
  e <- sapply(names(nets), function(nn) if (is.na(emp_methods[[mk]])) NA_real_ else emp[[paste(nn, emp_methods[[mk]])]])
  rows <- c(rows, sprintf("%s & %s & %s \\\\", main_methods[[mk]], paste(fmt(s), collapse = " & "), paste(fmt(e), collapse = " & ")))
}
hdr <- c("\\begin{tabular}{lccccccc}", "\\toprule",
         sprintf(" & \\multicolumn{%d}{c}{Simulated series} & \\multicolumn{4}{c}{Empirical networks} \\\\", length(ns)),
         sprintf("\\cmidrule(lr){2-%d} \\cmidrule(lr){%d-%d}", 1 + length(ns), 2 + length(ns), 5 + length(ns)),
         sprintf("Method & %s & %s \\\\", paste0("$n = ", ns, "$", collapse = " & "), paste(nets, collapse = " & ")),
         "\\midrule")
tex <- c(hdr, rows, "\\bottomrule", "\\end{tabular}")
writeLines(tex, file.path(outdir, "tab_app_runtime.tex"))
cat(tex, sep = "\n")
ratio <- sim$runtime_s[sim$method == "Dynamic SBM"] / sim$runtime_s[sim$method == "DynMux Jaccard"]
cat(sprintf("\nDynSBM / DynMux Jaccard runtime ratio by n: %s\n", paste(round(ratio, 1), collapse = ", ")))
cat("wrote tab_app_runtime.tex\n")
