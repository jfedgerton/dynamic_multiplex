#!/usr/bin/env Rscript
# 20_coverage_gate.R -- Reliability-gated coverage for co-membership CIs.
# Data: cov_task*.csv in manuscript/output/coverage3_grid (see DM_ROOT below).
#   Per-simulation rows: config cols (n, K, p_switch, p_in, p_out, density,
#   T_layers, weights, resample), spec cols (fit_type, algorithm),
#   cov_P_mean  = share of pairwise 95% co-membership CIs covering truth,
#   width_P_mean = mean interval width in that simulation.
# FINAL GATE RULE (cite this definition in the manuscript):
#   width_P_mean < 0.05  AND  n >= 100
# Calibration/validation split: unique configurations sorted, odd positions
#   -> calibration, even -> validation. Gate chosen on calibration; headline
#   number is the out-of-sample validation coverage.
# Outputs (tables are LaTeX only; nothing downstream consumed the CSVs):
#   manuscript/tables/tab_coverage.tex          (main text)
#   manuscript/tables/tab_coverage_spec.tex     (appendix per-spec)
#   manuscript/tables/tab_coverage_misspec.tex  (appendix, if data)
#   manuscript/tables/tab_coverage_valued.tex   (appendix, if data)
#   manuscript/figures/fig_coverage_curve.pdf/.png   (main text)
# The .tex files are bare booktabs tabulars with no caption or label; the
# manuscript wraps each in its own table float. Requires \usepackage{booktabs}.
set.seed(123)
suppressMessages({ library(ggplot2) })

# --- paths -------------------------------------------------------------
# DM_ROOT is the project root holding manuscript/output, manuscript/tables,
# and manuscript/figures. Set it before running, e.g.
#   export DM_ROOT=/path/to/dynamic_multiplex
# If unset, the current working directory is used, so the script also runs
# correctly when invoked from the project root.
ROOT <- Sys.getenv("DM_ROOT", unset = getwd())
OUT  <- file.path(ROOT, "manuscript", "output")
TAB  <- file.path(ROOT, "manuscript", "tables")
FIG  <- file.path(ROOT, "manuscript", "figures")
if (!dir.exists(OUT)) {
  stop("Simulation output not found at: ", OUT,
       "\n  Set DM_ROOT to the project root, or run from that directory.",
       call. = FALSE)
}
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
readdir <- function(sub) {
  fs <- list.files(file.path(OUT, sub), "^cov_task.*csv$", full.names = TRUE)
  if (!length(fs)) return(NULL)
  do.call(rbind, lapply(fs, read.csv))
}

# --- LaTeX table writer ------------------------------------------------
# Writes a bare booktabs tabular. No caption, no label, no table float:
# the manuscript supplies those, so recompiling the pipeline never
# overwrites caption text that lives in the paper.
write_tex <- function(header, body, align, path) {
  stopifnot(length(body) > 0, nchar(align) > 0)
  writeLines(c(sprintf("\\begin{tabular}{%s}", align),
               "\\toprule", header, "\\midrule",
               body,
               "\\bottomrule", "\\end{tabular}"), path)
  cat("wrote", basename(path), "(", length(body), "rows )\n")
}
commafmt <- function(x) format(x, big.mark = ",", trim = TRUE)

d <- readdir("coverage3_grid")
stopifnot(!is.null(d))
cat("rows:", nrow(d), "\n")

d$cfg  <- paste(d$n, d$K, d$p_switch, d$p_in, d$p_out, d$density,
                d$T_layers, d$weights, d$resample, sep = "|")
d$spec <- paste(d$fit_type, d$algorithm, sep = "/")
cfgs <- sort(unique(d$cfg))
cat("configs:", length(cfgs), " specs:", length(unique(d$spec)), "\n")
calib_cfg <- cfgs[seq(1, length(cfgs), by = 2)]
d$split <- ifelse(d$cfg %in% calib_cfg, "calibration", "validation")

gates <- list(
  "Ungated"                                   = rep(TRUE, nrow(d)),
  "Width $<$ 0.05"                            = d$width_P_mean < 0.05,
  "Width $<$ 0.05, $n \\geq 100$"             = d$width_P_mean < 0.05 & d$n >= 100,
  "Width $<$ 0.05, $n \\geq 100$, Louvain"    = d$width_P_mean < 0.05 & d$n >= 100 & d$algorithm == "louvain"
)

rows <- lapply(names(gates), function(g) {
  k <- gates[[g]]
  data.frame(gate = g,
    calib    = mean(d$cov_P_mean[k & d$split == "calibration"]),
    valid    = mean(d$cov_P_mean[k & d$split == "validation"]),
    retained = mean(k[d$split == "validation"]),
    n_sims   = sum(k))
})
tab <- do.call(rbind, rows)
print(tab)
write_tex(
  header = "Reliability gate & Calibration & Validation & Share retained & Simulations \\\\",
  body   = sprintf("%s & %.3f & %.3f & %.3f & %s \\\\",
                   tab$gate, tab$calib, tab$valid, tab$retained,
                   commafmt(tab$n_sims)),
  align  = "lcccc",
  path   = file.path(TAB, "tab_coverage.tex"))

# --- Per-spec coverage under the final gate, validation side (appendix) ---
k <- gates[[3]]
v <- d[k & d$split == "validation", ]
sp  <- aggregate(cov_P_mean ~ spec, data = v, FUN = mean)
spn <- aggregate(cbind(n_sims = cov_P_mean) ~ spec, data = v, FUN = length)
sp <- merge(sp, spn)
print(sp)

# Reshape to coupling x algorithm so the table shows that coverage is
# invariant to the coupling and varies only with the detection algorithm.
parts <- do.call(rbind, strsplit(sp$spec, "/", fixed = TRUE))
stopifnot(ncol(parts) == 2)
sp$coupling  <- parts[, 1]
sp$algorithm <- parts[, 2]

coup_lab <- c(identity         = "Identity",
              jaccard          = "Jaccard",
              overlap          = "Overlap",
              weighted_jaccard = "Weighted Jaccard",
              weighted_overlap = "Weighted overlap")
alg_lab  <- c(leiden = "Leiden", louvain = "Louvain")
stopifnot(all(sp$coupling  %in% names(coup_lab)))
stopifnot(all(sp$algorithm %in% names(alg_lab)))

cov_w <- tapply(sp$cov_P_mean, list(sp$coupling, sp$algorithm), identity)
n_w   <- tapply(sp$n_sims,     list(sp$coupling, sp$algorithm), identity)
stopifnot(!anyNA(cov_w))   # every coupling x algorithm cell must be present

coup_order <- names(coup_lab)[names(coup_lab) %in% rownames(cov_w)]
alg_order  <- names(alg_lab)[names(alg_lab)  %in% colnames(cov_w)]
stopifnot(length(coup_order) > 0, length(alg_order) > 0)

spec_body <- vapply(coup_order, function(cp) {
  cells <- vapply(alg_order, function(al)
    sprintf("%.3f (%s)", cov_w[cp, al], commafmt(n_w[cp, al])), character(1))
  sprintf("%s & %s \\\\", coup_lab[[cp]], paste(cells, collapse = " & "))
}, character(1))

write_tex(
  header = sprintf("Interlayer coupling & %s \\\\",
                   paste(unname(alg_lab[alg_order]), collapse = " & ")),
  body   = unname(spec_body),
  align  = paste0("l", strrep("c", length(alg_order))),
  path   = file.path(TAB, "tab_coverage_spec.tex"))

# --- Coverage curve: coverage vs binned interval width, by network size ---
d$wbin <- cut(d$width_P_mean, breaks = c(seq(0, 0.15, 0.01), Inf), right = FALSE)
agg <- aggregate(cov_P_mean ~ wbin + n, data = d, FUN = mean)
cnt <- aggregate(cbind(nsims = cov_P_mean) ~ wbin + n, data = d, FUN = length)
agg <- merge(agg, cnt)
agg <- agg[agg$nsims >= 200, ]  # drop unstable bins
mids <- seq(0.005, 0.155, 0.01)
agg$wmid <- mids[as.integer(agg$wbin)]
p <- ggplot(agg, aes(wmid, cov_P_mean, colour = factor(n),
                     linetype = factor(n), shape = factor(n))) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = 0.05, linetype = "dotted", colour = "grey60") +
  geom_line(linewidth = 0.6) + geom_point(size = 1.8) +
  scale_colour_brewer(palette = "Dark2", name = "Nodes (n)") +
  scale_linetype_manual(values = c("solid", "longdash", "dotdash", "twodash"),
                        name = "Nodes (n)") +
  scale_shape_manual(values = c(16, 17, 15, 18), name = "Nodes (n)") +
  labs(x = "Mean interval width", y = "Empirical coverage (nominal 0.95)") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(FIG, "fig_coverage_curve.pdf"), p, width = 6.5, height = 3.6)
ggsave(file.path(FIG, "fig_coverage_curve.png"), p, width = 6.5, height = 3.6, dpi = 300)

# --- Appendix robustness: same gate on misspecified and weighted grids ---
aux <- function(sub, stem) {
  x <- tryCatch(readdir(sub), error = function(e) NULL)
  if (is.null(x) || !all(c("cov_P_mean", "width_P_mean", "n") %in% names(x))) {
    cat("skip:", sub, "\n"); return(invisible(NULL))
  }
  k <- x$width_P_mean < 0.05 & x$n >= 100
  res <- data.frame(
    gate = c("Ungated", "Width $<$ 0.05, $n \\geq 100$"),
    coverage = c(mean(x$cov_P_mean), mean(x$cov_P_mean[k])),
    retained = c(1, mean(k)), n_sims = c(nrow(x), sum(k)))
  print(res)
  write_tex(
    header = "Reliability gate & Coverage & Share retained & Simulations \\\\",
    body   = sprintf("%s & %.3f & %.3f & %s \\\\",
                     res$gate, res$coverage, res$retained, commafmt(res$n_sims)),
    align  = "lccc",
    path   = file.path(TAB, paste0(stem, ".tex")))
}
aux("coverage3_misspec", "tab_coverage_misspec")
aux("coverage3_valued",  "tab_coverage_valued")
cat("done 20\n")
