#!/usr/bin/env Rscript
# 20_coverage_gate.R -- Reliability-gated coverage for co-membership CIs.
# Data: cov_task*.csv in jfe4_collab .../manuscript/output/coverage3_grid.
#   Per-simulation rows: config cols (n, K, p_switch, p_in, p_out, density,
#   T_layers, weights, resample), spec cols (fit_type, algorithm),
#   cov_P_mean  = share of pairwise 95% co-membership CIs covering truth,
#   width_P_mean = mean interval width in that simulation.
# FINAL GATE RULE (cite this definition in the manuscript):
#   width_P_mean < 0.05  AND  n >= 100
# Calibration/validation split: unique configurations sorted, odd positions
#   -> calibration, even -> validation. Gate chosen on calibration; headline
#   number is the out-of-sample validation coverage.
# Outputs:
#   manuscript/tables/tab_coverage.csv/.tex         (main text)
#   manuscript/tables/tab_coverage_spec.csv         (appendix per-spec)
#   manuscript/tables/tab_coverage_misspec.csv      (appendix, if data)
#   manuscript/tables/tab_coverage_valued.csv       (appendix, if data)
#   manuscript/figures/fig_coverage_curve.pdf/.png  (main text)
set.seed(123)
suppressMessages({ library(ggplot2) })

OUT <- "/storage/group/LiberalArts/default/jfe4_collab/dynamic_multiplex/manuscript/output"
readdir <- function(sub) {
  fs <- list.files(file.path(OUT, sub), "^cov_task.*csv$", full.names = TRUE)
  if (!length(fs)) return(NULL)
  do.call(rbind, lapply(fs, read.csv))
}

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
write.csv(tab, "manuscript/tables/tab_coverage.csv", row.names = FALSE)

tex <- c("\\begin{tabular}{lcccc}", "\\toprule",
  "Reliability gate & Calibration & Validation & Share retained & Simulations \\\\",
  "\\midrule",
  sprintf("%s & %.3f & %.3f & %.3f & %s \\\\", tab$gate, tab$calib, tab$valid,
          tab$retained, format(tab$n_sims, big.mark = ",", trim = TRUE)),
  "\\bottomrule", "\\end{tabular}")
writeLines(tex, "manuscript/tables/tab_coverage.tex")

# Per-spec coverage under the final gate, validation side (appendix)
k <- gates[[3]]
v <- d[k & d$split == "validation", ]
sp  <- aggregate(cov_P_mean ~ spec, data = v, FUN = mean)
spn <- aggregate(cbind(n_sims = cov_P_mean) ~ spec, data = v, FUN = length)
sp <- merge(sp, spn)
print(sp)
write.csv(sp, "manuscript/tables/tab_coverage_spec.csv", row.names = FALSE)

# Coverage curve: empirical coverage vs binned interval width, by network size
d$wbin <- cut(d$width_P_mean, breaks = c(seq(0, 0.15, 0.01), Inf), right = FALSE)
agg <- aggregate(cov_P_mean ~ wbin + n, data = d, FUN = mean)
cnt <- aggregate(cbind(nsims = cov_P_mean) ~ wbin + n, data = d, FUN = length)
agg <- merge(agg, cnt)
agg <- agg[agg$nsims >= 200, ]  # drop unstable bins
mids <- seq(0.005, 0.155, 0.01)
agg$wmid <- mids[as.integer(agg$wbin)]
p <- ggplot(agg, aes(wmid, cov_P_mean, colour = factor(n))) +
  geom_hline(yintercept = 0.95, linetype = 2, colour = "grey40") +
  geom_vline(xintercept = 0.05, linetype = 3, colour = "grey60") +
  geom_line() + geom_point(size = 1.3) +
  labs(x = "Mean interval width", y = "Empirical coverage (nominal 0.95)",
       colour = "Nodes (n)") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave("manuscript/figures/fig_coverage_curve.pdf", p, width = 6.5, height = 3.6)
ggsave("manuscript/figures/fig_coverage_curve.png", p, width = 6.5, height = 3.6, dpi = 300)

# Appendix robustness: same gate on misspecified and weighted grids, if present
aux <- function(sub, out) {
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
  write.csv(res, out, row.names = FALSE)
}
aux("coverage3_misspec", "manuscript/tables/tab_coverage_misspec.csv")
aux("coverage3_valued",  "manuscript/tables/tab_coverage_valued.csv")
cat("done 20\n")
