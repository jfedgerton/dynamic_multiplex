#!/usr/bin/env Rscript
# =============================================================================
# replication/post/12_appendix_coverage.R
# Reliability-gated coverage of the Wilson co-membership intervals: the gate
# ladder, the width x n gate grid, per-spec coverage, the coverage-vs-width
# curve, coverage by width band inside the gate, and the weighted and
# degree-corrected arms.
#
# Data (one row per simulation, written by sim/02, sim/03, sim/04):
#   $DM_ROOT/output/coverage_grid/cov_task*.csv     binary SBM  (594 cfg x 6 spec)
#   $DM_ROOT/output/coverage_valued/cov_task*.csv   weighted    (720 cfg)
#   $DM_ROOT/output/coverage_misspec/cov_task*.csv  DC-SBM      (432 cfg)
#   Columns used: n, K, p_switch, p_in, p_out, density, T_layers, weights,
#   resample, fit_type, algorithm, cov_P_mean, width_P_mean
#   (+ hetero, balance in the misspec arm).
#
# FINAL GATE RULE (cite this definition in the manuscript):
#   width_P_mean < 0.05  AND  n >= 100
# Calibration/validation split: unique configurations sorted, odd positions
#   -> calibration, even -> validation. Gate chosen on calibration; headline
#   number is the out-of-sample validation coverage.
#
# Outputs (bare booktabs tabulars, no caption/label; figures PDF + PNG):
#   tables/tab_coverage.tex                  gate ladder (main text / appendix)
#   tables/tab_coverage_gate_grid.tex        width threshold x minimum n
#   tables/tab_coverage_spec.tex             coupling x algorithm under the gate
#   tables/tab_coverage_width_bands.tex      coverage by width band, n >= 100
#   tables/tab_coverage_valued.tex           weighted arm, ungated vs gated
#   tables/tab_coverage_misspec.tex          DC-SBM arm, ungated vs gated
#   tables/tab_coverage_misspec_hetero.tex   DC-SBM arm by degree heterogeneity
#   tables/tab_coverage_misspec_balance.tex  DC-SBM arm by community balance
#   figures/fig_coverage_curve.pdf           coverage vs width, by n
#   figures/fig_misspec_curve.pdf            coverage vs width, by heterogeneity
#
# Formerly paper_scripts/20_coverage_gate.R and 34_misspec_breakdown.R.
# Requires \usepackage{booktabs}.
#
# Usage:  DM_ROOT=/path/to/dynamic_multiplex Rscript replication/post/12_appendix_coverage.R
# =============================================================================
set.seed(123)
suppressMessages({ library(ggplot2) })

ROOT <- Sys.getenv("DM_ROOT", unset = getwd())
OUT  <- file.path(ROOT, "output")
TAB  <- file.path(ROOT, "manuscript", "tables")
FIG  <- file.path(ROOT, "manuscript", "figures")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

readdir <- function(sub) {
  fs <- list.files(file.path(OUT, sub), "^cov_task.*csv$", full.names = TRUE)
  if (!length(fs)) return(NULL)
  cat(sub, ": files", length(fs), "\n")
  do.call(rbind, lapply(fs, read.csv))
}
write_tex <- function(header, body, align, path) {
  stopifnot(length(body) > 0, nchar(align) > 0)
  writeLines(c(sprintf("\\begin{tabular}{%s}", align),
               "\\toprule", header, "\\midrule",
               body,
               "\\bottomrule", "\\end{tabular}"), path)
  cat("wrote", basename(path), "(", length(body), "rows )\n")
}
commafmt <- function(x) format(x, big.mark = ",", trim = TRUE)
W_GATE <- 0.05
N_GATE <- 100

# =============================================================================
# Binary arm
# =============================================================================
d <- readdir("coverage_grid")
if (is.null(d)) stop("No coverage output in ", file.path(OUT, "coverage_grid"),
                     " -- run sim/02 first.", call. = FALSE)
stopifnot(all(c("n", "K", "p_switch", "p_in", "p_out", "density", "T_layers", "weights",
                "resample", "fit_type", "algorithm", "cov_P_mean", "width_P_mean") %in% names(d)))
cat("rows:", nrow(d), "\n")

d$cfg  <- paste(d$n, d$K, d$p_switch, d$p_in, d$p_out, d$density,
                d$T_layers, d$weights, d$resample, sep = "|")
d$spec <- paste(d$fit_type, d$algorithm, sep = "/")
cfgs <- sort(unique(d$cfg))
cat("configs:", length(cfgs), " specs:", length(unique(d$spec)), "\n")
calib_cfg <- cfgs[seq(1, length(cfgs), by = 2)]
d$split <- ifelse(d$cfg %in% calib_cfg, "calibration", "validation")

# --- gate ladder ----------------------------------------------------------
gates <- list(
  "Ungated"                                = rep(TRUE, nrow(d)),
  "Width $<$ 0.05"                         = d$width_P_mean < W_GATE,
  "Width $<$ 0.05, $n \\geq 100$"          = d$width_P_mean < W_GATE & d$n >= N_GATE,
  "Width $<$ 0.05, $n \\geq 100$, Louvain" = d$width_P_mean < W_GATE & d$n >= N_GATE & d$algorithm == "louvain"
)
tab <- do.call(rbind, lapply(names(gates), function(g) {
  k <- gates[[g]]
  data.frame(gate = g,
    calib    = mean(d$cov_P_mean[k & d$split == "calibration"]),
    valid    = mean(d$cov_P_mean[k & d$split == "validation"]),
    retained = mean(k),
    n_sims   = sum(k))
}))
print(tab)
write_tex(
  header = "Reliability gate & Calibration & Validation & Share retained & Simulations \\\\",
  body   = sprintf("%s & %.3f & %.3f & %.3f & %s \\\\",
                   tab$gate, tab$calib, tab$valid, tab$retained, commafmt(tab$n_sims)),
  align  = "lcccc",
  path   = file.path(TAB, "tab_coverage.tex"))

# --- gate grid: width threshold x minimum network size -------------------
# Width thresholds are round hundredths; n thresholds are the network sizes
# in the design. Neither contains a value chosen after inspecting results.
W_GRID <- c(0.03, 0.04, 0.05, 0.06)
N_GRID <- c(50, 100, 200, 400)
stopifnot(all(N_GRID %in% unique(d$n)))
grid <- expand.grid(w = W_GRID, nmin = N_GRID)
grid$calib <- NA_real_; grid$valid <- NA_real_; grid$retained <- NA_real_; grid$n_sims <- NA_integer_
for (i in seq_len(nrow(grid))) {
  k <- d$width_P_mean < grid$w[i] & d$n >= grid$nmin[i]
  grid$n_sims[i]   <- sum(k)
  grid$retained[i] <- mean(k)
  if (any(k & d$split == "calibration") && any(k & d$split == "validation")) {
    grid$calib[i] <- mean(d$cov_P_mean[k & d$split == "calibration"])
    grid$valid[i] <- mean(d$cov_P_mean[k & d$split == "validation"])
  }
}
cat("\n--- gate grid (width x minimum n) ---\n"); print(grid)
gap <- abs(grid$calib - grid$valid)
cat("max |calibration - validation| across grid:", sprintf("%.4f", max(gap, na.rm = TRUE)), "\n")
if (any(gap >= 0.02, na.rm = TRUE))               # split-specific noise would show here
  warning("CHECK: calibration and validation coverage differ by >= 0.02 in ",
          sum(gap >= 0.02, na.rm = TRUE), " gate-grid cell(s); the manuscript claims agreement.")
grid_body <- vapply(W_GRID, function(w) {
  cells <- vapply(N_GRID, function(nm) {
    j <- which(grid$w == w & grid$nmin == nm)
    if (is.na(grid$valid[j]) || grid$n_sims[j] == 0) "---"
    else sprintf("%.3f (%.1f\\%%)", grid$valid[j], 100 * grid$retained[j])
  }, character(1))
  sprintf("$<$ %.2f & %s \\\\", w, paste(cells, collapse = " & "))
}, character(1))
write_tex(
  header = sprintf("Maximum interval width & %s \\\\",
                   paste(sprintf("$n \\geq %d$", N_GRID), collapse = " & ")),
  body   = grid_body,
  align  = paste0("l", strrep("c", length(N_GRID))),
  path   = file.path(TAB, "tab_coverage_gate_grid.tex"))

# Resolution: simulations are clustered within configurations, so the
# effective sample size is the number of gated validation configurations.
k_adopted <- d$width_P_mean < W_GATE & d$n >= N_GATE
v_adopted <- d[k_adopted & d$split == "validation", ]
cfg_cov <- tapply(v_adopted$cov_P_mean, v_adopted$cfg, mean)
cat("gated validation configurations:", length(cfg_cov),
    " between-config SD:", sprintf("%.4f", sd(cfg_cov)),
    " implied SE:", sprintf("%.4f", sd(cfg_cov) / sqrt(length(cfg_cov))), "\n")

# --- per-spec coverage under the gate, validation side -------------------
v   <- d[k_adopted & d$split == "validation", ]
sp  <- aggregate(cov_P_mean ~ spec, data = v, FUN = mean)
spn <- aggregate(cbind(n_sims = cov_P_mean) ~ spec, data = v, FUN = length)
sp  <- merge(sp, spn)
parts <- do.call(rbind, strsplit(sp$spec, "/", fixed = TRUE))
sp$coupling <- parts[, 1]; sp$algorithm <- parts[, 2]
coup_lab <- c(identity = "Identity", jaccard = "Jaccard", overlap = "Overlap",
              weighted_jaccard = "Weighted Jaccard", weighted_overlap = "Weighted overlap")
alg_lab  <- c(leiden = "Leiden", louvain = "Louvain")
stopifnot(all(sp$coupling %in% names(coup_lab)), all(sp$algorithm %in% names(alg_lab)))
cov_w <- tapply(sp$cov_P_mean, list(sp$coupling, sp$algorithm), identity)
n_w   <- tapply(sp$n_sims,     list(sp$coupling, sp$algorithm), identity)
stopifnot(!anyNA(cov_w))
coup_order <- names(coup_lab)[names(coup_lab) %in% rownames(cov_w)]
alg_order  <- names(alg_lab)[names(alg_lab)  %in% colnames(cov_w)]
spec_body <- vapply(coup_order, function(cp) {
  cells <- vapply(alg_order, function(al)
    sprintf("%.3f (%s)", cov_w[cp, al], commafmt(n_w[cp, al])), character(1))
  sprintf("%s & %s \\\\", coup_lab[[cp]], paste(cells, collapse = " & "))
}, character(1))
write_tex(
  header = sprintf("Interlayer coupling & %s \\\\", paste(unname(alg_lab[alg_order]), collapse = " & ")),
  body   = unname(spec_body),
  align  = paste0("l", strrep("c", length(alg_order))),
  path   = file.path(TAB, "tab_coverage_spec.tex"))

# --- coverage by width band inside n >= 100 (conditional, not marginal) --
# The gate's 0.960 is a marginal average over bands with very different
# coverage; this table shows the composition.
bands <- c(0, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10, Inf)
dn <- d[d$n >= N_GATE, ]
dn$band <- cut(dn$width_P_mean, breaks = bands, right = TRUE, include.lowest = TRUE)
wb <- aggregate(cov_P_mean ~ band, data = dn, FUN = mean)
wbn <- aggregate(cbind(n_sims = cov_P_mean) ~ band, data = dn, FUN = length)
wb <- merge(wb, wbn)
wb <- wb[order(as.integer(wb$band)), ]
wb$share <- wb$n_sims / nrow(dn)
print(wb)
band_lab <- function(b) {
  lo <- bands[as.integer(b)]; hi <- bands[as.integer(b) + 1]
  ifelse(is.infinite(hi), sprintf("$>$ %.2f", lo),
         ifelse(lo == 0, sprintf("$\\leq$ %.2f", hi), sprintf("(%.2f, %.2f]", lo, hi)))
}
write_tex(
  header = "Interval width band & Coverage & Share of $n \\geq 100$ simulations & Simulations \\\\",
  body   = sprintf("%s & %.3f & %.1f\\%% & %s \\\\",
                   band_lab(wb$band), wb$cov_P_mean, 100 * wb$share, commafmt(wb$n_sims)),
  align  = "lccr",
  path   = file.path(TAB, "tab_coverage_width_bands.tex"))

# --- coverage curve: coverage vs binned width, by n ----------------------
d$wbin <- cut(d$width_P_mean, breaks = c(seq(0, 0.15, 0.01), Inf), right = FALSE)
agg <- aggregate(cov_P_mean ~ wbin + n, data = d, FUN = mean)
cnt <- aggregate(cbind(nsims = cov_P_mean) ~ wbin + n, data = d, FUN = length)
agg <- merge(agg, cnt)
agg <- agg[agg$nsims >= 200, ]                    # drop unstable bins
mids <- seq(0.005, 0.155, 0.01)
agg$wmid <- mids[as.integer(agg$wbin)]
p <- ggplot(agg, aes(wmid, cov_P_mean, colour = factor(n), linetype = factor(n), shape = factor(n))) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = W_GATE, linetype = "dotted", colour = "grey60") +
  geom_line(linewidth = 0.6) + geom_point(size = 1.8) +
  scale_colour_brewer(palette = "Dark2", name = "Nodes (n)") +
  scale_linetype_manual(values = c("solid", "longdash", "dotdash", "twodash"), name = "Nodes (n)") +
  scale_shape_manual(values = c(16, 17, 15, 18), name = "Nodes (n)") +
  labs(x = "Mean interval width", y = "Empirical coverage (nominal 0.95)") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(FIG, "fig_coverage_curve.pdf"), p, width = 6.5, height = 3.6)
ggsave(file.path(FIG, "fig_coverage_curve.png"), p, width = 6.5, height = 3.6, dpi = 300)
cat("wrote fig_coverage_curve.pdf/.png\n")

# =============================================================================
# Weighted and degree-corrected arms: same gate, ungated vs gated
# =============================================================================
arm_table <- function(x, stem) {
  k <- x$width_P_mean < W_GATE & x$n >= N_GATE
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

xv <- readdir("coverage_valued")
if (is.null(xv)) {
  cat("skip: coverage_valued (run sim/03)\n")
} else {
  arm_table(xv, "tab_coverage_valued")
}

xm <- readdir("coverage_misspec")
if (is.null(xm)) {
  cat("skip: coverage_misspec (run sim/04)\n")
} else {
  stopifnot(all(c("hetero", "balance") %in% names(xm)))
  arm_table(xm, "tab_coverage_misspec")

  # Breakdown by the two misspecification dials. hetero = "none" sets
  # theta == 1 and is the plain SBM: it doubles as a regression check
  # against the main sweep's gated coverage.
  xm$gated <- xm$width_P_mean < W_GATE & xm$n >= N_GATE
  het_levels <- c("none", "moderate", "severe")
  bal_levels <- c("balanced", "skewed")
  stopifnot(setequal(unique(xm$hetero), het_levels), setequal(unique(xm$balance), bal_levels))
  breakdown <- function(col, lv) {
    do.call(rbind, lapply(lv, function(g) {
      y <- xm[xm[[col]] == g, ]
      stopifnot(nrow(y) > 0)
      data.frame(level = g, ungated = mean(y$cov_P_mean), gated = mean(y$cov_P_mean[y$gated]),
                 retained = mean(y$gated), n_sims = nrow(y), stringsAsFactors = FALSE)
    }))
  }
  het <- breakdown("hetero",  het_levels)
  bal <- breakdown("balance", bal_levels)
  cat("\n--- by degree heterogeneity ---\n"); print(het)
  cat("\n--- by community balance ---\n");    print(bal)
  het_lab <- c(none = "None ($\\theta \\equiv 1$)", moderate = "Moderate (lognormal)",
               severe = "Severe (Pareto hubs)")
  bal_lab <- c(balanced = "Balanced", skewed = "Skewed (largest group $\\approx$ 40\\%)")
  write_tex(
    header = "Degree heterogeneity & Ungated & Gated & Share retained & Simulations \\\\",
    body   = sprintf("%s & %.3f & %.3f & %.3f & %s \\\\",
                     het_lab[het$level], het$ungated, het$gated, het$retained, commafmt(het$n_sims)),
    align  = "lcccc",
    path   = file.path(TAB, "tab_coverage_misspec_hetero.tex"))
  write_tex(
    header = "Community sizes & Ungated & Gated & Share retained & Simulations \\\\",
    body   = sprintf("%s & %.3f & %.3f & %.3f & %s \\\\",
                     bal_lab[bal$level], bal$ungated, bal$gated, bal$retained, commafmt(bal$n_sims)),
    align  = "lcccc",
    path   = file.path(TAB, "tab_coverage_misspec_balance.tex"))

  # coverage vs width by heterogeneity, on the same axes as fig_coverage_curve
  xm$wbin <- cut(xm$width_P_mean, breaks = c(seq(0, 0.15, 0.01), Inf), right = FALSE)
  am <- aggregate(cov_P_mean ~ wbin + hetero, data = xm, FUN = mean)
  cm <- aggregate(cbind(nsims = cov_P_mean) ~ wbin + hetero, data = xm, FUN = length)
  am <- merge(am, cm)
  am <- am[am$nsims >= 200, ]
  if (!nrow(am)) stop("fig_misspec_curve: no width bin has >= 200 simulations", call. = FALSE)
  am$wmid   <- mids[as.integer(am$wbin)]
  am$hetero <- factor(am$hetero, levels = het_levels, labels = c("None", "Moderate", "Severe"))
  pm <- ggplot(am, aes(wmid, cov_P_mean, colour = hetero, linetype = hetero, shape = hetero)) +
    geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey40") +
    geom_vline(xintercept = W_GATE, linetype = "dotted", colour = "grey60") +
    geom_line(linewidth = 0.6) + geom_point(size = 1.8) +
    scale_colour_brewer(palette = "Dark2", name = "Degree heterogeneity") +
    scale_linetype_manual(values = c("solid", "longdash", "dotdash"), name = "Degree heterogeneity") +
    scale_shape_manual(values = c(16, 17, 15), name = "Degree heterogeneity") +
    labs(x = "Mean interval width", y = "Empirical coverage (nominal 0.95)") +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
  ggsave(file.path(FIG, "fig_misspec_curve.pdf"), pm, width = 6.5, height = 3.6)
  ggsave(file.path(FIG, "fig_misspec_curve.png"), pm, width = 6.5, height = 3.6, dpi = 300)
  cat("wrote fig_misspec_curve.pdf/.png\n")
}

cat("done 12_appendix_coverage\n")
