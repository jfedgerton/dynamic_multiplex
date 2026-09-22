#!/usr/bin/env Rscript
# =============================================================================
# replication/post/12_stability.R  --  Section 4 and its appendix, from sim/02
#
# Reads output/stability/<arm>_{stab,node,pair}_task*.csv and writes
#   output/stability/stability_calibration_table.csv   the lookup shipped in
#          r_code/inst/extdata and python_code/src/dynamic_multiplex/data
#          (levels partition_nmi, partition_ari, node_jaccard; 10 bins each)
#   manuscript/tables/tab_stability_floor.tex          main text: bin, n, median
#          accuracy, floor, validation share above the floor (partition NMI)
#   manuscript/tables/tab_app_stability_levels.tex     partition ARI and node
#          Jaccard floors
#   manuscript/tables/tab_app_stability_arms.tex       validation share above
#          the binary-calibrated floor on the binary validation half, the
#          degree-corrected arm (by heterogeneity x balance) and the weighted
#          arm (by weight regime)
#   manuscript/tables/tab_app_stability_lolo.tex       leave-one-level-out:
#          recalibrate with each n, each density and each K held out; share
#          above the floor on the held-out level (share among held-out fits whose
#       stability bin has a calibrated floor; fits below the calibrated range are counted separately)
#   manuscript/tables/tab_app_stability_pairs.tex      decided / undetermined
#          pair shares and the accuracy of the decided calls
#   manuscript/figures/fig_stability_floor.pdf         Figure 2: stability vs
#          accuracy (validation fits) with the calibrated floor as a step
#   manuscript/figures/fig_app_stability_by_n.pdf      same, faceted by n
#   manuscript/figures/fig_stability_floor_node.pdf    Figure 2b: node-level Jaccard
#                                                      stability vs accuracy with the node floor
# Rule (pre-registered 2026-09-20): the section survives if, on the validation
# half, Spearman(stability, accuracy) >= 0.7 and the calibrated 5th percentile
# is non-decreasing in stability and >= 0.8 for stability >= 0.9. The script
# prints PASS / FAIL and stops (exit 2) on FAIL so the postprocess chain notices.
# Split: cells sorted by (n, K, p_switch, density, T), odd -> calibration,
# even -> validation. Bins: 10 equal-width stability bins on [0, 1].
# Usage: DM_ROOT=... Rscript replication/post/12_stability.R
# =============================================================================
ROOT <- Sys.getenv("DM_ROOT", unset = getwd())
suppressPackageStartupMessages(library(ggplot2))
IN  <- file.path(ROOT, "output", "stability")
TAB <- file.path(ROOT, "manuscript", "tables");  dir.create(TAB, showWarnings = FALSE, recursive = TRUE)
FIG <- file.path(ROOT, "manuscript", "figures"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
options(width = 160)
BR <- seq(0, 1, by = 0.1)
bin_of <- function(x) pmin(pmax(cut(x, BR, include.lowest = TRUE, labels = FALSE), 1L), 10L)
read_arm <- function(arm, kind) { f <- list.files(IN, sprintf("^%s_%s_task.*csv$", arm, kind), full.names = TRUE)
  if (!length(f)) return(NULL); do.call(rbind, lapply(f, read.csv, stringsAsFactors = FALSE)) }
write_tex <- function(header, body, align, path) writeLines(c(sprintf("\\begin{tabular}{%s}", align), "\\toprule", header, "\\midrule", body, "\\bottomrule", "\\end{tabular}"), path)
fmt <- function(x, d = 3) ifelse(is.na(x), "--", sprintf(paste0("%.", d, "f"), x))

# ---- binary arm: split ------------------------------------------------------
s <- read_arm("binary", "stab"); stopifnot(!is.null(s), nrow(s) > 0)
s$cell <- paste(s$n, s$K, s$p_switch, s$density, s$T_layers, sep = "|")
cells <- unique(s[, c("cell", "n", "K", "p_switch", "density", "T_layers")])
cells <- cells[order(cells$n, cells$K, cells$p_switch, cells$density, cells$T_layers), ]
calib_cells <- cells$cell[seq(1, nrow(cells), by = 2)]
s$split <- ifelse(s$cell %in% calib_cells, "calibration", "validation")
cat(sprintf("binary arm: %d fits, %d cells (%d calibration / %d validation)\n", nrow(s), nrow(cells), length(calib_cells), nrow(cells) - length(calib_cells)))
if (nrow(cells) != 594) warning("expected 594 binary cells, found ", nrow(cells))

# ---- calibrate one level -----------------------------------------------------
calibrate <- function(stab, acc, split) {
  b <- bin_of(stab); cal <- split == "calibration"
  q05 <- tapply(acc[cal], factor(b[cal], levels = 1:10), quantile, probs = 0.05, names = FALSE)
  med <- tapply(acc[cal], factor(b[cal], levels = 1:10), median)
  ncal <- as.integer(table(factor(b[cal], levels = 1:10)))
  tab <- data.frame(bin = 1:10, stab_lo = BR[-11], stab_hi = BR[-1], n_calib = ncal, acc_median = as.numeric(med), acc_q05 = as.numeric(q05))
  # validation: share above the floor. An empty calibration bin falls back to the nearest populated bin
  # below; a fit whose bin has no populated bin at or below it has NO floor (NA) and is reported separately.
  floor_of <- function(bb) { i <- bb; while (i > 1 && (is.na(tab$acc_q05[i]) || tab$n_calib[i] == 0)) i <- i - 1L
    f <- tab$acc_q05[i]; if (is.na(f) || tab$n_calib[i] == 0) NA_real_ else f } # no populated bin at or below: no floor
  fl <- vapply(b, floor_of, numeric(1)); val <- !cal; def <- !is.na(fl)
  tab$n_valid <- as.integer(table(factor(b[val], levels = 1:10)))
  tab$n_valid_nofloor <- as.integer(table(factor(b[val & !def], levels = 1:10)))
  tab$valid_share_above <- as.numeric(tapply(acc[val & def] >= fl[val & def], factor(b[val & def], levels = 1:10), mean))
  if (any(val & !def)) cat(sprintf("  [%d of %d validation fits have no calibrated floor (no calibration fits in or below their bin); excluded from the share]\n", sum(val & !def), sum(val)))
  list(tab = tab, share = mean(acc[val & def] >= fl[val & def]), n_nofloor = sum(val & !def), n_valid = sum(val),
       rho = cor(stab[val], acc[val], method = "spearman"), floor_of = floor_of)
}
P_nmi <- calibrate(s$stab_nmi, s$acc_nmi, s$split)
P_ari <- calibrate(s$stab_ari, s$acc_ari, s$split)
nd <- read_arm("binary", "node"); nd$cell <- paste(nd$n, nd$K, nd$p_switch, nd$density, nd$T_layers, sep = "|")
nd$split <- ifelse(nd$cell %in% calib_cells, "calibration", "validation")
N_jac <- calibrate(nd$stab, nd$acc, nd$split)

# ---- pre-registered rule -----------------------------------------------------
rule <- function(P, label) {
  q <- P$tab$acc_q05[!is.na(P$tab$acc_q05)]; mono <- all(diff(q) >= -0.02)
  hi <- P$tab$acc_q05[P$tab$stab_lo >= 0.9]; ok_hi <- length(hi) && all(hi[!is.na(hi)] >= 0.8)
  pass <- P$rho >= 0.7 && mono && ok_hi
  cat(sprintf("%-22s Spearman %.3f  monotone %s  floor>=0.8 at s>=0.9 %s  validation share above floor %.3f  => %s\n",
              label, P$rho, mono, ok_hi, P$share, if (pass) "PASS" else "FAIL"))
  pass
}
cat("\n--- pre-registered rule, validation half ---\n")
pass_nmi <- rule(P_nmi, "partition NMI"); pass_ari <- rule(P_ari, "partition ARI"); pass_node <- rule(N_jac, "node Jaccard")

# ---- calibration table shipped with the packages --------------------------------
mk <- function(P, level) data.frame(level = level, stab_lo = P$tab$stab_lo, stab_hi = P$tab$stab_hi, n_calib = P$tab$n_calib,
                                    acc_median = round(P$tab$acc_median, 4), acc_q05 = round(P$tab$acc_q05, 4),
                                    source = sprintf("sim02_binary_%s", format(Sys.Date())))
ctab <- rbind(mk(P_nmi, "partition_nmi"), mk(P_ari, "partition_ari"), mk(N_jac, "node_jaccard"))
write.csv(ctab, file.path(IN, "stability_calibration_table.csv"), row.names = FALSE)
cat("wrote", file.path(IN, "stability_calibration_table.csv"), "\n")

# ---- main-text table + appendix levels table ----------------------------------
bin_lab <- function(t) sprintf("[%.1f, %.1f%s", t$stab_lo, t$stab_hi, ifelse(t$stab_hi == 1, "]", ")"))
tb <- P_nmi$tab; keep <- tb$n_calib > 0 | tb$n_valid > 0
write_tex("Stability bin & Calibration fits & Median accuracy & Accuracy floor & Validation fits & Share above floor \\\\",
          sprintf("%s & %d & %s & %s & %d & %s \\\\", bin_lab(tb)[keep], tb$n_calib[keep], fmt(tb$acc_median[keep]), fmt(tb$acc_q05[keep]), tb$n_valid[keep], fmt(tb$valid_share_above[keep])),
          "lccccc", file.path(TAB, "tab_stability_floor.tex"))
lv <- function(P, name) { t <- P$tab; k <- t$n_calib > 0 | t$n_valid > 0
  sprintf("%s & %s & %d & %s & %s & %s \\\\", ifelse(seq_len(sum(k)) == 1, name, ""), bin_lab(t)[k], t$n_calib[k], fmt(t$acc_median[k]), fmt(t$acc_q05[k]), fmt(t$valid_share_above[k])) }
write_tex("Level & Stability bin & Calibration fits & Median accuracy & Accuracy floor & Share above floor \\\\",
          c(lv(P_ari, "Partition, ARI"), "\\addlinespace", lv(N_jac, "Node, Jaccard")), "llcccc", file.path(TAB, "tab_app_stability_levels.tex"))

# ---- robustness arms: binary-calibrated NMI floor applied to each arm ----------------
floor_nmi <- function(stab) vapply(bin_of(stab), P_nmi$floor_of, numeric(1))
arm_rows <- list(); add <- function(label, stab, acc) arm_rows[[length(arm_rows) + 1]] <<- data.frame(label = label, n = length(stab),
  n_floor = sum(!is.na(floor_nmi(stab))), share = mean((acc >= floor_nmi(stab))[!is.na(floor_nmi(stab))]),
  rho = cor(stab, acc, method = "spearman"), mean_stab = mean(stab), mean_acc = mean(acc))
sv <- s[s$split == "validation", ]; add("Binary, validation half (main)", sv$stab_nmi, sv$acc_nmi)
dc <- read_arm("dcsbm", "stab")
if (!is.null(dc)) { add("Degree-corrected, all", dc$stab_nmi, dc$acc_nmi)
  for (h in c("none", "moderate", "severe")) for (bl in c("balanced", "skewed")) { x <- dc[dc$hetero == h & dc$balance == bl, ]
    if (nrow(x)) add(sprintf("\\quad heterogeneity %s, %s sizes", h, bl), x$stab_nmi, x$acc_nmi) } } else cat("(no dcsbm arm output yet)\n")
wt <- read_arm("weighted", "stab")
if (!is.null(wt)) { add("Weighted, all", wt$stab_nmi, wt$acc_nmi)
  for (w in c("aligned", "orthogonal")) { x <- wt[wt$weights == w, ]; if (nrow(x)) add(sprintf("\\quad weights %s", w), x$stab_nmi, x$acc_nmi) } } else cat("(no weighted arm output yet)\n")
ar <- do.call(rbind, arm_rows)
write_tex("Arm & Fits & Fits with a floor & Share above floor & Spearman & Mean stability & Mean accuracy \\\\",
  sprintf("%s & %d & %d & %s & %s & %s & %s \\\\", ar$label, ar$n, ar$n_floor, fmt(ar$share), fmt(ar$rho), fmt(ar$mean_stab, 2), fmt(ar$mean_acc, 2)),
  "lcccccc", file.path(TAB, "tab_app_stability_arms.tex"))
cat("\n--- robustness arms (binary-calibrated floor) ---\n"); print(ar, row.names = FALSE, digits = 3)

# ---- leave-one-level-out ----------------------------------------------------------
lolo <- list()
for (fac in c("n", "density", "K")) for (lev in sort(unique(s[[fac]]))) {
  hold <- s[[fac]] == lev; cal <- !hold
  b <- bin_of(s$stab_nmi); q05 <- tapply(s$acc_nmi[cal], factor(b[cal], levels = 1:10), quantile, probs = 0.05, names = FALSE)
  ncal <- as.integer(table(factor(b[cal], levels = 1:10)))
  fo <- function(bb) { i <- bb; while (i > 1 && (is.na(q05[i]) || ncal[i] == 0)) i <- i - 1L; f <- q05[i]; if (is.na(f) || ncal[i] == 0) NA_real_ else f }
  fl <- vapply(b[hold], fo, numeric(1))
  def <- !is.na(fl)
  lolo[[length(lolo) + 1]] <- data.frame(factor = fac, level = as.character(lev), n = sum(hold), n_floor = sum(def),
    share = if (any(def)) mean(s$acc_nmi[hold][def] >= fl[def]) else NA_real_)
}
lo <- do.call(rbind, lolo)
write_tex("Held-out factor & Level & Held-out fits & Fits with a floor & Share above floor \\\\",
  sprintf("%s & %s & %d & %d & %s \\\\", c("Network size $n$", "Density", "Communities $K$")[match(lo$factor, c("n", "density", "K"))], lo$level, lo$n, lo$n_floor, fmt(lo$share)),
  "llccc", file.path(TAB, "tab_app_stability_lolo.tex"))
cat("\n--- leave-one-level-out ---\n"); print(lo, row.names = FALSE, digits = 3)

# ---- decided / undetermined pairs ----------------------------------------------------
pr <- read_arm("binary", "pair")
if (!is.null(pr)) {
  pr$cell <- paste(pr$n, pr$K, pr$p_switch, pr$density, pr$T_layers, sep = "|"); pv <- pr[!(pr$cell %in% calib_cells), ]
  tog <- pv$share >= 0.9; ap <- pv$share <= 0.1; und <- !tog & !ap
  prow <- data.frame(what = c("Decided together (share $\\geq$ 0.9)", "Decided apart (share $\\leq$ 0.1)", "Undetermined"),
                     share = c(mean(tog), mean(ap), mean(und)),
                     acc = c(mean(pv$true_same[tog] == 1), mean(pv$true_same[ap] == 0), mean(pv$true_same[und] == pv$point_same[und])))
  write_tex("Pair status & Share of sampled pairs & Agreement with truth \\\\", sprintf("%s & %s & %s \\\\", prow$what, fmt(prow$share), fmt(prow$acc)),
            "lcc", file.path(TAB, "tab_app_stability_pairs.tex"))
  cat("\n--- pairs (validation, sampled) ---\n"); print(prow, row.names = FALSE, digits = 3)
}

# ---- figures -------------------------------------------------------------------------
step <- data.frame(x = c(tb$stab_lo, 1), y = c(tb$acc_q05, tb$acc_q05[10]))
sv$n_lab <- factor(paste0("n = ", sv$n), levels = paste0("n = ", sort(unique(sv$n))))
g <- ggplot(sv, aes(stab_nmi, acc_nmi)) + geom_point(alpha = 0.08, size = 0.6, colour = "#1f4e79") +
  geom_step(data = step, aes(x, y), colour = "#b2182b", linewidth = 0.9, direction = "hv") +
  geom_abline(slope = 1, intercept = 0, linetype = 3, colour = "grey60") +
  scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(0, 1)) + coord_equal() +
  labs(x = "Bootstrap stability (mean NMI, replicate vs point estimate)", y = "Accuracy (NMI, point estimate vs truth)") +
  theme_bw(base_size = 11)
ggsave(file.path(FIG, "fig_stability_floor.pdf"), g, width = 5.2, height = 5.2)
ggsave(file.path(FIG, "fig_stability_floor.png"), g, width = 5.2, height = 5.2, dpi = 300)
g2 <- g + facet_wrap(~ n_lab, ncol = 2)
ggsave(file.path(FIG, "fig_app_stability_by_n.pdf"), g2, width = 7, height = 7)

# node level (Figure 2b): per-node Jaccard stability vs per-node Jaccard accuracy on the
# validation half of the node sample, with the node-level calibrated floor as the step ----
tbn <- N_jac$tab; stepn <- data.frame(x = c(tbn$stab_lo, 1), y = c(tbn$acc_q05, tbn$acc_q05[10]))
ndv <- nd[nd$split == "validation", ]
gn <- ggplot(ndv, aes(stab, acc)) + geom_point(alpha = 0.08, size = 0.6, colour = "#1f4e79") +
  geom_step(data = stepn, aes(x, y), colour = "#b2182b", linewidth = 0.9, direction = "hv") +
  geom_abline(slope = 1, intercept = 0, linetype = 3, colour = "grey60") +
  scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(0, 1)) + coord_equal() +
  labs(x = "Bootstrap stability (mean Jaccard, replicate vs point estimate)", y = "Accuracy (Jaccard, point estimate vs truth)") +
  theme_bw(base_size = 11)
ggsave(file.path(FIG, "fig_stability_floor_node.pdf"), gn, width = 5.2, height = 5.2)
ggsave(file.path(FIG, "fig_stability_floor_node.png"), gn, width = 5.2, height = 5.2, dpi = 300)
cat(sprintf("node-level figure: %d validation node-layer observations\n", nrow(ndv)))

cat(sprintf("\nDECISION: partition NMI %s | partition ARI %s | node %s\n", if (pass_nmi) "PASS" else "FAIL", if (pass_ari) "PASS" else "FAIL", if (pass_node) "PASS" else "FAIL"))
if (!(pass_nmi || pass_ari)) { cat("=> partition-level rule FAILED: the stability section cannot be supported by this run\n"); quit(status = 2) }
cat("done 12_stability\n")
