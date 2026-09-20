#!/usr/bin/env Rscript
# =============================================================================
# replication/post/16_stability.R   (EXPLORATORY: applies the pre-registered
# decision rule to sim/06 output; writes nothing to the manuscript)
#
# Rule (fixed 2026-09-20, before the data existed):
#   PARTITION level passes if, on VALIDATION configurations,
#     (a) Spearman(stab_nmi, acc_nmi) >= 0.7 across simulations, and
#     (b) the calibrated 5th percentile of acc_nmi given stab_nmi (fitted on
#         calibration configurations, 10 stability bins) is non-decreasing
#         and >= 0.8 for stability >= 0.9.
#   NODE level passes under the same two conditions for (stab, acc) per node.
#   K level is reported (coverage of the community-count percentile interval)
#   but is not part of the rule.
# Split: configurations sorted, odd -> calibration, even -> validation.
# Usage: DM_ROOT=... Rscript replication/post/16_stability.R
# =============================================================================
ROOT <- Sys.getenv("DM_ROOT", unset = getwd())
OUT  <- file.path(ROOT, "output", "stability")
ALT  <- file.path(ROOT, "output", "alternatives"); dir.create(ALT, showWarnings = FALSE, recursive = TRUE)
options(width = 160)

s <- do.call(rbind, lapply(list.files(OUT, "^stab_task.*csv$", full.names = TRUE), read.csv))
cat("stability rows:", nrow(s), " tasks:", length(unique(s$task)), "\n")
s$cfg <- paste(s$n, s$K, s$p_switch, s$density, s$T_layers, sep = "|")
cfgs <- sort(unique(s$cfg)); calib <- cfgs[seq(1, length(cfgs), by = 2)]
s$split <- ifelse(s$cfg %in% calib, "calibration", "validation")

evaluate <- function(d, stab, acc, label) {
  cal <- d[d$split == "calibration", ]; val <- d[d$split == "validation", ]
  rho_all <- cor(val[[stab]], val[[acc]], method = "spearman")
  rho_n <- sapply(split(val, val$n), function(x) cor(x[[stab]], x[[acc]], method = "spearman"))
  # calibrated lower bound: 5th percentile of accuracy within stability bin
  br <- seq(0, 1, by = 0.1); cal$bin <- cut(cal[[stab]], br, include.lowest = TRUE, labels = FALSE)
  val$bin <- cut(val[[stab]], br, include.lowest = TRUE, labels = FALSE)
  q05 <- tapply(cal[[acc]], cal$bin, quantile, probs = 0.05); med <- tapply(cal[[acc]], cal$bin, median)
  ncal <- tapply(cal[[acc]], cal$bin, length)
  lb <- q05[as.character(val$bin)]; hit <- mean(val[[acc]] >= lb, na.rm = TRUE)     # validation: share above the bound
  tab <- data.frame(stab_bin = sprintf("[%.1f,%.1f]", br[-11], br[-1]),
                    n_calib = as.integer(ncal[as.character(1:10)]),
                    acc_median = as.numeric(med[as.character(1:10)]),
                    acc_q05 = as.numeric(q05[as.character(1:10)]),
                    n_valid = as.integer(table(factor(val$bin, levels = 1:10))))
  vb <- tapply(val[[acc]] >= lb, val$bin, mean); tab$valid_share_above_q05 <- as.numeric(vb[as.character(1:10)])
  hi <- which(br[-1] > 0.9)                                          # bins with stability >= 0.9
  q05_ok <- !is.na(tab$acc_q05[hi]); lb_hi <- tab$acc_q05[hi][q05_ok]
  mono <- { q <- tab$acc_q05[!is.na(tab$acc_q05)]; all(diff(q) >= -0.02) }   # allow 0.02 sampling noise
  pass_a <- rho_all >= 0.7; pass_b <- mono && length(lb_hi) > 0 && all(lb_hi >= 0.8)
  cat(sprintf("\n=== %s ===\nSpearman(stability, accuracy) validation: %.3f   by n: %s\n", label, rho_all,
              paste(sprintf("n=%s %.2f", names(rho_n), rho_n), collapse = ", ")))
  print(tab, row.names = FALSE, digits = 3)
  cat(sprintf("validation share of accuracy >= calibrated 5th percentile: %.3f (target 0.95)\n", hit))
  cat(sprintf("RULE (a) Spearman >= 0.7: %s   RULE (b) q05 monotone and >= 0.8 at stability >= 0.9: %s   => %s\n",
              pass_a, pass_b, if (pass_a && pass_b) "PASS" else "FAIL"))
  write.csv(tab, file.path(ALT, sprintf("stability_%s.csv", gsub("[^a-z]", "_", tolower(label)))), row.names = FALSE)
  invisible(pass_a && pass_b)
}

cat("\n--- raw summary by design (all sims) ---\n")
print(aggregate(cbind(acc_nmi, stab_nmi, K_cov) ~ n + density, data = s, FUN = function(x) round(mean(x), 3)))
p1 <- evaluate(s, "stab_nmi", "acc_nmi", "PARTITION level, NMI")
p2 <- evaluate(s, "stab_ari", "acc_ari", "PARTITION level, ARI")
cat(sprintf("\nK-level: community-count percentile interval coverage = %.3f (nominal 0.95); mean |K_hat - K_true| = %.3f\n",
            mean(s$K_cov), mean(abs(s$K_hat_mean - s$K_true_mean))))

nd <- do.call(rbind, lapply(list.files(OUT, "^node_task.*csv$", full.names = TRUE), read.csv))
nd$cfg <- paste(nd$n, nd$K, nd$p_switch, nd$density, nd$T_layers, sep = "|")
nd$split <- ifelse(nd$cfg %in% calib, "calibration", "validation")
cat("\nnode rows:", nrow(nd), "\n")
p3 <- evaluate(nd, "stab", "acc", "NODE level, Jaccard")

cat("\n================ DECISION ================\n")
cat("partition (NMI):", if (p1) "PASS" else "FAIL", " partition (ARI):", if (p2) "PASS" else "FAIL",
    " node:", if (p3) "PASS" else "FAIL", "\n")
cat(if (p1 || p2 || p3) "=> at least one level recovers: keep the uncertainty section, reframed around stability.\n"
    else "=> nothing recovers: cut the uncertainty section; keep the bootstrap as a stability diagnostic only.\n")
