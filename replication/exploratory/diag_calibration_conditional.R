#!/usr/bin/env Rscript
# =============================================================================
# exploratory/diag_calibration_conditional.R
# DIAGNOSTIC, not part of run_all.sh. Question: how much narrower does the
# calibrated co-assignment interval get when the lookup is conditioned on
# something the analyst observes -- network size n, or the reliability gate
# (mean Wilson width < 0.05 and n >= 100, i.e. a polarised co-assignment
# matrix) -- instead of being pooled over the whole design?
#
# For each conditioning subset: fit L(phat), U(phat) on the subset's
# calibration tasks, evaluate on its validation tasks. Prints, by phat
# decile: calibrated width, calibrated coverage, Wilson coverage, pairs.
#
# Usage: DM_ROOT=... Rscript replication/exploratory/diag_calibration_conditional.R
# =============================================================================
ROOT <- Sys.getenv("DM_ROOT", unset = getwd())
OUT  <- file.path(ROOT, "output")
ALPHA <- 0.05; Z <- qnorm(1 - ALPHA / 2); B_BOOT <- 100

cov <- do.call(rbind, lapply(list.files(file.path(OUT, "coverage_grid"), "^cov_task.*csv$", full.names = TRUE), read.csv))
cov$cfg <- paste(cov$n, cov$K, cov$p_switch, cov$p_in, cov$p_out, cov$density, cov$T_layers, cov$weights, cov$resample, sep = "|")
cfgs <- sort(unique(cov$cfg)); calib_cfg <- cfgs[seq(1, length(cfgs), by = 2)]
tk <- aggregate(cbind(width = width_P_mean, cov_P = cov_P_mean) ~ task + cfg + n + K + density + algorithm + fit_type, data = cov, FUN = mean)
tk$gated_share <- tapply(cov$width_P_mean < 0.05, cov$task, mean)[as.character(tk$task)]
tk$split <- ifelse(tk$cfg %in% calib_cfg, "calibration", "validation")
cat("tasks:", nrow(tk), "\n")

fs <- list.files(file.path(OUT, "coverage_grid"), "^joint_task.*csv$", full.names = TRUE)
J <- vector("list", length(fs)); tid <- integer(length(fs))
for (i in seq_along(fs)) {
  x <- read.csv(fs[i]); NB <- max(x$bin_phat)
  m <- matrix(0, NB, NB); m[cbind(x$bin_phat, x$bin_pstar)] <- x$count
  J[[i]] <- m; tid[i] <- x$task[1]
}
names(J) <- tid
edges <- seq(0, 1, length.out = NB + 1)
cb <- do.call(rbind, lapply(list.files(file.path(OUT, "coverage_grid"), "^calib_task.*csv$", full.names = TRUE), read.csv))

fit_lut <- function(M) {
  lo <- rep(NA_real_, NB); hi <- rep(NA_real_, NB)
  for (b in seq_len(NB)) { tot <- sum(M[b, ]); if (tot == 0) next
    cdf <- cumsum(M[b, ]) / tot
    lo[b] <- edges[which(cdf > ALPHA / 2)[1]]; hi[b] <- edges[which(cdf >= 1 - ALPHA / 2)[1] + 1] }
  thin <- is.na(lo) | rowSums(M) < 1000
  if (any(thin)) { ok <- which(!thin); for (b in which(thin)) { nb <- ok[which.min(abs(ok - b))]
    lo[b] <- min(lo[nb], edges[b]); hi[b] <- max(hi[nb], edges[b + 1]) } }
  lo <- cummax(lo); hi <- rev(cummin(rev(hi)))
  list(lo = lo, hi = hi)
}
wilson_width <- function(p, B) { den <- 1 + Z^2 / B; 2 * (Z / den) * sqrt(p * (1 - p) / B + Z^2 / (4 * B^2)) }

evaluate <- function(label, keep) {
  tc <- tk$task[keep & tk$split == "calibration"]; tv <- tk$task[keep & tk$split == "validation"]
  if (length(tc) < 10 || length(tv) < 10) { cat(sprintf("\n=== %s: too few tasks (%d/%d)\n", label, length(tc), length(tv))); return(invisible()) }
  M  <- Reduce(`+`, J[as.character(tc)]); MV <- Reduce(`+`, J[as.character(tv)])
  L  <- fit_lut(M)
  C  <- outer(seq_len(NB), seq_len(NB), function(b, j) edges[j] >= L$lo[b] & edges[j + 1] <= L$hi[b])
  wb <- cb[cb$task %in% tv, ]; wb <- aggregate(cbind(n_wilson_covered, n_pairs) ~ bin, data = wb, FUN = sum)
  mid <- (edges[-1] + edges[-(NB + 1)]) / 2; dec <- pmin(10L, 1L + floor(mid * 10))
  rows <- lapply(1:10, function(k) { b <- which(dec == k); tot <- rowSums(MV)[b]; if (sum(tot) == 0) return(NULL)
    w <- tot / sum(tot)
    data.frame(decile = sprintf("[%.1f,%.1f)", (k - 1) / 10, k / 10),
               cal_width = sum(w * (L$hi[b] - L$lo[b])), cal_cov = sum(rowSums((MV * C))[b]) / sum(tot),
               wil_width = sum(w * wilson_width(mid[b], B_BOOT)),
               wil_cov = sum(wb$n_wilson_covered[wb$bin %in% b]) / sum(wb$n_pairs[wb$bin %in% b]),
               pairs = sum(tot)) })
  r <- do.call(rbind, rows)
  cat(sprintf("\n=== %s: calibration tasks %d, validation tasks %d\n", label, length(tc), length(tv)))
  cat(sprintf("pooled over phat: calibrated cov %.4f  width %.4f | Wilson cov %.4f  width %.4f\n",
              sum(MV * C) / sum(MV), sum(rowSums(MV) * (L$hi - L$lo)) / sum(MV),
              sum(wb$n_wilson_covered) / sum(wb$n_pairs), sum(rowSums(MV) * wilson_width(mid, B_BOOT)) / sum(MV)))
  print(r, row.names = FALSE, digits = 3)
}

evaluate("ALL tasks (pooled)", rep(TRUE, nrow(tk)))
for (nn in c(50, 100, 200, 400)) evaluate(sprintf("n = %d", nn), tk$n == nn)
evaluate("GATED tasks: mean width < 0.05 & n >= 100", tk$width < 0.05 & tk$n >= 100)
evaluate("GATED, n = 100", tk$width < 0.05 & tk$n == 100)
evaluate("GATED, n >= 200", tk$width < 0.05 & tk$n >= 200)
for (dd in c("weak", "default", "strong")) evaluate(sprintf("separation = %s (all n)", dd), tk$density == dd)
evaluate("GATED & strong separation", tk$width < 0.05 & tk$n >= 100 & tk$density == "strong")
evaluate("GATED & weak separation", tk$width < 0.05 & tk$n >= 100 & tk$density == "weak")
cat("\ndone diag\n")
