#!/usr/bin/env Rscript
# =============================================================================
# replication/post/15_alternatives.R   (EXPLORATORY: decision support, no
# manuscript tables are written by this script)
#
# Compares ways of turning the bootstrap co-assignment share phat into a
# statement about the fresh-draw propensity p*. Two parts:
#
# PART A -- from the existing grid (sim/02), always runs:
#   A1 one pooled lookup vs a second table fitted inside the reliability gate
#      (mean Wilson width < 0.05 & n >= 100): width and coverage by phat decile
#   A2 calibration curve E[p* | phat] and P(same planted community | phat):
#      is phat informative as a corrected POINT estimate even where the
#      interval is wide?
#
# PART B -- from sim/05_alt_bootstrap (runs only if output/alt_bootstrap exists):
#   B1 three bootstrap schemes (sbm = current, dcsbm, rewire): Wilson coverage
#      and calibrated width by phat decile, same networks
#   B2 feature-conditioned calibration on the per-pair sample: for each
#      observable feature, conditional quantiles of p* within
#      (phat decile x feature bin), fitted on calibration tasks, evaluated on
#      validation tasks. Which observable, if any, narrows the middle?
#
# Everything prints to stdout; the summary CSVs go to output/alternatives/.
# Usage: DM_ROOT=... Rscript replication/post/15_alternatives.R
# =============================================================================
ROOT <- Sys.getenv("DM_ROOT", unset = getwd())
OUT  <- file.path(ROOT, "output")
ALT  <- file.path(OUT, "alternatives"); dir.create(ALT, showWarnings = FALSE, recursive = TRUE)
ALPHA <- 0.05; Z <- qnorm(1 - ALPHA / 2); B_BOOT <- 100
options(width = 160)

# --------------------------- shared helpers ----------------------------------
read_joint <- function(sub, pat = "^joint_task.*csv$", variant = NULL) {
  fs <- list.files(file.path(OUT, sub), pat, full.names = TRUE)
  if (!length(fs)) return(NULL)
  J <- list(); des <- list(); k <- 0
  for (f in fs) {
    x <- read.csv(f)
    if (!is.null(variant)) x <- x[x$variant == variant, ]
    NB <- max(x$bin_phat); m <- matrix(0, NB, NB); m[cbind(x$bin_phat, x$bin_pstar)] <- x$count
    k <- k + 1; J[[k]] <- m
    des[[k]] <- x[1, intersect(names(x), c("task", "n", "K", "p_switch", "density", "T_", "T_layers", "algorithm", "fit_type"))]
  }
  list(J = J, des = do.call(rbind, des), NB = NB)
}
fit_lut <- function(M, edges, NB) {
  lo <- rep(NA_real_, NB); hi <- rep(NA_real_, NB)
  for (b in seq_len(NB)) { tot <- sum(M[b, ]); if (tot == 0) next
    cdf <- cumsum(M[b, ]) / tot
    lo[b] <- edges[which(cdf > ALPHA / 2)[1]]; hi[b] <- edges[which(cdf >= 1 - ALPHA / 2)[1] + 1] }
  thin <- is.na(lo) | rowSums(M) < 1000
  if (any(thin)) { ok <- which(!thin); for (b in which(thin)) { nb <- ok[which.min(abs(ok - b))]
    lo[b] <- min(lo[nb], edges[b]); hi[b] <- max(hi[nb], edges[b + 1]) } }
  data.frame(bin = seq_len(NB), phat_lo = edges[-(NB + 1)], phat_hi = edges[-1],
             lower = cummax(lo), upper = rev(cummin(rev(hi))), n_pairs = rowSums(M))
}
cover_mask <- function(lut, edges, NB) outer(seq_len(NB), seq_len(NB),
  function(b, j) edges[j] >= lut$lower[b] & edges[j + 1] <= lut$upper[b])
wilson_width <- function(p, B) { den <- 1 + Z^2 / B; 2 * (Z / den) * sqrt(p * (1 - p) / B + Z^2 / (4 * B^2)) }
decile_table <- function(MV, lut, C, edges, NB, wb = NULL) {
  mid <- (edges[-1] + edges[-(NB + 1)]) / 2; dec <- pmin(10L, 1L + floor(mid * 10))
  do.call(rbind, lapply(1:10, function(k) { b <- which(dec == k); tot <- rowSums(MV)[b]
    if (sum(tot) == 0) return(NULL); w <- tot / sum(tot)
    r <- data.frame(decile = sprintf("[%.1f,%.1f)", (k - 1) / 10, k / 10),
               cal_width = sum(w * (lut$upper[b] - lut$lower[b])),
               cal_cov = sum(rowSums(MV * C)[b]) / sum(tot),
               wil_width = sum(w * wilson_width(mid[b], B_BOOT)), pairs = sum(tot))
    if (!is.null(wb)) r$wil_cov <- sum(wb$n_wilson_covered[wb$bin %in% b]) / sum(wb$n_pairs[wb$bin %in% b])
    r }))
}
split_tasks <- function(cov) {
  cov$cfg <- do.call(paste, c(cov[intersect(names(cov), c("n","K","p_switch","p_in","p_out","density","T_layers","weights","resample"))], sep = "|"))
  cfgs <- sort(unique(cov$cfg)); calib <- cfgs[seq(1, length(cfgs), by = 2)]
  ts <- unique(cov[, c("task", "cfg")]); ts$split <- ifelse(ts$cfg %in% calib, "calibration", "validation"); ts
}

# =============================================================================
# PART A: existing grid
# =============================================================================
cat("################ PART A: pooled vs gated lookup, calibration curve ################\n")
cov <- do.call(rbind, lapply(list.files(file.path(OUT, "coverage_grid"), "^cov_task.*csv$", full.names = TRUE), read.csv))
ts  <- split_tasks(cov)
tw  <- aggregate(width_P_mean ~ task + n, data = cov, FUN = mean)
tw$gated <- tw$width_P_mean < 0.05 & tw$n >= 100
G <- read_joint("coverage_grid"); NB <- G$NB; edges <- seq(0, 1, length.out = NB + 1)
G$des$split <- ts$split[match(G$des$task, ts$task)]
G$des$gated <- tw$gated[match(G$des$task, tw$task)]
cb <- do.call(rbind, lapply(list.files(file.path(OUT, "coverage_grid"), "^calib_task.*csv$", full.names = TRUE), read.csv))
cb$split <- ts$split[match(cb$task, ts$task)]; cb$gated <- tw$gated[match(cb$task, tw$task)]
sumJ <- function(idx, J) Reduce(`+`, J[idx])

for (nm in c("pooled", "gated")) {
  keep <- if (nm == "pooled") rep(TRUE, nrow(G$des)) else G$des$gated
  if (sum(keep & G$des$split == "calibration") < 5 || sum(keep & G$des$split == "validation") < 5) {
    cat(sprintf("\n--- A1 %s table: too few tasks, skipped\n", nm)); next }
  M  <- sumJ(which(keep & G$des$split == "calibration"), G$J)
  MV <- sumJ(which(keep & G$des$split == "validation"), G$J)
  lut <- fit_lut(M, edges, NB); C <- cover_mask(lut, edges, NB); lut$width <- lut$upper - lut$lower
  wb <- aggregate(cbind(n_wilson_covered, n_pairs) ~ bin, data = cb[keep[match(cb$task, G$des$task)] & cb$split == "validation", ], FUN = sum)
  cat(sprintf("\n--- A1 %s table: calibration tasks %d, validation tasks %d; validation coverage %.4f, mean width %.4f\n",
              nm, sum(keep & G$des$split == "calibration"), sum(keep & G$des$split == "validation"),
              sum(MV * C) / sum(MV), sum(rowSums(MV) * lut$width) / sum(MV)))
  print(decile_table(MV, lut, C, edges, NB, wb), row.names = FALSE, digits = 3)
  write.csv(lut, file.path(ALT, sprintf("lut_%s.csv", nm)), row.names = FALSE)
}

# A2: calibration curve on the validation split, pooled and gated
cat("\n--- A2 calibration curve: E[p* | phat] and P(same planted | phat), validation split\n")
for (nm in c("pooled", "gated")) {
  keep <- if (nm == "pooled") rep(TRUE, nrow(G$des)) else G$des$gated
  if (sum(keep & G$des$split == "validation") < 5) { cat(sprintf("\n[%s] too few tasks, skipped\n", nm)); next }
  MV <- sumJ(which(keep & G$des$split == "validation"), G$J)
  mid <- (edges[-1] + edges[-(NB + 1)]) / 2
  Ep <- as.numeric(MV %*% mid) / rowSums(MV)                     # E[p* | phat bin]
  sdp <- sqrt(pmax(0, as.numeric(MV %*% mid^2) / rowSums(MV) - Ep^2))
  cbv <- cb[keep[match(cb$task, G$des$task)] & cb$split == "validation", ]
  pt <- aggregate(cbind(n_true, n_pairs) ~ bin, data = cbv, FUN = sum)
  curve <- data.frame(bin = seq_len(NB), phat_mid = mid, E_pstar = Ep, sd_pstar = sdp,
                      p_same_planted = pt$n_true[match(seq_len(NB), pt$bin)] / pt$n_pairs[match(seq_len(NB), pt$bin)],
                      n_pairs = rowSums(MV))
  dec <- pmin(10L, 1L + floor(mid * 10))
  agg <- do.call(rbind, lapply(1:10, function(k) { b <- which(dec == k); w <- curve$n_pairs[b]; if (sum(w) == 0) return(NULL)
    data.frame(decile = sprintf("[%.1f,%.1f)", (k - 1) / 10, k / 10), E_pstar = sum(w * curve$E_pstar[b]) / sum(w),
               sd_pstar = sum(w * curve$sd_pstar[b]) / sum(w), p_same_planted = sum(w * curve$p_same_planted[b], na.rm = TRUE) / sum(w),
               pairs = sum(w)) }))
  cat(sprintf("\n[%s]\n", nm)); print(agg, row.names = FALSE, digits = 3)
  write.csv(curve, file.path(ALT, sprintf("calibration_curve_%s.csv", nm)), row.names = FALSE)
}

# =============================================================================
# PART B: alternative bootstraps and feature conditioning
# =============================================================================
if (!dir.exists(file.path(OUT, "alt_bootstrap")) ||
    !length(list.files(file.path(OUT, "alt_bootstrap"), "^alt_joint"))) {
  cat("\n(no output/alt_bootstrap: PART B skipped; run sim/05 first)\n"); quit(save = "no")
}
cat("\n################ PART B: bootstrap variants, feature conditioning ################\n")
acov <- do.call(rbind, lapply(list.files(file.path(OUT, "alt_bootstrap"), "^alt_cov_task.*csv$", full.names = TRUE), read.csv))
cat("alt tasks:", length(unique(acov$task)), " sims:", nrow(acov) / 3, "\n")
acov$cfg <- paste(acov$n, acov$K, acov$p_switch, acov$density, acov$T_layers, sep = "|")
cfgs <- sort(unique(acov$cfg)); calib <- cfgs[seq(1, length(cfgs), by = 2)]
ats <- unique(acov[, c("task", "cfg")]); ats$split <- ifelse(ats$cfg %in% calib, "calibration", "validation")
acal <- do.call(rbind, lapply(list.files(file.path(OUT, "alt_bootstrap"), "^alt_calib_task.*csv$", full.names = TRUE), read.csv))
acal$split <- ats$split[match(acal$task, ats$task)]

cat("\n--- B1 per-sim Wilson coverage and width by bootstrap variant (all tasks)\n")
print(aggregate(cbind(cov_P_mean, width_P_mean) ~ variant, data = acov, FUN = mean), digits = 4)
for (v in c("sbm", "dcsbm", "rewire")) {
  A <- read_joint("alt_bootstrap", "^alt_joint_task.*csv$", variant = v)
  A$des$split <- ats$split[match(A$des$task, ats$task)]
  M  <- sumJ(which(A$des$split == "calibration"), A$J); MV <- sumJ(which(A$des$split == "validation"), A$J)
  lut <- fit_lut(M, edges, NB); C <- cover_mask(lut, edges, NB)
  wb <- aggregate(cbind(n_wilson_covered, n_pairs) ~ bin, data = acal[acal$variant == v & acal$split == "validation", ], FUN = sum)
  cat(sprintf("\n[%s] validation: calibrated cov %.4f width %.4f | Wilson cov %.4f\n", v,
              sum(MV * C) / sum(MV), sum(rowSums(MV) * (lut$upper - lut$lower)) / sum(MV),
              sum(wb$n_wilson_covered) / sum(wb$n_pairs)))
  print(decile_table(MV, lut, C, edges, NB, wb), row.names = FALSE, digits = 3)
}

# B2 feature-conditioned calibration on the pair sample
cat("\n--- B2 feature-conditioned quantiles of p* (pair sample; mid-range oversampled, so read by decile)\n")
pr <- do.call(rbind, lapply(list.files(file.path(OUT, "alt_bootstrap"), "^alt_pairs_task.*csv$", full.names = TRUE), read.csv))
pr$split <- ats$split[match(pr$task, ats$task)]
pr$dec <- pmin(10L, 1L + floor(pr$phat * 10))
pr$ent_max <- pmax(pr$ent_i, pr$ent_j); pr$deg_min <- pmin(pr$deg_i, pr$deg_j)
pr$csize_min <- pmin(pr$csize_i, pr$csize_j); pr$deg_rel <- pr$deg_min / pmax(1, pr$n)
cat("pair rows:", nrow(pr), " calibration:", sum(pr$split == "calibration"), " validation:", sum(pr$split == "validation"), "\n")
# how well does each phat variant track p* on the same pairs?
cat("correlation with p* on sampled pairs: sbm", sprintf("%.3f", cor(pr$phat, pr$pstar)),
    " dcsbm", sprintf("%.3f", cor(pr$phat_dcsbm, pr$pstar)), " rewire", sprintf("%.3f", cor(pr$phat_rewire, pr$pstar)), "\n")

tercile <- function(x, ref) cut(x, breaks = unique(quantile(ref, c(0, 1/3, 2/3, 1))), include.lowest = TRUE, labels = FALSE)
features <- list(
  none          = function(d, ref) rep(1L, nrow(d)),
  gate          = function(d, ref) as.integer(d$mean_width < 0.05 & d$n >= 100),
  n             = function(d, ref) match(d$n, c(50, 100, 200, 400)),
  same_detected = function(d, ref) d$same_detected,
  entropy_max   = function(d, ref) tercile(d$ent_max, ref$ent_max),
  degree_min    = function(d, ref) tercile(d$deg_rel, ref$deg_rel),
  csize_min     = function(d, ref) tercile(d$csize_min, ref$csize_min),
  gate_x_same   = function(d, ref) 2L * as.integer(d$mean_width < 0.05 & d$n >= 100) + d$same_detected,
  gate_x_entropy= function(d, ref) 3L * as.integer(d$mean_width < 0.05 & d$n >= 100) + tercile(d$ent_max, ref$ent_max),
  same_x_entropy= function(d, ref) 3L * d$same_detected + tercile(d$ent_max, ref$ent_max))
cal <- pr[pr$split == "calibration", ]; val <- pr[pr$split == "validation", ]
res <- list()
for (fn in names(features)) {
  cal$fb <- features[[fn]](cal, cal); val$fb <- features[[fn]](val, cal)
  # conditional 2.5/97.5% quantiles of p* in each (decile, feature bin), min 200 calibration pairs
  q <- aggregate(pstar ~ dec + fb, data = cal, FUN = function(x) c(lo = unname(quantile(x, 0.025)), hi = unname(quantile(x, 0.975)), n = length(x)))
  q <- data.frame(dec = q$dec, fb = q$fb, lo = q$pstar[, "lo"], hi = q$pstar[, "hi"], n = q$pstar[, "n"])
  qd <- aggregate(pstar ~ dec, data = cal, FUN = function(x) c(lo = unname(quantile(x, 0.025)), hi = unname(quantile(x, 0.975))))
  qd <- data.frame(dec = qd$dec, lo0 = qd$pstar[, "lo"], hi0 = qd$pstar[, "hi"])
  key <- paste(val$dec, val$fb); qk <- paste(q$dec, q$fb)
  lo <- q$lo[match(key, qk)]; hi <- q$hi[match(key, qk)]; nq <- q$n[match(key, qk)]
  fallback <- is.na(lo) | nq < 200                         # thin cell: fall back to decile-only quantiles
  lo[fallback] <- qd$lo0[match(val$dec[fallback], qd$dec)]; hi[fallback] <- qd$hi0[match(val$dec[fallback], qd$dec)]
  val$covered <- val$pstar >= lo & val$pstar <= hi; val$width <- hi - lo
  s <- aggregate(cbind(width, covered) ~ dec, data = val, FUN = mean)
  s$feature <- fn; s$n_val <- as.integer(table(val$dec)[as.character(s$dec)])
  res[[fn]] <- s
}
R <- do.call(rbind, res)
wide_w <- reshape(R[, c("feature", "dec", "width")], idvar = "feature", timevar = "dec", direction = "wide")
wide_c <- reshape(R[, c("feature", "dec", "covered")], idvar = "feature", timevar = "dec", direction = "wide")
names(wide_w) <- sub("width.", "d", names(wide_w)); names(wide_c) <- sub("covered.", "d", names(wide_c))
cat("\nvalidation WIDTH of the 95% interval by phat decile (d1 = [0,0.1) ... d10 = [0.9,1]):\n"); print(wide_w, row.names = FALSE, digits = 2)
cat("\nvalidation COVERAGE by phat decile:\n"); print(wide_c, row.names = FALSE, digits = 3)
mid_w <- aggregate(width ~ feature, data = R[R$dec %in% 3:8, ], FUN = mean)
mid_c <- aggregate(covered ~ feature, data = R[R$dec %in% 3:8, ], FUN = mean)
summ <- merge(mid_w, mid_c); names(summ) <- c("feature", "mid_width", "mid_coverage")
summ <- summ[order(summ$mid_width), ]
cat("\nmid-range (phat 0.2-0.8) mean width and coverage by conditioning feature, best first:\n"); print(summ, row.names = FALSE, digits = 3)
write.csv(R, file.path(ALT, "feature_conditioning.csv"), row.names = FALSE)
cat("\ndone 15_alternatives\n")
