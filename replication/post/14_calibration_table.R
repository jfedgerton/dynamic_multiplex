#!/usr/bin/env Rscript
# =============================================================================
# replication/post/14_calibration_table.R
# Calibrated co-assignment interval: fit the lookup table, validate it out of
# sample, and write the appendix tables that document it.
#
# Why. The Wilson interval treats the B bootstrap replicates as a binomial
# sample, so its width scales as 1/sqrt(B): it measures Monte Carlo error in
# the bootstrap, not uncertainty about the estimand p* (the fresh-draw
# co-assignment propensity). Its conditional coverage is therefore uneven
# across the range of phat, and raising B would make it worse, not better.
#
# What. From the simulations' joint counts of (phat bin, p* bin), pooled over
# the CALIBRATION configurations, take for every phat bin the conditional
# 2.5% and 97.5% quantiles of p*. That pair [L(phat), U(phat)] is the
# interval. Its width does not depend on B and its conditional coverage
# given phat is >= 0.95 on the calibration split by construction; the
# VALIDATION split is the out-of-sample check reported in the manuscript.
#
# Split: identical to post/12 (configurations sorted, odd -> calibration,
# even -> validation), carried over by task id from the cov_task files.
#
# Inputs:  $DM_ROOT/output/coverage_grid/{cov,calib,joint}_task*.csv  (sim/02)
#          $DM_ROOT/output/coverage_valued/joint_task*.csv            (sim/03, optional)
#          $DM_ROOT/output/coverage_misspec/joint_task*.csv           (sim/04, optional)
# Outputs:
#   output/calibration/coassign_calibration_table.csv
#       one row per phat bin: bin, phat_lo, phat_hi, lower, upper, n_pairs
#       -> bundled into r_code (inst/extdata) and python_code (package data)
#   tables/tab_calibrated_interval.tex   by phat decile: Wilson vs calibrated
#                                        width and validation coverage
#   tables/tab_calibrated_by_design.tex  validation coverage by n, K, density
#   tables/tab_calibrated_arms.tex       binary lookup applied to the weighted
#                                        and DC-SBM arms
#   figures/fig_calibrated_interval.pdf  L(phat), U(phat) and the Wilson band
#
# Usage:  DM_ROOT=/path/to/dynamic_multiplex Rscript replication/post/14_calibration_table.R
# =============================================================================
set.seed(123)
suppressMessages({ library(ggplot2) })

ROOT <- Sys.getenv("DM_ROOT", unset = getwd())
OUT  <- file.path(ROOT, "output")
CAL  <- file.path(OUT, "calibration")
TAB  <- file.path(ROOT, "manuscript", "tables")
FIG  <- file.path(ROOT, "manuscript", "figures")
dir.create(CAL, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
ALPHA   <- 0.05
B_BOOT  <- 100L                                   # as in sim/02-04
Z       <- qnorm(1 - ALPHA / 2)

write_tex <- function(header, body, align, path) {
  stopifnot(length(body) > 0, nchar(align) > 0)
  writeLines(c(sprintf("\\begin{tabular}{%s}", align),
               "\\toprule", header, "\\midrule",
               body,
               "\\bottomrule", "\\end{tabular}"), path)
  cat("wrote", basename(path), "(", length(body), "rows )\n")
}
commafmt <- function(x) format(x, big.mark = ",", trim = TRUE)
readall <- function(sub, pat) {
  fs <- list.files(file.path(OUT, sub), pat, full.names = TRUE)
  if (!length(fs)) return(NULL)
  cat(sub, pat, ": files", length(fs), "\n")
  do.call(rbind, lapply(fs, read.csv))
}

# =============================================================================
# 1. Read, and carry the calibration/validation split over by task id
# =============================================================================
cov <- readall("coverage_grid", "^cov_task.*csv$")
if (is.null(cov)) stop("No coverage output -- run sim/02 first.", call. = FALSE)
cov$cfg <- paste(cov$n, cov$K, cov$p_switch, cov$p_in, cov$p_out, cov$density,
                 cov$T_layers, cov$weights, cov$resample, sep = "|")
cfgs <- sort(unique(cov$cfg))
calib_cfg <- cfgs[seq(1, length(cfgs), by = 2)]
task_split <- unique(cov[, c("task", "cfg")])
stopifnot(!anyDuplicated(task_split$task))          # one configuration per task
task_split$split <- ifelse(task_split$cfg %in% calib_cfg, "calibration", "validation")
cat("tasks:", nrow(task_split), " calibration:", sum(task_split$split == "calibration"),
    " validation:", sum(task_split$split == "validation"), "\n")

# Joint files are 2,500 rows per task (3,564 tasks for the grid): reading them
# into one data frame needs several GB and gets the process killed on a login
# node. Instead every file is reduced on read to a 50 x 50 count matrix and one
# row of design columns; all downstream sums are matrix operations.
read_joint <- function(sub) {
  fs <- list.files(file.path(OUT, sub), "^joint_task.*csv$", full.names = TRUE)
  if (!length(fs)) return(NULL)
  cat(sub, "joint files:", length(fs), "\n")
  des <- vector("list", length(fs)); J <- vector("list", length(fs))
  for (i in seq_along(fs)) {
    x <- read.csv(fs[i])
    if (i == 1) {
      stopifnot(all(c("task", "n", "K", "density", "algorithm", "bin_phat", "bin_pstar", "count") %in% names(x)))
      NB <- max(x$bin_phat); stopifnot(NB == max(x$bin_pstar), nrow(x) == NB * NB)
    }
    stopifnot(length(unique(x$task)) == 1, nrow(x) == NB * NB)
    m <- matrix(0, NB, NB); m[cbind(x$bin_phat, x$bin_pstar)] <- x$count
    J[[i]]   <- m
    des[[i]] <- x[1, setdiff(names(x), c("bin_phat", "bin_pstar", "count"))]
  }
  list(des = do.call(rbind, des), J = J, NB = NB)
}
G <- read_joint("coverage_grid")
if (is.null(G)) stop("No joint files -- run sim/02 first.", call. = FALSE)
N_BINS <- G$NB; stopifnot(N_BINS >= 20)
G$des$split <- task_split$split[match(G$des$task, task_split$task)]
stopifnot(!anyNA(G$des$split))
cat("bins:", N_BINS, " total pairs:", commafmt(sum(vapply(G$J, sum, numeric(1)))), "\n")
sumJ <- function(idx, JJ = G$J) Reduce(`+`, JJ[idx], accumulate = FALSE)

cb <- readall("coverage_grid", "^calib_task.*csv$")      # Wilson coverage by phat bin (50 rows/task)
stopifnot(all(c("task", "bin", "n_pairs", "n_wilson_covered") %in% names(cb)))
cb$split <- task_split$split[match(cb$task, task_split$task)]
stopifnot(!anyNA(cb$split))

# =============================================================================
# 2. Fit the lookup on the calibration split
# Joint counts M[b, j] = # pairs with phat in bin b and p* in bin j.
# For each b: L = lower edge of the p* bin at which the conditional CDF first
# exceeds alpha/2 (mass strictly below L <= alpha/2); U = upper edge of the
# bin at which it first reaches 1 - alpha/2. Mass in [L, U] >= 1 - alpha.
# =============================================================================
M <- sumJ(which(G$des$split == "calibration"))
edges <- seq(0, 1, length.out = N_BINS + 1)

lut <- data.frame(bin = seq_len(N_BINS), phat_lo = edges[-(N_BINS + 1)], phat_hi = edges[-1],
                  lower = NA_real_, upper = NA_real_, n_pairs = rowSums(M))
for (b in seq_len(N_BINS)) {
  tot <- sum(M[b, ])
  if (tot == 0) next                                 # empty phat bin: filled below
  cdf <- cumsum(M[b, ]) / tot
  jL  <- which(cdf > ALPHA / 2)[1]                   # first bin with mass above alpha/2
  jU  <- which(cdf >= 1 - ALPHA / 2)[1]
  lut$lower[b] <- edges[jL]
  lut$upper[b] <- edges[jU + 1]
}
# Empty or near-empty phat bins (fewer than 1,000 calibration pairs) borrow
# the nearest populated bin's bounds, extended to be monotone so the interval
# never shrinks as phat moves away from the data.
thin <- is.na(lut$lower) | lut$n_pairs < 1000
if (any(thin)) {
  ok <- which(!thin)
  for (b in which(thin)) {
    nb <- ok[which.min(abs(ok - b))]
    lut$lower[b] <- min(lut$lower[nb], lut$phat_lo[b])
    lut$upper[b] <- max(lut$upper[nb], lut$phat_hi[b])
  }
  cat("thin phat bins filled from neighbours:", paste(which(thin), collapse = ","), "\n")
}
# Enforce monotone bounds in phat (isotonic, pool-adjacent-violators) so the
# table is a proper interval function of phat.
lut$lower <- cummax(lut$lower)
lut$upper <- rev(cummin(rev(lut$upper)))
stopifnot(all(lut$lower <= lut$upper), all(lut$lower >= 0), all(lut$upper <= 1))
lut$width <- lut$upper - lut$lower
write.csv(lut, file.path(CAL, "coassign_calibration_table.csv"), row.names = FALSE)
cat("wrote coassign_calibration_table.csv\n"); print(lut)

# =============================================================================
# 3. Validate: conditional coverage by phat bin on each split
# A pair (phat bin b, p* bin j) is covered iff edges[j] >= lower[b] and
# edges[j+1] <= upper[b]; i.e. the whole p* bin lies inside [L, U]. This is
# the conservative reading of the binned data. COVER is that N x N mask.
# =============================================================================
COVER <- outer(seq_len(N_BINS), seq_len(N_BINS),
               function(b, j) edges[j] >= lut$lower[b] & edges[j + 1] <= lut$upper[b])
# per-task covered / total pair counts -> a small data frame with design columns
G$des$cov <- vapply(G$J, function(m) sum(m * COVER), numeric(1))
G$des$tot <- vapply(G$J, sum, numeric(1))
cov_by <- function(x, by) {
  key <- do.call(paste, c(x[by], sep = "\r"))
  cv  <- tapply(x$cov, key, sum); tt <- tapply(x$tot, key, sum)
  a   <- unique(x[by]); rownames(a) <- NULL
  k   <- do.call(paste, c(a, sep = "\r"))
  a$cov <- as.numeric(cv[k]); a$tot <- as.numeric(tt[k])
  a <- a[a$tot > 0, , drop = FALSE]
  a$coverage <- a$cov / a$tot
  a <- a[do.call(order, unname(a[by])), , drop = FALSE]
  rownames(a) <- NULL
  a
}
pool <- cov_by(G$des, "split")
print(pool)
stopifnot(pool$coverage[pool$split == "calibration"] >= 1 - ALPHA - 1e-9)   # by construction

# Wilson conditional coverage by phat bin, from calib files (same split)
wb <- aggregate(cbind(n_wilson_covered, n_pairs) ~ bin + split, data = cb, FUN = sum)
wb <- wb[wb$n_pairs > 0, ]
wb$wilson_cov <- wb$n_wilson_covered / wb$n_pairs
# Wilson width at the bin midpoint with B = 100 (the width is a function of
# phat and B only)
wilson_width <- function(p, B) {
  den <- 1 + Z^2 / B
  2 * (Z / den) * sqrt(p * (1 - p) / B + Z^2 / (4 * B^2))
}
# calibrated coverage by phat bin on the validation split
MV <- sumJ(which(G$des$split == "validation"))
cal_bin <- data.frame(bin = seq_len(N_BINS), split = "validation",
                      cov = rowSums(MV * COVER), tot = rowSums(MV))
cal_bin <- cal_bin[cal_bin$tot > 0, ]
cal_bin$coverage <- cal_bin$cov / cal_bin$tot
v <- merge(cal_bin, wb[wb$split == "validation", ], by = c("bin", "split"))
v <- v[order(v$bin), ]
v$mid <- (edges[v$bin] + edges[v$bin + 1]) / 2
v$wilson_width <- wilson_width(v$mid, B_BOOT)
v$cal_width    <- lut$width[v$bin]

# --- appendix table by phat decile (pairs-weighted means of bin rows) -----
v$dec <- pmin(10L, 1L + floor(v$mid * 10))
dec <- do.call(rbind, lapply(sort(unique(v$dec)), function(k) {
  s <- v[v$dec == k, ]
  w <- s$tot / sum(s$tot)
  data.frame(dec = k, lo = (k - 1) / 10, hi = k / 10,
             wilson_width = sum(w * s$wilson_width), wilson_cov = sum(w * s$wilson_cov),
             cal_width = sum(w * s$cal_width), cal_cov = sum(w * s$coverage),
             n_pairs = sum(s$tot))
}))
print(dec)
write_tex(
  header = paste("$\\hat p$ range & \\multicolumn{2}{c}{Wilson ($B = 100$)} &",
                 "\\multicolumn{2}{c}{Calibrated} & Pairs \\\\",
                 "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}",
                 "& Width & Coverage & Width & Coverage & \\\\"),
  body   = sprintf("[%.1f, %.1f%s & %.3f & %.3f & %.3f & %.3f & %s \\\\",
                   dec$lo, dec$hi, ifelse(dec$dec == 10, "]", ")"),
                   dec$wilson_width, dec$wilson_cov, dec$cal_width, dec$cal_cov,
                   commafmt(dec$n_pairs)),
  align  = "lccccr",
  path   = file.path(TAB, "tab_calibrated_interval.tex"))
cat("validation, pooled over phat: Wilson",
    sprintf("%.4f", sum(v$wilson_cov * v$tot) / sum(v$tot)),
    " calibrated", sprintf("%.4f", sum(v$coverage * v$tot) / sum(v$tot)),
    " mean width Wilson", sprintf("%.4f", sum(v$wilson_width * v$tot) / sum(v$tot)),
    " calibrated", sprintf("%.4f", sum(v$cal_width * v$tot) / sum(v$tot)), "\n")

# --- validation coverage by design cell -----------------------------------
jv <- G$des[G$des$split == "validation", ]
by_n   <- cov_by(jv, "n");        by_n$dim   <- "$n$";     by_n$level   <- as.character(by_n$n)
by_K   <- cov_by(jv, "K");        by_K$dim   <- "$K$";     by_K$level   <- as.character(by_K$K)
by_d   <- cov_by(jv, "density");  by_d$dim   <- "Separation"
by_d$level <- c(weak = "Weak", default = "Moderate", strong = "Strong")[by_d$density]
by_d   <- by_d[order(factor(by_d$density, levels = c("weak", "default", "strong"))), ]
by_a   <- cov_by(jv, "algorithm"); by_a$dim  <- "Algorithm"
by_a$level <- c(leiden = "Leiden", louvain = "Louvain")[by_a$algorithm]
des <- rbind(by_n[, c("dim", "level", "coverage", "tot")], by_K[, c("dim", "level", "coverage", "tot")],
             by_d[, c("dim", "level", "coverage", "tot")], by_a[, c("dim", "level", "coverage", "tot")])
print(des)
# worst design cell: the number to quote as the floor of conditional coverage
cell <- cov_by(jv, c("n", "K", "density", "algorithm"))
cat("validation design cells:", nrow(cell), " min coverage:", sprintf("%.4f", min(cell$coverage)),
    " max:", sprintf("%.4f", max(cell$coverage)), "\n")
show <- c(TRUE, des$dim[-1] != des$dim[-nrow(des)])
write_tex(
  header = "Design dimension & Level & Calibrated coverage & Pairs \\\\",
  body   = sprintf("%s & %s & %.3f & %s \\\\", ifelse(show, des$dim, ""), des$level,
                   des$coverage, commafmt(des$tot)),
  align  = "llcr",
  path   = file.path(TAB, "tab_calibrated_by_design.tex"))

# =============================================================================
# 4. Robustness: apply the binary-arm lookup to the other arms, unchanged
# =============================================================================
arms <- list()
for (arm in c("coverage_valued", "coverage_misspec")) {
  A <- read_joint(arm)
  if (is.null(A)) { cat("skip:", arm, "\n"); next }
  stopifnot(A$NB == N_BINS)
  MA <- sumJ(seq_along(A$J), A$J)
  arms[[arm]] <- data.frame(arm = c(coverage_valued = "Weighted networks",
                                    coverage_misspec = "Degree-corrected generating model")[arm],
                            coverage = sum(MA * COVER) / sum(MA), n_pairs = sum(MA))
}
if (length(arms)) {
  arms <- rbind(data.frame(arm = "Binary (validation split)",
                           coverage = pool$coverage[pool$split == "validation"],
                           n_pairs = pool$tot[pool$split == "validation"]),
                do.call(rbind, arms))
  print(arms)
  write_tex(
    header = "Arm & Calibrated coverage & Pairs \\\\",
    body   = sprintf("%s & %.3f & %s \\\\", arms$arm, arms$coverage, commafmt(arms$n_pairs)),
    align  = "lcr",
    path   = file.path(TAB, "tab_calibrated_arms.tex"))
}

# =============================================================================
# 5. Figure: the interval as a function of phat, against the Wilson band
# =============================================================================
pf <- data.frame(mid = (lut$phat_lo + lut$phat_hi) / 2, lower = lut$lower, upper = lut$upper)
pf$w_lo <- pmax(0, pf$mid - wilson_width(pf$mid, B_BOOT) / 2)
pf$w_hi <- pmin(1, pf$mid + wilson_width(pf$mid, B_BOOT) / 2)
p <- ggplot(pf, aes(mid)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = "Calibrated"), alpha = 0.35) +
  geom_ribbon(aes(ymin = w_lo, ymax = w_hi, fill = "Wilson (B = 100)"), alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = c("Calibrated" = "#1b9e77", "Wilson (B = 100)" = "#d95f02"), name = NULL) +
  scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(0, 1)) +
  coord_equal() +
  labs(x = expression("Bootstrap co-assignment share " * hat(p)),
       y = expression("Interval for the co-assignment propensity " * p^"*")) +
  theme_bw(base_size = 9) + theme(legend.position = "bottom")
ggsave(file.path(FIG, "fig_calibrated_interval.pdf"), p, width = 4.6, height = 4.9)
ggsave(file.path(FIG, "fig_calibrated_interval.png"), p, width = 4.6, height = 4.9, dpi = 300)
cat("wrote fig_calibrated_interval.pdf/.png\n")

cat("done 14_calibration_table\n")
