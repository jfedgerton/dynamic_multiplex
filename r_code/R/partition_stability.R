#' @title Stability score and calibrated accuracy floor for a tracked partition
#'
#' @description Summarises a \code{\link{bootstrap_multilayer}} result into
#' one reliability report for the recovered (tracked) partition:
#' \enumerate{
#'   \item \strong{Stability score} \eqn{s}: the mean agreement (NMI by
#'     default, or ARI) between each bootstrap replicate's meta-partition and
#'     the point-estimate meta-partition, averaged over layers.
#'   \item \strong{Accuracy floor}: a calibrated lower bound on the accuracy
#'     of the point-estimate partition against the (unknown) truth. In the
#'     package's simulation study, fits were binned by \eqn{s} and the 5th
#'     percentile of accuracy (NMI or ARI to the planted partition) recorded
#'     in each bin on a calibration half of the configurations. On the
#'     held-out half, accuracy exceeded the floor in 95 percent of fits, in
#'     every bin. The floor is therefore read as "with this stability,
#'     partition accuracy was at least this value 95 percent of the time in
#'     calibration".
#'   \item \strong{Node-level stability}: per node, the mean Jaccard overlap
#'     between its replicate community and its point-estimate community, with
#'     a calibrated floor on node-level accuracy that is reported only above
#'     stability 0.9 (below that the floor is uninformative).
#'   \item \strong{Decided pairs}: per layer, node pairs whose bootstrap
#'     co-assignment share is below \code{decided[1]} (decidedly apart) or
#'     above \code{decided[2]} (decidedly together); everything in between is
#'     undetermined. No interval is attached to pairs: pair-level intervals
#'     were found not to be calibratable for ambiguous pairs.
#' }
#' The calibration is within the simulated family (planted-partition networks
#' with node switching, sizes 50-400, 3-10 communities, 5-15 layers). On
#' networks whose structure the planted-partition model does not describe, a
#' stable but wrong partition could receive a floor it does not deserve; the
#' accompanying paper reports the degree-corrected and weighted robustness
#' arms.
#'
#' @param boot_result Output from \code{\link{bootstrap_multilayer}} (version
#'   1.3.0 or later, which records \code{stability_samples}).
#' @param metric \code{"nmi"} (default) or \code{"ari"}: the agreement measure
#'   for the partition-level score and floor.
#' @param decided Length-2 numeric: co-assignment shares at or below the first
#'   value are "apart", at or above the second are "together".
#' @param calibration_table Optional path or data frame overriding the bundled
#'   table (\code{inst/extdata/stability_calibration_table.csv}; columns
#'   \code{level, stab_lo, stab_hi, n_calib, acc_median, acc_q05}).
#'
#' @return A list of class \code{"partition_stability"}:
#'   \describe{
#'     \item{stability}{Partition stability score \eqn{s} (mean over
#'       replicates and layers).}
#'     \item{stability_by_layer}{Per-layer mean agreement.}
#'     \item{stability_mc_se}{Monte Carlo standard error of \eqn{s} from the
#'       replicate spread (about 0.02 at 100 replicates).}
#'     \item{floor}{Calibrated 5th-percentile accuracy for the bin containing
#'       \eqn{s}; \code{floor_median} is the median accuracy in that bin.}
#'     \item{bin}{The stability bin used, as "[lo, hi)".}
#'     \item{node}{Per-layer data frames with node, stability (Jaccard) and
#'       floor (NA below 0.9).}
#'     \item{pairs}{Per-layer integer matrices: 1 = decidedly together, -1 =
#'       decidedly apart, 0 = undetermined; \code{pair_summary} gives the
#'       shares.}
#'     \item{report}{One-sentence plain-language summary.}
#'   }
#'
#' @examples
#' set.seed(123)
#' layers <- lapply(1:3, function(i) {
#'   m <- matrix(rbinom(400, 1, 0.3), nrow = 20)
#'   m <- pmax(m, t(m)); diag(m) <- 0; m
#' })
#' boot <- bootstrap_multilayer(layers, fit_type = "jaccard",
#'                              algorithm = "leiden", n_boot = 10, seed = 123)
#' ps <- partition_stability(boot)
#' ps$stability; ps$floor; ps$report
#'
#' @seealso \code{\link{bootstrap_multilayer}}, \code{\link{community_est}}
#' @export
partition_stability <- function(boot_result, metric = c("nmi", "ari"),
                                decided = c(0.1, 0.9), calibration_table = NULL) {
  metric <- match.arg(metric)
  if (is.null(boot_result$stability_samples)) {
    stop("`boot_result` has no stability_samples; rerun bootstrap_multilayer() ",
         "with dynamicmultiplex >= 1.3.0.", call. = FALSE)
  }
  if (boot_result$n_boot < 2) stop("At least two completed bootstrap replicates are needed.", call. = FALSE)
  S <- boot_result$stability_samples[[metric]]
  n_layers <- ncol(S)

  # partition-level score ----
  by_layer <- colMeans(S)
  per_rep <- rowMeans(S)
  s <- mean(per_rep)
  mc_se <- stats::sd(per_rep) / sqrt(length(per_rep))

  # calibrated floor ----
  tab <- .load_stability_table(calibration_table)
  lvl <- if (metric == "nmi") "partition_nmi" else "partition_ari"
  fl <- .floor_lookup(s, tab[tab$level == lvl, ])

  # node level ----
  node_tab <- tab[tab$level == "node_jaccard", ]
  node <- lapply(seq_len(n_layers), function(t) {
    st <- boot_result$node_jaccard_stability[[t]]
    fr <- vapply(st, function(v) if (v >= 0.9) .floor_lookup(v, node_tab)$floor else NA_real_, numeric(1))
    data.frame(node = seq_along(st), stability = st, floor = fr)
  })

  # decided / undetermined pairs ----
  pairs <- lapply(boot_result$co_assignment, function(P) {
    M <- matrix(0L, nrow(P), ncol(P)); M[P >= decided[2]] <- 1L; M[P <= decided[1]] <- -1L
    diag(M) <- 1L; dimnames(M) <- dimnames(P); M
  })
  up <- lapply(pairs, function(M) M[upper.tri(M)])
  allp <- unlist(up)
  pair_summary <- c(together = mean(allp == 1L), apart = mean(allp == -1L), undetermined = mean(allp == 0L))

  report <- sprintf(paste0(
    "Partition stability %.2f (%s over %d bootstrap replicates and %d layers; MC s.e. %.3f). ",
    "In calibration, fits with stability in %s had accuracy of at least %.2f in 95%% of cases (median %.2f). ",
    "%.0f%% of node pairs are decided (%.0f%% together, %.0f%% apart); %.0f%% are undetermined."),
    s, toupper(metric), boot_result$n_boot, n_layers, mc_se, fl$bin, fl$floor, fl$median,
    100 * (1 - pair_summary[["undetermined"]]), 100 * pair_summary[["together"]],
    100 * pair_summary[["apart"]], 100 * pair_summary[["undetermined"]])

  structure(list(stability = s, stability_by_layer = by_layer, stability_mc_se = mc_se,
                 metric = metric, floor = fl$floor, floor_median = fl$median, bin = fl$bin,
                 node = node, pairs = pairs, pair_summary = pair_summary, report = report,
                 calibration_source = attr(tab, "source")),
            class = "partition_stability")
}

#' @export
print.partition_stability <- function(x, ...) {
  cat(strwrap(x$report, width = 80), sep = "\n")
  invisible(x)
}

# Look up the calibration row whose [stab_lo, stab_hi) contains s.
.floor_lookup <- function(s, rows) {
  if (!nrow(rows)) return(list(floor = NA_real_, median = NA_real_, bin = NA_character_))
  rows <- rows[order(rows$stab_lo), ]
  idx <- findInterval(s, c(rows$stab_lo, 1), rightmost.closed = TRUE)
  idx <- min(max(idx, 1L), nrow(rows))
  # empty bins (no calibration fits) fall back to the nearest populated bin below
  while (idx > 1 && (is.na(rows$acc_q05[idx]) || rows$n_calib[idx] == 0)) idx <- idx - 1L
  list(floor = rows$acc_q05[idx], median = rows$acc_median[idx],
       bin = sprintf("[%.1f, %.1f)", rows$stab_lo[idx], rows$stab_hi[idx]))
}

# Bundled default: inst/extdata/stability_calibration_table.csv, written by
# replication/post/12_stability.R from the stability simulations.
.load_stability_table <- function(calibration_table = NULL) {
  if (is.null(calibration_table)) {
    path <- system.file("extdata", "stability_calibration_table.csv", package = "dynamicmultiplex")
    if (!nzchar(path)) stop("The bundled stability calibration table is missing; pass `calibration_table`.", call. = FALSE)
    tab <- utils::read.csv(path, stringsAsFactors = FALSE); src <- path
  } else if (is.character(calibration_table)) {
    tab <- utils::read.csv(calibration_table, stringsAsFactors = FALSE); src <- calibration_table
  } else {
    tab <- as.data.frame(calibration_table); src <- "user-supplied"
  }
  need <- c("level", "stab_lo", "stab_hi", "n_calib", "acc_median", "acc_q05")
  if (!all(need %in% names(tab))) stop("calibration table needs columns: ", paste(need, collapse = ", "), call. = FALSE)
  attr(tab, "source") <- src
  tab
}
