#' @title Bootstrap stability score of a tracked partition
#'
#' @description Summarises a \code{\link{bootstrap_multilayer}} result into the
#' stability score \eqn{s}: the mean agreement (NMI by default, or ARI) between
#' each bootstrap replicate's meta-partition and the point-estimate
#' meta-partition, averaged over layers and replicates. \eqn{s} is a point
#' estimate with Monte Carlo standard error \code{stability_mc_se}. It has no
#' validated mapping to accuracy: use it only to compare two fits of the same
#' network (for example DynMux against multislice), never as an estimate of
#' how accurate either fit is. Since 1.4.0 the package reports no interval or
#' accuracy floor.
#'
#' @param boot_result Output from \code{\link{bootstrap_multilayer}} (version
#'   1.3.0 or later, which records \code{stability_samples}).
#' @param metric \code{"nmi"} (default) or \code{"ari"}.
#'
#' @return A list of class \code{"partition_stability"}:
#'   \describe{
#'     \item{stability}{Stability score \eqn{s} (mean over replicates and layers).}
#'     \item{stability_by_layer}{Per-layer mean agreement.}
#'     \item{stability_mc_se}{Monte Carlo standard error of \eqn{s} from the
#'       replicate spread.}
#'     \item{metric}{The agreement measure used.}
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
#' ps$stability; ps$report
#'
#' @seealso \code{\link{bootstrap_multilayer}}
#' @export
partition_stability <- function(boot_result, metric = c("nmi", "ari")) {
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

  report <- sprintf(paste0(
    "Partition stability %.2f (%s over %d bootstrap replicates and %d layers; MC s.e. %.3f). ",
    "Use it only to compare fits of the same network."),
    s, toupper(metric), boot_result$n_boot, n_layers, mc_se)

  structure(list(stability = s, stability_by_layer = by_layer, stability_mc_se = mc_se,
                 metric = metric, report = report),
            class = "partition_stability")
}

#' @export
print.partition_stability <- function(x, ...) {
  cat(strwrap(x$report, width = 80), sep = "\n")
  invisible(x)
}
