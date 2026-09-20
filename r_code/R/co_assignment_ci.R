#' @title Descriptive Monte Carlo interval for node-pair co-assignment
#'
#' @description For every pair of nodes in every layer, returns the bootstrap
#' co-assignment share (the fraction of replicates in which the pair landed in
#' the same meta-community) with a Wilson score interval treating the
#' replicates as binomial draws.
#'
#' \strong{Read this interval as a diagnostic, not as a calibrated confidence
#' interval.} In the package's simulation study its width shrinks with the
#' number of replicates while its coverage of the fresh-data co-assignment
#' propensity was near nominal only for pairs whose share is close to 0 or 1
#' and fell to 0.03-0.15 for ambiguous pairs (share between 0.3 and 0.8).
#' No conditioning or alternative bootstrap repaired this, so versions
#' from 1.3.0 no longer offer a "calibrated" pair interval. The validated
#' reliability product is \code{\link{partition_stability}}: a stability score
#' for the tracked partition with a calibrated lower bound on its accuracy,
#' plus a decided / undetermined flag for each pair.
#'
#' @param boot_result Output from \code{\link{bootstrap_multilayer}}.
#'
#' @param alpha Significance level; \code{alpha = 0.05} gives 95 percent
#'   Wilson intervals.
#'
#' @param method \code{"wilson"}. \code{"calibrated"} is accepted only to
#'   raise an informative error pointing to \code{partition_stability()}.
#'
#' @param calibration_table Ignored; retained for backward compatibility.
#'
#' @return A list with one element per layer. Each element is a list with
#'   components:
#'   \describe{
#'     \item{estimate}{n x n matrix of co-assignment shares.}
#'     \item{lower}{n x n matrix of interval lower bounds.}
#'     \item{upper}{n x n matrix of interval upper bounds.}
#'   }
#'   Diagonals are 1 by construction. The list carries the attribute
#'   \code{"method"}.
#'
#' @examples
#' set.seed(123)
#' layers <- lapply(1:3, function(i) {
#'   m <- matrix(rbinom(64, 1, 0.35), nrow = 8)
#'   m <- pmax(m, t(m))
#'   diag(m) <- 0
#'   m
#' })
#' boot <- bootstrap_multilayer(
#'   layers,
#'   fit_type = "jaccard",
#'   algorithm = "louvain",
#'   n_boot = 5,
#'   seed = 123
#' )
#' pci <- co_assignment_ci(boot)
#' pci[[1]]$estimate[1:4, 1:4]
#'
#' @seealso \code{\link{partition_stability}} for the calibrated reliability
#' score; \code{\link{community_est}} for community-count summaries.
#'
#' @export
co_assignment_ci <- function(boot_result, alpha = 0.05,
                             method = c("wilson", "calibrated"),
                             calibration_table = NULL) {

  method <- match.arg(method)
  if (method == "calibrated") {
    stop("method = \"calibrated\" was removed in dynamicmultiplex 1.3.0: pair-level ",
         "intervals could not be calibrated for ambiguous pairs. Use ",
         "partition_stability() for the calibrated reliability score and the ",
         "decided / undetermined pair flags, or method = \"wilson\" for the ",
         "descriptive Monte Carlo interval.", call. = FALSE)
  }

  # check that completed bootstrap replicates are present ----
  if (boot_result$n_boot == 0) {
    stop("No completed bootstrap replicates.", call. = FALSE)
  }

  # Wilson score interval ingredients ----
  b <- boot_result$n_boot
  z <- stats::qnorm(1 - alpha / 2)
  z2 <- z^2

  layer_cis <- lapply(boot_result$co_assignment, function(phat) {
    denom  <- 1 + z2 / b
    center <- (phat + z2 / (2 * b)) / denom
    half   <- z * sqrt(phat * (1 - phat) / b + z2 / (4 * b^2)) / denom
    lower  <- pmax(center - half, 0)
    upper  <- pmin(center + half, 1)
    diag(lower) <- 1
    diag(upper) <- 1
    list(estimate = phat, lower = lower, upper = upper)
  })

  attr(layer_cis, "method") <- method
  return(layer_cis)
}
