#' @title Confidence intervals for node-pair co-assignment
#'
#' @description For every pair of nodes in every layer, computes an interval
#' for the co-clustering propensity: the probability that the fitted
#' community detection procedure places the two nodes in the same community
#' when the data are perturbed. The point estimate is the co-assignment
#' probability from \code{\link{bootstrap_multilayer}} (the share of
#' bootstrap replicates in which the pair was co-assigned).
#'
#' Two interval constructions are available.
#' \describe{
#'   \item{\code{"calibrated"} (default)}{The interval is read from a lookup
#'   table fitted in the package's simulation study: for each bin of the
#'   bootstrap co-assignment share \eqn{\hat p}, the table stores the 2.5 and
#'   97.5 percent conditional quantiles of the fresh-data co-assignment
#'   propensity \eqn{p^*}. Its width does not depend on \code{n_boot} and its
#'   conditional coverage given \eqn{\hat p} was validated out of sample on
#'   the simulation design (binary and weighted stochastic block models,
#'   50--400 nodes, 3--10 communities, 5--15 layers).}
#'   \item{\code{"wilson"}}{The Wilson score interval that treats the
#'   \code{n_boot} replicates as binomial draws. Its width scales as
#'   \eqn{1/\sqrt{n_{boot}}}, so it measures Monte Carlo error in the bootstrap
#'   rather than uncertainty about \eqn{p^*}; it is retained for comparison
#'   and for the coverage study reported in the paper.}
#' }
#'
#' Because co-assignment is label-invariant (it never compares community
#' labels across replicates, only whether two nodes sit together), it avoids
#' the label-switching problem that makes per-node membership intervals
#' ill-defined.
#'
#' @section Warning:
#' These intervals quantify the stability of the detection procedure, not
#' the probability that two nodes truly share a community. The calibrated
#' table was fitted on simulated networks; networks far outside the
#' simulation design (very small, very sparse, or with extreme degree
#' heterogeneity) are extrapolations. Interpret cautiously on networks with
#' fewer than 100 nodes, where community detection itself is unstable.
#'
#' @param boot_result Output from \code{\link{bootstrap_multilayer}}.
#'
#' @param alpha Significance level (default 0.05 for 95 percent intervals).
#'   Only \code{0.05} is available for \code{method = "calibrated"}, which is
#'   the level the bundled lookup table was fitted at.
#'
#' @param method Either \code{"calibrated"} (default) or \code{"wilson"}.
#'
#' @param calibration_table Optional replacement lookup table for
#'   \code{method = "calibrated"}: a data frame with columns \code{phat_lo},
#'   \code{phat_hi}, \code{lower}, \code{upper} (bins must partition
#'   \eqn{[0, 1]}), or a path to a CSV with those columns. Defaults to the
#'   table shipped in \code{inst/extdata/coassign_calibration_table.csv},
#'   produced by \code{replication/post/14_calibration_table.R}.
#'
#' @return A list with one element per layer. Each element is a list with
#'   components:
#'   \describe{
#'     \item{estimate}{n x n matrix of co-assignment probabilities.}
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
#' pci <- co_assignment_ci(boot, method = "wilson")
#' pci[[1]]$estimate[1:4, 1:4]
#' pci[[1]]$lower[1:4, 1:4]
#'
#' @seealso \code{\link{community_est}} for community-count point estimates and
#' node stability summaries.
#'
#' @export

co_assignment_ci <- function(boot_result, alpha = 0.05,
                             method = c("calibrated", "wilson"),
                             calibration_table = NULL) {

  method <- match.arg(method)

  # check that completed bootstrap replicates are present ----
  if (boot_result$n_boot == 0) {
    stop("No completed bootstrap replicates.", call. = FALSE)
  }

  if (method == "wilson") {

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

  } else {

    # calibrated interval: conditional quantiles of p* given phat ----
    if (!isTRUE(all.equal(alpha, 0.05))) {
      stop("method = \"calibrated\" is only available for alpha = 0.05, ",
           "the level the bundled lookup table was fitted at.", call. = FALSE)
    }
    tab <- .load_calibration_table(calibration_table)

    layer_cis <- lapply(boot_result$co_assignment, function(phat) {
      # bin index: phat in [phat_lo, phat_hi); phat == 1 falls in the last bin
      idx <- findInterval(phat, c(tab$phat_lo, 1), rightmost.closed = TRUE)
      idx <- pmin(pmax(idx, 1L), nrow(tab))
      lower <- matrix(tab$lower[idx], nrow(phat), ncol(phat))
      upper <- matrix(tab$upper[idx], nrow(phat), ncol(phat))
      dimnames(lower) <- dimnames(phat)
      dimnames(upper) <- dimnames(phat)
      diag(lower) <- 1
      diag(upper) <- 1
      list(estimate = phat, lower = lower, upper = upper)
    })
  }

  attr(layer_cis, "method") <- method
  return(layer_cis)
}

# Read and validate the calibrated-interval lookup table.
# Bundled default: inst/extdata/coassign_calibration_table.csv, written by
# replication/post/14_calibration_table.R from the coverage simulations.
.load_calibration_table <- function(calibration_table = NULL) {
  if (is.null(calibration_table)) {
    path <- system.file("extdata", "coassign_calibration_table.csv",
                        package = "dynamicmultiplex")
    if (!nzchar(path)) {
      stop("The bundled calibration table is missing. Run ",
           "replication/post/14_calibration_table.R and copy ",
           "output/calibration/coassign_calibration_table.csv to ",
           "r_code/inst/extdata/, or pass `calibration_table`, or use ",
           "method = \"wilson\".", call. = FALSE)
    }
    tab <- utils::read.csv(path, stringsAsFactors = FALSE)
  } else if (is.character(calibration_table)) {
    tab <- utils::read.csv(calibration_table, stringsAsFactors = FALSE)
  } else {
    tab <- as.data.frame(calibration_table)
  }
  need <- c("phat_lo", "phat_hi", "lower", "upper")
  if (!all(need %in% names(tab))) {
    stop("calibration_table needs columns: ", paste(need, collapse = ", "), call. = FALSE)
  }
  tab <- tab[order(tab$phat_lo), ]
  if (!isTRUE(all.equal(tab$phat_lo[1], 0)) ||
      !isTRUE(all.equal(tab$phat_hi[nrow(tab)], 1)) ||
      !isTRUE(all.equal(tab$phat_hi[-nrow(tab)], tab$phat_lo[-1]))) {
    stop("calibration_table bins must partition [0, 1].", call. = FALSE)
  }
  if (any(tab$lower > tab$upper) || any(tab$lower < 0) || any(tab$upper > 1)) {
    stop("calibration_table bounds must satisfy 0 <= lower <= upper <= 1.", call. = FALSE)
  }
  tab
}
