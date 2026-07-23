#' @title Wilson confidence intervals for node-pair co-assignment
#'
#' @description For every pair of nodes in every layer, computes a Wilson
#' score interval for the co-clustering propensity: the probability that the
#' fitted community detection procedure places the two nodes in the same
#' community when the data are perturbed. The point estimate is the
#' co-assignment probability from \code{\link{bootstrap_multilayer}} (the
#' share of bootstrap replicates in which the pair was co-assigned), and the
#' interval treats the \code{n_boot} replicates as binomial draws.
#'
#' Because co-assignment is label-invariant (it never compares community
#' labels across replicates, only whether two nodes sit together), it avoids
#' the label-switching problem that makes per-node membership intervals
#' ill-defined.
#'
#' @section Warning:
#' These intervals quantify the stability of the detection procedure, not
#' the probability that two nodes truly share a community. The Wilson
#' interval is exact for the binomial sampling of bootstrap replicates;
#' its coverage of the fresh-data co-clustering propensity was evaluated in
#' the package's simulation study (see the package NEWS for the release in
#' which validation results were incorporated). Interpret cautiously on
#' networks with fewer than 100 nodes, where community detection itself is
#' unstable.
#'
#' @param boot_result Output from \code{\link{bootstrap_multilayer}}.
#'
#' @param alpha Significance level (default 0.05 for 95 percent intervals).
#'
#' @return A list with one element per layer. Each element is a list with
#'   components:
#'   \describe{
#'     \item{estimate}{n x n matrix of co-assignment probabilities.}
#'     \item{lower}{n x n matrix of Wilson lower bounds.}
#'     \item{upper}{n x n matrix of Wilson upper bounds.}
#'   }
#'   Diagonals are 1 by construction.
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
#' pci <- co_assignment_ci(boot, alpha = 0.05)
#' pci[[1]]$estimate[1:4, 1:4]
#' pci[[1]]$lower[1:4, 1:4]
#'
#' @seealso \code{\link{community_est}} for community-count point estimates and
#' node stability summaries.
#'
#' @export

co_assignment_ci <- function(boot_result, alpha = 0.05) {

  # check that completed bootstrap replicates are present ----
  if (boot_result$n_boot == 0) {
    stop("No completed bootstrap replicates.", call. = FALSE)
  }

  # Wilson score interval ingredients ----
  b <- boot_result$n_boot
  z <- stats::qnorm(1 - alpha / 2)
  z2 <- z^2

  # compute per-layer Wilson bounds from the co-assignment matrices ----
  layer_cis <- lapply(boot_result$co_assignment, function(phat) {

    ## Wilson center and half-width ----
    denom <- 1 + z2 / b
    center <- (phat + z2 / (2 * b)) / denom
    half <- z * sqrt(phat * (1 - phat) / b + z2 / (4 * b^2)) / denom

    ## clamp to the unit interval (matrix first so dims are preserved) ----
    lower <- pmax(center - half, 0)
    upper <- pmin(center + half, 1)

    ## diagonals are degenerate: a node is always co-assigned with itself ----
    diag(lower) <- 1
    diag(upper) <- 1

    ## compile layer result ----
    list(estimate = phat, lower = lower, upper = upper)
  })

  # return per-layer confidence intervals ----
  return(layer_cis)
}
