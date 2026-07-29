#' @title Extract cross-layer meta-community membership
#'
#' @description Returns the tracked (cross-layer) community assignment produced
#' by the second-stage detection: for each layer, a vector giving every node's
#' meta-community. Unlike the per-layer \code{layer_communities} (detected
#' independently, ignoring the coupling), the meta-communities are the partition
#' that reflects the interlayer ties and any custom \code{layer_links}. This is
#' the membership that \code{\link{bootstrap_multilayer}} validates.
#'
#' @param fit A fit object from one of the \code{fit_multilayer_*} functions.
#'
#' @return A list with one integer vector per layer giving each node's
#'   meta-community assignment (node order).
#'
#' @examples
#' set.seed(123)
#' layers <- lapply(1:3, function(i) {
#'   m <- matrix(rbinom(64, 1, 0.35), nrow = 8)
#'   m <- pmax(m, t(m)); diag(m) <- 0; m
#' })
#' fit <- fit_multilayer_jaccard(layers, algorithm = "leiden")
#' extract_meta_membership(fit)
#'
#' @export
extract_meta_membership <- function(fit) {
  if (is.null(fit$meta_communities)) {
    stop("`fit` has no meta_communities; refit with a current ",
         "fit_multilayer_* function.", call. = FALSE)
  }
  fit$meta_communities
}
