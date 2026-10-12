#' @title Fit multilayer communities with Hungarian snapshot matching
#'
#' @description Detects communities independently in each layer and then tracks
#' them across time by matching community labels between consecutive layers with
#' the Hungarian (optimal linear-sum assignment) algorithm applied to the
#' community overlap (contingency) matrix. This is a two-stage
#' \emph{snapshot-and-match} tracker: unlike the coupling-based fits
#' (\code{\link{fit_multilayer_jaccard}}, \code{\link{fit_multilayer_overlap}},
#' \code{\link{fit_multilayer_identity_ties}}), communities are found per layer
#' and aligned \emph{post hoc}, not jointly optimised on a coupled supra-graph.
#' It provides the common "independent detection + optimal label matching"
#' baseline (e.g. per-layer Leiden matched with the Hungarian algorithm) as a
#' first-class method with the same interface as the other fits.
#'
#' @param layers List of \code{igraph} objects or square adjacency matrices.
#'
#' @param algorithm Community algorithm: \code{"leiden"} (default) or \code{"louvain"}.
#'
#' @param resolution_parameter Leiden resolution parameter.
#'
#' @param directed Logical; if \code{TRUE}, build directed graphs from adjacency
#' matrices (collapsed to undirected for detection, as elsewhere in the
#' package).
#'
#' @param objective Objective for the per-layer detection: \code{"modularity"}
#' or \code{"cpm"} (Leiden only). Defaults to \code{"cpm"} for directed layers
#' and \code{"modularity"} otherwise.
#'
#' @param seed Integer seed for reproducible community detection (default 123, as in
#'   the Python package; pass NULL for unseeded detection).
#' When supplied, the global RNG state is saved, the RNG is seeded for the
#' duration of the call, and the previous state is restored on exit, so the
#' caller's random number stream (e.g. bootstrap resampling) is unaffected.
#' Defaults to `NULL` (detection inherits the caller's RNG stream, matching
#' previous behavior; call `set.seed()` beforehand for reproducibility).
#'
#' @param allow_unequal_nodes Logical; if `TRUE`, layers may contain
#' different node sets. Layers must then carry vertex names (or row names on
#' the adjacency matrices); community overlap is counted on the nodes present
#' in both layers of each consecutive pair. Defaults to `FALSE`.
#'
#' @return A list of class \code{"multilayer_community_fit"} with components:
#'   \describe{
#'     \item{layer_communities}{Per-layer community detection (each with
#'       \code{membership} and \code{communities}), detected independently.}
#'     \item{meta_communities}{The tracked partition: one integer vector per
#'       layer where labels are aligned across consecutive layers by Hungarian
#'       matching on community overlap. A community with no positive-overlap
#'       match is given a new label (a birth); labels not carried forward are
#'       treated as deaths. See \code{\link{extract_meta_membership}}.}
#'     \item{interlayer_ties}{\code{NULL} - no coupled supra-graph is built.}
#'     \item{method}{\code{"hungarian"}.}
#'     \item{layer_links}{Sequential adjacent-layer links (the chain the
#'       matcher walks).}
#'   }
#'
#' @details Matching is sequential: layer \eqn{t} is aligned to layer
#' \eqn{t-1}. Nodes are matched across layers by vertex name when names are
#' present (by position otherwise), so the order in which vertices are stored
#' does not affect the result. For each adjacent pair the community overlap matrix is completed
#' to a square cost matrix and \code{\link[clue]{solve_LSAP}} maximises total
#' overlap; unmatched current communities receive fresh labels. The result is a
#' consistently-labelled sequence suitable for \code{\link{extract_meta_membership}},
#' \code{\link{bootstrap_multilayer}}, and the plotting helpers.
#'
#' @seealso \code{\link{fit_multilayer_jaccard}},
#'   \code{\link{fit_multilayer_overlap}},
#'   \code{\link{fit_multilayer_identity_ties}} for coupling-based tracking.
#'
#' @examples
#' set.seed(123)
#' layers <- lapply(1:3, function(i) {
#'   m <- matrix(rbinom(64, 1, 0.35), nrow = 8)
#'   m <- pmax(m, t(m)); diag(m) <- 0; m
#' })
#' fit <- fit_multilayer_hungarian(layers, algorithm = "leiden")
#' extract_meta_membership(fit)
#'
#' @export

fit_multilayer_hungarian <- function(
    layers,
    algorithm = c("leiden", "louvain"),
    resolution_parameter = 1,
    directed = FALSE,
    objective = NULL,
    seed = 123,
    allow_unequal_nodes = FALSE
  ) {

  algorithm <- match.arg(algorithm)

  # Scoped seed for reproducible detection (restores caller RNG state) ----
  if (!is.null(seed)) {
    rng_state <- save_rng_state()
    on.exit(restore_rng_state(rng_state), add = TRUE)
    set.seed(seed)
  }

  # Prepare graph layers ----
  graph_layers <- prepare_multilayer_graphs(
    layers,
    directed = directed,
    require_same_nodes = !allow_unequal_nodes
  )
  links <- make_layer_links(length(graph_layers), NULL)

  # Stage 1: independent per-layer community detection ----
  fit <- fit_layer_communities(
    graph_layers,
    algorithm = algorithm,
    resolution_parameter = resolution_parameter,
    directed = directed,
    objective = objective
  )

  # Stage 2: Hungarian label matching across consecutive layers ----
  # memberships are keyed by vertex name (position when unnamed), so layers
  # whose vertices are stored in different orders, or that hold different
  # node sets, are matched actor by actor ----
  memberships <- lapply(seq_along(fit), function(t) {
    m <- as.integer(fit[[t]]$membership)
    nm <- igraph::V(graph_layers[[t]])$name
    names(m) <- if (is.null(nm)) as.character(seq_along(m)) else as.character(nm)
    m
  })
  meta <- match_communities_hungarian(memberships)

  # Compile multilayer community fit object ----
  multilayer_hungarian <- structure(
    list(
      algorithm = algorithm,
      layer_communities = fit,
      meta_communities = meta,
      meta_ids = sort(unique(unlist(meta))),
      layer_links = links,
      interlayer_ties = NULL,
      method = "hungarian",
      directed = directed
    ),
    class = "multilayer_community_fit"
  )

  return(multilayer_hungarian)
}


# Internal: align a sequence of per-layer memberships by Hungarian matching of
# consecutive layers on the community overlap matrix (maximise total overlap).
# Each membership is an integer vector named by node id; overlap is counted on
# the node ids present in both layers. Returns unnamed vectors in each layer's
# own node order.
match_communities_hungarian <- function(memberships) {

  keyed <- lapply(memberships, function(m) {
    k <- as.integer(m)
    names(k) <- if (is.null(names(m))) as.character(seq_along(m)) else names(m)
    k
  })

  out <- vector("list", length(keyed))
  out[[1L]] <- keyed[[1L]]
  next_free <- max(c(out[[1L]], 0L)) + 1L

  if (length(keyed) >= 2L) {
    for (t in 2:length(keyed)) {

      prev <- out[[t - 1L]]
      cur  <- keyed[[t]]
      cl   <- sort(unique(cur))
      pl   <- sort(unique(prev))

      # community overlap (contingency) matrix on shared node ids ----
      shared <- intersect(names(cur), names(prev))
      overlap <- matrix(0L, nrow = length(cl), ncol = length(pl))
      if (length(shared) > 0L) {
        cs <- cur[shared]; ps <- prev[shared]
        for (i in seq_along(cl)) {
          for (j in seq_along(pl)) {
            overlap[i, j] <- sum(cs == cl[i] & ps == pl[j])
          }
        }
      }

      # pad to a square cost matrix and solve the assignment ----
      d <- max(nrow(overlap), ncol(overlap))
      square <- matrix(0L, d, d)
      square[seq_len(nrow(overlap)), seq_len(ncol(overlap))] <- overlap
      assignment <- as.integer(clue::solve_LSAP(max(square) - square))

      # map current labels to previous labels; births get fresh labels ----
      relabel <- rep(NA_integer_, length(cl))
      for (i in seq_along(cl)) {
        col <- assignment[i]
        if (col <= length(pl) && overlap[i, col] > 0L) relabel[i] <- pl[col]
      }
      for (i in seq_along(cl)) {
        if (is.na(relabel[i])) {
          relabel[i]  <- next_free
          next_free   <- next_free + 1L
        }
      }
      next_free <- max(next_free, max(relabel) + 1L)

      names(relabel) <- as.character(cl)
      out[[t]] <- stats::setNames(unname(relabel[as.character(cur)]), names(cur))
    }
  }

  lapply(out, function(x) as.integer(unname(x)))
}
