#' @title Fit multilayer communities and interlayer node-strength weighted
#' Jaccard ties
#'
#' @description Runs Louvain or Leiden community detection for each layer and
#' creates interlayer ties between communities in selected layer pairs using
#' node-strength weighted Jaccard similarity.
#'
#' @param layers List of `igraph` objects or square adjacency matrices.
#'
#' @param algorithm Community algorithm: `"louvain"` or `"leiden"`.
#'
#' @param layer_links Optional data.frame defining which layers to connect,
#' with columns `from`, `to`, and optional `weight`. If `NULL`, adjacent layers
#' are connected in sequence.
#'
#' @param min_similarity Minimum weighted similarity required to keep an
#' interlayer tie.
#'
#' @param resolution_parameter Leiden resolution parameter.
#'
#' @param directed Logical; if `TRUE`, build directed graphs from adjacency
#' matrices.
#'
#' @param add_self_loops Logical; if `TRUE`, add intra-layer community
#' self-loop ties.
#'
#' @param self_loop_multiplier Numeric multiplier applied to self-loop weighted
#' ties.
#'
#' @param objective One of "cpm" or "modularity" for directed networks only
#'
#' @param seed Optional integer seed for reproducible community detection.
#' When supplied, the global RNG state is saved, the RNG is seeded for the
#' duration of the call, and the previous state is restored on exit, so the
#' caller's random number stream (e.g. bootstrap resampling) is unaffected.
#' Defaults to `NULL` (detection inherits the caller's RNG stream, matching
#' previous behavior; call `set.seed()` beforehand for reproducibility).
#'
#' @return A list of class \code{"multilayer_community_fit"} with components:
#'   \describe{
#'     \item{layer_communities}{Per-layer community detection (each with
#'       \code{membership} and \code{communities}). Detected independently
#'       per layer.}
#'     \item{meta_communities}{The cross-layer tracked partition from the
#'       second-stage detection: one integer vector per layer giving each
#'       node's meta-community. This is the membership that reflects the
#'       interlayer ties and any custom \code{layer_links}, and the one
#'       validated by \code{\link{bootstrap_multilayer}}. See
#'       \code{\link{extract_meta_membership}}.}
#'     \item{interlayer_ties}{Interlayer similarity edges between communities
#'       (plus self-loops).}
#'     \item{layer_links}{The layer connectivity used.}
#'   }
#'
#' @section Directed networks:
#' Directed layers are stored as directed graphs and the interlayer self-loop
#' weighting is directed-aware, but community \emph{detection} collapses
#' directed layers to weighted undirected graphs on both Louvain and Leiden
#' (igraph's detectors are undirected-only). For detection that respects edge
#' direction, use the Python package with \code{algorithm = "leiden"}.
#'
#' @examples
#' set.seed(123)
#' layers <- lapply(1:3, function(i) {
#'   w <- matrix(runif(64), 8, 8) * matrix(rbinom(64, 1, 0.35), 8, 8)
#'   w <- (w + t(w)) / 2
#'   diag(w) <- 0
#'   w
#' })
#' fit <- fit_multilayer_weighted_jaccard(layers, algorithm = "louvain")
#' names(fit)
#'
#' @export

fit_multilayer_weighted_jaccard <- function(
    layers,
    algorithm = c("louvain", "leiden"),
    layer_links = NULL,
    min_similarity = 0,
    resolution_parameter = 1,
    directed = FALSE,
    add_self_loops = TRUE,
    self_loop_multiplier = 1,
    objective = NULL,
    seed = NULL
  ) {

  # Check arguments ----
  algorithm <- match.arg(algorithm)

  # Scoped seed for reproducible detection (restores caller RNG state) ----
  if (!is.null(seed)) {
    rng_state <- save_rng_state()
    on.exit(restore_rng_state(rng_state), add = TRUE)
    set.seed(seed)
  }

  # Prepare graph layers and layer links ----
  graph_layers <- prepare_multilayer_graphs(layers, directed = directed)
  links <- make_layer_links(length(graph_layers), layer_links)

  # Fit layer communities ----
  fit <- fit_layer_communities(
    graph_layers,
    algorithm = algorithm,
    resolution_parameter = resolution_parameter,
    directed = directed,
    objective = objective
  )

  # Assign node weights ----
  node_weights <- layer_node_strengths(graph_layers, directed = directed)

  # Calculate interlayer ties ----
  interlayer_ties <- community_overlap_edges(
    fit = fit,
    layer_links = links,
    graph_layers = graph_layers,
    metric = "jaccard",
    min_similarity = min_similarity,
    node_weights_by_layer = node_weights
  )

  # Add community self loops ----
  if (isTRUE(add_self_loops)) {
    interlayer_ties <- add_community_self_loops(
      edge_df = interlayer_ties,
      fit = fit,
      layer_links = links,
      self_loop_multiplier = self_loop_multiplier,
      min_similarity = min_similarity,
      directed = directed
    )
  }

  # Compile multilayer community fit object ----
  # second-stage detection: group per-layer communities into cross-layer
  # meta-communities from the interlayer ties (the tracked partition) ----
  meta <- detect_interlayer_communities(
    layer_communities = fit,
    interlayer_ties = interlayer_ties,
    algorithm = algorithm,
    resolution_parameter = resolution_parameter
  )

  multilayer_weighted_jaccard <- structure(
    list(
      algorithm = algorithm,
      layer_communities = fit,
      meta_communities = meta$membership,
      meta_ids = meta$meta_ids,
      layer_links = links,
      interlayer_ties = interlayer_ties,
      directed = directed,
      weighting = "node_strength"
    ),
    class = "multilayer_community_fit"
  )

  # Return multilayer community fit object ----
  return(multilayer_weighted_jaccard)
}
