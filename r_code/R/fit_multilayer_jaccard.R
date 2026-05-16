#' @title Fit multilayer communities and interlayer weighted Jaccard ties
#'
#' @description Runs Louvain or Leiden community detection for each layer and
#' creates interlayer ties between communities in selected layer pairs using
#' weighted Jaccard similarity.
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
#' matrices. For `algorithm = "louvain"`, directed layers are collapsed to
#' undirected weighted graphs before community detection.
#'
#' @param add_self_loops Logical; if `TRUE`, add intra-layer community
#' self-loop ties after similarity computation to mimic aggregation behavior.
#'
#' @param self_loop_multiplier Numeric multiplier applied to self-loop weighted
#' ties. The default `1` yields self-loop weighted similarity of
#' `2 * layer_weight` for undirected or `1 * layer_weight` for directed
#' networks, following the Louvain/Leiden community aggregation convention.
#'
#' @return A list with detected communities per layer and interlayer ties.
#'
#' @export

fit_multilayer_jaccard <- function(
    layers,
    algorithm = c("louvain", "leiden"),
    layer_links = NULL,
    min_similarity = 0,
    resolution_parameter = 1,
    directed = FALSE,
    add_self_loops = TRUE,
    self_loop_multiplier = 1,
    objective = NULL
  ) {

  # Check arguments ----
  algorithm <- match.arg(algorithm)

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

  # Calculate interlayer ties ----
  interlayer_ties <- community_overlap_edges(
    fit = fit,
    layer_links = links,
    metric = "jaccard",
    min_similarity = min_similarity
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
  multilayer_jaccard <- structure(
    list(
      algorithm = algorithm,
      layer_communities = fit,
      layer_links = links,
      interlayer_ties = interlayer_ties,
      directed = directed
    ),
    class = "multilayer_community_fit"
  )

  # Return multilayer community fit object ----
  return(multilayer_jaccard)
}
