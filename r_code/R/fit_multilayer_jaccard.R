#' Fit multilayer communities and interlayer weighted Jaccard ties
#'
#' Runs Louvain or Leiden community detection for each layer and creates
#' interlayer ties between communities in selected layer pairs using weighted
#' Jaccard similarity.
#'
#' @param layers List of `igraph` objects or square adjacency matrices.
#' @param algorithm Community algorithm: `"louvain"` or `"leiden"`.
#' @param layer_links Optional data.frame defining which layers to connect,
#' with columns `from`, `to`, and optional `weight`. If `NULL`, adjacent layers
#' are connected in sequence.
#' @param min_similarity Minimum weighted similarity required to keep an
#' interlayer tie.
#' @param resolution_parameter Leiden resolution parameter.
#' @param directed Logical; if `TRUE`, build directed graphs from adjacency matrices.
#'   For `algorithm = "louvain"`, directed layers are collapsed to undirected
#'   weighted graphs before community detection.
#'
#' @return A list with detected communities per layer and interlayer ties.
#' @export
fit_multilayer_jaccard <- function(layers,
                                   algorithm = c("louvain", "leiden"),
                                   layer_links = NULL,
                                   min_similarity = 0,
                                   resolution_parameter = 1,
                                   directed = FALSE) {
  algorithm <- match.arg(algorithm)

  graph_layers <- prepare_multilayer_graphs(layers, directed = directed)
  links <- make_layer_links(length(graph_layers), layer_links)
  fit <- fit_layer_communities(
    graph_layers,
    algorithm = algorithm,
    resolution_parameter = resolution_parameter,
    directed = directed
  )

  interlayer_ties <- community_overlap_edges(
    fit = fit,
    layer_links = links,
    metric = "jaccard",
    min_similarity = min_similarity
  )

  structure(
    list(
      algorithm = algorithm,
      layer_communities = fit,
      layer_links = links,
      interlayer_ties = interlayer_ties,
      directed = directed
    ),
    class = "multilayer_community_fit"
  )
}
