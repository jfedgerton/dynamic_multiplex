#' Fit multilayer communities with identity interlayer ties
#'
#' Runs Louvain or Leiden community detection for each layer and creates
#' interlayer ties between the same node in selected adjacent layers.
#'
#' @param layers List of `igraph` objects or square adjacency matrices.
#' @param algorithm Community algorithm: `"louvain"` or `"leiden"`.
#' @param layer_links Optional data.frame defining which layers to connect,
#' with columns `from`, `to`, and optional `weight`. If `NULL`, adjacent layers
#' are connected in sequence.
#' @param resolution_parameter Leiden resolution parameter.
#' @param directed Logical; if `TRUE`, build directed graphs from adjacency matrices.
#'   For `algorithm = "louvain"`, directed layers are collapsed to undirected
#'   weighted graphs before community detection.
#'
#' @return A list with detected communities per layer and node-level interlayer ties.
#' @export
fit_multilayer_identity_ties <- function(layers,
                                         algorithm = c("louvain", "leiden"),
                                         layer_links = NULL,
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

  n_nodes <- igraph::vcount(graph_layers[[1]])
  if (!all(vapply(graph_layers, igraph::vcount, integer(1)) == n_nodes)) {
    stop("All layers must have the same number of nodes for identity ties.", call. = FALSE)
  }

  ties <- do.call(
    rbind,
    lapply(seq_len(nrow(links)), function(i) {
      data.frame(
        from_layer = links$from[i],
        to_layer = links$to[i],
        node = seq_len(n_nodes),
        layer_weight = links$weight[i],
        stringsAsFactors = FALSE
      )
    })
  )

  structure(
    list(
      algorithm = algorithm,
      layer_communities = fit,
      layer_links = links,
      interlayer_ties = ties,
      directed = directed
    ),
    class = "multilayer_identity_fit"
  )
}
