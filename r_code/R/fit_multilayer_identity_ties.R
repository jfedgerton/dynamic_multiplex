#' Fit multilayer communities with identity interlayer ties
#'
#' Runs Louvain or Leiden community detection for each layer and creates
#' interlayer ties between the same node in selected adjacent layers.
#' Layers may contain different node sets; only nodes present in both layers
#' of a linked pair receive interlayer edges. Nodes are matched by
#' \code{V(g)$name} when available, or by vertex index otherwise.
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

  ties <- do.call(
    rbind,
    lapply(seq_len(nrow(links)), function(i) {
      g_from <- graph_layers[[links$from[i]]]
      g_to <- graph_layers[[links$to[i]]]

      names_from <- igraph::V(g_from)$name
      names_to <- igraph::V(g_to)$name

      if (is.null(names_from)) names_from <- seq_len(igraph::vcount(g_from))
      if (is.null(names_to)) names_to <- seq_len(igraph::vcount(g_to))

      shared <- intersect(names_from, names_to)

      if (length(shared) == 0) return(NULL)

      data.frame(
        from_layer = links$from[i],
        to_layer = links$to[i],
        node = shared,
        layer_weight = links$weight[i],
        stringsAsFactors = FALSE
      )
    })
  )

  if (is.null(ties)) {
    ties <- data.frame(
      from_layer = integer(0),
      to_layer = integer(0),
      node = character(0),
      layer_weight = numeric(0),
      stringsAsFactors = FALSE
    )
  }

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
