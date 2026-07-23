#' @title Fit multilayer communities with identity interlayer ties
#'
#' @description Runs Louvain or Leiden community detection for each layer and
#' creates interlayer ties between the same node in selected adjacent layers.
#' Layers may contain different node sets; only nodes present in both layers of
#' a linked pair receive interlayer edges. Nodes are matched by
#' \code{V(g)$name} when available, or by vertex index otherwise.
#'
#' For this specification the cross-layer \code{meta_communities} are obtained
#' by \strong{Mucha (2010) multislice modularity}: the layers are stacked into a
#' single supra-graph (intra-layer edges are each layer's own adjacency) with
#' interlayer identity edges joining each node to its copies in the coupled
#' layers (weighted by \code{layer_links}), and one community detection is run
#' on the whole supra-graph. The coupling strength is the layer-link weight.
#'
#' @param layers List of `igraph` objects or square adjacency matrices.
#'
#' @param algorithm Community algorithm: `"louvain"` or `"leiden"`.
#'
#' @param layer_links Optional data.frame defining which layers to connect,
#' with columns `from`, `to`, and optional `weight`. If `NULL`, adjacent layers
#' are connected in sequence.
#'
#' @param resolution_parameter Leiden resolution parameter.
#'
#' @param directed Logical; if `TRUE`, build directed graphs from adjacency
#' matrices. Directed layers are collapsed to undirected weighted graphs before
#' community detection (igraph supports Leiden on undirected graphs only;
#' a warning is issued for `algorithm = "leiden"`).
#'
#' @param objective One of "cpm" or "modularity" for directed networks only
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
#'   m <- matrix(rbinom(64, 1, 0.35), nrow = 8)
#'   m <- pmax(m, t(m))
#'   diag(m) <- 0
#'   m
#' })
#' fit <- fit_multilayer_identity_ties(layers, algorithm = "louvain")
#' names(fit)
#'
#' @export

fit_multilayer_identity_ties <- function(
    layers,
    algorithm = c("louvain", "leiden"),
    layer_links = NULL,
    resolution_parameter = 1,
    directed = FALSE,
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

  # Create specified layer ties ----
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

  # Create null ties data frame if necessary ----
  if (is.null(ties)) {
    ties <- data.frame(
      from_layer = integer(0),
      to_layer = integer(0),
      node = character(0),
      layer_weight = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  # Compile multilayer identity fit object ----
  # second-stage detection: group per-layer communities into cross-layer
  # meta-communities from the interlayer ties (the tracked partition) ----
  # node-level Mucha multislice: stack layers (intra = original adjacency) with
  # identity interlayer ties, single detection on the supra-graph ----
  meta_membership <- detect_multislice_communities(
    graph_layers = graph_layers,
    interlayer_ties = ties,
    algorithm = algorithm
  )

  multilayer_identity_ties <- structure(
    list(
      algorithm = algorithm,
      layer_communities = fit,
      meta_communities = meta_membership,
      meta_ids = NULL,
      layer_links = links,
      interlayer_ties = ties,
      directed = directed
    ),
    class = "multilayer_identity_fit"
  )

  # Return multilayer identity fit object ----
  return(multilayer_identity_ties)
}
