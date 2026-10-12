#' @title Fit multilayer communities and interlayer node-strength weighted
#' overlap ties
#'
#' @description Runs Louvain or Leiden community detection for each layer and
#' creates interlayer ties between communities in selected layer pairs using
#' node-strength weighted overlap coefficient.
#'
#' @param layers List of `igraph` objects or square adjacency matrices.
#'
#' @param algorithm Community algorithm: `"leiden"` (default) or `"louvain"`.
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
#' @param objective Objective for the first-stage (per-layer) detection:
#' `"modularity"` or `"cpm"` (CPM requires `algorithm = "leiden"`). Defaults to
#' `"cpm"` for directed layers and `"modularity"` otherwise. It does not
#' affect the second stage, which always maximizes modularity on the
#' community graph.
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
#' different node sets (nodes entering or exiting the system). Layers must
#' then carry vertex names (or row names on the adjacency matrices) so nodes
#' can be matched across layers; only nodes present in both layers of a
#' linked pair contribute to the interlayer similarity. Defaults to `FALSE`:
#' layers must share the same node universe.
#'
#' @param tracking \code{"smoothed"} (default, the 1.3.2 behaviour): stage 2
#'   is fitted on all layers at once, so a layer's meta ids can depend on later
#'   layers. \code{"filtered"}: stage 2 runs forward in time and the meta ids of
#'   layers 1..t depend only on layers 1..t (no future information). Filtered
#'   mode assumes the layers are in time order and needs a non-NULL
#'   \code{seed}.
#' @param repair_chains Logical. If \code{TRUE}, rejoin meta-communities that
#'   stage 2 cut along time: a group that ends at layer s and one that starts at
#'   a linked layer t are joined when their communities are mutual best matches
#'   with similarity at least \code{repair_min_similarity}. Default \code{FALSE}
#'   (1.3.2 behaviour).
#' @param repair_min_similarity Similarity floor for chain repair, between 0 and 1.
#'   Default 0.2.
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
#' fit <- fit_multilayer_weighted_overlap(layers, algorithm = "louvain")
#' names(fit)
#'
#' @export

fit_multilayer_weighted_overlap <- function(
    layers,
    algorithm = c("leiden", "louvain"),
    layer_links = NULL,
    min_similarity = 0,
    resolution_parameter = 1,
    directed = FALSE,
    add_self_loops = TRUE,
    self_loop_multiplier = 1,
    objective = NULL,
    seed = 123,
    allow_unequal_nodes = FALSE,
    tracking = "smoothed",
    repair_chains = FALSE,
    repair_min_similarity = 0.2
  ) {

  # Check arguments ----
  algorithm <- match.arg(algorithm)
  .check_repair_args(tracking, repair_chains, repair_min_similarity, seed)

  # Scoped seed for reproducible detection (restores caller RNG state) ----
  if (!is.null(seed)) {
    rng_state <- save_rng_state()
    on.exit(restore_rng_state(rng_state), add = TRUE)
    set.seed(seed)
  }

  # Prepare graph layers and layer links ----
  graph_layers <- prepare_multilayer_graphs(
    layers,
    directed = directed,
    require_same_nodes = !allow_unequal_nodes
  )
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
  pair_ties <- community_overlap_edges(
    fit = fit,
    layer_links = links,
    graph_layers = graph_layers,
    metric = "overlap",
    min_similarity = min_similarity,
    node_weights_by_layer = node_weights
  )

  # Compile multilayer community fit object ----
  # second-stage detection: group per-layer communities into cross-layer
  # meta-communities from the interlayer ties (the tracked partition) ----
  meta <- run_two_stage(
    fit, pair_ties, links,
    add_self_loops = add_self_loops,
    self_loop_multiplier = self_loop_multiplier,
    min_similarity = min_similarity,
    directed = directed,
    algorithm = algorithm,
    resolution_parameter = resolution_parameter,
    seed = seed,
    tracking = tracking,
    repair_chains = repair_chains,
    repair_min_similarity = repair_min_similarity
  )
  interlayer_ties <- meta$interlayer_ties

  multilayer_weighted_overlap <- structure(
    list(
      algorithm = algorithm,
      layer_communities = fit,
      meta_communities = meta$membership,
      meta_ids = meta$meta_ids,
      layer_links = links,
      interlayer_ties = interlayer_ties,
      directed = directed,
      weighting = "node_strength",
      tracking = meta$tracking,
      repair_chains = meta$repair_chains,
      repair_min_similarity = meta$repair_min_similarity,
      chain_repairs = meta$chain_repairs
    ),
    class = "multilayer_community_fit"
  )

  # Return multilayer community fit object ----
  return(multilayer_weighted_overlap)
}
