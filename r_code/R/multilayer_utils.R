#' @title Prepare Multilayer Graphs
#'
#' @noRd

prepare_multilayer_graphs <- function(layers, directed = FALSE) {

  # check for required package ----
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required.", call. = FALSE)
  }

  # check layers argument ----
  if (!is.list(layers) || length(layers) < 2) {
    stop("`layers` must be a list with at least two network layers.", call. = FALSE)
  }

  # create graph for each layer ----
  graph_layers <- lapply(seq_along(layers), function(i) {
    layer <- layers[[i]]
    if (inherits(layer, "igraph")) {

      ## return layer as-is if already an igraph object ----
      layer_graph <- layer
    } else {

      ## check that layer is a matrix ----
      if (!is.matrix(layer)) {
        stop(sprintf("Layer %d is not an igraph object or adjacency matrix.", i), call. = FALSE)
      }

      ## check that layer matrix is square ----
      if (nrow(layer) != ncol(layer)) {
        stop(sprintf("Layer %d adjacency matrix must be square.", i), call. = FALSE)
      }

      ## create layer graph from adjacency matrix ----
      layer_graph <- igraph::graph_from_adjacency_matrix(
        adjmatrix = layer,
        mode = if (directed) "directed" else "undirected",
        weighted = TRUE,
        diag = FALSE
      )
    }

    # return layer graph ----
    return(layer_graph)
  })

  # assign layer names when no names are present ----
  if (is.null(names(graph_layers))) {
    names(graph_layers) <- paste0("layer_", seq_along(graph_layers))
  }

  # return compiled graph layers ----
  return(graph_layers)
}

#' @title Make Layer Links
#'
#' @noRd

make_layer_links <- function(n_layers, layer_links = NULL) {

  # construct or check appropriate layer links ----
  if (is.null(layer_links)) {

    ## assign equally-weighted links if no links are provided ----
    layer_links <- data.frame(
      from = seq_len(n_layers - 1),
      to = seq(2, n_layers),
      weight = 1,
      stringsAsFactors = FALSE
    )
  } else {

    ## check that layer links consists of a data frame with the appropriate columns ----
    if (!is.data.frame(layer_links) || !all(c("from", "to") %in% names(layer_links))) {
      stop("`layer_links` must be a data.frame with columns `from` and `to`.", call. = FALSE)
    }

    ## assign equal weights if no weight column is present ----
    if (!"weight" %in% names(layer_links)) {
      layer_links$weight <- 1
    }

    ## check that layer link indices are within bounds ----
    if (any(layer_links$from < 1 | layer_links$to < 1 | layer_links$from > n_layers | layer_links$to > n_layers)) {
      stop("`layer_links` indices must be between 1 and number of layers.", call. = FALSE)
    }
  }

  # return layer links ----
  return(layer_links)
}



#' @title Fit Layer Communities
#'
#' @noRd

fit_layer_communities <- function(
    graph_layers,
    algorithm = c("louvain", "leiden"),
    resolution_parameter = 1,
    directed = FALSE,
    objective = NULL
  ) {

  # check algorithm argument ----
  algorithm <- match.arg(algorithm)

  # check objective argument ----
  if (!is.null(objective)) {
    objective <- tolower(objective)
    if (!objective %in% c("modularity", "cpm")) {
      stop("`objective` must be one of 'modularity' or 'cpm'.", call. = FALSE)
    }
    if (algorithm == "louvain" && objective == "cpm") {
      stop("Louvain does not support the CPM objective. Use algorithm = 'leiden'.", call. = FALSE)
    }
  }

  # resolve effective objective between explicit choice or default based on direction ----
  effective_objective <- if (!is.null(objective)) objective else if (directed) "cpm" else "modularity"

  # check for negative edge weights for modularity objectives ----
  if (effective_objective == "modularity") {
    for (i in seq_along(graph_layers)) {
      w <- igraph::E(graph_layers[[i]])$weight
      if (!is.null(w) && any(w < 0)) {
        stop(
          sprintf(
            fmt = paste0(
              "Layer %d contains negative edge weights. ",
              "Modularity-based methods do not support negative weights. ",
              "Use objective = \"cpm\" to select the CPM objective, ",
              "which handles negative weights correctly."
            ),
            i
          ),
          call. = FALSE
        )
      }
    }
  }

  # warn once for directed leiden (collapsed per layer below) ----
  if (algorithm == "leiden" && directed &&
      any(vapply(graph_layers, igraph::is_directed, logical(1)))) {
    warning(
      "Leiden in igraph does not support directed graphs; ",
      "collapsing directed layers to weighted undirected graphs. ",
      "For directed-aware detection use the Python package (leidenalg).",
      call. = FALSE
    )
  }

  # fit communities for each layer ----
  layer_communities <- lapply(graph_layers, function(g) {
    g_input <- g

    ## identify community clusters ----
    if (algorithm == "louvain") {

      ### convert to undirected graph ----
      if (directed && igraph::is_directed(g_input)) {
        g <- igraph::as_undirected(
          graph = g_input,
          mode = "collapse",
          edge.attr.comb = list(weight = "sum")
        )
      }

      ### find louvain community clusters ----
      cl <- igraph::cluster_louvain(
        graph = g,
        weights = igraph::E(g)$weight,
        resolution = resolution_parameter
      )
    } else {

      ### convert to undirected graph ----
      ### (igraph's cluster_leiden supports undirected graphs only; see
      ### community/leiden.c. Mirror the louvain path: collapse directed
      ### layers to weighted undirected graphs; warned once above.)
      if (directed && igraph::is_directed(g_input)) {
        g <- igraph::as_undirected(
          graph = g_input,
          mode = "collapse",
          edge.attr.comb = list(weight = "sum")
        )
      }

      ### find leiden community clusters ----
      cl <- igraph::cluster_leiden(
        graph = g,
        objective_function = if (effective_objective == "cpm") "CPM" else "modularity",
        resolution = resolution_parameter,
        weights = igraph::E(g)$weight
      )
    }

    ## compile communities information for a single layer ----
    communities <- list(
      membership = igraph::membership(cl),
      modularity = if (igraph::is_directed(g_input) || effective_objective == "cpm") {
        NA_real_
      } else {
        igraph::modularity(g, igraph::membership(cl), weights = igraph::E(g)$weight)
      },
      communities = split(seq_along(igraph::membership(cl)), igraph::membership(cl))
    )

    ## return communities information for a single layer ----
    return(communities)
  })

  # return layer communities ----
  return(layer_communities)
}



#' @title Weighted Jaccard
#'
#' @noRd

weighted_jaccard <- function(a, b) {
  inter <- length(intersect(a, b))
  union <- length(union(a, b))
  if (union == 0) {
    jaccard <- 0
  } else {
    jaccard <- inter / union
  }
  return(jaccard)
}



#' @title Weighted Overlap
#'
#' @noRd

weighted_overlap <- function(a, b) {
  inter <- length(intersect(a, b))
  min_size <- min(length(a), length(b))
  if (min_size == 0) {
    overlap <- 0
  } else {
    overlap <- inter / min_size
  }
  return(overlap)
}



#' @title Weighted Jaccard Similarity
#'
#' @noRd

weighted_jaccard_similarity <- function(a, b, weights_a, weights_b) {

  # assign nodes ----
  nodes <- union(a, b)

  # assign weighted similarity ----
  if (length(nodes) == 0) {
    weighted_similarity <- 0
  } else {
    node_keys <- as.character(nodes)
    wa <- weights_a[as.character(nodes)]
    wb <- weights_b[as.character(nodes)]
    wa[is.na(wa)] <- 0
    wb[is.na(wb)] <- 0
    inter_weight <- sum(pmin(wa, wb))
    union_weight <- sum(pmax(wa, wb))
    if (union_weight == 0) {
      weighted_similarity <- 0
    } else {
      weighted_similarity <- inter_weight / union_weight
    }
  }

  # return weighted similarity ----
  return(weighted_similarity)
}



#' @title Weighted Overlap Similarity
#'
#' @noRd

weighted_overlap_similarity <- function(a, b, weights_a, weights_b) {

  # assign weighted similarity ----
  inter <- intersect(a, b)
  inter_weight <- sum(vapply(inter, function(node) {
    min(weights_a[[as.character(node)]], weights_b[[as.character(node)]])
  }, numeric(1)))
  a_weight <- sum(vapply(a, function(node) weights_a[[as.character(node)]], numeric(1)))
  b_weight <- sum(vapply(b, function(node) weights_b[[as.character(node)]], numeric(1)))
  min_weight <- min(a_weight, b_weight)
  if (min_weight == 0) {
    weighted_similarity <- 0
  } else {
    weighted_similarity <- inter_weight / min_weight
  }

  # return weighted similarity ----
  return(weighted_similarity)
}



#' @title Layer Node Strengths
#'
#' @noRd

layer_node_strengths <- function(graph_layers, directed = FALSE) {

  # assign node strengths ----
  node_strengths <- lapply(graph_layers, function(g) {
    strength_vals <- igraph::strength(g, mode = "all", loops = FALSE, weights = igraph::E(g)$weight)
    names(strength_vals) <- as.character(seq_along(strength_vals))
    return(strength_vals)
  })

  # return node strengths ----
  return(node_strengths)
}



#' @title Community Overlap Edges
#'
#' @noRd

community_overlap_edges <- function(
    fit,
    layer_links,
    metric = c("jaccard", "overlap"),
    min_similarity = 0,
    node_weights_by_layer = NULL
  ) {

  # check metric argument ----
  metric <- match.arg(metric)

  # assign appropriate simulation functions ----
  sim_fun <- if (metric == "jaccard") weighted_jaccard else weighted_overlap
  weighted_sim_fun <- if (metric == "jaccard") weighted_jaccard_similarity else weighted_overlap_similarity

  # build overlap edges for each layer link ----
  edge_rows <- list()
  for (i in seq_len(nrow(layer_links))) {

    ## extract layer link information ----
    from_idx <- layer_links$from[i]
    to_idx <- layer_links$to[i]
    layer_weight <- layer_links$weight[i]
    from_comms <- fit[[from_idx]]$communities
    to_comms <- fit[[to_idx]]$communities
    from_ids <- as.integer(names(from_comms))
    to_ids <- as.integer(names(to_comms))

    ## calculate overlap for each community edge ----
    for (from_c in from_ids) {
      for (to_c in to_ids) {
        from_c_str <- as.character(from_c)
        to_c_str <- as.character(to_c)

        ### calculate similarity ----
        if (is.null(node_weights_by_layer)) {
          sim <- sim_fun(from_comms[[from_c_str]], to_comms[[to_c_str]])
        } else {
          sim <- weighted_sim_fun(
            from_comms[[from_c_str]],
            to_comms[[to_c_str]],
            node_weights_by_layer[[from_idx]],
            node_weights_by_layer[[to_idx]]
          )
        }

        ### calculate weighted similarity ----
        weighted_sim <- sim * layer_weight

        ### compile edge parameters ----
        if (weighted_sim >= min_similarity) {
          edge_rows[[length(edge_rows) + 1]] <- data.frame(
            from_layer = from_idx,
            to_layer = to_idx,
            from_community = from_c,
            to_community = to_c,
            similarity = sim,
            layer_weight = layer_weight,
            weighted_similarity = weighted_sim,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  # compile overlap edges data ----
  if (length(edge_rows) == 0) {
    overlap_edges <- data.frame(
      from_layer = integer(0),
      to_layer = integer(0),
      from_community = integer(0),
      to_community = integer(0),
      similarity = numeric(0),
      layer_weight = numeric(0),
      weighted_similarity = numeric(0)
    )
  } else {
    overlap_edges <- do.call(rbind, edge_rows)
  }

  # return overlap edges ----
  return(overlap_edges)
}



#' @title Add Community Self Loops
#'
#' @noRd

add_community_self_loops <- function(
    edge_df,
    fit,
    layer_links,
    self_loop_multiplier = 1,
    min_similarity = 0,
    directed = FALSE
  ) {

  # initialize loop rows ----
  loop_rows <- list()

  # Undirected self-loops count each internal edge twice (A[i,j] + A[j,i])
  self_sim <- if (directed) 1 else 2

  # compute maximum layer weight for each unique layer across all layer links ----
  all_layers <- sort(unique(c(layer_links$from, layer_links$to)))
  layer_weights <- vapply(
    all_layers,
    function(idx) {
      max(
        c(
          layer_links$weight[layer_links$from == idx],
          layer_links$weight[layer_links$to == idx]
        )
      )
    },
    numeric(1)
  )

  # assign names to layer weights ----
  names(layer_weights) <- as.character(all_layers)

  # assemble loop data for each layer and community ----
  for (layer_idx in all_layers) {
    layer_weight <- layer_weights[as.character(layer_idx)]
    comms <- fit[[layer_idx]]$communities
    for (comm_idx in as.integer(names(comms))) {
      weighted_sim <- self_sim * layer_weight * self_loop_multiplier
      if (weighted_sim >= min_similarity) {
        loop_rows[[length(loop_rows) + 1]] <- data.frame(
          from_layer = layer_idx,
          to_layer = layer_idx,
          from_community = comm_idx,
          to_community = comm_idx,
          similarity = self_sim,
          layer_weight = layer_weight,
          weighted_similarity = weighted_sim,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  # compile community self loops ----
  if (length(loop_rows) == 0) {
    community_self_loops <- edge_df
  } else {
    community_self_loops <- rbind(edge_df, do.call(rbind, loop_rows))
  }

  # return community self loops ----
  return(community_self_loops)
}


#' @title Detect cross-layer (meta) communities from interlayer ties
#'
#' @description Second-stage community detection. Treats each per-layer
#' community as a node in a community graph whose edges are the interlayer
#' similarity ties (plus the community self-loops), then runs community
#' detection on that graph to group per-layer communities into cross-layer
#' \emph{meta-communities}. The result is the tracked partition: which
#' communities persist, merge, or split across layers. This is the step that
#' makes custom \code{layer_links} and the interlayer coupling actually affect
#' the membership that is returned and validated.
#'
#' @param layer_communities The \code{layer_communities} element of a fit
#'   (per-layer detection results, each with \code{membership} and
#'   \code{communities}).
#'
#' @param interlayer_ties The \code{interlayer_ties} data frame of a fit
#'   (columns \code{from_layer}, \code{to_layer}, \code{from_community},
#'   \code{to_community}, \code{weighted_similarity}), including self-loops.
#'
#' @param algorithm Community algorithm for the second stage: \code{"louvain"}
#'   or \code{"leiden"} (match the per-layer algorithm).
#'
#' @param resolution_parameter Resolution for the second-stage detection.
#'
#' @return A list with \describe{
#'   \item{meta_ids}{Named integer vector mapping each per-layer community
#'     (key \code{"L<layer>C<community>"}) to its meta-community id.}
#'   \item{membership}{List with one integer vector per layer giving each
#'     node's meta-community assignment (node order).}
#' }
#'
#' @noRd
detect_interlayer_communities <- function(
    layer_communities,
    interlayer_ties,
    algorithm = c("louvain", "leiden"),
    resolution_parameter = 1
  ) {

  algorithm <- match.arg(algorithm)
  key <- function(l, c) paste0("L", l, "C", c)
  n_layers <- length(layer_communities)

  # The second stage needs community-level ties. Node-level identity ties are
  # not handled here (their multislice form is a separate node-level build);
  # fall back to per-layer communities made globally distinct, i.e. no
  # cross-layer merging. ----
  has_comm_ties <- !is.null(interlayer_ties) && nrow(interlayer_ties) > 0 &&
    all(c("from_community", "to_community", "weighted_similarity") %in%
          names(interlayer_ties))
  if (!has_comm_ties) {
    offset <- 0L
    membership <- vector("list", n_layers)
    for (t in seq_len(n_layers)) {
      mem <- as.integer(layer_communities[[t]]$membership)
      membership[[t]] <- mem + offset
      offset <- offset + max(mem)
    }
    return(list(meta_ids = NULL, membership = membership))
  }

  # enumerate every per-layer community as a super-node ----
  supernodes <- unique(unlist(lapply(seq_len(n_layers), function(t) {
    key(t, as.integer(names(layer_communities[[t]]$communities)))
  })))

  # build the community graph from interlayer ties (incl. self-loops) ----
  if (!is.null(interlayer_ties) && nrow(interlayer_ties) > 0) {
    edge_df <- data.frame(
      from = key(interlayer_ties$from_layer, interlayer_ties$from_community),
      to = key(interlayer_ties$to_layer, interlayer_ties$to_community),
      weight = interlayer_ties$weighted_similarity,
      stringsAsFactors = FALSE
    )
  } else {
    edge_df <- data.frame(from = character(0), to = character(0),
                          weight = numeric(0), stringsAsFactors = FALSE)
  }

  g <- igraph::graph_from_data_frame(
    d = edge_df,
    directed = FALSE,
    vertices = data.frame(name = supernodes, stringsAsFactors = FALSE)
  )

  # second-stage detection on the community graph ----
  if (igraph::ecount(g) == 0) {
    # no ties: every per-layer community is its own meta-community
    meta_ids <- stats::setNames(seq_along(supernodes), supernodes)
  } else {
    if (algorithm == "louvain") {
      cl <- igraph::cluster_louvain(
        graph = g, weights = igraph::E(g)$weight,
        resolution = resolution_parameter
      )
    } else {
      cl <- igraph::cluster_leiden(
        graph = g, objective_function = "modularity",
        weights = igraph::E(g)$weight,
        resolution = resolution_parameter, n_iterations = 3
      )
    }
    meta_ids <- as.integer(igraph::membership(cl))
    names(meta_ids) <- igraph::V(g)$name
  }

  # map each node to its meta-community, per layer ----
  membership <- lapply(seq_len(n_layers), function(t) {
    mem <- layer_communities[[t]]$membership
    as.integer(meta_ids[key(t, as.integer(mem))])
  })

  list(meta_ids = meta_ids, membership = membership)
}
