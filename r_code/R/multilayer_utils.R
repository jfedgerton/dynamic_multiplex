#' @keywords internal
prepare_multilayer_graphs <- function(layers, directed = FALSE) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required.", call. = FALSE)
  }

  if (!is.list(layers) || length(layers) < 2) {
    stop("`layers` must be a list with at least two network layers.", call. = FALSE)
  }

  graph_layers <- lapply(seq_along(layers), function(i) {
    layer <- layers[[i]]
    if (inherits(layer, "igraph")) {
      return(layer)
    }

    if (!is.matrix(layer)) {
      stop(sprintf("Layer %d is not an igraph object or adjacency matrix.", i), call. = FALSE)
    }

    if (nrow(layer) != ncol(layer)) {
      stop(sprintf("Layer %d adjacency matrix must be square.", i), call. = FALSE)
    }

    igraph::graph_from_adjacency_matrix(
      adjmatrix = layer,
      mode = if (directed) "directed" else "undirected",
      weighted = TRUE,
      diag = FALSE
    )
  })

  if (is.null(names(graph_layers))) {
    names(graph_layers) <- paste0("layer_", seq_along(graph_layers))
  }

  graph_layers
}

#' @keywords internal
make_layer_links <- function(n_layers, layer_links = NULL) {
  if (is.null(layer_links)) {
    return(data.frame(
      from = seq_len(n_layers - 1),
      to = seq(2, n_layers),
      weight = 1,
      stringsAsFactors = FALSE
    ))
  }

  if (!is.data.frame(layer_links) || !all(c("from", "to") %in% names(layer_links))) {
    stop("`layer_links` must be a data.frame with columns `from` and `to`.", call. = FALSE)
  }

  if (!"weight" %in% names(layer_links)) {
    layer_links$weight <- 1
  }

  if (any(layer_links$from < 1 | layer_links$to < 1 |
          layer_links$from > n_layers | layer_links$to > n_layers)) {
    stop("`layer_links` indices must be between 1 and number of layers.", call. = FALSE)
  }

  layer_links
}

#' @keywords internal
fit_layer_communities <- function(graph_layers, algorithm = c("louvain", "leiden"),
                                  resolution_parameter = 1,
                                  directed = FALSE,
                                  objective = NULL) {
  algorithm <- match.arg(algorithm)

  if (!is.null(objective)) {
    objective <- tolower(objective)
    if (!objective %in% c("modularity", "cpm")) {
      stop("`objective` must be one of 'modularity' or 'cpm'.", call. = FALSE)
    }
    if (algorithm == "louvain" && objective == "cpm") {
      stop("Louvain does not support the CPM objective. Use algorithm = 'leiden'.", call. = FALSE)
    }
  }

  # Resolve effective objective: explicit choice, or default based on direction
  effective_objective <- if (!is.null(objective)) objective else if (directed) "cpm" else "modularity"

  if (effective_objective == "modularity") {
    for (i in seq_along(graph_layers)) {
      w <- igraph::E(graph_layers[[i]])$weight
      if (!is.null(w) && any(w < 0)) {
        stop(
          sprintf(
            "Layer %d contains negative edge weights. Modularity-based methods do not support negative weights. Use objective = \"cpm\" to select the CPM objective, which handles negative weights correctly.",
            i
          ),
          call. = FALSE
        )
      }
    }
  }

  lapply(graph_layers, function(g) {
    g_input <- g

    if (algorithm == "louvain") {
      if (directed && igraph::is_directed(g_input)) {
        g <- igraph::as.undirected(g_input, mode = "collapse", edge.attr.comb = list(weight = "sum"))
      }
      cl <- igraph::cluster_louvain(g, weights = igraph::E(g)$weight,
                                       resolution = resolution_parameter)
    } else {
      cl <- igraph::cluster_leiden(
        g,
        objective_function = if (effective_objective == "cpm") "CPM" else "modularity",
        resolution_parameter = resolution_parameter,
        weights = igraph::E(g)$weight
      )
    }

    list(
      membership = igraph::membership(cl),
      modularity = if (igraph::is_directed(g_input) || effective_objective == "cpm") NA_real_ else igraph::modularity(g, igraph::membership(cl),
                                      weights = igraph::E(g)$weight),
      communities = split(seq_along(igraph::membership(cl)), igraph::membership(cl))
    )
  })
}

#' @keywords internal
weighted_jaccard <- function(a, b) {
  inter <- length(intersect(a, b))
  union <- length(union(a, b))
  if (union == 0) {
    return(0)
  }
  inter / union
}

#' @keywords internal
weighted_overlap <- function(a, b) {
  inter <- length(intersect(a, b))
  min_size <- min(length(a), length(b))
  if (min_size == 0) {
    return(0)
  }
  inter / min_size
}


#' @keywords internal
weighted_jaccard_similarity <- function(a, b, weights_a, weights_b) {
  nodes <- union(a, b)
  if (length(nodes) == 0) {
    return(0)
  }

  node_keys <- as.character(nodes)
  wa <- weights_a[as.character(nodes)]
  wb <- weights_b[as.character(nodes)]
  wa[is.na(wa)] <- 0
  wb[is.na(wb)] <- 0

  inter_weight <- sum(pmin(wa, wb))
  union_weight <- sum(pmax(wa, wb))

  if (union_weight == 0) {
    return(0)
  }
  inter_weight / union_weight
}

#' @keywords internal
weighted_overlap_similarity <- function(a, b, weights_a, weights_b) {
  inter <- intersect(a, b)
  inter_weight <- sum(vapply(inter, function(node) {
    min(weights_a[[as.character(node)]], weights_b[[as.character(node)]])
  }, numeric(1)))

  a_weight <- sum(vapply(a, function(node) weights_a[[as.character(node)]], numeric(1)))
  b_weight <- sum(vapply(b, function(node) weights_b[[as.character(node)]], numeric(1)))
  min_weight <- min(a_weight, b_weight)

  if (min_weight == 0) {
    return(0)
  }
  inter_weight / min_weight
}

#' @keywords internal
layer_node_strengths <- function(graph_layers, directed = FALSE) {
  lapply(graph_layers, function(g) {
    strength_vals <- igraph::strength(g, mode = "all", loops = FALSE, weights = igraph::E(g)$weight)
    names(strength_vals) <- as.character(seq_along(strength_vals))
    strength_vals
  })
}

#' @keywords internal
community_overlap_edges <- function(fit, layer_links, metric = c("jaccard", "overlap"),
                                    min_similarity = 0,
                                    node_weights_by_layer = NULL) {
  metric <- match.arg(metric)
  sim_fun <- if (metric == "jaccard") weighted_jaccard else weighted_overlap
  weighted_sim_fun <- if (metric == "jaccard") weighted_jaccard_similarity else weighted_overlap_similarity

  edge_rows <- list()

  for (i in seq_len(nrow(layer_links))) {
    from_idx <- layer_links$from[i]
    to_idx <- layer_links$to[i]
    layer_weight <- layer_links$weight[i]

    from_comms <- fit[[from_idx]]$communities
    to_comms <- fit[[to_idx]]$communities

    from_ids <- as.integer(names(from_comms))
    to_ids <- as.integer(names(to_comms))

    for (from_c in from_ids) {
      for (to_c in to_ids) {
        from_c_str <- as.character(from_c)
        to_c_str <- as.character(to_c)
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
        weighted_sim <- sim * layer_weight

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

  if (length(edge_rows) == 0) {
    return(data.frame(
      from_layer = integer(0),
      to_layer = integer(0),
      from_community = integer(0),
      to_community = integer(0),
      similarity = numeric(0),
      layer_weight = numeric(0),
      weighted_similarity = numeric(0)
    ))
  }

  do.call(rbind, edge_rows)
}


#' @keywords internal
add_community_self_loops <- function(edge_df, fit, layer_links,
                                     self_loop_multiplier = 1,
                                     min_similarity = 0,
                                     directed = FALSE) {
  loop_rows <- list()

  # Undirected self-loops count each internal edge twice (A[i,j] + A[j,i])
  self_sim <- if (directed) 1 else 2

  # Compute max layer_weight for each unique layer across all links
  all_layers <- sort(unique(c(layer_links$from, layer_links$to)))
  layer_weights <- vapply(all_layers, function(idx) {
    max(c(layer_links$weight[layer_links$from == idx],
          layer_links$weight[layer_links$to == idx]))
  }, numeric(1))
  names(layer_weights) <- as.character(all_layers)

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

  if (length(loop_rows) == 0) {
    return(edge_df)
  }

  rbind(edge_df, do.call(rbind, loop_rows))
}
