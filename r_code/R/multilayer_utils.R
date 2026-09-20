#' @title Prepare Multilayer Graphs
#'
#' @noRd

prepare_multilayer_graphs <- function(layers, directed = FALSE, require_same_nodes = TRUE) {

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

  # validate shared node universe across all layers ----
  # Layers are compared on vertex names when present, on vertex indices
  # otherwise. Interlayer ties and meta-communities assume node i in one
  # layer is node i in every other layer, so unequal universes are rejected
  # unless the caller explicitly opts out (identity ties only).
  if (require_same_nodes) {
    node_ids <- lapply(graph_layers, function(g) {
      nm <- igraph::V(g)$name
      if (is.null(nm)) as.character(seq_len(igraph::vcount(g))) else as.character(nm)
    })
    ref <- sort(node_ids[[1]])
    same <- vapply(node_ids[-1], function(ids) identical(sort(ids), ref), logical(1))
    if (!all(same)) {
      stop(
        "All layers must share the same node set. Found layers with different ",
        "nodes; align the node universe (adding isolates where needed) before ",
        "fitting, or use allow_unequal_nodes = TRUE with named vertices.",
        call. = FALSE
      )
    }
  } else {
    ## unequal node sets are matched by vertex NAME, so every layer must carry
    ## names; positions would silently pair unrelated nodes ----
    named <- vapply(graph_layers, function(g) !is.null(igraph::V(g)$name), logical(1))
    sizes <- vapply(graph_layers, function(g) as.numeric(igraph::vcount(g)), numeric(1))
    if (!all(named) && length(unique(sizes)) > 1L) {
      stop(
        "allow_unequal_nodes = TRUE requires vertex names on every layer ",
        "(igraph vertex names, or dimnames on adjacency matrices) so that nodes ",
        "can be matched across layers.",
        call. = FALSE
      )
    }
  }

  # assign layer names when no names are present ----
  if (is.null(names(graph_layers))) {
    names(graph_layers) <- paste0("layer_", seq_along(graph_layers))
  }

  # return compiled graph layers ----
  return(graph_layers)
}


#' @title Save and Restore the Global RNG State
#'
#' @description Used by the `seed` argument of the fit functions: the caller
#' saves the state, seeds the RNG for reproducible community detection, and
#' restores the state on exit so seeded detection never disturbs the caller's
#' random number stream (e.g. the bootstrap's resampling draws).
#'
#' @noRd

save_rng_state <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else {
    NULL
  }
}

#' @noRd
restore_rng_state <- function(state) {
  if (is.null(state)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  } else {
    assign(".Random.seed", state, envir = globalenv())
  }
  invisible(NULL)
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
  # weights_a / weights_b are LAYER-WIDE node strengths, so a node's weight
  # counts toward community a only when the node is a MEMBER of a (and
  # likewise for b). Without this restriction two disjoint communities scored
  # 1.0 whenever their members had non-zero strength in both layers.
  if (length(nodes) == 0) {
    weighted_similarity <- 0
  } else {
    node_keys <- as.character(nodes)
    wa <- weights_a[node_keys]
    wb <- weights_b[node_keys]
    wa[is.na(wa)] <- 0
    wb[is.na(wb)] <- 0
    wa[!(nodes %in% a)] <- 0
    wb[!(nodes %in% b)] <- 0
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
  # Nodes with no recorded strength (absent from a layer) contribute 0,
  # mirroring weighted_jaccard_similarity() and the Python implementation.
  inter <- intersect(a, b)
  wa_inter <- weights_a[as.character(inter)]
  wb_inter <- weights_b[as.character(inter)]
  wa_inter[is.na(wa_inter)] <- 0
  wb_inter[is.na(wb_inter)] <- 0
  inter_weight <- sum(pmin(wa_inter, wb_inter))
  wa <- weights_a[as.character(a)]
  wb <- weights_b[as.character(b)]
  wa[is.na(wa)] <- 0
  wb[is.na(wb)] <- 0
  a_weight <- sum(wa)
  b_weight <- sum(wb)
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
    ## key by vertex name when present, so lookups match the node ids that
    ## community_overlap_edges() uses (names for named graphs, positions
    ## otherwise) ----
    nm <- igraph::V(g)$name
    names(strength_vals) <- if (is.null(nm)) as.character(seq_along(strength_vals)) else as.character(nm)
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
    node_weights_by_layer = NULL,
    graph_layers = NULL
  ) {

  # check metric argument ----
  metric <- match.arg(metric)

  # assign appropriate simulation functions ----
  sim_fun <- if (metric == "jaccard") weighted_jaccard else weighted_overlap
  weighted_sim_fun <- if (metric == "jaccard") weighted_jaccard_similarity else weighted_overlap_similarity

  # per-layer node identity: community index-sets are positions in each layer's
  # own graph, so across layers with different node sets they must be translated
  # to node NAMES before overlap is computed. Falls back to positions when the
  # graphs carry no names (e.g. fixed-node simulations, where positions align). ----
  node_ids <- function(l) {
    if (!is.null(graph_layers)) {
      nm <- igraph::V(graph_layers[[l]])$name
      if (!is.null(nm)) return(as.character(nm))
    }
    NULL
  }

  # build overlap edges for each layer link ----
  edge_rows <- vector("list", nrow(layer_links))
  for (i in seq_len(nrow(layer_links))) {

    ## extract layer link information ----
    from_idx <- layer_links$from[i]
    to_idx <- layer_links$to[i]
    layer_weight <- layer_links$weight[i]
    from_comms <- fit[[from_idx]]$communities
    to_comms <- fit[[to_idx]]$communities
    fnm <- node_ids(from_idx); tnm <- node_ids(to_idx)

    ## translate community position-sets to node ids (names when available) ----
    fset <- if (!is.null(fnm)) lapply(from_comms, function(p) fnm[p]) else lapply(from_comms, as.integer)
    tset <- if (!is.null(tnm)) lapply(to_comms,   function(p) tnm[p]) else lapply(to_comms,   as.integer)

    if (is.null(node_weights_by_layer)) {
      ## FAST unweighted path: only co-occurring community pairs contribute a
      ## non-zero overlap. Map each shared node to its (from_comm, to_comm) and
      ## count co-occurrences = intersection sizes -- O(n), not O(C x C). ----
      fsize <- lengths(fset); tsize <- lengths(tset)
      node2f <- stats::setNames(rep(names(fset), fsize), unlist(fset, use.names = FALSE))
      node2t <- stats::setNames(rep(names(tset), tsize), unlist(tset, use.names = FALSE))
      shared <- intersect(unlist(fset, use.names = FALSE), unlist(tset, use.names = FALSE))
      if (length(shared) == 0L) next
      pair <- paste(node2f[shared], node2t[shared], sep = "\r")
      ct <- table(pair)
      pk <- do.call(rbind, strsplit(names(ct), "\r", fixed = TRUE))
      fc_id <- as.integer(pk[, 1]); tc_id <- as.integer(pk[, 2]); inter <- as.integer(ct)
      sz_f <- fsize[as.character(fc_id)]; sz_t <- tsize[as.character(tc_id)]
      sim <- if (metric == "jaccard") inter / (sz_f + sz_t - inter) else inter / pmin(sz_f, sz_t)
      weighted_sim <- sim * layer_weight
      keep <- weighted_sim >= min_similarity & inter > 0L
      if (any(keep)) {
        edge_rows[[i]] <- data.frame(
          from_layer = from_idx, to_layer = to_idx,
          from_community = fc_id[keep], to_community = tc_id[keep],
          similarity = sim[keep], layer_weight = layer_weight,
          weighted_similarity = weighted_sim[keep], stringsAsFactors = FALSE)
      }
    } else {
      ## weighted path: only iterate co-occurring pairs (name-based), call the
      ## weighted similarity on the node-id sets. ----
      wf <- node_weights_by_layer[[from_idx]]; wt <- node_weights_by_layer[[to_idx]]
      rows <- list()
      for (fc in names(fset)) {
        cand <- names(tset)[vapply(tset, function(s) length(intersect(fset[[fc]], s)) > 0, logical(1))]
        for (tc in cand) {
          sim <- weighted_sim_fun(fset[[fc]], tset[[tc]], wf, wt)
          weighted_sim <- sim * layer_weight
          if (weighted_sim >= min_similarity)
            rows[[length(rows) + 1]] <- data.frame(from_layer = from_idx, to_layer = to_idx,
              from_community = as.integer(fc), to_community = as.integer(tc),
              similarity = sim, layer_weight = layer_weight,
              weighted_similarity = weighted_sim, stringsAsFactors = FALSE)
        }
      }
      if (length(rows)) edge_rows[[i]] <- do.call(rbind, rows)
    }
  }
  edge_rows <- edge_rows[!vapply(edge_rows, is.null, logical(1))]

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


#' @title Multislice (Mucha) meta-communities for node-identity coupling
#'
#' @description Node-level second stage for the identity specification:
#' optimises Mucha et al. (2010) multislice modularity
#' \deqn{Q = \frac{1}{2\mu}\sum_{ijsr}\Big[\big(A_{ijs} - \gamma\,\tfrac{k_{is}k_{js}}{2m_s}\big)\delta_{sr} + \delta_{ij} C_{jsr}\Big]\delta(g_{is}, g_{jr})}
#' with a generalized Louvain (local moving + aggregation, as in GenLouvain).
#' Intra-slice edges use each slice's OWN configuration null model
#' \eqn{k_{is}k_{js}/2m_s}; interlayer identity ties \eqn{C_{jsr}} (weight =
#' layer-link weight times \code{omega}) carry no null term. Running plain
#' single-graph modularity on the stacked supra-graph is NOT equivalent: its
#' pooled null model \eqn{k_i k_j / 2M} is roughly \eqn{T} times too small and
#' penalises cross-slice co-membership, so at \eqn{\omega \le 1} it returns
#' each slice as one community. That was the behaviour of this function before
#' package version 1.2.1.
#'
#' @param graph_layers List of per-layer \code{igraph} objects.
#'
#' @param interlayer_ties Node-level identity ties (columns \code{from_layer},
#'   \code{to_layer}, \code{node}, \code{layer_weight}).
#'
#' @param algorithm Kept for API compatibility; the optimiser is always the
#'   generalized Louvain below (igraph's Louvain/Leiden cannot take a per-slice
#'   null model).
#'
#' @param omega Interlayer coupling strength (Mucha's omega). Multiplies the
#'   interlayer identity-edge weights. omega = 0 decouples the slices
#'   (independent per-slice modularity); large omega forces each node's copies
#'   into one community.
#'
#' @param resolution_parameter Mucha's gamma (resolution of the per-slice null
#'   model). Larger values yield more, smaller communities.
#'
#' @param seed Optional integer; the node visiting order is randomised.
#'
#' @param max_levels,max_passes Safety caps on aggregation levels and local
#'   moving passes per level.
#'
#' @return List with one integer vector per layer giving each node's
#'   meta-community assignment (node order).
#'
#' @references Mucha, P. J., Richardson, T., Macon, K., Porter, M. A., and
#'   Onnela, J.-P. (2010). Community structure in time-dependent, multiscale,
#'   and multiplex networks. \emph{Science} 328(5980): 876-878.
#'
#' @noRd
detect_multislice_communities <- function(
    graph_layers,
    interlayer_ties,
    algorithm = c("louvain", "leiden"),
    omega = 1,
    resolution_parameter = 1,
    seed = NULL,
    max_levels = 20L,
    max_passes = 50L
  ) {

  algorithm <- match.arg(algorithm)
  vkey <- function(l, node) paste0("L", l, "N", node)
  n_layers <- length(graph_layers)

  layer_nodes <- lapply(graph_layers, function(g) {
    nm <- igraph::V(g)$name
    if (is.null(nm)) as.character(seq_len(igraph::vcount(g))) else as.character(nm)
  })

  # supra vertices: one per (layer, node); slice index per vertex ----
  verts <- unlist(lapply(seq_len(n_layers), function(t) vkey(t, layer_nodes[[t]])))
  slice <- unlist(lapply(seq_len(n_layers), function(t) rep(t, length(layer_nodes[[t]]))))
  N <- length(verts); vid <- stats::setNames(seq_len(N), verts)

  # intra-slice edges (each slice's own adjacency) ----
  intra <- do.call(rbind, lapply(seq_len(n_layers), function(t) {
    g <- graph_layers[[t]]
    el <- igraph::as_edgelist(g, names = TRUE)
    if (nrow(el) == 0) return(NULL)
    w <- igraph::E(g)$weight
    if (is.null(w)) w <- rep(1, nrow(el))
    keep <- el[, 1] != el[, 2]                                   # drop self loops
    data.frame(from = vid[vkey(t, el[keep, 1])], to = vid[vkey(t, el[keep, 2])],
               w = w[keep], stringsAsFactors = FALSE)
  }))

  # per-slice degree k_is (intra edges only) and 2m_s ----
  k <- numeric(N)
  if (!is.null(intra) && nrow(intra) > 0) {
    kf <- tapply(intra$w, intra$from, sum); kt <- tapply(intra$w, intra$to, sum)
    k[as.integer(names(kf))] <- k[as.integer(names(kf))] + kf
    k[as.integer(names(kt))] <- k[as.integer(names(kt))] + kt
  }
  twom <- as.numeric(tapply(k, factor(slice, levels = seq_len(n_layers)), sum))
  twom[is.na(twom)] <- 0

  # interlayer identity ties: weight = layer weight * omega, no null term ----
  inter <- NULL
  if (!is.null(interlayer_ties) && nrow(interlayer_ties) > 0 && omega > 0) {
    w <- interlayer_ties$layer_weight
    if (is.null(w)) w <- rep(1, nrow(interlayer_ties))
    inter <- data.frame(from = vid[vkey(interlayer_ties$from_layer, interlayer_ties$node)],
                        to   = vid[vkey(interlayer_ties$to_layer,   interlayer_ties$node)],
                        w = w * omega, stringsAsFactors = FALSE)
    inter <- inter[!is.na(inter$from) & !is.na(inter$to), , drop = FALSE]
  }
  edges <- rbind(intra, inter)
  if (is.null(edges) || nrow(edges) == 0) {
    return(lapply(seq_len(n_layers), function(t) as.integer(vid[vkey(t, layer_nodes[[t]])])))
  }

  # K: N x S matrix of per-slice degree (node i in slice s has k_is in column s) ----
  K <- matrix(0, N, n_layers); K[cbind(seq_len(N), slice)] <- k

  if (!is.null(seed)) { rng <- save_rng_state(); on.exit(restore_rng_state(rng), add = TRUE); set.seed(seed) }
  membership <- genlouvain_multislice(edges, K, twom, gamma = resolution_parameter,
                                      max_levels = max_levels, max_passes = max_passes)
  names(membership) <- verts

  # map back to per-layer node order ----
  lapply(seq_len(n_layers), function(t) as.integer(membership[vkey(t, layer_nodes[[t]])]))
}


#' @title Generalized Louvain for multislice modularity
#'
#' @description Optimises \eqn{\sum_{ij} [w_{ij} - \gamma \sum_s K_{is} K_{js} / 2m_s]\,\delta(g_i, g_j)}
#' over a weighted undirected graph whose vertices carry a per-slice degree
#' vector (rows of \code{K}). Local moving (each vertex joins the neighbouring
#' community with the largest gain) alternates with aggregation (communities
#' become vertices whose degree vectors and edge weights are summed) until no
#' move improves the objective. Gains are computed exactly; the \eqn{1/2\mu}
#' normalisation is a constant and is dropped.
#'
#' @param edges data.frame with integer columns \code{from}, \code{to} and
#'   numeric \code{w}; each undirected edge once.
#' @param K numeric matrix (vertices x slices) of per-slice degrees.
#' @param twom numeric vector of per-slice total degree \eqn{2m_s}.
#' @param gamma resolution.
#' @param max_levels,max_passes safety caps.
#' @param n_starts number of random restarts (random vertex visiting order);
#'   the partition with the highest objective is returned.
#'
#' @return Integer membership vector for the original vertices (1..C).
#'
#' @noRd
genlouvain_multislice <- function(edges, K, twom, gamma = 1, max_levels = 20L, max_passes = 50L, n_starts = 3L) {
  inv2m <- ifelse(twom > 0, 1 / twom, 0)
  # objective (unnormalised): sum of within-community edge weight minus
  # gamma * sum_s sum_c Ktot[c, s]^2 / 2m_s ----
  quality <- function(memb) {
    within <- sum(edges$w[memb[edges$from] == memb[edges$to]])
    Ktot <- rowsum(K, memb)
    within - gamma * sum(sweep(Ktot^2, 2, inv2m, "*"))
  }
  best <- NULL; best_q <- -Inf
  for (start in seq_len(n_starts)) {
    memb <- genlouvain_multislice_once(edges, K, twom, gamma, max_levels, max_passes)
    q <- quality(memb)
    if (q > best_q) { best_q <- q; best <- memb }
  }
  best
}

genlouvain_multislice_once <- function(edges, K, twom, gamma, max_levels, max_passes) {
  inv2m <- ifelse(twom > 0, 1 / twom, 0)
  N0 <- nrow(K)
  assign0 <- seq_len(N0)                      # original vertex -> current level vertex
  ef <- as.integer(edges$from); et <- as.integer(edges$to); ew <- as.numeric(edges$w)
  Kcur <- K

  for (level in seq_len(max_levels)) {
    N <- nrow(Kcur)
    # symmetric adjacency lists ----
    nb_from <- c(ef, et); nb_to <- c(et, ef); nb_w <- c(ew, ew)
    ord <- order(nb_from); nb_from <- nb_from[ord]; nb_to <- nb_to[ord]; nb_w <- nb_w[ord]
    starts <- c(1L, cumsum(tabulate(nb_from, N)) + 1L)     # CSR offsets, length N+1
    comm <- seq_len(N)
    Ktot <- Kcur                                          # community x slice degree totals
    KS <- Kcur * matrix(inv2m, N, ncol(Kcur), byrow = TRUE) # k_is / 2m_s, precomputed
    moved_any <- FALSE
    for (pass in seq_len(max_passes)) {
      moved <- 0L
      for (v in sample.int(N)) {
        a <- starts[v]; b <- starts[v + 1L] - 1L
        if (b < a) next
        nbr <- nb_to[a:b]; wv <- nb_w[a:b]
        cv <- comm[v]
        # remove v from its community ----
        Ktot[cv, ] <- Ktot[cv, ] - Kcur[v, ]
        # weight from v to each neighbouring community ----
        nc <- comm[nbr]
        wc <- rowsum(wv, nc, reorder = FALSE)
        cand <- as.integer(rownames(wc))
        # gain of joining community c: w_vc - gamma * sum_s (k_vs / 2m_s) * Ktot[c, s] ----
        gain <- as.numeric(wc) - gamma * as.numeric(Ktot[cand, , drop = FALSE] %*% KS[v, ])
        ib <- which.max(gain); bg <- gain[ib]; best <- cand[ib]
        icv <- match(cv, cand)
        # staying put has gain 0 when the current community holds no neighbour,
        # else its own gain; move only on a strict improvement ----
        stay <- if (is.na(icv)) 0 else gain[icv]
        if (bg <= stay + 1e-12) best <- cv
        Ktot[best, ] <- Ktot[best, ] + Kcur[v, ]
        if (best != cv) { comm[v] <- best; moved <- moved + 1L }
      }
      if (moved > 0L) moved_any <- TRUE else break
    }
    # relabel 1..C ----
    comm <- match(comm, sort(unique(comm)))
    C <- max(comm)
    assign0 <- comm[assign0]
    if (!moved_any || C == N) break
    # aggregate: communities become vertices ----
    Kcur <- rowsum(Kcur, comm, reorder = TRUE)
    cf <- comm[ef]; ct <- comm[et]
    keep <- cf != ct                                      # internal weight is constant for later gains
    if (!any(keep)) break
    key <- ifelse(cf[keep] < ct[keep], paste(cf[keep], ct[keep]), paste(ct[keep], cf[keep]))
    agg <- tapply(ew[keep], key, sum)
    pr <- do.call(rbind, strsplit(names(agg), " ", fixed = TRUE))
    ef <- as.integer(pr[, 1]); et <- as.integer(pr[, 2]); ew <- as.numeric(agg)
  }
  as.integer(assign0)
}
