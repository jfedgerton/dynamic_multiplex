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

    # reject layers without vertices ----
    if (igraph::vcount(layer_graph) == 0) {
      stop(sprintf("Layer %d has no vertices. Drop empty layers before fitting.", i), call. = FALSE)
    }

    # input policy (same in R and Python): within-layer self-loops and
    # zero-weight edges are dropped, so an igraph input and the equivalent
    # adjacency matrix give the same graph ----
    loops <- which(igraph::which_loop(layer_graph))
    if (length(loops)) layer_graph <- igraph::delete_edges(layer_graph, loops)
    w <- igraph::E(layer_graph)$weight
    if (!is.null(w)) {
      zero <- which(!is.na(w) & w == 0)
      if (length(zero)) layer_graph <- igraph::delete_edges(layer_graph, zero)
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
    ## (1.4.0) checked whatever the layer sizes: nodes can enter and leave
    ## while the size stays the same ----
    named <- vapply(graph_layers, function(g) !is.null(igraph::V(g)$name), logical(1))
    if (!all(named)) {
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
    algorithm = c("leiden", "louvain"),
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

  # directed layers with directed = FALSE: collapsed per layer below (before
  # 1.3.2 igraph stopped with an error here) ----
  if (!directed && any(vapply(graph_layers, igraph::is_directed, logical(1)))) {
    warning(
      "Directed layers were passed with directed = FALSE; collapsing them to ",
      "weighted undirected graphs (reciprocal arc weights summed) for ",
      "per-layer detection.",
      call. = FALSE
    )
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
      if (igraph::is_directed(g_input)) {
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
      if (igraph::is_directed(g_input)) {
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
    algorithm = c("leiden", "louvain"),
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


# ---------------------------------------------------------------------------
# 1.4.0: chain repair and filtered (past-only) tracking for the two-stage fits
# (mirrors python_code/src/dynamic_multiplex/multilayer_utils.py)
# ---------------------------------------------------------------------------

.parse_key <- function(k) {
  m <- regmatches(k, regexec("^L([0-9]+)C([0-9]+)$", k))
  cbind(layer = as.integer(vapply(m, `[`, "", 2)), comm = as.integer(vapply(m, `[`, "", 3)))
}

# per-layer node membership (layer vertex order) from supernode meta ids
.membership_from_meta <- function(layer_communities, meta_ids) {
  lapply(seq_along(layer_communities), function(t) {
    mem <- layer_communities[[t]]$membership
    as.integer(meta_ids[paste0("L", t, "C", as.integer(mem))])
  })
}

# relabel meta ids 1..K in order of first appearance (layer, then community id)
.renumber_first_appearance <- function(meta_ids) {
  pk <- .parse_key(names(meta_ids))
  o <- order(pk[, "layer"], pk[, "comm"])
  ids <- meta_ids[o]
  out <- integer(length(meta_ids))
  out[o] <- match(ids, unique(ids))
  stats::setNames(out, names(meta_ids))
}

# best match of every supernode into every linked layer, by similarity:
# key "layer|comm|other_layer" -> c(similarity, other_comm); ties -> lower id
.best_matches <- function(pair_ties) {
  best <- new.env(hash = TRUE)
  if (is.null(pair_ties) || !nrow(pair_ties)) return(best)
  for (i in seq_len(nrow(pair_ties))) {
    al <- pair_ties$from_layer[i]; ac <- pair_ties$from_community[i]
    bl <- pair_ties$to_layer[i]; bc <- pair_ties$to_community[i]
    if (al == bl) next
    s <- pair_ties$similarity[i]
    for (d in 1:2) {
      xl <- if (d == 1) al else bl; xc <- if (d == 1) ac else bc
      yl <- if (d == 1) bl else al; yc <- if (d == 1) bc else ac
      key <- paste(xl, xc, yl, sep = "|")
      cur <- best[[key]]
      if (is.null(cur) || s > cur[1] || (s == cur[1] && yc < cur[2])) assign(key, c(s, yc), envir = best)
    }
  }
  best
}

.check_repair_args <- function(tracking, repair_chains, repair_min_similarity, seed = 123) {
  if (!is.character(tracking) || length(tracking) != 1L || !tracking %in% c("smoothed", "filtered"))
    stop('`tracking` must be "smoothed" or "filtered".', call. = FALSE)
  # Unseeded stage 1 on later layers advances the RNG that filtered stage 2
  # draws from, which would break append-invariance.
  if (identical(tracking, "filtered") && is.null(seed))
    stop('`tracking = "filtered"` needs a seed; `seed = NULL` is not supported.', call. = FALSE)
  if (!is.logical(repair_chains) || length(repair_chains) != 1L || is.na(repair_chains))
    stop("`repair_chains` must be TRUE or FALSE.", call. = FALSE)
  if (!is.numeric(repair_min_similarity) || length(repair_min_similarity) != 1L ||
      is.na(repair_min_similarity) || repair_min_similarity < 0 || repair_min_similarity > 1)
    stop("`repair_min_similarity` must be between 0 and 1.", call. = FALSE)
  invisible(TRUE)
}

.empty_repair_log <- function() {
  data.frame(from_layer = integer(0), from_community = integer(0), to_layer = integer(0),
             to_community = integer(0), similarity = numeric(0))
}

#' @title Rejoin meta-communities cut along time (chain repair)
#'
#' @description Stage-2 modularity tends to cut a long-lived community into
#' time segments (a resolution-limit effect on chains of communities). For
#' each layer link (s, t) with s < t, in increasing order, and each community
#' a in layer s (by community id), the meta-community of community b in layer
#' t is joined to that of a when (i) a and b are each other's best match
#' across (s, t) by the fit's own similarity (mutual best match, as in
#' Kawadia and Sreenivasan 2012), (ii) their similarity is at least
#' \code{min_similarity} (a matching floor, as in Greene et al. 2010), and
#' (iii) a's meta-community has no community after s and b's none before t,
#' so groups that coexist are never merged.
#'
#' @param meta_ids Named integer vector of meta ids, names "L<layer>C<community>".
#' @param pair_ties Interlayer ties without self-loops (columns from_layer,
#'   to_layer, from_community, to_community, similarity).
#' @param layer_links Layer links data frame (from, to).
#' @param min_similarity Similarity floor in [0, 1].
#' @return A list with \code{meta_ids} (repaired, renumbered 1..K by first
#'   appearance) and \code{log} (one row per join).
#' @noRd
repair_meta_chains <- function(meta_ids, pair_ties, layer_links, min_similarity = 0.2) {
  ids <- unique(as.integer(meta_ids))
  root <- stats::setNames(ids, as.character(ids))
  pk <- .parse_key(names(meta_ids))
  layers_of <- lapply(split(pk[, "layer"], as.integer(meta_ids)), unique)
  find <- function(m) {
    m <- as.character(m)
    while (root[[m]] != as.integer(m)) m <- as.character(root[[m]])
    as.integer(m)
  }
  best <- .best_matches(pair_ties)
  lp <- unique(t(apply(cbind(as.integer(layer_links$from), as.integer(layer_links$to)), 1, sort)))
  lp <- lp[lp[, 1] != lp[, 2], , drop = FALSE]
  lp <- lp[order(lp[, 1], lp[, 2]), , drop = FALSE]
  log <- list()
  for (r in seq_len(nrow(lp))) {
    s <- lp[r, 1]; t <- lp[r, 2]
    for (ca in sort(pk[pk[, "layer"] == s, "comm"])) {
      fwd <- best[[paste(s, ca, t, sep = "|")]]
      if (is.null(fwd)) next
      sim <- fwd[1]; cb <- fwd[2]
      bwd <- best[[paste(t, cb, s, sep = "|")]]
      if (is.null(bwd) || bwd[2] != ca || sim < min_similarity) next
      ra <- find(meta_ids[[paste0("L", s, "C", ca)]]); rb <- find(meta_ids[[paste0("L", t, "C", cb)]])
      if (ra == rb) next
      if (max(layers_of[[as.character(ra)]]) == s && min(layers_of[[as.character(rb)]]) == t) {
        root[[as.character(rb)]] <- ra
        layers_of[[as.character(ra)]] <- union(layers_of[[as.character(ra)]], layers_of[[as.character(rb)]])
        layers_of[[as.character(rb)]] <- NULL
        log[[length(log) + 1L]] <- data.frame(from_layer = as.integer(s), from_community = as.integer(ca),
                                              to_layer = as.integer(t), to_community = as.integer(cb),
                                              similarity = sim)
      }
    }
  }
  repaired <- stats::setNames(vapply(as.integer(meta_ids), find, integer(1)), names(meta_ids))
  list(meta_ids = .renumber_first_appearance(repaired),
       log = if (length(log)) do.call(rbind, log) else .empty_repair_log())
}

# stage-2 solver call on a weighted graph (same igraph calls as
# detect_interlayer_communities); nodes = character names
.cluster_weighted <- function(nodes, edge_df, algorithm, resolution_parameter) {
  g <- igraph::graph_from_data_frame(d = edge_df, directed = FALSE,
                                     vertices = data.frame(name = nodes, stringsAsFactors = FALSE))
  cl <- if (algorithm == "louvain") {
    igraph::cluster_louvain(graph = g, weights = igraph::E(g)$weight, resolution = resolution_parameter)
  } else {
    igraph::cluster_leiden(graph = g, objective_function = "modularity", weights = igraph::E(g)$weight,
                           resolution = resolution_parameter, n_iterations = 3)
  }
  stats::setNames(as.integer(igraph::membership(cl)), igraph::V(g)$name)
}

#' @title Past-only (filtered) second stage
#'
#' @description The meta id of a layer-t community depends only on layers
#' 1..t. At each layer t the package's community graph is built on layers 1..t
#' from the links whose ends are both at or before t; every meta-community
#' already assigned in layers before t is contracted into one read-only anchor
#' (contraction keeps degrees and internal weights) and the graph is clustered
#' with the same Leiden or Louvain call as the smoothed stage 2, reseeding
#' with \code{seed} at every layer. A layer-t community takes the id of the
#' anchor in its cluster with the strongest direct tie to it (then the
#' strongest tie to the cluster's other layer-t communities, then the lower
#' id); a cluster without an anchor starts a new id. Anchors are never
#' relabelled. With \code{repair_chains}, a layer-t community that started a
#' new id is joined to a group that ended at a linked earlier layer s under
#' the rule of \code{repair_meta_chains}, using layers s and t only.
#' @noRd
detect_interlayer_communities_filtered <- function(
    layer_communities, pair_ties, layer_links, algorithm = c("leiden", "louvain"),
    resolution_parameter = 1, seed = 123, add_self_loops = TRUE, self_loop_multiplier = 1,
    min_similarity = 0, directed = FALSE, repair_chains = FALSE, repair_min_similarity = 0.2) {

  algorithm <- match.arg(algorithm)
  n_layers <- length(layer_communities)
  comms1 <- sort(as.integer(names(layer_communities[[1]]$communities)))
  meta <- stats::setNames(seq_along(comms1), paste0("L1C", comms1))
  next_id <- length(meta) + 1L
  last_layer <- stats::setNames(rep(1L, length(meta)), as.character(meta))
  log <- list()
  has_ties <- !is.null(pair_ties) && nrow(pair_ties) > 0
  best <- if (repair_chains && has_ties) .best_matches(pair_ties) else new.env()

  for (t in 2:n_layers) {
    if (n_layers < 2) break
    lk <- layer_links[layer_links$from <= t & layer_links$to <= t, , drop = FALSE]
    pt <- if (has_ties) pair_ties[pair_ties$from_layer <= t & pair_ties$to_layer <= t, , drop = FALSE] else pair_ties
    ties <- pt
    if (isTRUE(add_self_loops) && nrow(lk) > 0) {
      ties <- add_community_self_loops(edge_df = pt, fit = layer_communities[seq_len(t)], layer_links = lk,
                                       self_loop_multiplier = self_loop_multiplier,
                                       min_similarity = min_similarity, directed = directed)
    }
    cur_c <- sort(as.integer(names(layer_communities[[t]]$communities)))
    cur <- paste0("L", t, "C", cur_c)
    start_id <- next_id

    node_of <- function(k) {
      lay <- .parse_key(k)[, "layer"]
      ifelse(lay < t, paste0("A", meta[k]), k)
    }
    anchors <- paste0("A", sort(unique(as.integer(meta))))
    nodes <- c(anchors, cur)
    ord <- stats::setNames(seq_along(nodes), nodes)
    edges <- NULL
    if (!is.null(ties) && nrow(ties) > 0) {
      u <- node_of(paste0("L", ties$from_layer, "C", ties$from_community))
      v <- node_of(paste0("L", ties$to_layer, "C", ties$to_community))
      swap <- ord[u] > ord[v]
      a <- ifelse(swap, v, u); b <- ifelse(swap, u, v)
      agg <- stats::aggregate(w ~ a + b, data = data.frame(a = a, b = b, w = ties$weighted_similarity,
                                                           stringsAsFactors = FALSE), FUN = sum)
      agg <- agg[order(ord[agg$a], ord[agg$b]), , drop = FALSE]
      edges <- data.frame(from = agg$a, to = agg$b, weight = agg$w, stringsAsFactors = FALSE)
    }

    if (is.null(edges) || nrow(edges) == 0) {
      for (k in cur) { meta[[k]] <- next_id; next_id <- next_id + 1L }
    } else {
      if (!is.null(seed)) set.seed(seed)
      mem <- .cluster_weighted(nodes, edges, algorithm, resolution_parameter)
      wkey <- paste(edges$from, edges$to)
      w_between <- function(x, y) {
        k <- if (ord[[x]] <= ord[[y]]) paste(x, y) else paste(y, x)
        i <- match(k, wkey); if (is.na(i)) 0 else edges$weight[i]
      }
      clusters <- split(nodes, mem[nodes])
      groups <- list()
      for (members in clusters) {
        frees <- members[!startsWith(members, "A")]
        if (!length(frees)) next
        groups[[length(groups) + 1L]] <- list(first = min(.parse_key(frees)[, "comm"]),
                                              anc = members[startsWith(members, "A")], frees = frees)
      }
      groups <- groups[order(vapply(groups, `[[`, numeric(1), "first"))]
      for (gr in groups) {
        if (!length(gr$anc)) {
          for (n in gr$frees) meta[[n]] <- next_id
          next_id <- next_id + 1L
          next
        }
        for (n in gr$frees) {
          sc <- t(vapply(gr$anc, function(a) c(w_between(a, n), sum(vapply(gr$frees, function(f) w_between(a, f), numeric(1))),
                                               -as.numeric(sub("^A", "", a))), numeric(3)))
          o <- order(-sc[, 1], -sc[, 2], -sc[, 3])
          meta[[n]] <- as.integer(sub("^A", "", gr$anc[o[1]]))
        }
      }
    }

    if (repair_chains && nrow(lk) > 0) {
      lp <- unique(t(apply(cbind(as.integer(lk$from), as.integer(lk$to)), 1, sort)))
      lp <- lp[lp[, 1] != lp[, 2] & lp[, 2] == t, , drop = FALSE]
      for (s in sort(lp[, 1])) {
        for (k in cur) {
          if (meta[[k]] < start_id) next                       # only groups that start at t
          cb <- .parse_key(k)[, "comm"]
          bwd <- best[[paste(t, cb, s, sep = "|")]]
          if (is.null(bwd)) next
          sim <- bwd[1]; ca <- bwd[2]
          fwd <- best[[paste(s, ca, t, sep = "|")]]
          if (is.null(fwd) || fwd[2] != cb || sim < repair_min_similarity) next
          ma <- meta[[paste0("L", s, "C", ca)]]
          # a's group must have ended at s: last seen at s and absent from t
          if (!identical(unname(last_layer[as.character(ma)]), as.integer(s)) || any(meta[cur] == ma)) next
          mb <- meta[[k]]
          meta[cur][meta[cur] == mb] <- ma
          log[[length(log) + 1L]] <- data.frame(from_layer = as.integer(s), from_community = as.integer(ca),
                                                to_layer = as.integer(t), to_community = as.integer(cb),
                                                similarity = sim)
        }
      }
    }
    last_layer[as.character(meta[cur])] <- as.integer(t)
  }

  meta <- .renumber_first_appearance(meta)
  list(meta_ids = meta, membership = .membership_from_meta(layer_communities, meta),
       chain_repairs = if (length(log)) do.call(rbind, log) else .empty_repair_log())
}

# Second stage shared by the four two-stage fitters. tracking = "smoothed" with
# repair_chains = FALSE is the 1.3.2 code path, unchanged.
run_two_stage <- function(fit, pair_ties, links, add_self_loops, self_loop_multiplier, min_similarity,
                          directed, algorithm, resolution_parameter, seed, tracking, repair_chains,
                          repair_min_similarity) {
  .check_repair_args(tracking, repair_chains, repair_min_similarity, seed)
  interlayer_ties <- pair_ties
  if (isTRUE(add_self_loops)) {
    interlayer_ties <- add_community_self_loops(edge_df = pair_ties, fit = fit, layer_links = links,
                                                self_loop_multiplier = self_loop_multiplier,
                                                min_similarity = min_similarity, directed = directed)
  }
  if (tracking == "smoothed") {
    meta <- detect_interlayer_communities(layer_communities = fit, interlayer_ties = interlayer_ties,
                                          algorithm = algorithm, resolution_parameter = resolution_parameter)
    repairs <- .empty_repair_log()
    if (repair_chains && !is.null(meta$meta_ids)) {
      rep <- repair_meta_chains(meta$meta_ids, pair_ties, links, repair_min_similarity)
      meta <- list(meta_ids = rep$meta_ids, membership = .membership_from_meta(fit, rep$meta_ids))
      repairs <- rep$log
    }
  } else {
    meta <- detect_interlayer_communities_filtered(
      fit, pair_ties, links, algorithm = algorithm, resolution_parameter = resolution_parameter,
      seed = seed, add_self_loops = add_self_loops, self_loop_multiplier = self_loop_multiplier,
      min_similarity = min_similarity, directed = directed, repair_chains = repair_chains,
      repair_min_similarity = repair_min_similarity)
    repairs <- meta$chain_repairs
  }
  list(meta_ids = meta$meta_ids, membership = meta$membership, interlayer_ties = interlayer_ties,
       chain_repairs = repairs, tracking = tracking, repair_chains = isTRUE(repair_chains),
       repair_min_similarity = as.numeric(repair_min_similarity))
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
#' @param directed Logical. If \code{TRUE}, every layer must be a directed
#'   graph and each slice keeps its arc directions: the slice term uses Leicht
#'   and Newman's directed null model \eqn{k^{out}_{is} k^{in}_{js} / m_s}, so a
#'   directed network in which every arc is reciprocated gives the same fit as
#'   its undirected twin at the same \code{omega}. The interlayer identity
#'   ties stay undirected (the coupling enters only through
#'   \eqn{C_{jsr} + C_{jrs}}). If \code{FALSE} (default), directed layers are
#'   symmetrised (\eqn{A + A^T}, a reciprocated pair counts twice) with a
#'   warning, as before 1.3.2.
#'
#' @param max_levels,max_passes Safety caps on aggregation levels and local
#'   moving passes per level.
#'
#' @return List with one integer vector per layer giving each node's
#'   meta-community assignment (node order).
#'
#' @references Leicht, E. A., and Newman, M. E. J. (2008). Community structure
#'   in directed networks. \emph{Physical Review Letters} 100: 118703.
#'
#'   Mucha, P. J., Richardson, T., Macon, K., Porter, M. A., and
#'   Onnela, J.-P. (2010). Community structure in time-dependent, multiscale,
#'   and multiplex networks. \emph{Science} 328(5980): 876-878.
#'
#' @noRd
detect_multislice_communities <- function(
    graph_layers,
    interlayer_ties,
    algorithm = c("leiden", "louvain"),
    omega = 1,
    resolution_parameter = 1,
    seed = NULL,
    max_levels = 20L,
    max_passes = 50L,
    directed = FALSE
  ) {

  algorithm <- match.arg(algorithm)
  vkey <- function(l, node) paste0("L", l, "N", node)
  n_layers <- length(graph_layers)

  # directedness: all slices directed (directed null model) or none ----
  is_dir <- vapply(graph_layers, igraph::is_directed, logical(1))
  if (isTRUE(directed) && !all(is_dir)) {
    stop("directed = TRUE requires every layer to be a directed graph; layers ",
         paste(which(!is_dir), collapse = ", "), " are undirected. Pass directed ",
         "graphs or adjacency matrices with directed = TRUE, or use directed = FALSE.",
         call. = FALSE)
  }
  if (!isTRUE(directed) && any(is_dir)) {
    warning("Directed layers were passed with directed = FALSE: the multislice fit ",
            "symmetrised them (A + A^T, a reciprocated pair counts twice). Use ",
            "directed = TRUE for the directed null model.", call. = FALSE)
  }

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

  # directed slices: out- and in-degree per slice and m_s = total arc weight.
  # The optimiser works on half the ordered-pair objective, where an arc
  # counts once and an undirected tie (interlayer identity tie) twice, so each
  # arc enters the edge list with weight w / 2 ----
  if (isTRUE(directed)) {
    k_out <- numeric(N); k_in <- numeric(N)
    if (!is.null(intra) && nrow(intra) > 0) {
      kf <- tapply(intra$w, intra$from, sum); kt <- tapply(intra$w, intra$to, sum)
      k_out[as.integer(names(kf))] <- kf
      k_in[as.integer(names(kt))] <- kt
      intra$w <- intra$w / 2
    }
    m_dir <- as.numeric(tapply(k_out, factor(slice, levels = seq_len(n_layers)), sum))
    m_dir[is.na(m_dir)] <- 0
    stopifnot(isTRUE(all.equal(sum(k_out), sum(k_in))), all(m_dir >= 0),
              isTRUE(all.equal(sum(m_dir), sum(k_out))))
  }

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
  if (isTRUE(directed)) {
    K_out <- matrix(0, N, n_layers); K_out[cbind(seq_len(N), slice)] <- k_out
    K_in <- matrix(0, N, n_layers); K_in[cbind(seq_len(N), slice)] <- k_in
    membership <- genlouvain_multislice(edges, K_out, m_dir, gamma = resolution_parameter,
                                        max_levels = max_levels, max_passes = max_passes,
                                        K_in = K_in)
  } else {
    membership <- genlouvain_multislice(edges, K, twom, gamma = resolution_parameter,
                                        max_levels = max_levels, max_passes = max_passes)
  }
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
#'   the partition with the highest objective is returned. Default 10 (was 3
#'   in 1.3.0): with 3 restarts the joint NMI of a fit on small networks
#'   (n <= 100) at an abrupt change-point varied by up to 0.2 across runs.
#'
#' @param K_in Optional in-degree matrix (vertices x slices). If supplied the
#'   slices are directed: \code{K} holds out-degrees, \code{twom} holds the
#'   per-slice total arc weight \eqn{m_s}, the null term is
#'   \eqn{\gamma \sum_s K^{out}_{is} K^{in}_{js} / m_s}, and each arc enters
#'   \code{edges} with half its weight.
#'
#' @return Integer membership vector for the original vertices (1..C).
#'
#' @noRd
genlouvain_multislice <- function(edges, K, twom, gamma = 1, max_levels = 20L, max_passes = 50L, n_starts = 10L,
                                  K_in = NULL) {
  if (!is.null(K_in)) {
    stopifnot(identical(dim(K_in), dim(K)), all(is.finite(K_in)), all(is.finite(K)))
    inv_m <- ifelse(twom > 0, 1 / twom, 0)
    # objective (unnormalised directed multislice modularity, ordered pairs):
    # stored edge weights are half the ordered-pair weight, so double them ----
    quality_dir <- function(memb) {
      within <- sum(edges$w[memb[edges$from] == memb[edges$to]])
      2 * within - gamma * sum(sweep(rowsum(K, memb) * rowsum(K_in, memb), 2, inv_m, "*"))
    }
    best <- NULL; best_q <- -Inf
    for (start in seq_len(n_starts)) {
      memb <- genlouvain_multislice_once(edges, K, twom, gamma, max_levels, max_passes, K_in = K_in)
      q <- quality_dir(memb)
      if (q > best_q) { best_q <- q; best <- memb }
    }
    return(best)
  }
  inv2m <- ifelse(twom > 0, 1 / twom, 0)
  # objective (unnormalised multislice modularity, 2 mu Q): every undirected
  # edge is stored once but counts twice in sum_ij A_ij, so the within-community
  # weight is doubled before the null term gamma * sum_s sum_c Ktot[c, s]^2 / 2m_s
  # is subtracted. (Before 1.3.2 the factor 2 was missing, which ranked the
  # restarts as if gamma were 2; the local moves were unaffected.) ----
  quality <- function(memb) {
    within <- sum(edges$w[memb[edges$from] == memb[edges$to]])
    Ktot <- rowsum(K, memb)
    2 * within - gamma * sum(sweep(Ktot^2, 2, inv2m, "*"))
  }
  best <- NULL; best_q <- -Inf
  for (start in seq_len(n_starts)) {
    memb <- genlouvain_multislice_once(edges, K, twom, gamma, max_levels, max_passes)
    q <- quality(memb)
    if (q > best_q) { best_q <- q; best <- memb }
  }
  best
}

genlouvain_multislice_once <- function(edges, K, twom, gamma, max_levels, max_passes, K_in = NULL) {
  if (!is.null(K_in)) return(genlouvain_multislice_once_directed(edges, K, K_in, twom, gamma, max_levels, max_passes))
  inv2m <- ifelse(twom > 0, 1 / twom, 0)
  N0 <- nrow(K)
  ef0 <- as.integer(edges$from); et0 <- as.integer(edges$to); ew0 <- as.numeric(edges$w)
  relabel <- function(x) match(x, sort(unique(x)))
  aggregate_edges <- function(ef, et, ew, comm) {           # community graph, internal weight dropped
    cf <- comm[ef]; ct <- comm[et]; keep <- cf != ct
    if (!any(keep)) return(list(ef = integer(0), et = integer(0), ew = numeric(0)))
    a <- pmin(cf[keep], ct[keep]); b <- pmax(cf[keep], ct[keep])
    key <- a * (max(comm) + 1) + b
    agg <- rowsum(ew[keep], key, reorder = TRUE); k <- as.numeric(rownames(agg))
    list(ef = as.integer(k %/% (max(comm) + 1)), et = as.integer(k %% (max(comm) + 1)), ew = as.numeric(agg))
  }
  comm <- seq_len(N0)
  for (outer in seq_len(max_levels)) {
    # node-level local moving from the current partition (this is the refinement
    # step: a node whose copies were merged into one community at a coarser
    # level can leave it once the slices are looked at one node at a time) ----
    comm <- relabel(genlouvain_local_moving(ef0, et0, ew0, K, inv2m, gamma, comm, max_passes))
    # multilevel aggregation from that partition ----
    assign <- comm; Kc <- rowsum(K, comm, reorder = TRUE); E <- aggregate_edges(ef0, et0, ew0, comm)
    for (level in seq_len(max_levels)) {
      if (!length(E$ef)) break
      c2 <- relabel(genlouvain_local_moving(E$ef, E$et, E$ew, Kc, inv2m, gamma, seq_len(nrow(Kc)), max_passes))
      if (max(c2) == nrow(Kc)) break                        # no community merged: stable
      assign <- c2[assign]; Kc <- rowsum(Kc, c2, reorder = TRUE); E <- aggregate_edges(E$ef, E$et, E$ew, c2)
    }
    if (identical(assign, comm)) break                     # refinement moved nothing and nothing merged
    comm <- assign
  }
  as.integer(comm)
}

# One round of local moving: every vertex joins the neighbouring community
# with the largest gain, repeated until no vertex moves (or max_passes).
# K is vertices x slices; comm_init the starting partition.
genlouvain_local_moving <- function(ef, et, ew, Kcur, inv2m, gamma, comm_init, max_passes) {
  N <- nrow(Kcur)
  if (!length(ef)) return(comm_init)
  nb_from <- c(ef, et); nb_to <- c(et, ef); nb_w <- c(ew, ew)
  ord <- order(nb_from); nb_from <- nb_from[ord]; nb_to <- nb_to[ord]; nb_w <- nb_w[ord]
  starts <- c(1L, cumsum(tabulate(nb_from, N)) + 1L)       # CSR offsets, length N+1
  comm <- as.integer(comm_init)
  Ktot <- rowsum(Kcur, comm, reorder = TRUE); Ktot <- Ktot[as.character(seq_len(max(comm))), , drop = FALSE]
  Ktot[is.na(Ktot)] <- 0
  KS <- Kcur * matrix(inv2m, N, ncol(Kcur), byrow = TRUE)   # k_is / 2m_s, precomputed
  for (pass in seq_len(max_passes)) {
    moved <- 0L
    for (v in sample.int(N)) {
      a <- starts[v]; b <- starts[v + 1L] - 1L
      if (b < a) next
      nbr <- nb_to[a:b]; wv <- nb_w[a:b]
      cv <- comm[v]
      Ktot[cv, ] <- Ktot[cv, ] - Kcur[v, ]                 # remove v from its community
      nc <- comm[nbr]
      wc <- rowsum(wv, nc, reorder = FALSE)
      cand <- as.integer(rownames(wc))
      # gain of joining c: w_vc - gamma * sum_s (k_vs / 2m_s) * Ktot[c, s] ----
      gain <- as.numeric(wc) - gamma * as.numeric(Ktot[cand, , drop = FALSE] %*% KS[v, ])
      ib <- which.max(gain); bg <- gain[ib]; best <- cand[ib]
      icv <- match(cv, cand)
      # (1.4.0) exact gain of staying in cv even when v has no tie into cv;
      # 0 was too high whenever cv still had other members, so the mover could
      # refuse an improving move (same fix as the directed mover in 1.3.2) ----
      stay <- if (!is.na(icv)) gain[icv] else -gamma * sum(Ktot[cv, ] * KS[v, ])
      if (bg <= stay + 1e-12) best <- cv
      Ktot[best, ] <- Ktot[best, ] + Kcur[v, ]
      if (best != cv) { comm[v] <- best; moved <- moved + 1L }
    }
    if (moved == 0L) break
  }
  comm
}

# Directed-slice variants of the two functions above. K_out / K_in are
# vertices x slices, inv_m = 1 / m_s (total arc weight per slice). Edge
# weights are half the ordered-pair weight (arcs w / 2, identity ties w), so
# the gain of joining community c is the half-scale change in
# sum_ij [A_ijs - gamma k_out_is k_in_js / m_s] + coupling:
#   w_vc - (gamma / 2) sum_s (k_out_vs Kin_cs + k_in_vs Kout_cs) / m_s.
# With every arc reciprocated (k_out = k_in = k, m_s = 2m_s) this equals the
# undirected gain exactly.
genlouvain_multislice_once_directed <- function(edges, K_out, K_in, m_s, gamma, max_levels, max_passes) {
  inv_m <- ifelse(m_s > 0, 1 / m_s, 0)
  N0 <- nrow(K_out)
  ef0 <- as.integer(edges$from); et0 <- as.integer(edges$to); ew0 <- as.numeric(edges$w)
  relabel <- function(x) match(x, sort(unique(x)))
  aggregate_edges <- function(ef, et, ew, comm) {
    cf <- comm[ef]; ct <- comm[et]; keep <- cf != ct
    if (!any(keep)) return(list(ef = integer(0), et = integer(0), ew = numeric(0)))
    a <- pmin(cf[keep], ct[keep]); b <- pmax(cf[keep], ct[keep])
    key <- a * (max(comm) + 1) + b
    agg <- rowsum(ew[keep], key, reorder = TRUE); k <- as.numeric(rownames(agg))
    list(ef = as.integer(k %/% (max(comm) + 1)), et = as.integer(k %% (max(comm) + 1)), ew = as.numeric(agg))
  }
  comm <- seq_len(N0)
  for (outer in seq_len(max_levels)) {
    comm <- relabel(genlouvain_local_moving_directed(ef0, et0, ew0, K_out, K_in, inv_m, gamma, comm, max_passes))
    assign <- comm
    Ko <- rowsum(K_out, comm, reorder = TRUE); Ki <- rowsum(K_in, comm, reorder = TRUE)
    E <- aggregate_edges(ef0, et0, ew0, comm)
    for (level in seq_len(max_levels)) {
      if (!length(E$ef)) break
      c2 <- relabel(genlouvain_local_moving_directed(E$ef, E$et, E$ew, Ko, Ki, inv_m, gamma, seq_len(nrow(Ko)), max_passes))
      if (max(c2) == nrow(Ko)) break
      assign <- c2[assign]
      Ko <- rowsum(Ko, c2, reorder = TRUE); Ki <- rowsum(Ki, c2, reorder = TRUE)
      E <- aggregate_edges(E$ef, E$et, E$ew, c2)
    }
    if (identical(assign, comm)) break
    comm <- assign
  }
  as.integer(comm)
}

genlouvain_local_moving_directed <- function(ef, et, ew, Kout, Kin, inv_m, gamma, comm_init, max_passes) {
  N <- nrow(Kout)
  if (!length(ef)) return(comm_init)
  nb_from <- c(ef, et); nb_to <- c(et, ef); nb_w <- c(ew, ew)
  ord <- order(nb_from); nb_from <- nb_from[ord]; nb_to <- nb_to[ord]; nb_w <- nb_w[ord]
  starts <- c(1L, cumsum(tabulate(nb_from, N)) + 1L)
  comm <- as.integer(comm_init)
  tot <- function(K) {
    Kt <- rowsum(K, comm, reorder = TRUE); Kt <- Kt[as.character(seq_len(max(comm))), , drop = FALSE]
    Kt[is.na(Kt)] <- 0
    Kt
  }
  Ktot_out <- tot(Kout); Ktot_in <- tot(Kin)
  KSo <- Kout * matrix(inv_m, N, ncol(Kout), byrow = TRUE)   # k_out_vs / m_s
  KSi <- Kin * matrix(inv_m, N, ncol(Kin), byrow = TRUE)     # k_in_vs / m_s
  for (pass in seq_len(max_passes)) {
    moved <- 0L
    for (v in sample.int(N)) {
      a <- starts[v]; b <- starts[v + 1L] - 1L
      if (b < a) next
      nbr <- nb_to[a:b]; wv <- nb_w[a:b]
      cv <- comm[v]
      Ktot_out[cv, ] <- Ktot_out[cv, ] - Kout[v, ]
      Ktot_in[cv, ] <- Ktot_in[cv, ] - Kin[v, ]
      nc <- comm[nbr]
      wc <- rowsum(wv, nc, reorder = FALSE)
      cand <- as.integer(rownames(wc))
      null <- as.numeric(Ktot_in[cand, , drop = FALSE] %*% KSo[v, ]) +
              as.numeric(Ktot_out[cand, , drop = FALSE] %*% KSi[v, ])
      gain <- as.numeric(wc) - (gamma / 2) * null
      ib <- which.max(gain); bg <- gain[ib]; best <- cand[ib]
      icv <- match(cv, cand)
      # staying: v's weight to its own community (0 if no edge) minus the null
      # term against the rest of that community; zero only if v is alone ----
      stay <- if (!is.na(icv)) gain[icv] else
        -(gamma / 2) * (sum(Ktot_in[cv, ] * KSo[v, ]) + sum(Ktot_out[cv, ] * KSi[v, ]))
      if (bg <= stay + 1e-12) best <- cv
      Ktot_out[best, ] <- Ktot_out[best, ] + Kout[v, ]
      Ktot_in[best, ] <- Ktot_in[best, ] + Kin[v, ]
      if (best != cv) { comm[v] <- best; moved <- moved + 1L }
    }
    if (moved == 0L) break
  }
  comm
}
