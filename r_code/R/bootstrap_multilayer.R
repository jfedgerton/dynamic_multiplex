#' @title Parametric network bootstrap for the stability score
#'
#' @description Quantifies the uncertainty of multilayer community detection
#' by refitting communities on \code{n_boot} resampled networks, using a
#' parametric network bootstrap: within- and between-community edge
#' probabilities (and, for weighted networks, edge-weight pools) are
#' estimated from the observed network using the point-estimate partition,
#' and each replicate redraws the full edge set from those estimates. This
#' reproduces the variability of fresh data, including which edges exist.
#'
#' Since 1.4.0 the package uses the bootstrap only for the stability score
#' \eqn{s} reported by \code{\link{partition_stability}}, which compares fits
#' of the same network. The interval outputs of earlier versions
#' (\code{co_assignment_ci()}, \code{community_est()}, the calibrated accuracy
#' floor) were removed because their coverage was not reliable.
#'
#' Agreement is measured on the cross-layer \emph{meta-communities}
#' (the tracked partition from the second-stage detection), not the
#' independently-detected per-layer communities. Co-assignment therefore
#' answers "do these two nodes belong to the same persistent community,"
#' and the community count is the number of meta-communities per layer.
#'
#' Versions before 1.1.0 instead used a Bayesian bootstrap on edge weights
#' (Exponential(1) multipliers on a fixed topology). That scheme was
#' removed: because it never varies which edges exist, it understates the
#' variability of fresh data, and in simulation studies confidence
#' intervals built from it undercovered substantially (pairwise
#' co-assignment intervals covered ~45-48 percent at a nominal 95 percent).
#'
#' @param layers List of adjacency matrices or igraph objects.
#'
#' @param fit_type One of \code{"jaccard"}, \code{"overlap"},
#' \code{"weighted_jaccard"}, \code{"weighted_overlap"}, \code{"identity"}.
#'
#' @param algorithm Community detection algorithm: \code{"leiden"} (default)
#' or \code{"louvain"}.
#'
#' @param n_boot Number of bootstrap replicates.
#'
#' @param layer_links Optional data.frame defining layer connectivity.
#'
#' @param min_similarity Minimum weighted similarity for interlayer ties.
#'
#' @param resolution_parameter Resolution parameter for community detection.
#'
#' @param directed Logical; if \code{TRUE}, treat networks as directed.
#'
#' @param seed Optional random seed for reproducibility.
#'
#' @param tracking,repair_chains,repair_min_similarity Passed to the two-stage
#' fit functions (see \code{\link{fit_multilayer_jaccard}}) so that every
#' refit uses the same settings as the point estimate. Not available for
#' \code{fit_type = "identity"}.
#'
#' @param allow_unequal_nodes Logical. If \code{TRUE}, layers may hold
#' different node sets (nodes entering and leaving, with or without a change
#' in layer size). Every layer must then carry vertex names; each layer is
#' redrawn on its own nodes and every refit matches nodes by name. Default
#' \code{FALSE}: all layers must share one node set.
#'
#' @param objective Objective for the first-stage (per-layer) detection,
#' passed to the fit function: `"modularity"` or `"cpm"` (Leiden only).
#' Defaults to `"cpm"` for directed layers and `"modularity"` otherwise.
#'
#' @return A list of class \code{"multilayer_bootstrap"} with components:
#'   \describe{
#'     \item{n_boot}{Number of completed bootstrap replicates.}
#'     \item{co_assignment}{Per-layer co-assignment probability matrices on
#'       the meta-communities (probability two nodes share a persistent
#'       community).}
#'     \item{node_stability}{Per-layer vectors giving the fraction of
#'       replicates in which each node was assigned to its modal community.}
#'     \item{modularity_samples}{Per-layer vectors of bootstrap modularity
#'       values.}
#'     \item{community_count_reproducibility}{Per-layer numeric vector: the
#'       share of completed replicates whose meta-community count equals
#'       the observed-network count. A descriptive stability measure. The raw
#'       per-replicate community counts are intentionally not returned.}
#'     \item{stability_samples}{List with matrices \code{nmi} and \code{ari}
#'       (completed replicates x layers): agreement of each replicate's
#'       meta-partition with the point-estimate partition. Summarised by
#'       \code{\link{partition_stability}}.}
#'     \item{node_jaccard_stability}{Per-layer vectors: for each node, the
#'       mean Jaccard overlap between its replicate community and its
#'       point-estimate community.}
#'     \item{point_estimate}{The fit result from the original data.}
#'   }
#'
#' @examples
#' set.seed(123)
#' layers <- lapply(1:3, function(i) {
#'   m <- matrix(rbinom(64, 1, 0.35), nrow = 8)
#'   m <- pmax(m, t(m))
#'   diag(m) <- 0
#'   m
#' })
#' boot <- bootstrap_multilayer(
#'   layers,
#'   fit_type = "jaccard",
#'   algorithm = "louvain",
#'   n_boot = 5,
#'   seed = 123
#' )
#' names(boot)
#'
#' @export

bootstrap_multilayer <- function(
    layers,
    fit_type = c(
      "jaccard",
      "overlap",
      "weighted_jaccard",
      "weighted_overlap",
      "identity"
    ),
    algorithm = c("leiden", "louvain"),
    n_boot = 100,
    layer_links = NULL,
    min_similarity = 0,
    resolution_parameter = 1,
    directed = FALSE,
    seed = NULL,
    objective = NULL,
    tracking = "smoothed",
    repair_chains = FALSE,
    repair_min_similarity = 0.2,
    allow_unequal_nodes = FALSE
  ) {

  # check arguments ----
  fit_type <- match.arg(fit_type)
  algorithm <- match.arg(algorithm)
  .check_repair_args(tracking, repair_chains, repair_min_similarity, seed)

  # assign seed unless a seed is provided ----
  if (!is.null(seed)) set.seed(seed)

  # select fitting function ----
  fit_fn <- switch(
    fit_type,
    jaccard = fit_multilayer_jaccard,
    overlap = fit_multilayer_overlap,
    weighted_jaccard = fit_multilayer_weighted_jaccard,
    weighted_overlap = fit_multilayer_weighted_overlap,
    identity = fit_multilayer_identity_ties
  )

  # convert layers to matrices for resampling ----
  mat_layers <- lapply(layers, function(layer) {
    if (inherits(layer, "igraph")) {
      as.matrix(igraph::as_adjacency_matrix(layer, attr = "weight", sparse = FALSE))
    } else {
      as.matrix(layer)
    }
  })

  # Assign layer and node counts ----
  n_layers <- length(mat_layers)
  sizes <- vapply(mat_layers, nrow, integer(1))

  # node identity is the vertex name, never the position (1.4.0): nodes can
  # enter and leave with or without a change in layer size ----
  layer_names <- lapply(mat_layers, rownames)
  named <- !vapply(layer_names, is.null, logical(1))
  same_sets <- all(named) && all(vapply(layer_names, function(n) setequal(n, layer_names[[1]]), logical(1)))
  unequal <- FALSE
  if (isTRUE(allow_unequal_nodes)) {
    if (!all(named)) {
      stop("allow_unequal_nodes = TRUE requires vertex names (igraph names or dimnames) on every layer.",
           call. = FALSE)
    }
    unequal <- !same_sets
  } else if (all(named) && !same_sets) {
    stop("Layers have different node sets; pass allow_unequal_nodes = TRUE (with named vertices) ",
         "to match nodes by name.", call. = FALSE)
  }

  # assign fit arguments ----
  fit_args <- list(
    algorithm = algorithm,
    layer_links = layer_links,
    directed = directed,
    objective = objective
  )

  # assign additional arguments for multilayer identity ties fit ---
  if (fit_type != "identity") {
    fit_args$min_similarity <- min_similarity
    fit_args$resolution_parameter <- resolution_parameter
    fit_args$tracking <- tracking
    fit_args$repair_chains <- repair_chains
    fit_args$repair_min_similarity <- repair_min_similarity
  } else if (!identical(tracking, "smoothed") || repair_chains) {
    stop("`tracking` and `repair_chains` apply to the two-stage fit types, not 'identity'.", call. = FALSE)
  }
  if (unequal) fit_args$allow_unequal_nodes <- TRUE

  # execute multilayer fit for point estimate assignment ----
  point_estimate <- do.call(fit_fn, c(list(layers = mat_layers), fit_args))

  # precompute per-layer edge models for the parametric network bootstrap ----
  # (estimated once from the observed network + point-estimate partition;
  # each replicate then redraws the full edge set from these estimates)
  edge_models <- lapply(seq_len(n_layers), function(layer_idx) {

      ## observed adjacency, detected partition, dyad selector ----
      A <- mat_layers[[layer_idx]]
      n_nodes <- nrow(A)
      mem <- point_estimate$layer_communities[[layer_idx]]$membership
      same <- outer(mem, mem, "==")
      sel <- if (directed) {
        !diag(TRUE, n_nodes)
      } else {
        upper.tri(A)
      }
      E <- A > 0

      ## within/between edge probabilities (fallback: overall density) ----
      in_dyads <- sel & same
      out_dyads <- sel & !same
      p_all <- mean(E[sel])
      p_in <- if (any(in_dyads)) mean(E[in_dyads]) else p_all
      p_out <- if (any(out_dyads)) mean(E[out_dyads]) else p_all

      ## edge-weight pools (resampled with replacement for new edges) ----
      w_in <- A[in_dyads & E]
      w_out <- A[out_dyads & E]
      w_all <- A[sel & E]
      if (length(w_all) == 0) w_all <- 1
      if (length(w_in) == 0) w_in <- w_all
      if (length(w_out) == 0) w_out <- w_all

      list(same = same, sel = sel, p_in = p_in, p_out = p_out,
           w_in = w_in, w_out = w_out)
  })

  # initialize accumulator matrices ----
  co_assign_accum <- lapply(seq_len(n_layers), function(i) {
    matrix(0, nrow = sizes[i], ncol = sizes[i])
  })

  # initialize empty vectors of community assignments ----
  membership_records <- lapply(seq_len(n_layers), function(i) {
    lapply(seq_len(sizes[i]), function(j) integer(0))
  })

  # initialize additional variables ----
  mod_samples <- lapply(seq_len(n_layers), function(i) numeric(0))
  count_samples <- lapply(seq_len(n_layers), function(i) integer(0))
  n_completed <- 0L

  # partition-level agreement of every replicate with the point estimate
  # (rows = completed replicates, columns = layers); node-level Jaccard
  # agreement accumulated per node. These feed partition_stability(). ----
  point_meta <- lapply(point_estimate$meta_communities, as.integer)
  nmi_samples <- matrix(numeric(0), nrow = 0, ncol = n_layers)
  ari_samples <- matrix(numeric(0), nrow = 0, ncol = n_layers)
  node_jaccard_accum <- lapply(seq_len(n_layers), function(i) numeric(sizes[i]))
  node_jaccard_one <- function(mem_a, mem_b) {
    same_a <- outer(mem_a, mem_a, "=="); same_b <- outer(mem_b, mem_b, "==")
    rowSums(same_a & same_b) / rowSums(same_a | same_b)
  }

  # build multilayer fits for all bootstrap replicates ----
  for (b in seq_len(n_boot)) {

    ## parametric network bootstrap: redraw the full edge set per layer ----
    perturbed <- lapply(seq_len(n_layers), function(layer_idx) {
      em <- edge_models[[layer_idx]]
      idx <- which(em$sel)
      probs <- ifelse(em$same[idx], em$p_in, em$p_out)
      on <- idx[stats::rbinom(length(idx), 1, probs) == 1]
      n_nodes <- sizes[layer_idx]
      M <- matrix(0, nrow = n_nodes, ncol = n_nodes)
      if (length(on) > 0) {
        same_on <- em$same[on]
        w <- numeric(length(on))
        if (any(same_on)) {
          w[same_on] <- sample(em$w_in, sum(same_on), replace = TRUE)
        }
        if (any(!same_on)) {
          w[!same_on] <- sample(em$w_out, sum(!same_on), replace = TRUE)
        }
        M[on] <- w
      }
      if (!directed) M <- M + t(M)
      diag(M) <- 0
      if (unequal) dimnames(M) <- dimnames(mat_layers[[layer_idx]])
      return(M)
    })

    ## create multilayer fit ----
    boot_fit <- tryCatch(
      do.call(fit_fn, c(list(layers = perturbed), fit_args)),
      error = function(e) NULL
    )

    ## update completed bootstrap count ----
    if (is.null(boot_fit)) next
    n_completed <- n_completed + 1L

    ## agreement with the point estimate, per layer ----
    boot_meta <- lapply(boot_fit$meta_communities, as.integer)
    nmi_row <- vapply(seq_len(n_layers), function(t)
      igraph::compare(boot_meta[[t]], point_meta[[t]], method = "nmi"), numeric(1))
    ari_row <- vapply(seq_len(n_layers), function(t)
      igraph::compare(boot_meta[[t]], point_meta[[t]], method = "adjusted.rand"), numeric(1))
    nmi_samples <- rbind(nmi_samples, nmi_row)
    ari_samples <- rbind(ari_samples, ari_row)
    for (t in seq_len(n_layers))
      node_jaccard_accum[[t]] <- node_jaccard_accum[[t]] + node_jaccard_one(boot_meta[[t]], point_meta[[t]])

    ## define outputs for each layer ----
    for (layer_idx in seq_len(n_layers)) {
      lc <- boot_fit$layer_communities[[layer_idx]]

      ### meta (cross-layer) membership is the validated partition; the
      ### co-assignment and community counts below are computed on it ----
      mem <- boot_fit$meta_communities[[layer_idx]]
      comms <- split(seq_along(mem), mem)

      ### define co-assignment ----
      for (comm_nodes in comms) {
        if (length(comm_nodes) >= 2) {
          pairs <- utils::combn(comm_nodes, 2)
          for (p in seq_len(ncol(pairs))) {
            ni <- pairs[1, p]
            nj <- pairs[2, p]
            co_assign_accum[[layer_idx]][ni, nj] <- co_assign_accum[[layer_idx]][ni, nj] + 1
            co_assign_accum[[layer_idx]][nj, ni] <- co_assign_accum[[layer_idx]][nj, ni] + 1
          }
        }
      }

      ### define membership records ----
      for (node_id in seq_along(mem)) {
        membership_records[[layer_idx]][[node_id]] <- c(
          membership_records[[layer_idx]][[node_id]],
          mem[node_id]
        )
      }

      ### define modularity ----
      mod_val <- lc$modularity
      mod_samples[[layer_idx]] <- c(
        mod_samples[[layer_idx]],
        if (is.null(mod_val) || is.na(mod_val)) NA_real_ else mod_val
      )

      ### define community count ----
      count_samples[[layer_idx]] <- c(
        count_samples[[layer_idx]],
        length(comms)
      )
    }
  }

  # normalize co-assignment ----
  if (n_completed > 0) {
    co_assignment <- lapply(co_assign_accum, function(m) {
      m <- m / n_completed
      diag(m) <- 1
      return(m)
    })
  } else {
    co_assignment <- co_assign_accum
  }

  # calculate node stabilities ----
  node_stability <- lapply(seq_len(n_layers), function(layer_idx) {
    vapply(seq_len(sizes[layer_idx]), function(node_idx) {
      records <- membership_records[[layer_idx]][[node_idx]]
      if (length(records) == 0) return(0)
      return(max(tabulate(records)) / length(records))
    }, numeric(1))
  })

  # per-layer bootstrap reproducibility of the community count: the share of
  # completed replicates whose community count equals the observed-network
  # (point-estimate) count. Computed here, from the raw per-replicate counts,
  # so those raw counts are never exposed on the returned object. ----
  community_count_reproducibility <- vapply(seq_len(n_layers), function(layer_idx) {
    est <- length(unique(point_estimate$meta_communities[[layer_idx]]))
    s <- count_samples[[layer_idx]]
    if (length(s) == 0) return(NA_real_)
    mean(s == est)
  }, numeric(1))

  # compile multilayer bootstrap ----
  multilayer_bootstrap <- structure(
    list(
      n_boot = n_completed,
      co_assignment = co_assignment,
      node_stability = node_stability,
      modularity_samples = mod_samples,
      community_count_reproducibility = community_count_reproducibility,
      stability_samples = list(nmi = unname(nmi_samples), ari = unname(ari_samples)),
      node_jaccard_stability = lapply(node_jaccard_accum, function(v)
        if (n_completed > 0) v / n_completed else v),
      point_estimate = point_estimate
    ),
    class = "multilayer_bootstrap"
  )

  # return multilayer bootstrap ----
  return(multilayer_bootstrap)
}
