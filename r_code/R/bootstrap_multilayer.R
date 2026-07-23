#' @title Bootstrap confidence intervals for multilayer community detection
#'
#' @description Quantifies the uncertainty of multilayer community detection
#' by refitting communities on \code{n_boot} resampled networks, using a
#' parametric network bootstrap: within- and between-community edge
#' probabilities (and, for weighted networks, edge-weight pools) are
#' estimated from the observed network using the point-estimate partition,
#' and each replicate redraws the full edge set from those estimates. This
#' reproduces the variability of fresh data, including which edges exist.
#'
#' Uncertainty is quantified on the cross-layer \emph{meta-communities}
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
#' @param algorithm Community detection algorithm: \code{"louvain"} or
#' \code{"leiden"}.
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
#' @param objective One of "cpm" or "modularity" for directed networks only
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
    algorithm = c("louvain", "leiden"),
    n_boot = 100,
    layer_links = NULL,
    min_similarity = 0,
    resolution_parameter = 1,
    directed = FALSE,
    seed = NULL,
    objective = NULL
  ) {

  # check arguments ----
  fit_type <- match.arg(fit_type)
  algorithm <- match.arg(algorithm)

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
  n_nodes <- nrow(mat_layers[[1]])

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
  }

  # execute multilayer fit for point estimate assignment ----
  point_estimate <- do.call(fit_fn, c(list(layers = mat_layers), fit_args))

  # precompute per-layer edge models for the parametric network bootstrap ----
  # (estimated once from the observed network + point-estimate partition;
  # each replicate then redraws the full edge set from these estimates)
  edge_models <- lapply(seq_len(n_layers), function(layer_idx) {

      ## observed adjacency, detected partition, dyad selector ----
      A <- mat_layers[[layer_idx]]
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
    matrix(0, nrow = n_nodes, ncol = n_nodes)
  })

  # initialize empty vectors of community assignments ----
  membership_records <- lapply(seq_len(n_layers), function(i) {
    lapply(seq_len(n_nodes), function(j) integer(0))
  })

  # initialize additional variables ----
  mod_samples <- lapply(seq_len(n_layers), function(i) numeric(0))
  count_samples <- lapply(seq_len(n_layers), function(i) integer(0))
  n_completed <- 0L

  # build multilayer fits for all bootstrap replicates ----
  for (b in seq_len(n_boot)) {

    ## parametric network bootstrap: redraw the full edge set per layer ----
    perturbed <- lapply(seq_len(n_layers), function(layer_idx) {
      em <- edge_models[[layer_idx]]
      idx <- which(em$sel)
      probs <- ifelse(em$same[idx], em$p_in, em$p_out)
      on <- idx[stats::rbinom(length(idx), 1, probs) == 1]
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
    vapply(seq_len(n_nodes), function(node_idx) {
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
      point_estimate = point_estimate
    ),
    class = "multilayer_bootstrap"
  )

  # return multilayer bootstrap ----
  return(multilayer_bootstrap)
}


#' @title Summarize bootstrap community-count reproducibility (meta-communities)
#'
#' @description Reports, for each layer, the community count from the
#' observed network together with its \emph{bootstrap reproducibility}: the
#' proportion of bootstrap replicates in which the fitted number of
#' communities equals the observed-network count. It also returns mean node
#' stability, per-node stability, and the co-assignment matrices.
#'
#' @section Why this is not a confidence interval:
#' Earlier versions returned a percentile \code{community_count_ci}. It was
#' replaced in version 1.1.0 with a reproducibility summary because the
#' interval's coverage is not robust to model misspecification. In a large
#' simulation study the nominal 95 percent community-count interval covered
#' the truth at or above the nominal level on well-specified
#' planted-partition networks (about 0.99 for n >= 100 nodes), but coverage
#' collapsed to about 0.62 when community sizes were strongly skewed, and no
#' observable diagnostic reliably separated the trustworthy cases from the
#' rest. Rather than ship an interval that silently undercovers, the
#' function now reports how often the community count reproduces under
#' resampling. This is a descriptive stability measure, not a calibrated
#' interval: it makes no claim about the probability that any range contains
#' the true community count. For a validated interval, use
#' \code{\link{co_assignment_ci}}, whose node-pair coverage held across the
#' same misspecification stress tests. The raw per-replicate community counts
#' are intentionally not exposed anywhere in the package output; only this
#' reproducibility summary is returned.
#'
#' @param boot_result Output from \code{\link{bootstrap_multilayer}}.
#'
#' @return A list with components:
#'   \describe{
#'     \item{community_count}{Data frame with columns layer, estimate (the
#'       observed-network meta-community count), and reproducibility (share of
#'       bootstrap replicates whose community count equals estimate, in
#'       [0, 1]).}
#'     \item{report}{Character vector, one plain-language sentence per layer.}
#'     \item{mean_node_stability}{Data frame with columns layer, mean_stability.}
#'     \item{node_stability}{Per-layer stability vectors.}
#'     \item{co_assignment}{Per-layer co-assignment matrices.}
#'   }
#'
#' @seealso \code{\link{co_assignment_ci}} for calibrated node-pair
#' co-assignment intervals.
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
#' est <- community_est(boot)
#' est$community_count
#' est$report
#'
#' @export

community_est <- function(boot_result) {

  # check that completed bootstrap replicates are present ----
  if (boot_result$n_boot == 0) {
    stop("No completed bootstrap replicates.", call. = FALSE)
  }

  # gather point estimate and layer count ----
  point <- boot_result$point_estimate
  n_layers <- length(boot_result$modularity_samples)

  # per-layer community count and bootstrap reproducibility ----
  count_rows <- lapply(seq_len(n_layers), function(i) {

    ## observed-network community count for this layer (the anchor) ----
    est <- length(unique(point$meta_communities[[i]]))

    ## reproducibility = share of bootstrap replicates whose community count
    ## equals the observed-network count. Descriptive stability, not a
    ## coverage claim. Precomputed in bootstrap_multilayer() from the raw
    ## per-replicate counts, which are deliberately not exposed. ----
    reproducibility <- boot_result$community_count_reproducibility[[i]]

    ## compile row ----
    data.frame(
      layer = i,
      estimate = est,
      reproducibility = reproducibility,
      stringsAsFactors = FALSE
    )
  })
  community_count <- do.call(rbind, count_rows)

  # one plain-language sentence per layer ----
  report <- sprintf(
    "Layer %d: community count (K = %d) reproduced in %d%% of bootstrap resamples.",
    community_count$layer,
    community_count$estimate,
    round(100 * community_count$reproducibility)
  )

  # calculate mean node stability for each layer ----
  stab_rows <- lapply(seq_len(n_layers), function(i) {
    data.frame(
      layer = i,
      mean_stability = mean(boot_result$node_stability[[i]]),
      stringsAsFactors = FALSE
    )
  })

  # compile results ----
  est_result <- list(
    community_count = community_count,
    report = report,
    mean_node_stability = do.call(rbind, stab_rows),
    node_stability = boot_result$node_stability,
    co_assignment = boot_result$co_assignment
  )

  # return results ----
  return(est_result)
}
