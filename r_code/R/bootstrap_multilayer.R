#' Bootstrap confidence intervals for multilayer community detection
#'
#' Uses a Bayesian bootstrap on edge weights: each bootstrap replicate
#' multiplies every edge weight by an independent Exponential(1) draw,
#' preserving graph topology while perturbing the weighting that drives
#' community detection.
#'
#' @param layers List of adjacency matrices or igraph objects.
#' @param fit_type One of \code{"jaccard"}, \code{"overlap"},
#'   \code{"weighted_jaccard"}, \code{"weighted_overlap"}, \code{"identity"}.
#' @param algorithm Community detection algorithm: \code{"louvain"} or
#'   \code{"leiden"}.
#' @param n_boot Number of bootstrap replicates.
#' @param layer_links Optional data.frame defining layer connectivity.
#' @param min_similarity Minimum weighted similarity for interlayer ties.
#' @param resolution_parameter Resolution parameter for community detection.
#' @param directed Logical; if \code{TRUE}, treat networks as directed.
#' @param seed Optional random seed for reproducibility.
#'
#' @return A list of class \code{"multilayer_bootstrap"} with components:
#'   \describe{
#'     \item{n_boot}{Number of completed bootstrap replicates.}
#'     \item{co_assignment}{Per-layer co-assignment probability matrices.}
#'     \item{node_stability}{Per-layer vectors giving the fraction of
#'       replicates in which each node was assigned to its modal community.}
#'     \item{modularity_samples}{Per-layer vectors of bootstrap modularity
#'       values.}
#'     \item{community_count_samples}{Per-layer vectors of bootstrap community
#'       counts.}
#'     \item{point_estimate}{The fit result from the original data.}
#'   }
#' @export
bootstrap_multilayer <- function(layers,
                                  fit_type = c("jaccard", "overlap",
                                               "weighted_jaccard",
                                               "weighted_overlap",
                                               "identity"),
                                  algorithm = c("louvain", "leiden"),
                                  n_boot = 100,
                                  layer_links = NULL,
                                  min_similarity = 0,
                                  resolution_parameter = 1,
                                  directed = FALSE,
                                  seed = NULL) {
  fit_type <- match.arg(fit_type)
  algorithm <- match.arg(algorithm)

  if (!is.null(seed)) set.seed(seed)

  # Select fitting function
  fit_fn <- switch(
    fit_type,
    jaccard          = fit_multilayer_jaccard,
    overlap          = fit_multilayer_overlap,
    weighted_jaccard = fit_multilayer_weighted_jaccard,
    weighted_overlap = fit_multilayer_weighted_overlap,
    identity         = fit_multilayer_identity_ties
  )

  # Convert layers to matrices for resampling
  mat_layers <- lapply(layers, function(layer) {
    if (inherits(layer, "igraph")) {
      as.matrix(igraph::as_adjacency_matrix(layer, attr = "weight", sparse = FALSE))
    } else {
      as.matrix(layer)
    }
  })

  n_layers <- length(mat_layers)
  n_nodes <- nrow(mat_layers[[1]])

  # Build common args
  fit_args <- list(
    algorithm = algorithm,
    layer_links = layer_links,
    directed = directed
  )
  if (fit_type != "identity") {
    fit_args$min_similarity <- min_similarity
    fit_args$resolution_parameter <- resolution_parameter
  }

  # Point estimate
  point_estimate <- do.call(fit_fn, c(list(layers = mat_layers), fit_args))

  # Accumulators
  co_assign_accum <- lapply(seq_len(n_layers), function(i) {
    matrix(0, nrow = n_nodes, ncol = n_nodes)
  })
  # membership_records[[layer]][[node]] = vector of community assignments
  membership_records <- lapply(seq_len(n_layers), function(i) {
    lapply(seq_len(n_nodes), function(j) integer(0))
  })
  mod_samples <- lapply(seq_len(n_layers), function(i) numeric(0))
  count_samples <- lapply(seq_len(n_layers), function(i) integer(0))

  n_completed <- 0L

 for (b in seq_len(n_boot)) {
    # Bayesian bootstrap: multiply edge weights by Exp(1) draws
    perturbed <- lapply(mat_layers, function(mat) {
      noise <- matrix(rexp(n_nodes * n_nodes, rate = 1),
                       nrow = n_nodes, ncol = n_nodes)
      if (!directed) {
        noise <- (noise + t(noise)) / 2
      }
      p_mat <- mat * noise
      diag(p_mat) <- 0
      p_mat
    })

    boot_fit <- tryCatch(
      do.call(fit_fn, c(list(layers = perturbed), fit_args)),
      error = function(e) NULL
    )

    if (is.null(boot_fit)) next
    n_completed <- n_completed + 1L

    for (layer_idx in seq_len(n_layers)) {
      lc <- boot_fit$layer_communities[[layer_idx]]
      mem <- lc$membership
      comms <- lc$communities

      # Co-assignment
      for (comm_nodes in comms) {
        if (length(comm_nodes) >= 2) {
          pairs <- utils::combn(comm_nodes, 2)
          for (p in seq_len(ncol(pairs))) {
            ni <- pairs[1, p]
            nj <- pairs[2, p]
            co_assign_accum[[layer_idx]][ni, nj] <-
              co_assign_accum[[layer_idx]][ni, nj] + 1
            co_assign_accum[[layer_idx]][nj, ni] <-
              co_assign_accum[[layer_idx]][nj, ni] + 1
          }
        }
      }

      # Membership records
      for (node_id in seq_along(mem)) {
        membership_records[[layer_idx]][[node_id]] <-
          c(membership_records[[layer_idx]][[node_id]], mem[node_id])
      }

      # Modularity
      mod_val <- lc$modularity
      mod_samples[[layer_idx]] <- c(
        mod_samples[[layer_idx]],
        if (is.null(mod_val) || is.na(mod_val)) NA_real_ else mod_val
      )

      # Community count
      count_samples[[layer_idx]] <- c(
        count_samples[[layer_idx]],
        length(comms)
      )
    }
  }

  # Normalize co-assignment
  if (n_completed > 0) {
    co_assignment <- lapply(co_assign_accum, function(m) {
      m <- m / n_completed
      diag(m) <- 1
      m
    })
  } else {
    co_assignment <- co_assign_accum
  }

  # Node stability
  node_stability <- lapply(seq_len(n_layers), function(layer_idx) {
    vapply(seq_len(n_nodes), function(node_idx) {
      records <- membership_records[[layer_idx]][[node_idx]]
      if (length(records) == 0) return(0)
      max(tabulate(records)) / length(records)
    }, numeric(1))
  })

  structure(
    list(
      n_boot = n_completed,
      co_assignment = co_assignment,
      node_stability = node_stability,
      modularity_samples = mod_samples,
      community_count_samples = count_samples,
      point_estimate = point_estimate
    ),
    class = "multilayer_bootstrap"
  )
}


#' Summarize bootstrap results into confidence intervals
#'
#' @param boot_result Output from \code{\link{bootstrap_multilayer}}.
#' @param alpha Significance level (default 0.05 for 95\% CIs).
#'
#' @return A list with components:
#'   \describe{
#'     \item{modularity_ci}{Data frame with columns layer, estimate, lower,
#'       upper.}
#'     \item{community_count_ci}{Data frame with columns layer, estimate,
#'       lower, upper.}
#'     \item{mean_node_stability}{Data frame with columns layer,
#'       mean_stability.}
#'     \item{node_stability}{Per-layer stability vectors.}
#'     \item{co_assignment}{Per-layer co-assignment matrices.}
#'   }
#' @export
community_ci <- function(boot_result, alpha = 0.05) {
  if (boot_result$n_boot == 0) {
    stop("No completed bootstrap replicates.", call. = FALSE)
  }

  lower_q <- alpha / 2
  upper_q <- 1 - alpha / 2

  point <- boot_result$point_estimate
  n_layers <- length(boot_result$modularity_samples)

  # Modularity CIs
  mod_rows <- lapply(seq_len(n_layers), function(i) {
    lc <- point$layer_communities[[i]]
    est <- lc$modularity
    if (is.null(est) || is.na(est)) est <- NA_real_
    samples <- boot_result$modularity_samples[[i]]
    valid <- samples[!is.na(samples)]
    if (length(valid) > 0) {
      qs <- quantile(valid, probs = c(lower_q, upper_q), names = FALSE)
      lo <- qs[1]
      hi <- qs[2]
    } else {
      lo <- NA_real_
      hi <- NA_real_
    }
    data.frame(layer = i, estimate = est, lower = lo, upper = hi,
               stringsAsFactors = FALSE)
  })

  # Community count CIs
  count_rows <- lapply(seq_len(n_layers), function(i) {
    lc <- point$layer_communities[[i]]
    est <- length(lc$communities)
    samples <- boot_result$community_count_samples[[i]]
    qs <- quantile(samples, probs = c(lower_q, upper_q), names = FALSE)
    data.frame(layer = i, estimate = est, lower = qs[1], upper = qs[2],
               stringsAsFactors = FALSE)
  })

  # Mean node stability
  stab_rows <- lapply(seq_len(n_layers), function(i) {
    data.frame(layer = i,
               mean_stability = mean(boot_result$node_stability[[i]]),
               stringsAsFactors = FALSE)
  })

  list(
    modularity_ci = do.call(rbind, mod_rows),
    community_count_ci = do.call(rbind, count_rows),
    mean_node_stability = do.call(rbind, stab_rows),
    node_stability = boot_result$node_stability,
    co_assignment = boot_result$co_assignment
  )
}
