#' @title Bootstrap confidence intervals for multilayer community detection
#'
#' @description Quantifies the uncertainty of multilayer community detection
#' by refitting communities on \code{n_boot} resampled versions of the data.
#' Two resampling schemes are available via \code{resample}:
#'
#' \describe{
#'   \item{\code{"edges"} (default)}{Parametric network bootstrap. Within-
#'   and between-community edge probabilities (and, for weighted networks,
#'   edge-weight pools) are estimated from the observed network using the
#'   point-estimate partition, and each replicate redraws the full edge set
#'   from those estimates. This reproduces the variability of fresh data,
#'   including which edges exist.}
#'   \item{\code{"weights"} (legacy)}{Bayesian bootstrap on edge weights:
#'   each replicate multiplies every observed edge weight by an independent
#'   Exponential(1) draw. Topology is held fixed, so replicates see less
#'   variability than fresh data. In simulation studies, confidence
#'   intervals built from this scheme undercovered substantially (pairwise
#'   co-assignment intervals covered ~45-48 percent at a nominal 95
#'   percent); it is retained for backward compatibility and comparison
#'   only.}
#' }
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
#' @param resample Resampling scheme: \code{"edges"} (parametric network
#' bootstrap, default) or \code{"weights"} (legacy Bayesian weight
#' bootstrap). See Description.
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
    objective = NULL,
    resample = c("edges", "weights")
  ) {

  # check arguments ----
  fit_type <- match.arg(fit_type)
  algorithm <- match.arg(algorithm)
  resample <- match.arg(resample)

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
  if (resample == "edges") {
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
  }

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

    ## build one resampled network per layer ----
    if (resample == "edges") {

      ### parametric network bootstrap: redraw the full edge set ----
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
    } else {

      ### legacy bayesian bootstrap: multiply edge weights by exp(1) draws ----
      perturbed <- lapply(mat_layers, function(mat) {
        noise <- matrix(
          data = stats::rexp(n_nodes * n_nodes, rate = 1),
          nrow = n_nodes,
          ncol = n_nodes
        )
        if (!directed) noise <- (noise + t(noise)) / 2
        p_mat <- mat * noise
        diag(p_mat) <- 0
        return(p_mat)
      })
    }

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
      mem <- lc$membership
      comms <- lc$communities

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

  # compile multilayer bootstrap ----
  multilayer_bootstrap <- structure(
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

  # return multilayer bootstrap ----
  return(multilayer_bootstrap)
}


#' @title Summarize bootstrap results into confidence intervals
#'
#' @description Compiles bootstrapping results into community count
#' confidence intervals and reports mean node stability, node stability,
#' and co-assignment metrics.
#'
#' @section Warning:
#' \strong{Community count intervals undercover badly on small networks.}
#' In a large simulation study on planted-partition multilayer networks
#' (n = 50 to 400 nodes; 3 to 10 communities; 5 to 15 layers; varying
#' switching rates and densities), the nominal 95 percent
#' \code{community_count_ci} contained the true community count in only
#' about 40 percent of simulations at n = 50 nodes. Coverage recovers to
#' at or above the nominal level for networks of roughly 100 nodes or
#' more. A runtime warning is issued when layers have fewer than 100
#' nodes; on such networks, treat the intervals as descriptive stability
#' summaries, not calibrated confidence intervals.
#'
#' Earlier versions of this package also returned \code{modularity_ci}.
#' It was removed in version 1.1.0: the same simulation study showed its
#' empirical coverage is never close to the nominal level at any network
#' size (about 0.40 at n = 50, 0.00 at n = 100, and vacuously 1.00 at
#' n = 200), because community detection maximizes modularity and the
#' bootstrap interval concentrates around that optimized, upwardly
#' biased value. Raw bootstrap modularity draws remain available as
#' \code{modularity_samples} in the \code{\link{bootstrap_multilayer}}
#' output for descriptive use.
#'
#' @param boot_result Output from \code{\link{bootstrap_multilayer}}.
#'
#' @param alpha Significance level (default 0.05 for 95 percent CIs).
#'
#' @return A list with components:
#'   \describe{
#'     \item{community_count_ci}{Data frame with columns layer, estimate, lower, upper.}
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
#' ci <- community_ci(boot, alpha = 0.05)
#' str(ci, max.level = 1)
#'
#' @export

community_ci <- function(boot_result, alpha = 0.05) {

  # check that completed bootstrap replicates are present ----
  if (boot_result$n_boot == 0) {
    stop("No completed bootstrap replicates.", call. = FALSE)
  }

  # assign lower and upper quantiles ----
  lower_q <- alpha / 2
  upper_q <- 1 - alpha / 2

  # gather point estimate and layer count ----
  point <- boot_result$point_estimate
  n_layers <- length(boot_result$modularity_samples)

  # warn loudly on small networks: coverage study showed ~40 percent
  # empirical coverage for the nominal 95 percent community count CI at
  # n = 50 nodes, recovering to nominal at n >= 100 ----
  n_nodes <- length(boot_result$node_stability[[1]])
  if (n_nodes < 100) {
    warning(
      "Layers have ", n_nodes, " nodes (< 100). In simulation studies the ",
      "nominal 95% community_count_ci covered the true community count in ",
      "only ~40% of small-network simulations (n = 50). Treat these ",
      "intervals as descriptive stability summaries, not calibrated ",
      "confidence intervals. See ?community_ci section 'Warning'.",
      call. = FALSE
    )
  }

  # calculate community count confidence intervals for each layer ----
  count_rows <- lapply(seq_len(n_layers), function(i) {

    ## calculate confidence interval bounds ----
    lc <- point$layer_communities[[i]]
    est <- length(lc$communities)
    samples <- boot_result$community_count_samples[[i]]
    qs <- stats::quantile(samples, probs = c(lower_q, upper_q), names = FALSE)

    ## compile community count confidence interval ----
    community_count_ci <- data.frame(
      layer = i,
      estimate = est,
      lower = qs[1],
      upper = qs[2],
      stringsAsFactors = FALSE
    )

    ## return community count confidence interval ----
    return(community_count_ci)
  })

  # calculate mean node stability for each layer ----
  stab_rows <- lapply(seq_len(n_layers), function(i) {

    ## calculate mean node stability ----
    mean_node_stability <- data.frame(
      layer = i,
      mean_stability = mean(boot_result$node_stability[[i]]),
      stringsAsFactors = FALSE
    )

    ## return mean node stability ----
    return(mean_node_stability)
  })

  # compile confidence interval results ----
  ci_result <- list(
    community_count_ci = do.call(rbind, count_rows),
    mean_node_stability = do.call(rbind, stab_rows),
    node_stability = boot_result$node_stability,
    co_assignment = boot_result$co_assignment
  )

  # return confidence interval results ----
  return(ci_result)
}
