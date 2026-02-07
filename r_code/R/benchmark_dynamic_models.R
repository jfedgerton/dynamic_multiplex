#' Simulate evolving multiplex networks with ground-truth communities
#'
#' Simulates a temporal multiplex system where node communities can persist or
#' switch over time, and generates one network layer per time point.
#'
#' @param n_nodes Number of nodes per layer.
#' @param n_layers Number of temporal layers.
#' @param n_communities Number of communities.
#' @param persistence Probability a node keeps its previous community at the
#'   next layer.
#' @param p_in Tie probability for within-community dyads.
#' @param p_out Tie probability for between-community dyads.
#' @param directed Logical; if `TRUE`, generate directed layers.
#' @param seed Optional random seed.
#'
#' @return List with `layers` (adjacency matrices) and `membership` matrix
#'   (rows = nodes, columns = layers).
#' @export
simulate_evolving_multilayer <- function(n_nodes = 150,
                                         n_layers = 6,
                                         n_communities = 4,
                                         persistence = 0.85,
                                         p_in = 0.20,
                                         p_out = 0.03,
                                         directed = FALSE,
                                         seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  membership <- matrix(NA_integer_, nrow = n_nodes, ncol = n_layers)
  membership[, 1] <- sample(seq_len(n_communities), n_nodes, replace = TRUE)

  for (t in seq(2, n_layers)) {
    keep <- stats::rbinom(n_nodes, 1, persistence) == 1
    membership[, t] <- membership[, t - 1]
    membership[!keep, t] <- sample(seq_len(n_communities), sum(!keep), replace = TRUE)
  }

  layers <- lapply(seq_len(n_layers), function(t) {
    mat <- matrix(0, n_nodes, n_nodes)
    if (directed) {
      for (i in seq_len(n_nodes)) {
        for (j in seq_len(n_nodes)) {
          if (i == j) next
          p <- if (membership[i, t] == membership[j, t]) p_in else p_out
          mat[i, j] <- stats::rbinom(1, 1, p)
        }
      }
    } else {
      for (i in seq_len(n_nodes - 1)) {
        for (j in seq((i + 1), n_nodes)) {
          p <- if (membership[i, t] == membership[j, t]) p_in else p_out
          tie <- stats::rbinom(1, 1, p)
          mat[i, j] <- tie
          mat[j, i] <- tie
        }
      }
    }
    mat
  })

  list(layers = layers, membership = membership)
}

#' @keywords internal
community_pairs_from_membership <- function(membership_vec) {
  split(seq_along(membership_vec), membership_vec)
}

#' @keywords internal
construct_truth_interlayer_ties <- function(membership, layer_links,
                                            metric = c("jaccard", "overlap"),
                                            min_similarity = 0) {
  metric <- match.arg(metric)
  sim_fun <- if (metric == "jaccard") weighted_jaccard else weighted_overlap

  rows <- list()
  for (i in seq_len(nrow(layer_links))) {
    from <- layer_links$from[i]
    to <- layer_links$to[i]
    w <- layer_links$weight[i]

    comm_from <- community_pairs_from_membership(membership[, from])
    comm_to <- community_pairs_from_membership(membership[, to])

    for (a in seq_along(comm_from)) {
      for (b in seq_along(comm_to)) {
        sim <- sim_fun(comm_from[[a]], comm_to[[b]])
        ws <- sim * w
        if (ws >= min_similarity) {
          rows[[length(rows) + 1]] <- data.frame(
            from_layer = from,
            to_layer = to,
            from_community = as.integer(names(comm_from)[a]),
            to_community = as.integer(names(comm_to)[b]),
            weighted_similarity = ws,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  if (length(rows) == 0) {
    return(data.frame(
      from_layer = integer(0), to_layer = integer(0),
      from_community = integer(0), to_community = integer(0),
      weighted_similarity = numeric(0)
    ))
  }

  do.call(rbind, rows)
}

#' @keywords internal
score_tie_recovery <- function(pred, truth, threshold = 0.05) {
  key <- function(df) paste(df$from_layer, df$to_layer, df$from_community, df$to_community, sep = "::")

  pred_pos <- pred[pred$weighted_similarity >= threshold, , drop = FALSE]
  truth_pos <- truth[truth$weighted_similarity >= threshold, , drop = FALSE]

  pred_keys <- unique(key(pred_pos))
  truth_keys <- unique(key(truth_pos))

  tp <- length(intersect(pred_keys, truth_keys))
  fp <- length(setdiff(pred_keys, truth_keys))
  fn <- length(setdiff(truth_keys, pred_keys))

  precision <- if ((tp + fp) == 0) 0 else tp / (tp + fp)
  recall <- if ((tp + fn) == 0) 0 else tp / (tp + fn)
  f1 <- if ((precision + recall) == 0) 0 else 2 * precision * recall / (precision + recall)

  data.frame(tp = tp, fp = fp, fn = fn,
             precision = precision, recall = recall, f1 = f1,
             predicted_ties = length(pred_keys),
             truth_ties = length(truth_keys),
             stringsAsFactors = FALSE)
}

#' @keywords internal
adjusted_rand_index <- function(x, y) {
  tab <- table(x, y)
  n <- sum(tab)
  if (n <= 1) return(NA_real_)

  sum_choose <- function(v) sum(v * (v - 1) / 2)

  a <- sum_choose(tab)
  b <- sum_choose(rowSums(tab))
  c <- sum_choose(colSums(tab))
  d <- n * (n - 1) / 2

  expected <- (b * c) / d
  max_idx <- (b + c) / 2
  denom <- max_idx - expected

  if (denom == 0) return(0)
  (a - expected) / denom
}

#' @keywords internal
mean_layer_ari <- function(fit, truth_membership) {
  aris <- vapply(seq_len(ncol(truth_membership)), function(t) {
    adjusted_rand_index(fit$layer_communities[[t]]$membership, truth_membership[, t])
  }, numeric(1))
  mean(aris, na.rm = TRUE)
}

#' @keywords internal
mean_layer_modularity <- function(fit) {
  mods <- vapply(fit$layer_communities, function(x) x$modularity, numeric(1))
  mean(mods, na.rm = TRUE)
}

#' Benchmark interlayer tie recovery against multinet-style all-to-all baselines
#'
#' Compares the new bayesnet interlayer strategies with all-to-all coupling
#' baselines that mirror common multilayer defaults.
#'
#' @param n_reps Number of simulation replications.
#' @param simulation_args Named list passed to `simulate_evolving_multilayer()`.
#' @param algorithm Community detection algorithm (`"louvain"` or `"leiden"`).
#' @param min_similarity Minimum weighted similarity to keep ties.
#' @param threshold Threshold used to score predicted ties as positive.
#' @param seed Optional random seed.
#'
#' @return List with `raw_results` and aggregated `summary`.
#' @export
benchmark_interlayer_tie_recovery <- function(n_reps = 50,
                                              simulation_args = list(),
                                              algorithm = c("louvain", "leiden"),
                                              min_similarity = 0.05,
                                              threshold = 0.05,
                                              seed = NULL) {
  algorithm <- match.arg(algorithm)
  if (!is.null(seed)) set.seed(seed)

  sim_defaults <- list(
    n_nodes = 150,
    n_layers = 6,
    n_communities = 4,
    persistence = 0.85,
    p_in = 0.20,
    p_out = 0.03,
    directed = FALSE
  )
  sim_par <- utils::modifyList(sim_defaults, simulation_args)

  one_step_links <- make_layer_links(sim_par$n_layers)
  all_to_all_links <- expand.grid(
    from = seq_len(sim_par$n_layers - 1),
    to = seq(2, sim_par$n_layers),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  all_to_all_links <- all_to_all_links[all_to_all_links$from < all_to_all_links$to, , drop = FALSE]
  all_to_all_links$weight <- 1

  out <- vector("list", n_reps * 5)
  idx <- 1

  for (rep_id in seq_len(n_reps)) {
    sim <- do.call(simulate_evolving_multilayer, sim_par)

    truth <- construct_truth_interlayer_ties(
      membership = sim$membership,
      layer_links = one_step_links,
      metric = "jaccard",
      min_similarity = min_similarity
    )

    t_j <- system.time({
      fit_j <- fit_multilayer_jaccard(
        sim$layers,
        algorithm = algorithm,
        layer_links = one_step_links,
        min_similarity = min_similarity,
        directed = sim_par$directed
      )
    })["elapsed"]

    t_o <- system.time({
      fit_o <- fit_multilayer_overlap(
        sim$layers,
        algorithm = algorithm,
        layer_links = one_step_links,
        min_similarity = min_similarity,
        directed = sim_par$directed
      )
    })["elapsed"]

    t_i <- system.time({
      fit_i <- fit_multilayer_identity_ties(
        sim$layers,
        algorithm = algorithm,
        layer_links = one_step_links,
        directed = sim_par$directed
      )
    })["elapsed"]

    t_all <- system.time({
      fit_all <- fit_multilayer_jaccard(
        sim$layers,
        algorithm = algorithm,
        layer_links = all_to_all_links,
        min_similarity = min_similarity,
        directed = sim_par$directed
      )
    })["elapsed"]

    identity_edges <- fit_i$interlayer_ties
    identity_edges$from_community <- mapply(
      function(layer_id, node_id) fit_i$layer_communities[[layer_id]]$membership[node_id],
      identity_edges$from_layer,
      identity_edges$node
    )
    identity_edges$to_community <- mapply(
      function(layer_id, node_id) fit_i$layer_communities[[layer_id]]$membership[node_id],
      identity_edges$to_layer,
      identity_edges$node
    )
    identity_edges$weighted_similarity <- identity_edges$layer_weight
    identity_edges <- unique(identity_edges[, c("from_layer", "to_layer", "from_community", "to_community", "weighted_similarity")])

    model_preds <- list(
      bayesnet_jaccard = fit_j$interlayer_ties[, c("from_layer", "to_layer", "from_community", "to_community", "weighted_similarity")],
      bayesnet_overlap = fit_o$interlayer_ties[, c("from_layer", "to_layer", "from_community", "to_community", "weighted_similarity")],
      bayesnet_identity = identity_edges,
      multinet_all_to_all = fit_all$interlayer_ties[, c("from_layer", "to_layer", "from_community", "to_community", "weighted_similarity")],
      random_baseline = {
        rb <- fit_j$interlayer_ties[, c("from_layer", "to_layer", "from_community", "to_community")]
        rb$weighted_similarity <- stats::runif(nrow(rb))
        rb
      }
    )

    model_meta <- list(
      bayesnet_jaccard = list(elapsed_sec = as.numeric(t_j), mean_layer_ari = mean_layer_ari(fit_j, sim$membership),
                              mean_layer_modularity = mean_layer_modularity(fit_j)),
      bayesnet_overlap = list(elapsed_sec = as.numeric(t_o), mean_layer_ari = mean_layer_ari(fit_o, sim$membership),
                              mean_layer_modularity = mean_layer_modularity(fit_o)),
      bayesnet_identity = list(elapsed_sec = as.numeric(t_i), mean_layer_ari = mean_layer_ari(fit_i, sim$membership),
                               mean_layer_modularity = mean_layer_modularity(fit_i)),
      multinet_all_to_all = list(elapsed_sec = as.numeric(t_all), mean_layer_ari = mean_layer_ari(fit_all, sim$membership),
                                 mean_layer_modularity = mean_layer_modularity(fit_all)),
      random_baseline = list(elapsed_sec = 0, mean_layer_ari = 0, mean_layer_modularity = NA_real_)
    )

    for (m in names(model_preds)) {
      sc <- score_tie_recovery(model_preds[[m]], truth, threshold = threshold)
      sc$elapsed_sec <- model_meta[[m]]$elapsed_sec
      sc$mean_layer_ari <- model_meta[[m]]$mean_layer_ari
      sc$mean_layer_modularity <- model_meta[[m]]$mean_layer_modularity
      sc$f1_per_second <- if (sc$elapsed_sec > 0) sc$f1 / sc$elapsed_sec else NA_real_
      sc$replicate <- rep_id
      sc$model <- m
      out[[idx]] <- sc
      idx <- idx + 1
    }
  }

  raw <- do.call(rbind, out)
  summary <- stats::aggregate(
    cbind(precision, recall, f1, elapsed_sec, f1_per_second,
          mean_layer_ari, mean_layer_modularity, predicted_ties, truth_ties) ~ model,
    data = raw,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  list(raw_results = raw, summary = summary)
}
