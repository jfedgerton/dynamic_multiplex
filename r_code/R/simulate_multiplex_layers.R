#' Simulate multiplex layers and fit interlayer models
#'
#' Generates synthetic multiplex network layers using a planted partition model,
#' then fits one of the provided interlayer tie strategies.
#'
#' @param n_nodes Number of nodes per layer.
#' @param n_layers Number of temporal layers.
#' @param n_communities Number of latent communities.
#' @param p_in Probability of an in-community edge.
#' @param p_out Probability of an out-community edge.
#' @param fit_type One of `"jaccard"`, `"overlap"`, or `"identity"`.
#' @param algorithm Community algorithm for fitting: `"louvain"` or `"leiden"`.
#' @param layer_links Optional layer connectivity specification.
#' @param min_similarity Minimum similarity threshold for overlap-based methods.
#' @param seed Optional random seed.
#' @param directed Logical; if `TRUE`, simulate directed adjacency matrices.
#'
#' @return A list containing simulated layers, true memberships, and fit results.
#' @export
simulate_and_fit_multilayer <- function(n_nodes = 100,
                                        n_layers = 4,
                                        n_communities = 4,
                                        p_in = 0.2,
                                        p_out = 0.05,
                                        fit_type = c("jaccard", "overlap", "identity"),
                                        algorithm = c("louvain", "leiden"),
                                        layer_links = NULL,
                                        min_similarity = 0,
                                        seed = NULL,
                                        directed = FALSE) {
  fit_type <- match.arg(fit_type)
  algorithm <- match.arg(algorithm)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  memberships <- sample(seq_len(n_communities), size = n_nodes, replace = TRUE)

  layers <- lapply(seq_len(n_layers), function(layer_id) {
    mat <- matrix(0, nrow = n_nodes, ncol = n_nodes)

    if (directed) {
      for (i in seq_len(n_nodes)) {
        for (j in seq_len(n_nodes)) {
          if (i == j) {
            next
          }
          prob <- if (memberships[i] == memberships[j]) p_in else p_out
          mat[i, j] <- stats::rbinom(1, 1, prob)
        }
      }
    } else {
      for (i in seq_len(n_nodes - 1)) {
        for (j in seq((i + 1), n_nodes)) {
          prob <- if (memberships[i] == memberships[j]) p_in else p_out
          tie <- stats::rbinom(1, 1, prob)
          mat[i, j] <- tie
          mat[j, i] <- tie
        }
      }
    }

    mat
  })

  fit <- switch(
    fit_type,
    jaccard = fit_multilayer_jaccard(
      layers,
      algorithm = algorithm,
      layer_links = layer_links,
      min_similarity = min_similarity,
      directed = directed
    ),
    overlap = fit_multilayer_overlap(
      layers,
      algorithm = algorithm,
      layer_links = layer_links,
      min_similarity = min_similarity,
      directed = directed
    ),
    identity = fit_multilayer_identity_ties(
      layers,
      algorithm = algorithm,
      layer_links = layer_links,
      directed = directed
    )
  )

  list(
    layers = layers,
    true_membership = memberships,
    fit = fit,
    directed = directed
  )
}
