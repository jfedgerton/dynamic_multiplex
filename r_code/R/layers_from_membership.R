#' @title Build one weighted layer per time slice from group memberships
#'
#' @description Membership panels (states in IGOs, legislators in caucuses) are
#' usually projected to a node-by-node network whose tie weight counts shared
#' groups. \code{layers_from_membership()} builds that projection, or one of
#' two size-aware alternatives, slice by slice, and returns layers that the
#' \code{fit_multilayer_*} functions accept.
#'
#' Let \eqn{B} be the \eqn{N \times M} 0/1 membership matrix of a slice and
#' \eqn{s_j} the size of group \eqn{j}.
#' \describe{
#'   \item{\code{"count"}}{\eqn{W = B B^\top} (number of shared groups), diagonal 0.}
#'   \item{\code{"size"}}{\eqn{W_{ik} = \sum_j B_{ij} B_{kj} / (s_j - 1)}, so each
#'     group adds total weight proportional to its size rather than its size
#'     squared.}
#'   \item{\code{"kumar"}}{Kumar et al.'s hypergraph reweighting (as in
#'     hypernetx 2.4.3). Group weights \eqn{w_j} start at 1; the slice is
#'     projected as in \code{"size"} with weights \eqn{w_j} and clustered with
#'     \code{algorithm}; each group gets the target
#'     \eqn{r_j = [\sum_c 1/(1 + n_{jc})] (s_j + C) / M} and
#'     \eqn{w_j \leftarrow (w_j + r_j)/2}, until
#'     \eqn{\max_j |w_j - r_j|/2 \le} \code{tol}. The projection with the final
#'     \eqn{w_j} is returned (the reference returns only the partition). The
#'     clustering calls the same per-layer routine as the fitters, so
#'     \code{algorithm} and \code{resolution_parameter} must be the ones you
#'     pass to the fitter.}
#' }
#'
#' Every slice is built on its own and reseeded with \code{seed}: the weights
#' of slice \eqn{t} never depend on any other slice, so no later information
#' enters an earlier layer.
#'
#' No weighting is safe in general. In small Python simulations (hypotheses,
#' not established results): the count projection failed badly when a few
#' very large groups were present, where size-weighted and Kumar both
#' recovered the planted blocs; Kumar beat size-weighting only under noise
#' groups of the same size as real ones; under small groups bridging blocs,
#' count was best and both size-aware weightings failed; and per-layer
#' Louvain reproduced Kumar's own partition less reliably than Leiden.
#'
#' @param memberships List with one entry per slice: a 0/1 matrix with node
#'   labels as row names (columns are groups), or a data frame with exactly
#'   the columns \code{node} and \code{group}.
#' @param weighting \code{"count"}, \code{"size"} or \code{"kumar"}. Required;
#'   there is no default.
#' @param algorithm \code{"leiden"} or \code{"louvain"}; required for
#'   \code{"kumar"} and not accepted otherwise.
#' @param resolution_parameter Resolution for Kumar's clustering (> 0).
#' @param seed Seed for Kumar's clustering, reused for every slice (default
#'   123). The caller's RNG state is restored on exit.
#' @param max_iter Maximum number of Kumar reweighting passes (default 50).
#'   Failing to converge is an error.
#' @param tol Convergence tolerance on the largest half-change in a group
#'   weight.
#' @param min_group_size Groups with fewer members (default 2) are dropped,
#'   with a warning.
#' @param rescale \code{"none"} (default) or \code{"mean"} (each slice divided
#'   by its mean nonzero tie). Matters only where tie scale matters, such as
#'   the multislice coupling \code{omega}; modularity-based two-stage fits are
#'   scale-invariant.
#'
#' @return A list of weighted undirected igraph graphs, one per slice, with
#'   vertex names. Pass \code{allow_unequal_nodes = TRUE} to the fitter when
#'   node sets differ across slices. The graph attribute \code{dm_membership}
#'   holds \code{weighting}, \code{group_weights} (final \eqn{w_j} by group),
#'   \code{iterations}, \code{final_diff}, \code{groups_dropped},
#'   \code{rescale} and \code{partition} (Kumar's final partition, a named
#'   integer vector; \code{NULL} otherwise).
#'
#' @details Projected ties are not independent across pairs (one group
#'   creates many ties), so the bootstrap stability score has not been checked
#'   for layers built this way.
#'
#' @examples
#' m1 <- data.frame(node = c("a", "b", "c", "c", "d"), group = c("g1", "g1", "g1", "g2", "g2"))
#' m2 <- data.frame(node = c("a", "b", "b", "d"), group = c("g1", "g1", "g3", "g3"))
#' layers <- layers_from_membership(list(m1, m2), weighting = "size")
#' igraph::E(layers[[1]])$weight
#'
#' @export
layers_from_membership <- function(memberships, weighting, algorithm = NULL, resolution_parameter = 1,
                                   seed = 123, max_iter = 50L, tol = 0.01, min_group_size = 2L,
                                   rescale = "none") {
  # check arguments ----
  if (!is.list(memberships) || is.data.frame(memberships)) {
    stop("`memberships` must be a list with one entry per time slice.", call. = FALSE)
  }
  if (length(memberships) == 0L) stop("`memberships` must contain at least one slice.", call. = FALSE)
  if (missing(weighting) || !is.character(weighting) || length(weighting) != 1L ||
      !weighting %in% c("count", "size", "kumar")) {
    stop('`weighting` must be "count", "size" or "kumar"; there is no default.', call. = FALSE)
  }
  if (weighting == "kumar") {
    if (is.null(algorithm) || !algorithm %in% c("leiden", "louvain")) {
      stop('weighting = "kumar" needs `algorithm` ("leiden" or "louvain"), the one you pass to the fitter.',
           call. = FALSE)
    }
  } else if (!is.null(algorithm)) {
    stop('`algorithm` is used only by weighting = "kumar"; leave it as NULL.', call. = FALSE)
  }
  if (!is.numeric(resolution_parameter) || length(resolution_parameter) != 1L || !(resolution_parameter > 0)) {
    stop("`resolution_parameter` must be a number > 0.", call. = FALSE)
  }
  if (!is.numeric(max_iter) || length(max_iter) != 1L || max_iter < 1 || max_iter != round(max_iter)) {
    stop("`max_iter` must be an integer >= 1.", call. = FALSE)
  }
  if (!is.numeric(tol) || length(tol) != 1L || !(tol > 0)) stop("`tol` must be a number > 0.", call. = FALSE)
  if (!is.numeric(min_group_size) || length(min_group_size) != 1L || min_group_size < 2 ||
      min_group_size != round(min_group_size)) {
    stop("`min_group_size` must be an integer >= 2.", call. = FALSE)
  }
  if (!is.character(rescale) || length(rescale) != 1L || !rescale %in% c("none", "mean")) {
    stop('`rescale` must be "none" or "mean".', call. = FALSE)
  }

  # scoped RNG: restore the caller's stream on exit ----
  if (weighting == "kumar" && !is.null(seed)) {
    rng_state <- save_rng_state()
    on.exit(restore_rng_state(rng_state), add = TRUE)
  }

  out <- vector("list", length(memberships))
  for (t in seq_along(memberships)) {
    B_all <- .membership_matrix(memberships[[t]], t)
    labels <- rownames(B_all)
    size_all <- colSums(B_all)
    keep <- size_all >= min_group_size
    n_drop <- sum(!keep)
    if (n_drop > 0) {
      warning(sprintf("In slice %d, dropped %d group(s) with fewer than %d members.", t, n_drop, min_group_size),
              call. = FALSE)
    }
    B <- B_all[, keep, drop = FALSE]
    size <- size_all[keep]
    orphan <- rowSums(B) == 0
    if (any(orphan)) {
      stop(sprintf(paste0("In slice %d, %d node(s) belong to no group of at least %d members, for example %s. ",
                          "Drop absent nodes from the slice before building layers."),
                   t, sum(orphan), min_group_size, paste0("'", utils::head(labels[orphan], 10), "'", collapse = ", ")),
           call. = FALSE)
    }

    w <- rep(1, ncol(B))
    iterations <- 0L; final_diff <- NULL; partition <- NULL
    if (weighting == "kumar") {
      comm <- .kumar_cluster(.membership_projection(B, w, size, "size"), algorithm, resolution_parameter, seed)
      repeat {
        r <- .kumar_targets(B, comm)
        final_diff <- max(0.5 * abs(w - r))
        w <- 0.5 * w + 0.5 * r
        iterations <- iterations + 1L
        comm <- .kumar_cluster(.membership_projection(B, w, size, "size"), algorithm, resolution_parameter, seed)
        if (final_diff <= tol) break
        if (iterations >= max_iter) {
          stop(sprintf(paste0("In slice %d, Kumar reweighting did not converge within max_iter = %d passes ",
                              "(largest half-change %.3g > tol = %g)."), t, as.integer(max_iter), final_diff, tol),
               call. = FALSE)
        }
      }
      partition <- stats::setNames(as.integer(comm), labels)
    }

    W <- .membership_projection(B, w, size, if (weighting == "count") "count" else "size")
    up <- which(upper.tri(W) & W != 0, arr.ind = TRUE)
    up <- up[order(up[, 1], up[, 2]), , drop = FALSE]
    vals <- W[up]
    if (!length(vals) || !(sum(vals) > 0)) stop(sprintf("In slice %d, the projection has zero total weight.", t), call. = FALSE)
    if (rescale == "mean") vals <- vals / mean(vals)

    g <- igraph::make_empty_graph(n = length(labels), directed = FALSE)
    igraph::V(g)$name <- labels
    if (nrow(up)) g <- igraph::add_edges(g, as.vector(t(up)), weight = vals)
    g <- igraph::set_graph_attr(g, "dm_membership", list(
      weighting = weighting, group_weights = stats::setNames(w, colnames(B)), iterations = iterations,
      final_diff = final_diff, groups_dropped = as.integer(n_drop), rescale = rescale, partition = partition))
    out[[t]] <- g
  }
  out
}

# validated 0/1 membership matrix (rows nodes, columns groups) for slice t
.membership_matrix <- function(entry, t) {
  if (is.data.frame(entry) && any(c("node", "group") %in% names(entry))) {
    if (!setequal(names(entry), c("node", "group")) || ncol(entry) != 2L) {
      stop(sprintf("In slice %d, a long table must have exactly the columns 'node' and 'group'.", t), call. = FALSE)
    }
    if (nrow(entry) == 0L) stop(sprintf("In slice %d, the membership table is empty.", t), call. = FALSE)
    if (anyNA(entry$node) || anyNA(entry$group)) {
      stop(sprintf("In slice %d, the long table has missing values.", t), call. = FALSE)
    }
    dup <- duplicated(entry[, c("node", "group")])
    if (any(dup)) stop(sprintf("In slice %d, the long table has %d duplicate (node, group) rows.", t, sum(dup)), call. = FALSE)
    nodes <- sort(unique(entry$node)); groups <- sort(unique(entry$group))
    B <- matrix(0, length(nodes), length(groups), dimnames = list(as.character(nodes), as.character(groups)))
    B[cbind(match(entry$node, nodes), match(entry$group, groups))] <- 1
    return(B)
  }
  if (is.data.frame(entry)) entry <- as.matrix(entry)
  if (!is.matrix(entry)) {
    stop(sprintf(paste0("In slice %d, the membership is not a matrix or data frame. Pass a 0/1 matrix with ",
                        "node labels as row names, or a data frame with columns 'node' and 'group'."), t), call. = FALSE)
  }
  if (nrow(entry) == 0L) stop(sprintf("In slice %d, the membership table is empty.", t), call. = FALSE)
  rn <- rownames(entry)
  if (is.null(rn) || anyNA(rn) || any(rn == "")) {
    stop(sprintf("In slice %d, every row needs a node label (row names).", t), call. = FALSE)
  }
  if (anyDuplicated(rn)) stop(sprintf("In slice %d, there are duplicate row labels: %s.", t,
                                      paste(utils::head(unique(rn[duplicated(rn)]), 10), collapse = ", ")), call. = FALSE)
  cn <- colnames(entry)
  if (!is.null(cn) && anyDuplicated(cn)) {
    stop(sprintf("In slice %d, there are duplicate column labels: %s.", t,
                 paste(utils::head(unique(cn[duplicated(cn)]), 10), collapse = ", ")), call. = FALSE)
  }
  if (anyNA(entry)) stop(sprintf("In slice %d, the membership table has missing values.", t), call. = FALSE)
  if (!is.numeric(entry) && !is.logical(entry)) stop(sprintf("In slice %d, membership entries must be 0/1.", t), call. = FALSE)
  if (!all(entry %in% c(0, 1))) stop(sprintf("In slice %d, membership entries must be 0/1.", t), call. = FALSE)
  if (is.null(cn)) colnames(entry) <- paste0("g", seq_len(ncol(entry)))
  storage.mode(entry) <- "double"
  entry
}

# Kumar reweighting targets r_j for a fixed node partition
.kumar_targets <- function(B, labels) {
  inv <- match(labels, unique(labels))
  C <- max(inv)
  H <- matrix(0, nrow(B), C); H[cbind(seq_len(nrow(B)), inv)] <- 1
  n <- crossprod(B, H)                                   # M x C member counts
  rowSums(1 / (1 + n)) * (colSums(B) + C) / ncol(B)
}

.membership_projection <- function(B, w, size, weighting) {
  W <- if (weighting == "count") tcrossprod(B) else B %*% (t(B) * (w / (size - 1)))   # rows of t(B) scaled by group
  diag(W) <- 0
  W
}

# same per-layer routine (and settings) the fitters use, on one weighted slice
.kumar_cluster <- function(W, algorithm, resolution_parameter, seed) {
  up <- which(upper.tri(W) & W != 0, arr.ind = TRUE)
  up <- up[order(up[, 1], up[, 2]), , drop = FALSE]
  g <- igraph::make_empty_graph(n = nrow(W), directed = FALSE)
  if (nrow(up)) g <- igraph::add_edges(g, as.vector(t(up)), weight = W[up])
  if (!is.null(seed)) set.seed(seed)
  fit <- fit_layer_communities(list(g), algorithm = algorithm, resolution_parameter = resolution_parameter,
                               directed = FALSE)[[1]]
  as.integer(fit$membership)
}
