#' @title Plot multiplex networks as wrapped panels
#'
#' @description Draws one panel per layer with nodes colored by community
#' assignment, using a faceted/wrapped layout.
#'
#' @param layers List of `igraph` objects or square adjacency matrices.
#'
#' @param fit Optional fit object from `fit_multilayer_jaccard()`,
#' `fit_multilayer_overlap()`, or `fit_multilayer_identity_ties()`.
#'
#' @param community_memberships Optional list of membership vectors
#' (one per layer). Ignored when `fit` is provided.
#'
#' @param directed Logical; if `TRUE`, adjacency matrices are treated as
#' directed.
#'
#' @param layout `igraph` layout function name. Defaults to `"layout_with_fr"`.
#'
#' @param ncol Number of wrap columns in the panel layout.
#'
#' @param palette ColorBrewer qualitative palette, one of `"Set2"` or `"Dark2"`.
#'
#' @return A `ggplot` object.
#'
#' @importFrom rlang .data
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(123)
#'   layers <- lapply(1:3, function(i) {
#'     m <- matrix(rbinom(64, 1, 0.35), nrow = 8)
#'     m <- pmax(m, t(m))
#'     diag(m) <- 0
#'     m
#'   })
#'   fit <- fit_multilayer_jaccard(layers, algorithm = "louvain")
#'   p <- plot_multilayer_series(layers, fit = fit)
#' }
#'
#' @export

plot_multilayer_series <- function(
    layers,
    fit = NULL,
    community_memberships = NULL,
    directed = FALSE,
    layout = "layout_with_fr",
    ncol = 3,
    palette = "Dark2"
  ) {

  # check for presence of suggested packages ----
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.", call. = FALSE)
  }

  # assign graph layers and memberships ----
  graph_layers <- prepare_multilayer_graphs(layers, directed = directed)
  memberships <- .resolve_memberships(graph_layers, fit, community_memberships)

  # format figure data for each graph layer ----
  plot_rows <- lapply(seq_along(graph_layers), function(i) {
    .build_layer_plot_data(graph_layers[[i]], memberships[[i]], i, layout)
  })

  # bind figure data for nodes and edges into singular data frames ----
  node_df <- do.call(what = rbind, args = lapply(plot_rows, `[[`, "nodes"))
  edge_df <- do.call(what = rbind, args = lapply(plot_rows, `[[`, "edges"))

  # create multilayer series figure ----
  multilayer_series <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edge_df,
      mapping = ggplot2::aes(
        x = .data[["x"]],
        y = .data[["y"]],
        xend = .data[["xend"]],
        yend = .data[["yend"]]),
      alpha = 0.3,
      linewidth = 0.4
    ) +
    ggplot2::geom_point(
      data = node_df,
      mapping = ggplot2::aes(
        x = .data[["x"]],
        y = .data[["y"]],
        color = .data[["community"]]
      ),
      size = 2.7
    ) +
    ggplot2::facet_wrap(facets = ~layer, ncol = ncol) +
    ggplot2::scale_color_manual(values = .brewer_community_palette(palette = palette)) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(color = "Community")

  # return multilayer series figure ----
  return(multilayer_series)
}

#' @title Animate multiplex networks as a GIF by layer
#'
#' @description Creates one frame per layer and colors nodes by community using
#' colorblind-friendly ColorBrewer palettes (`Set2` or `Dark2`).
#'
#' @param layers List of `igraph` objects or square adjacency matrices.
#'
#' @param fit Optional fit object from `fit_multilayer_jaccard()`,
#' `fit_multilayer_overlap()`, or `fit_multilayer_identity_ties()`.
#'
#' @param community_memberships Optional list of membership vectors
#' (one per layer). Ignored when `fit` is provided.
#'
#' @param output_file Output GIF file path. Required; the function writes the
#' animation to this location and nowhere else. Use
#' `tempfile(fileext = ".gif")` to write to the session temporary directory.
#'
#' @param directed Logical; if `TRUE`, adjacency matrices are treated as
#' directed.
#' @param fps Frames per second.
#'
#' @param width Width of GIF in pixels.
#'
#' @param height Height of GIF in pixels.
#'
#' @param layout `igraph` layout function name. Defaults to `"layout_with_fr"`.
#'
#' @param palette ColorBrewer qualitative palette, one of `"Set2"` or `"Dark2"`.
#'
#' @return The path to the created GIF.
#'
#' @importFrom rlang .data
#'
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("gganimate", quietly = TRUE) &&
#'     requireNamespace("gifski", quietly = TRUE)) {
#'   set.seed(123)
#'   layers <- lapply(1:3, function(i) {
#'     m <- matrix(rbinom(64, 1, 0.35), nrow = 8)
#'     m <- pmax(m, t(m))
#'     diag(m) <- 0
#'     m
#'   })
#'   fit <- fit_multilayer_jaccard(layers, algorithm = "louvain")
#'   gif_path <- animate_multilayer_gif(
#'     layers,
#'     fit = fit,
#'     output_file = tempfile(fileext = ".gif"),
#'     fps = 2,
#'     width = 200,
#'     height = 150
#'   )
#' }
#' }
#'
#' @export

animate_multilayer_gif <- function(
    layers,
    fit = NULL,
    community_memberships = NULL,
    output_file,
    directed = FALSE,
    fps = 2,
    width = 800,
    height = 600,
    layout = "layout_with_fr",
    palette = "Dark2"
  ) {

  # check for presence of suggested packages ----
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("gganimate", quietly = TRUE)) {
    stop("Packages 'ggplot2' and 'gganimate' are required.", call. = FALSE)
  }

  # require an explicit output path: the function must never write to the
  # user's filespace by default (CRAN policy) ----
  if (missing(output_file)) {
    stop(
      "`output_file` must be supplied; there is no default write location. ",
      'Use tempfile(fileext = ".gif") to write to the temporary directory.',
      call. = FALSE
    )
  }

  # assign graph layers and memberships ----
  graph_layers <- prepare_multilayer_graphs(layers, directed = directed)
  memberships <- .resolve_memberships(graph_layers, fit, community_memberships)

  # format figure data for each graph layer ----
  plot_rows <- lapply(seq_along(graph_layers), function(i) {
    .build_layer_plot_data(graph_layers[[i]], memberships[[i]], i, layout)
  })

  # bind figure data for nodes and edges into singular data frames ----
  node_df <- do.call(what = rbind, args = lapply(plot_rows, `[[`, "nodes"))
  edge_df <- do.call(what = rbind, args = lapply(plot_rows, `[[`, "edges"))

  # create multilayer animation ----
  multilayer_animation <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edge_df,
      mapping = ggplot2::aes(
        x = .data[["x"]],
        y = .data[["y"]],
        xend = .data[["xend"]],
        yend = .data[["yend"]]
      ),
      alpha = 0.3,
      linewidth = 0.4
    ) +
    ggplot2::geom_point(
      data = node_df,
      mapping = ggplot2::aes(
        x = .data[["x"]],
        y = .data[["y"]],
        color = .data[["community"]]
      ),
      size = 3
    ) +
    ggplot2::scale_color_manual(values = .brewer_community_palette(palette = palette)) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(title = "{closest_state}", color = "Community") +
    gganimate::transition_states(
      states = .data[["layer"]],
      state_length = 1,
      transition_length = 1
    ) +
    gganimate::ease_aes(default = "linear")

  # save multilayer animation ----
  gganimate::anim_save(
    filename = output_file,
    animation = multilayer_animation |>
      gganimate::animate(fps = fps, width = width, height = height)
  )

  # return saved animation file name ----
  return(output_file)
}

#' @title Plot an alluvial view of community transitions over time
#'
#' @param fit Fit object from `fit_multilayer_jaccard()`,
#' `fit_multilayer_overlap()`, or `fit_multilayer_identity_ties()`.
#'
#' @param max_nodes Optional cap on number of nodes to include for readability.
#'
#' @param palette ColorBrewer qualitative palette, one of `"Set2"` or `"Dark2"`.
#'
#' @return A `ggplot` object.
#'
#' @importFrom rlang .data
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("ggalluvial", quietly = TRUE)) {
#'   set.seed(123)
#'   layers <- lapply(1:3, function(i) {
#'     m <- matrix(rbinom(64, 1, 0.35), nrow = 8)
#'     m <- pmax(m, t(m))
#'     diag(m) <- 0
#'     m
#'   })
#'   fit <- fit_multilayer_jaccard(layers, algorithm = "louvain")
#'   p <- plot_multilayer_alluvial(fit)
#' }
#'
#' @export

plot_multilayer_alluvial <- function(fit, max_nodes = NULL, palette = "Dark2") {

  # check for presence of suggested packages ----
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("ggalluvial", quietly = TRUE)) {
    stop("Packages 'ggplot2' and 'ggalluvial' are required.", call. = FALSE)
  }

  # check that layer communities are present ----
  if (is.null(fit$layer_communities)) {
    stop("`fit` must contain `layer_communities`.", call. = FALSE)
  }

  # assign layers and check layer length ----
  layers <- fit$layer_communities
  if (length(layers) < 2) {
    stop("At least two layers are required for alluvial plotting.", call. = FALSE)
  }

  # bind formatted rows ----
  rows <- do.call(
    what = rbind,
    args = lapply(seq_along(layers), function(i) {
      mem <- layers[[i]]$membership
      node_ids <- suppressWarnings(as.integer(names(mem)))
      if (all(is.na(node_ids))) node_ids <- seq_along(mem)
      data.frame(
        node = node_ids,
        layer = paste0("Layer ", i),
        community = as.factor(mem),
        stringsAsFactors = FALSE
      )
    })
  )

  # filter to maximum node count if specified ----
  if (!is.null(max_nodes)) {
    keep <- unique(rows$node)[seq_len(min(max_nodes, length(unique(rows$node))))]
    rows <- rows[rows$node %in% keep, ]
  }

  # create multilayer alluvial figure ----
  multilayer_alluvial <- rows |>
    ggplot2::ggplot(
      mapping = ggplot2::aes(
        x = .data[["layer"]],
        stratum = .data[["community"]],
        alluvium = .data[["node"]],
        y = 1,
        fill = .data[["community"]]
      )
    ) +
    ggalluvial::geom_alluvium(alpha = 0.5) +
    ggalluvial::geom_stratum(alpha = 0.9) +
    ggplot2::scale_fill_manual(values = .brewer_community_palette(palette = palette)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Temporal layer", y = "Nodes", fill = "Community")

  # return multilayer alluvial figure ----
  return(multilayer_alluvial)
}

.resolve_memberships <- function(graph_layers, fit = NULL, community_memberships = NULL) {
  if (!is.null(fit) && !is.null(fit$layer_communities)) {
    return(lapply(fit$layer_communities, function(x) x$membership))
  }

  if (!is.null(community_memberships)) {
    return(community_memberships)
  }

  return(lapply(graph_layers, function(g) seq_len(igraph::vcount(g))))
}

.brewer_community_palette <- function(palette = "Dark2") {
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("Package 'RColorBrewer' is required for palette support.", call. = FALSE)
  }
  palette <- match.arg(palette, c("Set2", "Dark2"))
  return(RColorBrewer::brewer.pal(8, palette))
}

.build_layer_plot_data <- function(g, membership, layer_idx, layout) {
  layout_fn <- if (is.function(layout)) layout else get(layout, envir = asNamespace("igraph"))
  coords <- layout_fn(g)
  nodes <- data.frame(
    node = seq_len(igraph::vcount(g)),
    x = coords[, 1],
    y = coords[, 2],
    community = as.factor(membership),
    layer = paste0("Layer ", layer_idx),
    stringsAsFactors = FALSE
  )

  edge_df <- igraph::as_data_frame(g, what = "edges")
  if (nrow(edge_df) == 0) {
    edge_plot <- data.frame(x = numeric(0), y = numeric(0), xend = numeric(0), yend = numeric(0), layer = character(0))
  } else {
    from_idx <- edge_df$from
    to_idx <- edge_df$to
    edge_plot <- data.frame(
      x = nodes$x[from_idx],
      y = nodes$y[from_idx],
      xend = nodes$x[to_idx],
      yend = nodes$y[to_idx],
      layer = paste0("Layer ", layer_idx)
    )
  }

  return(list(nodes = nodes, edges = edge_plot))
}
