# Input policy, node alignment, defaults and plotting fixes (version 1.3.2)

test_that("default algorithm is Leiden", {
  A <- matrix(c(0, 1, 1, 0), 2)
  expect_identical(fit_multilayer_jaccard(list(A, A), seed = 123)$algorithm, "leiden")
  expect_identical(fit_multilayer_hungarian(list(A, A), seed = 123)$algorithm, "leiden")
  expect_identical(fit_multilayer_identity_ties(list(A, A), seed = 123)$algorithm, "leiden")
})

test_that("self-loops and zero-weight edges are dropped for matrix and igraph input", {
  L <- matrix(c(7, 1, 1, 0), 2)
  gl <- igraph::graph_from_adjacency_matrix(L, mode = "undirected", weighted = TRUE, diag = TRUE)
  for (layers in list(list(L, L), list(gl, gl))) {
    g <- dynamicmultiplex:::prepare_multilayer_graphs(layers)[[1]]
    expect_false(any(igraph::which_loop(g)))
    expect_equal(igraph::ecount(g), 1)
  }
  gz <- igraph::make_graph(c(1, 2, 2, 3), directed = FALSE)
  igraph::E(gz)$weight <- c(0, 1)
  g <- dynamicmultiplex:::prepare_multilayer_graphs(list(gz, gz))[[1]]
  expect_equal(igraph::ecount(g), 1)
  expect_equal(igraph::E(g)$weight, 1)
})

test_that("zero-vertex layers are rejected with a clear error", {
  expect_error(fit_multilayer_jaccard(list(matrix(0, 0, 0), matrix(0, 0, 0))), "no vertices")
})

test_that("Hungarian matching follows vertex names, not storage order", {
  g <- igraph::disjoint_union(igraph::make_full_graph(3), igraph::make_full_graph(3))
  igraph::V(g)$name <- letters[1:6]
  g2 <- igraph::permute(g, c(4:6, 1:3))
  h <- fit_multilayer_hungarian(list(g, g2), seed = 123)
  lab1 <- stats::setNames(h$meta_communities[[1]], igraph::V(g)$name)
  lab2 <- stats::setNames(h$meta_communities[[2]], igraph::V(g2)$name)
  expect_identical(lab1[letters[1:6]], lab2[letters[1:6]])
})

test_that("Hungarian matching supports unequal named node sets", {
  g1 <- igraph::disjoint_union(igraph::make_full_graph(4), igraph::make_full_graph(4))
  igraph::V(g1)$name <- letters[1:8]
  g2 <- igraph::delete_vertices(g1, "a")
  g2 <- igraph::add_vertices(g2, 1, name = "z")
  g2 <- igraph::add_edges(g2, c(rbind("z", c("e", "f", "g", "h"))))
  expect_error(fit_multilayer_hungarian(list(g1, g2)), "same node set")
  h <- fit_multilayer_hungarian(list(g1, g2), allow_unequal_nodes = TRUE, seed = 123)
  lab1 <- stats::setNames(h$meta_communities[[1]], igraph::V(g1)$name)
  lab2 <- stats::setNames(h$meta_communities[[2]], igraph::V(g2)$name)
  shared <- intersect(names(lab1), names(lab2))
  expect_identical(lab1[shared], lab2[shared])
  expect_identical(lab2[["z"]], lab2[["e"]])
})

test_that("series plot keeps edges of named graphs and handles > 8 communities", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("RColorBrewer")
  g <- igraph::make_ring(3)
  igraph::V(g)$name <- c("a", "b", "c")
  p <- plot_multilayer_series(list(g, g))
  e <- p$layers[[1]]$data
  expect_equal(nrow(e), 6)
  expect_false(anyNA(e[, c("x", "y", "xend", "yend")]))
  empty9 <- igraph::make_empty_graph(9, directed = FALSE)
  expect_no_error(ggplot2::ggplot_build(plot_multilayer_series(list(empty9, empty9))))
})

test_that("alluvial plot follows actors by name", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggalluvial")
  skip_if_not_installed("RColorBrewer")
  fit <- list(layer_communities = list(
    list(membership = c(alice = 1L, bob = 2L)),
    list(membership = c(bob = 2L, alice = 1L))))
  d <- plot_multilayer_alluvial(fit)$data
  expect_identical(as.integer(as.character(d$community[d$node == "alice"])), c(1L, 1L))
  expect_identical(as.integer(as.character(d$community[d$node == "bob"])), c(2L, 2L))
})
