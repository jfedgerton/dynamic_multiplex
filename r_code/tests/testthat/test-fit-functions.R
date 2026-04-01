test_that("default layer links are adjacent", {
  layers <- list(diag(4), diag(4), diag(4))
  fit <- fit_multilayer_jaccard(layers, algorithm = "louvain")

  expect_equal(nrow(fit$layer_links), 2)
  expect_equal(fit$layer_links$from, c(1, 2))
  expect_equal(fit$layer_links$to, c(2, 3))
})

test_that("identity ties contain one record per node per layer pair", {
  layers <- list(diag(5), diag(5), diag(5))
  links <- data.frame(from = c(1, 2), to = c(2, 3), weight = c(1, 0.5))
  fit <- fit_multilayer_identity_ties(layers, algorithm = "louvain", layer_links = links)

  expect_equal(nrow(fit$interlayer_ties), 10)
  expect_true(all(c("from_layer", "to_layer", "node", "layer_weight") %in% names(fit$interlayer_ties)))
})

test_that("identity ties handle variable node sets across layers", {
  g1 <- igraph::make_empty_graph(n = 0, directed = FALSE)
  g1 <- igraph::add_vertices(g1, 3, name = c("A", "B", "C"))
  g1 <- igraph::add_edges(g1, c("A", "B", "A", "C", "B", "C"))

  g2 <- igraph::make_empty_graph(n = 0, directed = FALSE)
  g2 <- igraph::add_vertices(g2, 3, name = c("B", "C", "D"))
  g2 <- igraph::add_edges(g2, c("B", "C", "B", "D", "C", "D"))

  g3 <- igraph::make_empty_graph(n = 0, directed = FALSE)
  g3 <- igraph::add_vertices(g3, 3, name = c("A", "C", "D"))
  g3 <- igraph::add_edges(g3, c("A", "C", "A", "D", "C", "D"))

  fit <- fit_multilayer_identity_ties(list(g1, g2, g3), algorithm = "louvain")

  # Layer 1->2: shared = B, C (2 ties)
  # Layer 2->3: shared = C, D (2 ties)
  # Total: 4 ties
  expect_equal(nrow(fit$interlayer_ties), 4)

  ties_12 <- fit$interlayer_ties[fit$interlayer_ties$from_layer == 1 &
                                   fit$interlayer_ties$to_layer == 2, ]
  expect_true(setequal(ties_12$node, c("B", "C")))

  ties_23 <- fit$interlayer_ties[fit$interlayer_ties$from_layer == 2 &
                                   fit$interlayer_ties$to_layer == 3, ]
  expect_true(setequal(ties_23$node, c("C", "D")))
})

test_that("jaccard undirected self loops use aggregation convention", {
  layers <- list(diag(6), diag(6), diag(6))
  fit <- fit_multilayer_jaccard(layers, algorithm = "louvain")

  loops <- fit$interlayer_ties[
    fit$interlayer_ties$from_layer == fit$interlayer_ties$to_layer &
      fit$interlayer_ties$from_community == fit$interlayer_ties$to_community,
  ]

  expect_gt(nrow(loops), 0)
  # Undirected self-loops use similarity=2 (Louvain aggregation convention)
  expect_true(all(loops$similarity == 2))
  expect_true(all(loops$weighted_similarity == 2 * loops$layer_weight))
})

test_that("overlap self loop multiplier scales weighted similarity", {
  layers <- list(diag(6), diag(6), diag(6))
  fit <- fit_multilayer_overlap(
    layers,
    algorithm = "louvain",
    self_loop_multiplier = 2
  )

  loops <- fit$interlayer_ties[
    fit$interlayer_ties$from_layer == fit$interlayer_ties$to_layer &
      fit$interlayer_ties$from_community == fit$interlayer_ties$to_community,
  ]

  expect_gt(nrow(loops), 0)
  # weighted_similarity = self_sim(2) * layer_weight * multiplier(2)
  expect_true(all(loops$weighted_similarity == 4 * loops$layer_weight))
})


test_that("weighted jaccard uses node strength weighting", {
  layer1 <- matrix(c(
    0, 10, 1,
    10, 0, 0,
    1, 0, 0
  ), nrow = 3, byrow = TRUE)
  layer2 <- matrix(c(
    0, 1, 10,
    1, 0, 0,
    10, 0, 0
  ), nrow = 3, byrow = TRUE)

  fit <- fit_multilayer_weighted_jaccard(
    list(layer1, layer2),
    algorithm = "louvain",
    add_self_loops = FALSE
  )

  expect_equal(fit$interlayer_ties$similarity[1], 13 / 31)
})

test_that("weighted overlap uses node strength weighting", {
  layer1 <- matrix(c(
    0, 10, 1,
    10, 0, 0,
    1, 0, 0
  ), nrow = 3, byrow = TRUE)
  layer2 <- matrix(c(
    0, 1, 10,
    1, 0, 0,
    10, 0, 0
  ), nrow = 3, byrow = TRUE)

  fit <- fit_multilayer_weighted_overlap(
    list(layer1, layer2),
    algorithm = "louvain",
    add_self_loops = FALSE
  )

  expect_equal(fit$interlayer_ties$similarity[1], 13 / 22)
})

test_that("negative weights rejected for louvain", {
  layer <- matrix(c(
    0, 1, -1,
    1, 0, 1,
    -1, 1, 0
  ), nrow = 3, byrow = TRUE)

  expect_error(
    fit_multilayer_jaccard(list(layer, layer), algorithm = "louvain"),
    "negative edge weights"
  )
})

test_that("negative weights rejected for undirected leiden", {
  layer <- matrix(c(
    0, 1, -1,
    1, 0, 1,
    -1, 1, 0
  ), nrow = 3, byrow = TRUE)

  expect_error(
    fit_multilayer_jaccard(list(layer, layer), algorithm = "leiden"),
    "negative edge weights"
  )
})

test_that("negative weights accepted with cpm objective", {
  layer <- matrix(c(
    0, 1, -1,
    1, 0, 1,
    -1, 1, 0
  ), nrow = 3, byrow = TRUE)

  fit <- fit_multilayer_jaccard(list(layer, layer), algorithm = "leiden",
                                 objective = "cpm")
  expect_length(fit$layer_communities, 2)
})

test_that("louvain with cpm objective raises error", {
  layer <- diag(4)

  expect_error(
    fit_multilayer_jaccard(list(layer, layer), algorithm = "louvain",
                            objective = "cpm"),
    "Louvain does not support the CPM objective"
  )
})
