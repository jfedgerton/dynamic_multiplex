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

test_that("jaccard adds self loops with default weighted value", {
  layers <- list(diag(6), diag(6), diag(6))
  fit <- fit_multilayer_jaccard(layers, algorithm = "louvain")

  loops <- fit$interlayer_ties[
    fit$interlayer_ties$from_layer == fit$interlayer_ties$to_layer &
      fit$interlayer_ties$from_community == fit$interlayer_ties$to_community,
  ]

  expect_gt(nrow(loops), 0)
  expect_true(all(loops$similarity == 1))
  expect_true(all(loops$weighted_similarity == loops$layer_weight))
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
  expect_true(all(loops$weighted_similarity == 2 * loops$layer_weight))
})
