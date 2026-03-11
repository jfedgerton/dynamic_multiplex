# Tests for node-universe alignment, and seed reproducibility
#
# Covers:
# 1. Unequal node-universe rejection across all fit methods
# 2. Seed reproducibility for community detection

# -----------------------------------------------------------------------
# 1. Unequal node-universe rejection
# -----------------------------------------------------------------------

test_that("different matrix sizes rejected by jaccard", {
  expect_error(
    fit_multilayer_jaccard(list(diag(4), diag(6)), algorithm = "louvain"),
    "same node"
  )
})

test_that("different matrix sizes rejected by overlap", {
  expect_error(
    fit_multilayer_overlap(list(diag(4), diag(6)), algorithm = "louvain"),
    "same node"
  )
})

test_that("different matrix sizes rejected by weighted_jaccard", {
  expect_error(
    fit_multilayer_weighted_jaccard(list(diag(4), diag(6)), algorithm = "louvain"),
    "same node"
  )
})

test_that("different matrix sizes rejected by weighted_overlap", {
  expect_error(
    fit_multilayer_weighted_overlap(list(diag(4), diag(6)), algorithm = "louvain"),
    "same node"
  )
})

test_that("different matrix sizes rejected by identity", {
  expect_error(
    fit_multilayer_identity_ties(list(diag(4), diag(6)), algorithm = "louvain"),
    "same node"
  )
})

test_that("three layers one different size rejected", {
  expect_error(
    fit_multilayer_jaccard(list(diag(5), diag(5), diag(4)), algorithm = "louvain"),
    "same node"
  )
})

test_that("same size matrices accepted", {
  fit <- fit_multilayer_jaccard(list(diag(5), diag(5)), algorithm = "louvain")
  expect_equal(length(fit$layer_communities), 2)
})


# -----------------------------------------------------------------------
# 2. Seed reproducibility
# -----------------------------------------------------------------------

make_test_layers <- function(seed = 42) {
  set.seed(seed)
  n <- 20
  memberships <- sample(1:3, n, replace = TRUE)
  layers <- lapply(1:3, function(l) {
    mat <- matrix(0, n, n)
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        prob <- if (memberships[i] == memberships[j]) 0.4 else 0.05
        tie <- rbinom(1, 1, prob)
        mat[i, j] <- tie
        mat[j, i] <- tie
      }
    }
    mat
  })
  layers
}

test_that("seed produces identical results for jaccard", {
  layers <- make_test_layers()
  fit1 <- fit_multilayer_jaccard(layers, algorithm = "louvain", seed = 123)
  fit2 <- fit_multilayer_jaccard(layers, algorithm = "louvain", seed = 123)

  for (i in seq_along(fit1$layer_communities)) {
    expect_equal(
      as.integer(fit1$layer_communities[[i]]$membership),
      as.integer(fit2$layer_communities[[i]]$membership)
    )
  }
})

test_that("seed produces identical results for overlap", {
  layers <- make_test_layers()
  fit1 <- fit_multilayer_overlap(layers, algorithm = "louvain", seed = 456)
  fit2 <- fit_multilayer_overlap(layers, algorithm = "louvain", seed = 456)

  for (i in seq_along(fit1$layer_communities)) {
    expect_equal(
      as.integer(fit1$layer_communities[[i]]$membership),
      as.integer(fit2$layer_communities[[i]]$membership)
    )
  }
})

test_that("seed produces identical results for identity", {
  layers <- make_test_layers()
  fit1 <- fit_multilayer_identity_ties(layers, algorithm = "louvain", seed = 789)
  fit2 <- fit_multilayer_identity_ties(layers, algorithm = "louvain", seed = 789)

  for (i in seq_along(fit1$layer_communities)) {
    expect_equal(
      as.integer(fit1$layer_communities[[i]]$membership),
      as.integer(fit2$layer_communities[[i]]$membership)
    )
  }
})

test_that("no seed runs without error", {
  layers <- make_test_layers()
  fit <- fit_multilayer_jaccard(layers, algorithm = "louvain", seed = NULL)
  expect_equal(length(fit$layer_communities), 3)
})

test_that("simulate_and_fit_multilayer with seed is reproducible", {
  out1 <- simulate_and_fit_multilayer(n_nodes = 15, n_layers = 2, seed = 77)
  out2 <- simulate_and_fit_multilayer(n_nodes = 15, n_layers = 2, seed = 77)

  for (i in seq_along(out1$fit$layer_communities)) {
    expect_equal(
      as.integer(out1$fit$layer_communities[[i]]$membership),
      as.integer(out2$fit$layer_communities[[i]]$membership)
    )
  }
})
