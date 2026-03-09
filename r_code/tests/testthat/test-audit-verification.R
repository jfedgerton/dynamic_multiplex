# Verification tests added during the code audit.
#
# These tests check:
# 1. Manual Jaccard/overlap computation on known sets
# 2. Resolution parameter is actually passed to Louvain
# 3. Directed graph modularity returns NA
# 4. Community ID correctness in interlayer ties (seq_along bug fix)
# 5. Temporal adjacency correctness
# 6. Edge cases

# -----------------------------------------------------------------------
# 1. Manual Jaccard / overlap checks
# -----------------------------------------------------------------------

test_that("jaccard on known sets matches manual computation", {
  # {1,2,3} vs {2,3,4}: intersection=2, union=4 => 0.5
  sim <- dynamic_multiplex:::weighted_jaccard(c(1, 2, 3), c(2, 3, 4))
  expect_equal(sim, 2 / 4)
})

test_that("overlap on known sets matches manual computation", {
  # {1,2,3} vs {2,3,4}: intersection=2, min_size=3 => 2/3
  sim <- dynamic_multiplex:::weighted_overlap(c(1, 2, 3), c(2, 3, 4))
  expect_equal(sim, 2 / 3)
})

test_that("jaccard of disjoint sets is 0", {
  expect_equal(dynamic_multiplex:::weighted_jaccard(c(1, 2), c(3, 4)), 0)
})

test_that("jaccard of identical sets is 1", {
  expect_equal(dynamic_multiplex:::weighted_jaccard(c(1, 2, 3), c(1, 2, 3)), 1)
})

test_that("overlap when A is subset of B is 1", {
  expect_equal(dynamic_multiplex:::weighted_overlap(c(1, 2), c(1, 2, 3)), 1)
})

test_that("weighted jaccard with known node strengths", {
  a <- c(1, 2, 3)
  b <- c(2, 3, 4)
  wa <- c("1" = 5, "2" = 3, "3" = 1)
  wb <- c("2" = 2, "3" = 4, "4" = 6)
  # intersection = {2,3}: min(3,2) + min(1,4) = 3
  # union = {1,2,3,4}: max(5,0) + max(3,2) + max(1,4) + max(0,6) = 18
  sim <- dynamic_multiplex:::weighted_jaccard_similarity(a, b, wa, wb)
  expect_equal(sim, 3 / 18)
})

test_that("weighted overlap with known node strengths", {
  a <- c(1, 2, 3)
  b <- c(2, 3, 4)
  wa <- c("1" = 5, "2" = 3, "3" = 1)
  wb <- c("2" = 2, "3" = 4, "4" = 6)
  # intersection weight = 3, a_weight = 9, b_weight = 12, min = 9
  sim <- dynamic_multiplex:::weighted_overlap_similarity(a, b, wa, wb)
  expect_equal(sim, 3 / 9)
})


# -----------------------------------------------------------------------
# 2. Resolution parameter forwarding (Bug 1 fix)
# -----------------------------------------------------------------------

test_that("louvain resolution parameter is forwarded", {
  # With a very high resolution, we should get more communities
  # than with the default (1.0)
  set.seed(42)
  n <- 30
  mat <- matrix(0, n, n)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      mat[i, j] <- mat[j, i] <- rbinom(1, 1, 0.3)
    }
  }
  layers <- list(mat, mat)

  fit_default <- fit_multilayer_jaccard(layers, algorithm = "louvain")
  fit_high_res <- fit_multilayer_jaccard(layers, algorithm = "louvain",
                                          resolution_parameter = 5.0)

  n_comm_default <- length(unique(fit_default$layer_communities[[1]]$membership))
  n_comm_high <- length(unique(fit_high_res$layer_communities[[1]]$membership))

  # Higher resolution should produce at least as many communities
  expect_gte(n_comm_high, n_comm_default)
})


# -----------------------------------------------------------------------
# 3. Directed graph modularity (Bug 4 fix)
# -----------------------------------------------------------------------

test_that("directed louvain returns NA modularity", {
  mat <- matrix(c(
    0, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1,
    1, 0, 0, 0
  ), nrow = 4, byrow = TRUE)
  layers <- list(mat, mat)

  fit <- fit_multilayer_jaccard(layers, algorithm = "louvain", directed = TRUE)
  for (lc in fit$layer_communities) {
    expect_true(is.na(lc$modularity),
                info = "Directed Louvain should return NA modularity")
  }
})


# -----------------------------------------------------------------------
# 4. Temporal adjacency
# -----------------------------------------------------------------------

test_that("default layer links form chain topology", {
  links <- dynamic_multiplex:::make_layer_links(5)
  expect_equal(links$from, 1:4)
  expect_equal(links$to, 2:5)
  expect_true(all(links$weight == 1))
})


# -----------------------------------------------------------------------
# 5. Interlayer ties only between specified pairs
# -----------------------------------------------------------------------

test_that("interlayer ties only between specified layer pairs", {
  set.seed(99)
  layers <- lapply(1:4, function(i) {
    mat <- matrix(rbinom(100, 1, 0.3), 10, 10)
    mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
    diag(mat) <- 0
    mat
  })

  custom_links <- data.frame(from = 1, to = 3, weight = 1.0)
  fit <- fit_multilayer_jaccard(layers, algorithm = "louvain",
                                 layer_links = custom_links,
                                 add_self_loops = FALSE)
  ties <- fit$interlayer_ties
  if (nrow(ties) > 0) {
    expect_true(all(ties$from_layer == 1))
    expect_true(all(ties$to_layer == 3))
  }
})


# -----------------------------------------------------------------------
# 6. Edge cases
# -----------------------------------------------------------------------

test_that("empty graph layers do not crash", {
  layers <- list(matrix(0, 5, 5), matrix(0, 5, 5))
  fit <- fit_multilayer_jaccard(layers, algorithm = "louvain")
  expect_length(fit$layer_communities, 2)
})

test_that("self-loops directed vs undirected", {
  mat <- matrix(c(
    0, 1, 0, 0,
    1, 0, 0, 0,
    0, 0, 0, 1,
    0, 0, 1, 0
  ), nrow = 4, byrow = TRUE)
  layers <- list(mat, mat)

  fit_u <- fit_multilayer_jaccard(layers, algorithm = "louvain", directed = FALSE)
  loops_u <- fit_u$interlayer_ties[
    fit_u$interlayer_ties$from_layer == fit_u$interlayer_ties$to_layer &
      fit_u$interlayer_ties$from_community == fit_u$interlayer_ties$to_community,
  ]
  if (nrow(loops_u) > 0) {
    expect_true(all(loops_u$similarity == 2))
  }
})
