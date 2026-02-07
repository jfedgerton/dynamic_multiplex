test_that("simulation returns expected structure", {
  sim <- simulate_and_fit_multilayer(
    n_nodes = 30,
    n_layers = 3,
    n_communities = 3,
    fit_type = "jaccard",
    algorithm = "louvain",
    seed = 123
  )

  expect_true(is.list(sim))
  expect_length(sim$layers, 3)
  expect_equal(length(sim$true_membership), 30)
  expect_true(is.list(sim$fit))
})

test_that("directed simulation keeps asymmetric option", {
  sim <- simulate_and_fit_multilayer(
    n_nodes = 20,
    n_layers = 3,
    fit_type = "identity",
    algorithm = "leiden",
    directed = TRUE,
    seed = 1
  )

  expect_true(isTRUE(sim$directed))
  expect_false(isSymmetric(sim$layers[[1]]))
})
