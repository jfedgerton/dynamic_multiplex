# 1.4.0: the bootstrap is kept only for the stability score s (manuscript
# decision tree). Interval and floor outputs are removed. Seed 123.

test_that("interval and floor functions are removed", {
  ns <- asNamespace("dynamicmultiplex")
  for (nm in c("co_assignment_ci", "community_est", ".load_stability_table", ".floor_lookup")) {
    expect_false(exists(nm, envir = ns, inherits = FALSE), info = nm)
  }
  expect_false("calibration_table" %in% names(formals(partition_stability)))
  expect_identical(system.file("extdata", "stability_calibration_table.csv", package = "dynamicmultiplex"), "")
})

test_that("partition_stability returns s only", {
  L <- lapply(sbm140_layers(), function(g) igraph::as_adjacency_matrix(g, attr = "weight", sparse = FALSE))
  p <- partition_stability(bootstrap_multilayer(L, n_boot = 3, seed = 123))
  expect_setequal(names(p), c("stability", "stability_by_layer", "stability_mc_se", "metric", "report"))
  expect_false(grepl("floor|calibrat", p$report))
})

test_that("stability score is reproducible with seed 123", {
  L <- lapply(sbm140_layers(), function(g) igraph::as_adjacency_matrix(g, attr = "weight", sparse = FALSE))
  for (ft in c("jaccard", "identity")) {
    a <- partition_stability(bootstrap_multilayer(L, fit_type = ft, n_boot = 6, seed = 123))
    b <- partition_stability(bootstrap_multilayer(L, fit_type = ft, n_boot = 6, seed = 123))
    expect_identical(a$stability, b$stability)
    expect_identical(a$stability_by_layer, b$stability_by_layer)
  }
})
