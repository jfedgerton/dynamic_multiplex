test_that("benchmark returns quality and speed metrics", {
  b <- benchmark_interlayer_tie_recovery(
    n_reps = 2,
    simulation_args = list(n_nodes = 40, n_layers = 4, n_communities = 3),
    algorithm = "louvain",
    seed = 42
  )

  expect_true(all(c("raw_results", "summary") %in% names(b)))
  expect_true(all(c(
    "precision", "recall", "f1", "elapsed_sec", "f1_per_second",
    "mean_layer_ari", "mean_layer_modularity", "predicted_ties", "truth_ties"
  ) %in% names(b$summary)))
})
