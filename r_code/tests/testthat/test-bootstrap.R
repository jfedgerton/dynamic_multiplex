make_planted_layers <- function(n_nodes = 30, n_layers = 3, seed = 42) {
  set.seed(seed)
  memberships <- sample(1:3, n_nodes, replace = TRUE)
  lapply(seq_len(n_layers), function(t) {
    mat <- matrix(0, n_nodes, n_nodes)
    for (i in seq_len(n_nodes - 1)) {
      for (j in (i + 1):n_nodes) {
        prob <- if (memberships[i] == memberships[j]) 0.4 else 0.05
        tie <- rbinom(1, 1, prob)
        mat[i, j] <- tie
        mat[j, i] <- tie
      }
    }
    mat
  })
}

test_that("bootstrap_multilayer smoke test", {
  layers <- make_planted_layers()
  result <- bootstrap_multilayer(layers, fit_type = "jaccard", n_boot = 5,
                                  seed = 1)

  expect_s3_class(result, "multilayer_bootstrap")
  expect_equal(result$n_boot, 5)
  expect_length(result$co_assignment, 3)
  expect_length(result$node_stability, 3)
  expect_length(result$modularity_samples, 3)
  expect_length(result$community_count_samples, 3)
  expect_true(!is.null(result$point_estimate))
})

test_that("co-assignment matrices have correct shape and range", {
  layers <- make_planted_layers(n_nodes = 20, n_layers = 2)
  result <- bootstrap_multilayer(layers, fit_type = "jaccard", n_boot = 10,
                                  seed = 2)

  for (co in result$co_assignment) {
    expect_equal(dim(co), c(20, 20))
    expect_true(all(co >= 0))
    expect_true(all(co <= 1))
    # Diagonal should be 1
    expect_equal(diag(co), rep(1, 20))
    # Symmetric
    expect_equal(co, t(co))
  }
})

test_that("node stability values are in [0, 1]", {
  layers <- make_planted_layers(n_nodes = 20, n_layers = 2)
  result <- bootstrap_multilayer(layers, fit_type = "overlap", n_boot = 10,
                                  seed = 3)

  for (stab in result$node_stability) {
    expect_length(stab, 20)
    expect_true(all(stab >= 0))
    expect_true(all(stab <= 1))
  }
})

test_that("modularity samples have correct length", {
  layers <- make_planted_layers(n_nodes = 20, n_layers = 2)
  result <- bootstrap_multilayer(layers, fit_type = "jaccard", n_boot = 8,
                                  seed = 4)

  for (mod_s in result$modularity_samples) {
    expect_length(mod_s, result$n_boot)
  }
})

test_that("community count samples are positive", {
  layers <- make_planted_layers(n_nodes = 20, n_layers = 2)
  result <- bootstrap_multilayer(layers, fit_type = "jaccard", n_boot = 5,
                                  seed = 5)

  for (count_s in result$community_count_samples) {
    expect_length(count_s, result$n_boot)
    expect_true(all(count_s >= 1))
  }
})

test_that("all fit types work with bootstrap", {
  layers <- make_planted_layers(n_nodes = 20, n_layers = 2)
  for (ft in c("jaccard", "overlap", "weighted_jaccard",
               "weighted_overlap", "identity")) {
    result <- bootstrap_multilayer(layers, fit_type = ft, n_boot = 3, seed = 6)
    expect_equal(result$n_boot, 3)
  }
})

test_that("custom layer links are respected", {
  layers <- make_planted_layers(n_nodes = 15, n_layers = 4)
  links <- data.frame(from = 1, to = 3, weight = 0.5)
  result <- bootstrap_multilayer(layers, fit_type = "jaccard", n_boot = 3,
                                  seed = 7, layer_links = links)
  expect_equal(result$n_boot, 3)
  expect_equal(nrow(result$point_estimate$layer_links), 1)
})

test_that("community_ci returns correct structure", {
  layers <- make_planted_layers(n_nodes = 20, n_layers = 3)
  boot <- bootstrap_multilayer(layers, fit_type = "jaccard", n_boot = 10,
                                seed = 10)
  expect_warning(ci <- community_ci(boot), "nodes")

  expect_false("modularity_ci" %in% names(ci))
  expect_true("community_count_ci" %in% names(ci))
  expect_true("mean_node_stability" %in% names(ci))
  expect_true("node_stability" %in% names(ci))
  expect_true("co_assignment" %in% names(ci))

  expect_equal(nrow(ci$community_count_ci), 3)
  expect_equal(nrow(ci$mean_node_stability), 3)
})

test_that("community_ci columns are correct", {
  layers <- make_planted_layers(n_nodes = 20, n_layers = 2)
  boot <- bootstrap_multilayer(layers, fit_type = "jaccard", n_boot = 10,
                                seed = 11)
  ci <- suppressWarnings(community_ci(boot))

  expected_cols <- c("layer", "estimate", "lower", "upper")
  expect_true(all(expected_cols %in% names(ci$community_count_ci)))
})

test_that("CI lower <= upper", {
  layers <- make_planted_layers(n_nodes = 20, n_layers = 2)
  boot <- bootstrap_multilayer(layers, fit_type = "jaccard", n_boot = 20,
                                seed = 12)
  ci <- suppressWarnings(community_ci(boot))

  count <- ci$community_count_ci
  expect_true(all(count$lower <= count$upper))
})

test_that("mean stability is in [0, 1]", {
  layers <- make_planted_layers(n_nodes = 20, n_layers = 2)
  boot <- bootstrap_multilayer(layers, fit_type = "jaccard", n_boot = 10,
                                seed = 14)
  ci <- suppressWarnings(community_ci(boot))

  expect_true(all(ci$mean_node_stability$mean_stability >= 0))
  expect_true(all(ci$mean_node_stability$mean_stability <= 1))
})


test_that("co_assignment_ci returns Wilson bounds with correct structure", {
  layers <- make_planted_layers(n_nodes = 20, n_layers = 2)
  boot <- bootstrap_multilayer(layers, fit_type = "jaccard", n_boot = 10,
                                seed = 21)
  pci <- co_assignment_ci(boot)

  expect_equal(length(pci), 2)
  for (layer_ci in pci) {
    expect_true(all(c("estimate", "lower", "upper") %in% names(layer_ci)))
    expect_equal(dim(layer_ci$estimate), c(20, 20))
    off <- !diag(20)
    expect_true(all(layer_ci$lower[off] <= layer_ci$estimate[off] + 1e-12))
    expect_true(all(layer_ci$estimate[off] <= layer_ci$upper[off] + 1e-12))
    expect_true(all(layer_ci$lower >= 0) && all(layer_ci$upper <= 1))
    expect_true(all(diag(layer_ci$lower) == 1))
    expect_true(all(diag(layer_ci$upper) == 1))
  }
})

test_that("co_assignment_ci intervals narrow as alpha grows", {
  layers <- make_planted_layers(n_nodes = 20, n_layers = 2)
  boot <- bootstrap_multilayer(layers, fit_type = "jaccard", n_boot = 20,
                                seed = 22)
  pci_95 <- co_assignment_ci(boot, alpha = 0.05)
  pci_50 <- co_assignment_ci(boot, alpha = 0.50)
  off <- !diag(20)
  w95 <- (pci_95[[1]]$upper - pci_95[[1]]$lower)[off]
  w50 <- (pci_50[[1]]$upper - pci_50[[1]]$lower)[off]
  expect_true(all(w50 <= w95 + 1e-12))
})

test_that("community_ci warns on small networks", {
  layers <- make_planted_layers(n_nodes = 20, n_layers = 2)
  boot <- bootstrap_multilayer(layers, fit_type = "jaccard", n_boot = 5,
                                seed = 23)
  expect_warning(community_ci(boot), "descriptive stability summaries")
})

test_that("edge-resampling bootstrap runs and differs from weights scheme", {
  layers <- make_planted_layers(n_nodes = 20, n_layers = 2)
  boot_e <- bootstrap_multilayer(layers, fit_type = "jaccard", n_boot = 8,
                                  seed = 31, resample = "edges")
  boot_w <- bootstrap_multilayer(layers, fit_type = "jaccard", n_boot = 8,
                                  seed = 31, resample = "weights")
  expect_equal(boot_e$n_boot, 8)
  expect_equal(boot_w$n_boot, 8)
  expect_equal(dim(boot_e$co_assignment[[1]]), c(20, 20))
  # both give probabilities in [0, 1]
  expect_true(all(boot_e$co_assignment[[1]] >= 0 & boot_e$co_assignment[[1]] <= 1))
  # edges scheme should generally produce different co-assignment than weights
  expect_false(identical(boot_e$co_assignment[[1]], boot_w$co_assignment[[1]]))
})

test_that("edges is the default resampling scheme", {
  layers <- make_planted_layers(n_nodes = 15, n_layers = 2)
  boot_default <- bootstrap_multilayer(layers, fit_type = "jaccard",
                                        n_boot = 5, seed = 32)
  boot_edges <- bootstrap_multilayer(layers, fit_type = "jaccard",
                                      n_boot = 5, seed = 32,
                                      resample = "edges")
  expect_equal(boot_default$co_assignment[[1]], boot_edges$co_assignment[[1]])
})
