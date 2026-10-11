# 1.4.0: nodes may enter and leave, with or without a change in layer size, in
# every function. Node identity is the vertex name, never the position. Seed 123.

blocks140 <- function(nodes, seed) {
  set.seed(seed)
  n <- length(nodes); A <- matrix(0, n, n)
  for (i in 1:(n - 1)) for (j in (i + 1):n) {
    p <- if ((nodes[i] - 1) %/% 5 == (nodes[j] - 1) %/% 5) 0.8 else 0.03
    if (stats::runif(1) < p) A[i, j] <- A[j, i] <- 1
  }
  g <- igraph::graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE)
  igraph::V(g)$name <- as.character(nodes)
  g
}
same_size140 <- function() list(blocks140(1:20, 1), blocks140(6:25, 2), blocks140(11:30, 3))
diff_size140 <- function() list(blocks140(1:20, 1), blocks140(1:25, 2), blocks140(6:23, 3))
all_fitters140 <- function() c(twostage_fitters(), list(
  identity = dynamicmultiplex::fit_multilayer_identity_ties,
  hungarian = dynamicmultiplex::fit_multilayer_hungarian))

test_that("every fitter rejects changed node sets by default and accepts them by name", {
  for (L in list(same_size140(), diff_size140())) {
    for (nm in names(all_fitters140())) {
      f <- all_fitters140()[[nm]]
      expect_error(suppressWarnings(f(L, seed = 123)), "same node set", info = nm)
      out <- suppressWarnings(f(L, seed = 123, allow_unequal_nodes = TRUE))
      expect_equal(lengths(out$meta_communities), vapply(L, igraph::vcount, numeric(1)), info = nm)
    }
  }
})

test_that("allow_unequal_nodes = TRUE needs names even when layer sizes are equal", {
  u <- lapply(same_size140(), function(g) igraph::delete_vertex_attr(g, "name"))
  expect_error(fit_multilayer_jaccard(u, allow_unequal_nodes = TRUE, seed = 123), "vertex names")
})

test_that("identity ties link a node only to itself", {
  out <- fit_multilayer_identity_ties(same_size140(), seed = 123, allow_unequal_nodes = TRUE)
  ti <- out$interlayer_ties
  got <- sort(paste(ti$from_layer, ti$to_layer, ti$node))
  exp <- sort(c(paste(1, 2, 6:20), paste(2, 3, 11:25)))
  expect_identical(got, exp)
})

test_that("bootstrap handles nodes entering and leaving", {
  for (L in list(same_size140(), diff_size140())) for (ft in c("jaccard", "identity")) {
    expect_error(bootstrap_multilayer(L, fit_type = ft, n_boot = 2, seed = 123), "different node sets")
    res <- suppressWarnings(bootstrap_multilayer(L, fit_type = ft, n_boot = 4, seed = 123, allow_unequal_nodes = TRUE))
    expect_equal(res$n_boot, 4)
    expect_equal(lengths(res$node_stability), vapply(L, igraph::vcount, numeric(1)))
    p <- partition_stability(res)
    expect_true(p$stability >= 0 && p$stability <= 1)
    expect_length(p$stability_by_layer, 3)
  }
})

test_that("bootstrap with equal node sets is unchanged by the flag", {
  L <- list(blocks140(1:20, 1), blocks140(1:20, 2))
  a <- partition_stability(bootstrap_multilayer(L, n_boot = 4, seed = 123))
  b <- partition_stability(bootstrap_multilayer(L, n_boot = 4, seed = 123, allow_unequal_nodes = TRUE))
  expect_identical(a$stability, b$stability)
})
