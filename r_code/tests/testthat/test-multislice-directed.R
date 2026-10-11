# Directed within-slice edges in the multislice (Mucha) fit.
#
# Mirrors python_code/tests/test_multislice_directed.py. Seed 123; the
# known-answer networks are hard coded.
#
# * Reciprocal equivalence: a directed network with every arc reciprocated
#   must give the same fit as its undirected twin at the same omega (before
#   1.3.2 its arcs were stacked as parallel undirected edges, which is the
#   undirected fit at omega / 2).
# * Directed null model: with directed = TRUE each slice uses the
#   k_out * k_in / m null model; a 2-slice, 4-node network is solved by brute
#   force over all 4,140 partitions of the 8 supra-nodes.
# * Regression: undirected inputs, and directed layers fitted with
#   directed = FALSE, reproduce the pre-change memberships exactly.

seed_dm <- 123

block_edges <- function(assign) {
  n <- length(assign)
  e <- NULL
  for (i in 1:(n - 1)) for (j in (i + 1):n) if (assign[i] == assign[j]) e <- rbind(e, c(i, j))
  rbind(e, c(1, 12))
}

planted12 <- function(directed) {
  base <- c(rep(1, 6), rep(2, 6))
  s3 <- base; s3[6] <- 2; s3[7] <- 1
  lapply(list(base, base, s3), function(a) {
    e <- block_edges(a)
    if (directed) e <- rbind(e, e[, 2:1])
    g <- igraph::graph_from_edgelist(e, directed = directed)
    igraph::V(g)$name <- as.character(1:12)
    igraph::E(g)$weight <- 1
    g
  })
}

# 2 slices x 4 nodes, one-way arcs (1-indexed); unique optimum at omega = 0.25
# splits {1, 4} from {2, 3} in both slices with Q = 46 / 15.
directed_arcs <- list(
  rbind(c(2, 1), c(2, 3), c(2, 4), c(3, 1), c(3, 2), c(4, 1)),
  rbind(c(1, 2), c(1, 3), c(1, 4), c(3, 2), c(4, 3))
)
directed_omega <- 0.25
directed_best_q <- 46 / 15
directed_best_partition <- c(1, 2, 2, 1, 1, 2, 2, 1)

directed4 <- function() {
  lapply(directed_arcs, function(e) {
    g <- igraph::make_empty_graph(n = 4, directed = TRUE)
    g <- igraph::add_edges(g, as.vector(t(e)))
    igraph::V(g)$name <- as.character(1:4)
    igraph::E(g)$weight <- 1
    g
  })
}

sbm_layers <- function(drop_node_from_layer2 = FALSE) {
  ed <- utils::read.csv(testthat::test_path("fixtures", "multislice_sbm_edges.csv"))
  out <- lapply(1:4, function(t) {
    e <- ed[ed$layer == t, ]
    g <- igraph::make_empty_graph(n = 30, directed = FALSE)
    igraph::V(g)$name <- as.character(1:30)
    g <- igraph::add_edges(g, as.vector(t(as.matrix(e[, c("from", "to")]))), weight = e$weight)
    g
  })
  if (drop_node_from_layer2) out[[2]] <- igraph::delete_vertices(out[[2]], "30")
  out
}

canon <- function(meta) {
  f <- unlist(meta)
  as.integer(factor(f, levels = unique(f)))
}

all_partitions <- function(n) {
  out <- list(); a <- integer(n)
  rec <- function(i, mx) {
    if (i > n) { out[[length(out) + 1L]] <<- a; return(invisible()) }
    for (v in 0:(mx + 1L)) { a[i] <<- v; rec(i + 1L, max(mx, v)) }
  }
  a[1] <- 0L
  rec(2L, 0L)
  out
}

directed_multislice_q <- function(memb, arcs_by_slice, n, omega, gamma = 1) {
  q <- 0
  for (s in seq_along(arcs_by_slice)) {
    a <- matrix(0, n, n)
    a[arcs_by_slice[[s]]] <- 1
    k_out <- rowSums(a); k_in <- colSums(a); m <- sum(a)
    idx <- (s - 1) * n + seq_len(n)
    same <- outer(memb[idx], memb[idx], "==")
    q <- q + sum((a - gamma * outer(k_out, k_in) / m)[same])
  }
  q + 2 * omega * sum(memb[1:n] == memb[n + 1:n])
}

regression_cases <- function() {
  sbm <- sbm_layers()
  pairs <- expand.grid(from = 1:4, to = 1:4)
  pairs <- pairs[pairs$from < pairs$to, ]
  pairs$weight <- 1 / (pairs$to - pairs$from)
  list(
    planted12_omega1 = list(layers = planted12(FALSE), omega = 1),
    planted12_omega6 = list(layers = planted12(FALSE), omega = 6),
    sbm_omega1 = list(layers = sbm, omega = 1),
    sbm_omega05_gamma12 = list(layers = sbm, omega = 0.5, resolution_parameter = 1.2),
    sbm_all_pairs = list(layers = sbm, omega = 1, layer_links = pairs),
    sbm_unequal_nodes = list(layers = sbm_layers(TRUE), omega = 1, allow_unequal_nodes = TRUE)
  )
}

run_case <- function(case) {
  fit <- suppressWarnings(do.call(fit_multilayer_identity_ties, c(case, list(seed = seed_dm))))
  lapply(fit$meta_communities, as.integer)
}

test_that("reciprocal directed layers give the undirected fit at the same omega", {
  for (om in c(0.5, 1, 3, 5, 6, 8, 10)) {
    und <- suppressWarnings(fit_multilayer_identity_ties(planted12(FALSE), omega = om, seed = seed_dm))
    dig <- suppressWarnings(fit_multilayer_identity_ties(planted12(TRUE), omega = om,
                                                         directed = TRUE, seed = seed_dm))
    expect_identical(canon(dig$meta_communities), canon(und$meta_communities), info = paste("omega", om))
  }
})

test_that("reciprocal directed fit no longer behaves like omega / 2", {
  und6 <- canon(fit_multilayer_identity_ties(planted12(FALSE), omega = 6, seed = seed_dm)$meta_communities)
  und3 <- canon(fit_multilayer_identity_ties(planted12(FALSE), omega = 3, seed = seed_dm)$meta_communities)
  expect_false(identical(und6, und3))
  dig6 <- suppressWarnings(fit_multilayer_identity_ties(planted12(TRUE), omega = 6,
                                                        directed = TRUE, seed = seed_dm))
  expect_identical(canon(dig6$meta_communities), und6)
})

test_that("directed null model matches the brute-force optimum", {
  parts <- all_partitions(8)
  expect_length(parts, 4140)
  qs <- vapply(parts, directed_multislice_q, numeric(1), arcs_by_slice = directed_arcs,
               n = 4, omega = directed_omega)
  expect_equal(max(qs), directed_best_q)
  expect_equal(sum(abs(qs - max(qs)) < 1e-9), 1)
  fit <- suppressWarnings(fit_multilayer_identity_ties(directed4(), omega = directed_omega,
                                                       directed = TRUE, seed = seed_dm))
  got <- canon(fit$meta_communities)
  expect_identical(got, as.integer(directed_best_partition))
  expect_equal(directed_multislice_q(got, directed_arcs, 4, directed_omega), directed_best_q)
})

test_that("detect_multislice_communities symmetrises directed layers when directed = FALSE", {
  # fit_multilayer_identity_ties() stops earlier (per-layer igraph Leiden is
  # undirected only), so call the multislice stage directly
  ties <- data.frame(from_layer = 1L, to_layer = 2L, node = as.character(1:4),
                     layer_weight = 1, stringsAsFactors = FALSE)
  set.seed(seed_dm)
  expect_warning(
    meta <- detect_multislice_communities(directed4(), ties, omega = directed_omega, directed = FALSE),
    "symmetrised"
  )
  expect_identical(canon(meta), rep(1L, 8))
})

test_that("directed = TRUE with undirected layers fails loudly", {
  layers <- planted12(TRUE); layers[[2]] <- igraph::as_undirected(layers[[2]])
  expect_error(
    detect_multislice_communities(layers, interlayer_ties = NULL, omega = 1, directed = TRUE),
    "requires every layer to be a directed"
  )
})

test_that("public fit symmetrises directed layers when directed = FALSE, with warnings", {
  # before 1.3.2 the per-layer igraph detection stopped with an error here
  w <- character(0)
  fit <- withCallingHandlers(
    fit_multilayer_identity_ties(directed4(), omega = directed_omega, directed = FALSE, seed = seed_dm),
    warning = function(cond) { w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning") }
  )
  expect_true(any(grepl("collapsing them to weighted undirected", w)))
  expect_true(any(grepl("symmetrised", w)))
  expect_identical(canon(fit$meta_communities), rep(1L, 8))
})

test_that("directed local moving leaves a community it has no ties to", {
  # arcs 1->2 and 3->4; start with 1 and 4 together and 2 and 3 together.
  # Moving vertex 1 to {2, 3} raises Q from -3 to -2 (gamma = 3); before the
  # fix the mover scored staying as 0 and kept the start partition.
  lm_dir <- dynamicmultiplex:::genlouvain_local_moving_directed
  ef <- c(1L, 3L); et <- c(2L, 4L); ew <- c(0.5, 0.5)
  k_out <- matrix(c(1, 0, 1, 0), 4, 1); k_in <- matrix(c(0, 1, 0, 1), 4, 1)
  q <- function(memb) {
    2 * sum(ew[memb[ef] == memb[et]]) - 3 * sum(rowsum(k_out, memb) * rowsum(k_in, memb)) / 2
  }
  set.seed(seed_dm)
  res <- lm_dir(ef, et, ew, k_out, k_in, 1 / 2, 3, c(1L, 2L, 2L, 1L), 50L)
  expect_gt(q(res), q(c(1L, 2L, 2L, 1L)))
})

test_that("a directed slice without arcs runs and returns finite memberships", {
  layers <- directed4()
  g3 <- igraph::make_empty_graph(n = 4, directed = TRUE); igraph::V(g3)$name <- as.character(1:4)
  layers[[3]] <- g3
  fit <- suppressWarnings(fit_multilayer_identity_ties(layers, omega = directed_omega,
                                                       directed = TRUE, seed = seed_dm))
  expect_length(fit$meta_communities, 3)
  expect_true(all(is.finite(unlist(fit$meta_communities))))
})

test_that("directed multislice fit is seed reproducible and keeps the fit fields", {
  a <- suppressWarnings(fit_multilayer_identity_ties(planted12(TRUE), omega = 6, directed = TRUE, seed = seed_dm))
  b <- suppressWarnings(fit_multilayer_identity_ties(planted12(TRUE), omega = 6, directed = TRUE, seed = seed_dm))
  expect_identical(a$meta_communities, b$meta_communities)
  expect_true(all(c("algorithm", "layer_communities", "layer_links", "interlayer_ties", "directed") %in% names(a)))
  expect_true(a$directed)
})

test_that("undirected inputs reproduce the pre-change memberships", {
  expected <- readRDS(testthat::test_path("fixtures", "multislice_regression_r.rds"))
  cases <- regression_cases()
  expect_setequal(names(expected), names(cases))
  for (nm in names(cases)) expect_identical(run_case(cases[[nm]]), expected[[nm]], info = nm)
})
