# 1.4.0: chain repair and filtered tracking (mirrors python_code/tests/test_140_tracking.py).
# Fixtures: twostage_regression_r.rds and other_regression_r.rds recorded from
# 1.3.2 (commit b55c66c); chain_panel_edges.csv shared with Python. Seed 123.

ids_of <- function(meta) sort(unique(unlist(meta)))
same_partition <- function(a, b) {
  x <- unlist(a); y <- unlist(b)
  identical(outer(x, x, "=="), outer(y, y, "=="))
}

test_that("defaults reproduce 1.3.2 exactly", {
  rec <- readRDS(testthat::test_path("fixtures", "twostage_regression_r.rds"))
  expect_identical(rec$recorded_from, "b55c66c")
  panels <- list(sbm = sbm140_layers(), chain = chain140_layers())
  for (key in names(rec$cases)) {
    p <- strsplit(key, "|", fixed = TRUE)[[1]]
    f <- twostage_fitters()[[p[2]]]
    got <- lapply(f(panels[[p[1]]], algorithm = p[3], seed = 123)$meta_communities, as.integer)
    expect_identical(got, rec$cases[[key]], info = key)
    got2 <- lapply(f(panels[[p[1]]], algorithm = p[3], seed = 123, repair_chains = FALSE,
                     tracking = "smoothed")$meta_communities, as.integer)
    expect_identical(got2, rec$cases[[key]], info = key)
  }
})

test_that("Hungarian output reproduces 1.3.2 exactly", {
  rec <- readRDS(testthat::test_path("fixtures", "other_regression_r.rds"))
  L <- chain140_layers()
  for (alg in c("leiden", "louvain")) {
    got <- lapply(fit_multilayer_hungarian(L, algorithm = alg, seed = 123)$meta_communities, as.integer)
    expect_identical(got, rec$cases[[paste0("chain|hungarian|", alg)]], info = alg)
  }
})

# hand-built case: the same case and answers as the Python tests
hand_meta <- c(L1C1 = 1L, L2C1 = 1L, L3C1 = 2L, L4C1 = 2L, L1C2 = 3L, L2C2 = 3L, L3C2 = 4L, L4C2 = 4L,
               L1C3 = 6L, L3C3 = 5L)
hand_ties <- data.frame(
  from_layer = c(1, 1, 1, 2, 2, 2, 3, 3, 3), to_layer = c(2, 2, 2, 3, 3, 3, 4, 4, 4),
  from_community = c(1, 2, 3, 1, 2, 2, 1, 2, 3), to_community = c(1, 2, 1, 1, 2, 3, 1, 2, 2),
  similarity = c(.8, .9, .9, .7, .15, .1, .9, .6, .3))
hand_links <- data.frame(from = 1:3, to = 2:4, weight = 1)

test_that("repair rule on the hand case, floor 0.2", {
  r <- repair_meta_chains(hand_meta, hand_ties, hand_links, min_similarity = 0.2)
  exp <- c(L1C1 = 1L, L2C1 = 1L, L3C1 = 1L, L4C1 = 1L, L1C2 = 2L, L2C2 = 2L, L1C3 = 3L,
           L3C2 = 4L, L4C2 = 4L, L3C3 = 5L)
  expect_identical(r$meta_ids[names(exp)], exp)
  expect_equal(nrow(r$log), 1L)
  expect_equal(unlist(r$log[1, ]), c(from_layer = 2, from_community = 1, to_layer = 3, to_community = 1,
                                     similarity = 0.7))
})

test_that("repair rule on the hand case, floor 0.1", {
  r <- repair_meta_chains(hand_meta, hand_ties, hand_links, min_similarity = 0.1)
  exp <- c(L1C1 = 1L, L2C1 = 1L, L3C1 = 1L, L4C1 = 1L, L1C2 = 2L, L2C2 = 2L, L3C2 = 2L, L4C2 = 2L,
           L1C3 = 3L, L3C3 = 4L)
  expect_identical(r$meta_ids[names(exp)], exp)
  expect_equal(nrow(r$log), 2L)
})

test_that("identical layers keep one id per block with repair", {
  base <- chain140_layers()[[1]]
  for (fn in names(twostage_fitters())) for (alg in c("leiden", "louvain")) for (tr in c("smoothed", "filtered")) {
    for (TT in c(4, 12, 20)) {
      out <- twostage_fitters()[[fn]](rep(list(base), TT), algorithm = alg, seed = 123,
                                      repair_chains = TRUE, tracking = tr)
      M <- out$meta_communities
      expect_length(ids_of(M), 3)
      expect_true(all(vapply(M, function(m) identical(m, M[[1]]), logical(1))), info = paste(fn, alg, tr, TT))
    }
  }
})

test_that("filtered output never depends on later layers", {
  L <- chain140_layers()
  for (fn in names(twostage_fitters())) for (alg in c("leiden", "louvain")) for (rep in c(FALSE, TRUE)) {
    f <- twostage_fitters()[[fn]]
    full <- f(L, algorithm = alg, seed = 123, tracking = "filtered", repair_chains = rep)$meta_communities
    for (t in c(2, 5, 7)) {
      trunc <- f(L[seq_len(t)], algorithm = alg, seed = 123, tracking = "filtered", repair_chains = rep)$meta_communities
      expect_identical(full[seq_len(t)], trunc, info = paste(fn, alg, rep, t))
    }
  }
})

test_that("filtered output ignores link weights that touch only later layers", {
  L <- chain140_layers()
  base <- data.frame(from = 1:9, to = 2:10, weight = 1)
  heavy <- base; heavy$weight[8] <- 5                       # link 8 -> 9
  for (rep in c(FALSE, TRUE)) {
    a <- fit_multilayer_jaccard(L, layer_links = base, seed = 123, tracking = "filtered", repair_chains = rep)
    b <- fit_multilayer_jaccard(L, layer_links = heavy, seed = 123, tracking = "filtered", repair_chains = rep)
    expect_identical(a$meta_communities[1:8], b$meta_communities[1:8], info = rep)
  }
})

test_that("repair only joins groups that never coexist", {
  for (L in list(chain140_layers(), sbm140_layers())) {
    for (fn in names(twostage_fitters())) for (alg in c("leiden", "louvain")) {
      f <- twostage_fitters()[[fn]]
      pre <- f(L, algorithm = alg, seed = 123)$meta_communities
      post <- f(L, algorithm = alg, seed = 123, repair_chains = TRUE)$meta_communities
      d <- data.frame(layer = rep(seq_along(pre), lengths(pre)), pre = unlist(pre), post = unlist(post))
      d <- unique(d)
      for (y in unique(d$post)) {
        sub <- d[d$post == y, ]
        spans <- split(sub$layer, sub$pre)
        if (length(spans) > 1) for (i in 1:(length(spans) - 1)) for (j in (i + 1):length(spans))
          expect_length(intersect(spans[[i]], spans[[j]]), 0)
      }
    }
  }
})

test_that("repair respects the floor and logs every join", {
  L <- chain140_layers()
  out <- fit_multilayer_jaccard(L, seed = 123, repair_chains = TRUE)
  log <- out$chain_repairs
  expect_setequal(names(log), c("from_layer", "from_community", "to_layer", "to_community", "similarity"))
  expect_gt(nrow(log), 0)
  expect_true(all(log$similarity >= 0.2))
  strict <- fit_multilayer_jaccard(L, seed = 123, repair_chains = TRUE, repair_min_similarity = 0.6)$chain_repairs
  expect_true(all(strict$similarity >= 0.6) && nrow(strict) <= nrow(log))
  before <- length(ids_of(fit_multilayer_jaccard(L, seed = 123)$meta_communities))
  expect_lt(length(ids_of(out$meta_communities)), before)
})

test_that("fit fields and argument checks", {
  L <- sbm140_layers()
  out <- fit_multilayer_jaccard(L, seed = 123)
  expect_identical(out$tracking, "smoothed")
  expect_false(out$repair_chains)
  expect_equal(out$repair_min_similarity, 0.2)
  expect_equal(nrow(out$chain_repairs), 0)
  expect_error(fit_multilayer_jaccard(L, tracking = "causal"), "tracking")
  expect_error(fit_multilayer_jaccard(L, repair_chains = TRUE, repair_min_similarity = -0.1), "repair_min_similarity")
  expect_error(fit_multilayer_jaccard(L, repair_chains = TRUE, repair_min_similarity = 1.5), "repair_min_similarity")
})

test_that("unequal node sets and gapped links work with repair and filtered tracking", {
  L <- chain140_layers()
  L[[4]] <- igraph::delete_vertices(L[[4]], "60")
  links <- data.frame(from = c(1, 2, 3, 4, 5, 6, 7, 8), to = c(3, 4, 5, 6, 7, 8, 9, 10), weight = 1)
  for (tr in c("smoothed", "filtered")) {
    out <- fit_multilayer_jaccard(L, layer_links = links, seed = 123, allow_unequal_nodes = TRUE,
                                  repair_chains = TRUE, tracking = tr)
    expect_equal(lengths(out$meta_communities), c(60, 60, 60, 59, 60, 60, 60, 60, 60, 60))
  }
  full <- fit_multilayer_jaccard(L, layer_links = links, seed = 123, allow_unequal_nodes = TRUE,
                                 repair_chains = TRUE, tracking = "filtered")$meta_communities
  trunc <- fit_multilayer_jaccard(L[1:6], layer_links = links[links$to <= 6, ], seed = 123,
                                  allow_unequal_nodes = TRUE, repair_chains = TRUE, tracking = "filtered")$meta_communities
  expect_identical(full[1:6], trunc)
})

test_that("bootstrap passes the new settings through to the refits", {
  L <- lapply(chain140_layers()[1:3], function(g) igraph::as_adjacency_matrix(g, attr = "weight", sparse = FALSE))
  res <- bootstrap_multilayer(L, n_boot = 2, seed = 123, repair_chains = TRUE, repair_min_similarity = 0.3)
  expect_true(res$point_estimate$repair_chains)
  expect_equal(res$point_estimate$repair_min_similarity, 0.3)
  res <- bootstrap_multilayer(L, n_boot = 2, seed = 123, tracking = "filtered")
  expect_identical(res$point_estimate$tracking, "filtered")
  expect_error(bootstrap_multilayer(L, fit_type = "identity", n_boot = 2, seed = 123, repair_chains = TRUE),
               "repair_chains")
})

test_that("undirected multislice local mover leaves a community it has no ties to", {
  # edges 1-2 and 3-4, start with {1, 4} and {2, 3}; gamma = 3. Moving vertex 1
  # to {2, 3} raises Q; before 1.4.0 the mover scored staying as 0 and could
  # keep the start partition.
  lm <- genlouvain_local_moving
  ef <- c(1L, 3L); et <- c(2L, 4L); ew <- c(1, 1)
  K <- matrix(c(1, 1, 1, 1), 4, 1)
  q <- function(memb) {
    2 * sum(ew[memb[ef] == memb[et]]) - 3 * sum(rowsum(K, memb)^2) / 4
  }
  set.seed(123)
  res <- lm(ef, et, ew, K, 1 / 4, 3, c(1L, 2L, 2L, 1L), 50L)
  expect_gt(q(res), q(c(1L, 2L, 2L, 1L)))
})

test_that("filtered tracking requires a seed (Codex review)", {
  # Unseeded stage 1 on later layers advances the RNG that filtered stage 2
  # draws from, so append-invariance needs a seed.
  L <- chain140_layers()
  expect_error(fit_multilayer_jaccard(L, seed = NULL, tracking = "filtered"), "seed")
  expect_error(bootstrap_multilayer(lapply(L, igraph::as_adjacency_matrix, attr = "weight", sparse = FALSE),
                                    n_boot = 2, seed = NULL, tracking = "filtered"), "seed")
  expect_no_error(fit_multilayer_jaccard(L, seed = NULL))
})

test_that("bootstrap forwards and validates the repair arguments (Codex review)", {
  L <- lapply(sbm140_layers(), igraph::as_adjacency_matrix, attr = "weight", sparse = FALSE)
  res <- bootstrap_multilayer(L, n_boot = 2, seed = 123, repair_min_similarity = 0.8)
  expect_identical(res$point_estimate$repair_min_similarity, 0.8)
  expect_false(res$point_estimate$repair_chains)
  expect_error(bootstrap_multilayer(L, n_boot = 2, seed = 123, repair_min_similarity = 1.5), "repair_min_similarity")
  expect_error(bootstrap_multilayer(L, n_boot = 2, seed = 123, repair_chains = 1), "repair_chains")
  expect_error(bootstrap_multilayer(L, n_boot = 2, seed = 123, tracking = "smooth"), "tracking")
})
