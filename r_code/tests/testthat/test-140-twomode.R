# 1.4.0: layers_from_membership() (mirrors python_code/tests/test_140_twomode.py).
# twomode_planted.csv and the reference CSVs were written by the Python package
# (same fixture), so the R-Python parity checks compare like with like. Seed 123.

hand_wide <- function() {
  m <- cbind(g1 = c(1, 1, 1, 0, 0), g2 = c(0, 0, 1, 1, 0), g3 = c(0, 1, 1, 1, 1))
  rownames(m) <- letters[1:5]
  m
}
edges140 <- function(g) {
  el <- igraph::as_edgelist(g, names = TRUE)
  w <- igraph::E(g)$weight
  a <- pmin(el[, 1], el[, 2]); b <- pmax(el[, 1], el[, 2])
  d <- data.frame(u = a, v = b, w = w, stringsAsFactors = FALSE)
  d[order(d$u, d$v), , drop = FALSE]
}
planted140 <- function() utils::read.csv(testthat::test_path("fixtures", "twomode_planted.csv"), stringsAsFactors = FALSE)

test_that("hand-computed count and size weights", {
  g <- layers_from_membership(list(hand_wide()), weighting = "count")[[1]]
  e <- edges140(g)
  expect_identical(paste(e$u, e$v), c("a b", "a c", "b c", "b d", "b e", "c d", "c e", "d e"))
  expect_equal(e$w, c(1, 1, 2, 1, 1, 2, 1, 1))
  g <- layers_from_membership(list(hand_wide()), weighting = "size")[[1]]
  expect_equal(edges140(g)$w, c(.5, .5, .5 + 1 / 3, 1 / 3, 1 / 3, 1 + 1 / 3, 1 / 3, 1 / 3), tolerance = 1e-12)
})

test_that("hand-computed Kumar targets", {
  r <- .kumar_targets(hand_wide(), c(1, 1, 1, 2, 2))
  expect_equal(unname(r), c(1.25 * 5 / 3, 4 / 3, 4 / 3), tolerance = 1e-12)
})

test_that("count equals the dyadic number of shared groups", {
  long <- planted140()
  m <- merge(long, long, by = "group")
  m <- m[m$node.x < m$node.y, ]
  tab <- stats::aggregate(group ~ node.x + node.y, data = m, FUN = length)
  g <- layers_from_membership(list(long), weighting = "count")[[1]]
  e <- edges140(g)
  got <- stats::setNames(e$w, paste(e$u, e$v))
  exp <- stats::setNames(as.numeric(tab$group), paste(pmin(as.character(tab$node.x), as.character(tab$node.y)),
                                                      pmax(as.character(tab$node.x), as.character(tab$node.y))))
  expect_equal(got[sort(names(got))], exp[sort(names(exp))])
})

test_that("R matches the Python package on the shared fixture", {
  long <- planted140()
  ref_e <- utils::read.csv(testthat::test_path("fixtures", "twomode_reference_edges.csv"), stringsAsFactors = FALSE)
  ref_w <- utils::read.csv(testthat::test_path("fixtures", "twomode_reference_weights.csv"), stringsAsFactors = FALSE)
  for (case in c("count|None", "size|None", "kumar|leiden", "kumar|louvain")) {
    p <- strsplit(case, "|", fixed = TRUE)[[1]]
    alg <- if (p[2] == "None") NULL else p[2]
    g <- layers_from_membership(list(long), weighting = p[1], algorithm = alg)[[1]]
    e <- edges140(g); e$u <- as.integer(e$u); e$v <- as.integer(e$v)
    a <- pmin(e$u, e$v); b <- pmax(e$u, e$v); e <- data.frame(u = a, v = b, w = e$w)[order(a, b), ]
    re <- ref_e[ref_e$case == case, ]
    expect_identical(paste(e$u, e$v), paste(re$u, re$v), info = case)
    expect_equal(e$w, re$weight, tolerance = 1e-10, info = case)
    rw <- ref_w[ref_w$case == case, ]
    att <- igraph::graph_attr(g, "dm_membership")
    expect_equal(unname(att$group_weights[as.character(rw$group)]), rw$weight, tolerance = 1e-10, info = case)
    expect_equal(att$iterations, rw$iterations[1], info = case)
  }
})

test_that("per-layer detection on Kumar's weights reproduces Kumar's partition", {
  long <- planted140()
  long <- long[!grepl("^gn", long$group), ]               # the clean home groups only
  for (alg in c("leiden", "louvain")) {
    g <- layers_from_membership(list(long), weighting = "kumar", algorithm = alg)[[1]]
    part <- igraph::graph_attr(g, "dm_membership")$partition
    set.seed(123)
    fit <- fit_layer_communities(list(g), algorithm = alg)[[1]]
    expect_equal(igraph::compare(as.integer(fit$membership), unname(part[igraph::V(g)$name]), method = "nmi"), 1,
                 info = alg)
  }
})

test_that("slice t depends only on slice t", {
  long <- planted140()
  S <- lapply(0:3, function(k) { d <- long; d$group <- paste0("t", k, "_", d$group)
                                  d[(seq_len(nrow(d)) + k) %% 7 != 0 | grepl("fix", d$group), ] })
  for (spec in list(list("count", NULL), list("size", NULL), list("kumar", "leiden"), list("kumar", "louvain"))) {
    full <- suppressWarnings(layers_from_membership(S, weighting = spec[[1]], algorithm = spec[[2]]))
    for (t in 1:4) {
      alone <- suppressWarnings(layers_from_membership(S[t], weighting = spec[[1]], algorithm = spec[[2]]))[[1]]
      expect_equal(edges140(alone), edges140(full[[t]]), info = paste(spec[[1]], t))
      expect_identical(igraph::graph_attr(alone, "dm_membership")$group_weights,
                       igraph::graph_attr(full[[t]], "dm_membership")$group_weights)
    }
  }
})

test_that("unequal node sets and group turnover are accepted by the fitters", {
  long <- planted140()
  s2 <- long; s2$node[s2$node == 39] <- 99
  s3 <- rbind(long, data.frame(node = 40:42, group = "new"))
  G <- layers_from_membership(list(long, s2, s3), weighting = "size")
  expect_equal(vapply(G, igraph::vcount, numeric(1)), c(40, 40, 43))
  out <- fit_multilayer_jaccard(G, allow_unequal_nodes = TRUE, seed = 123)
  expect_equal(lengths(out$meta_communities), c(40, 40, 43))
})

test_that("bad input fails loudly", {
  W <- hand_wide()
  bad <- W; bad[1, 1] <- 2
  expect_error(layers_from_membership(list(bad), weighting = "count"), "0/1")
  bad <- W; bad[1, 1] <- NA
  expect_error(layers_from_membership(list(bad), weighting = "count"), "missing")
  bad <- W; rownames(bad) <- c("a", "b", "c", "d", "a")
  expect_error(layers_from_membership(list(bad), weighting = "count"), "duplicate")
  bad <- W; colnames(bad) <- c("g1", "g1", "g3")
  expect_error(layers_from_membership(list(bad), weighting = "count"), "duplicate")
  long <- planted140()
  expect_error(layers_from_membership(list(rbind(long, long[1, ])), weighting = "count"), "duplicate")
  expect_error(layers_from_membership(list(data.frame(node = 1, org = "a", x = 1)), weighting = "count"), "node.*group")
  bad <- W; rownames(bad) <- NULL
  expect_error(layers_from_membership(list(bad), weighting = "count"), "row names")
  expect_error(layers_from_membership(W, weighting = "count"), "list")
  expect_error(layers_from_membership(list(), weighting = "count"), "at least one")
  expect_error(layers_from_membership(list(W[0, , drop = FALSE]), weighting = "count"), "empty")
  expect_error(layers_from_membership(list(W)), "no default")
  expect_error(layers_from_membership(list(W), weighting = "hyper"), "weighting")
  expect_error(layers_from_membership(list(W), weighting = "kumar"), "algorithm")
  expect_error(layers_from_membership(list(W), weighting = "count", algorithm = "leiden"), "algorithm")
  expect_error(layers_from_membership(list(W), weighting = "kumar", algorithm = "leiden", resolution_parameter = 0), "resolution_parameter")
  expect_error(layers_from_membership(list(W), weighting = "kumar", algorithm = "leiden", max_iter = 0), "max_iter")
  expect_error(layers_from_membership(list(W), weighting = "kumar", algorithm = "leiden", tol = 0), "tol")
  expect_error(layers_from_membership(list(W), weighting = "kumar", algorithm = "leiden", min_group_size = 1), "min_group_size")
  expect_error(layers_from_membership(list(W), weighting = "kumar", algorithm = "leiden", rescale = "max"), "rescale")
})

test_that("dropped groups warn and orphaned nodes fail", {
  W <- cbind(hand_wide(), solo = c(1, 0, 0, 0, 0))
  expect_warning(g <- layers_from_membership(list(W), weighting = "count")[[1]], "dropped 1 group")
  expect_equal(igraph::graph_attr(g, "dm_membership")$groups_dropped, 1L)
  O <- rbind(hand_wide(), f = 0); O <- cbind(O, solo = c(0, 0, 0, 0, 0, 1))
  expect_error(suppressWarnings(layers_from_membership(list(O), weighting = "size")), "slice 1.*'f'")
})

test_that("Kumar non-convergence fails loudly", {
  expect_error(layers_from_membership(list(planted140()), weighting = "kumar", algorithm = "leiden",
                                      max_iter = 1, tol = 1e-15), "did not converge")
})

test_that("attributes, seed and RNG restore", {
  long <- planted140()
  for (spec in list(list("count", NULL), list("size", NULL), list("kumar", "leiden"), list("kumar", "louvain"))) {
    set.seed(1); before <- stats::runif(1); set.seed(1)
    a <- layers_from_membership(list(long), weighting = spec[[1]], algorithm = spec[[2]])[[1]]
    expect_identical(stats::runif(1), before)                # caller's RNG stream untouched
    b <- layers_from_membership(list(long), weighting = spec[[1]], algorithm = spec[[2]])[[1]]
    expect_equal(edges140(a), edges140(b))
    att <- igraph::graph_attr(a, "dm_membership")
    expect_setequal(names(att), c("weighting", "group_weights", "iterations", "final_diff", "groups_dropped",
                                  "rescale", "partition"))
    expect_length(att$group_weights, length(unique(long$group)))
    if (spec[[1]] == "kumar") {
      expect_true(att$iterations >= 1 && att$iterations <= 50 && att$final_diff <= 0.01)
    } else {
      expect_equal(att$iterations, 0L); expect_null(att$partition); expect_true(all(att$group_weights == 1))
    }
  }
})

test_that("rescale = 'mean' gives mean nonzero tie 1", {
  long <- planted140()
  for (spec in list(list("count", NULL), list("size", NULL), list("kumar", "leiden"))) {
    raw <- edges140(layers_from_membership(list(long), weighting = spec[[1]], algorithm = spec[[2]])[[1]])$w
    sc <- edges140(layers_from_membership(list(long), weighting = spec[[1]], algorithm = spec[[2]], rescale = "mean")[[1]])$w
    expect_equal(mean(sc), 1, tolerance = 1e-12)
    expect_equal(sc, raw / mean(raw), tolerance = 1e-12)
  }
})
