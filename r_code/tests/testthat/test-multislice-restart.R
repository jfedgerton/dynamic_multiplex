test_that("multislice keeps the restart with the highest multislice modularity", {
  # two identical 8-node layers joined by identity ties (omega = 1)
  e1 <- data.frame(from = c(1, 2, 2, 3, 1, 2, 3, 4, 5, 2),
                   to   = c(2, 4, 5, 5, 6, 6, 6, 7, 7, 8), w = 1)
  k <- tabulate(c(e1$from, e1$to), 8)
  e <- rbind(e1, transform(e1, from = from + 8, to = to + 8),
             data.frame(from = 1:8, to = 9:16, w = 1))
  K <- rbind(cbind(k, 0), cbind(0, k)); tm <- colSums(K)
  mucha <- function(m) {
    2 * sum(e$w[m[e$from] == m[e$to]]) - sum(sweep(rowsum(K, m)^2, 2, 1 / tm, "*"))
  }
  set.seed(123)
  cands <- replicate(10, dynamicmultiplex:::genlouvain_multislice_once(e, K, tm, 1, 20L, 50L),
                     simplify = FALSE)
  set.seed(123)
  chosen <- dynamicmultiplex:::genlouvain_multislice(e, K, tm)
  expect_equal(mucha(chosen), max(vapply(cands, mucha, numeric(1))))
})
