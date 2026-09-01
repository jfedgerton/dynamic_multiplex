# Build DCA network (Kinne DCAD v1.0) for 08_alliance_dca_empirical.R
# Edge = any DCA, narrow coding (dcaAnyV1 == 1), undirected. Outputs dca_series.rds + dca_union.rds to EMP_DATA.
set.seed(123)
suppressMessages({library(igraph)})
DATA <- Sys.getenv("EMP_DATA", "/tmp")
d <- read.csv(file.path(DATA, "DCAD-v1.0-dyadic.csv"), stringsAsFactors = FALSE)
tie <- ifelse(is.na(d$dcaAnyV1), 0L, d$dcaAnyV1) == 1L
E <- d[tie, c("ccode1","ccode2","year")]
YEARS <- sort(unique(E$year))
ucodes <- sort(unique(c(E$ccode1, E$ccode2))); N <- length(ucodes)
idx <- setNames(seq_len(N), as.character(ucodes))
graph_layers <- vector("list", length(YEARS)); layers <- vector("list", length(YEARS)); active <- vector("list", length(YEARS))
for (k in seq_along(YEARS)) {
  e <- E[E$year == YEARS[k], ]
  A <- matrix(0L, N, N, dimnames = list(as.character(ucodes), as.character(ucodes)))
  if (nrow(e) > 0) { ii <- idx[as.character(e$ccode1)]; jj <- idx[as.character(e$ccode2)]
    A[cbind(ii,jj)] <- 1L; A[cbind(jj,ii)] <- 1L }
  layers[[k]] <- A; active[[k]] <- rowSums(A) > 0
  g <- graph_from_adjacency_matrix(A, mode = "undirected", diag = FALSE)
  igraph::V(g)$name <- as.character(ucodes); graph_layers[[k]] <- g
}
saveRDS(list(years = YEARS, graph_layers = graph_layers), file.path(DATA, "dca_series.rds"))
saveRDS(list(layers = layers, active = active, names = as.character(ucodes)), file.path(DATA, "dca_union.rds"))
cat("BUILT dca | years", length(YEARS), sprintf("(%d-%d)", min(YEARS), max(YEARS)), "| union", N,
    "| mean active", round(mean(sapply(active, sum)), 1), "| mean edges", round(mean(sapply(layers, function(m) sum(m)/2)), 1), "\n")
