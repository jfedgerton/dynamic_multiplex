# Build ATOP defense|offense dynamic alliance network for 06_alliance_dca_empirical.R
# Outputs (to EMP_DATA): atop_series.rds (native per-year graphs) + atop_union.rds (union set + active masks)
set.seed(123)
suppressMessages({library(peacesciencer); library(dplyr); library(igraph)})
YEARS <- eval(parse(text = Sys.getenv("ATOP_YEARS", "1816:2018")))
OUT   <- Sys.getenv("EMP_DATA", "/tmp")

dy <- create_dyadyears(subset_years = YEARS, directed = FALSE)
dy <- add_atop_alliance(dy)
def <- ifelse(is.na(dy$atop_defense), 0L, dy$atop_defense)
off <- ifelse(is.na(dy$atop_offense), 0L, dy$atop_offense)
dy$tie <- as.integer(def == 1L | off == 1L)
E <- dy[dy$tie == 1L, c("ccode1","ccode2","year")]

ucodes <- sort(unique(c(E$ccode1, E$ccode2))); N <- length(ucodes)
idx <- setNames(seq_len(N), as.character(ucodes))

graph_layers <- vector("list", length(YEARS))
layers <- vector("list", length(YEARS)); active <- vector("list", length(YEARS))
for (k in seq_along(YEARS)) {
  e <- E[E$year == YEARS[k], ]
  A <- matrix(0L, N, N, dimnames = list(as.character(ucodes), as.character(ucodes)))
  if (nrow(e) > 0) { ii <- idx[as.character(e$ccode1)]; jj <- idx[as.character(e$ccode2)]
    A[cbind(ii,jj)] <- 1L; A[cbind(jj,ii)] <- 1L }
  layers[[k]] <- A; active[[k]] <- rowSums(A) > 0
  g <- graph_from_adjacency_matrix(A, mode = "undirected", diag = FALSE)
  igraph::V(g)$name <- as.character(ucodes); graph_layers[[k]] <- g
}
saveRDS(list(years = YEARS, graph_layers = graph_layers), file.path(OUT, "atop_series.rds"))
saveRDS(list(layers = layers, active = active, names = as.character(ucodes)), file.path(OUT, "atop_union.rds"))
cat("BUILT atop | years", length(YEARS), "| union", N,
    "| mean active", round(mean(sapply(active, sum)), 1),
    "| mean edges", round(mean(sapply(layers, function(m) sum(m)/2)), 1), "\n")
