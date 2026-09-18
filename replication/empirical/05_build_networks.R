# =============================================================================
# 05_build_networks.R
# Build the four real dynamic international networks used in the empirical
# application:
#   atop  : ATOP defense|offense alliances, 1816-2018 (peacesciencer)
#   dca   : Kinne DCAD v1.0 defense cooperation agreements (dcaAnyV1 == 1)
#   igo   : COW/peacesciencer shared-IGO membership counts, 1816-2014
#   trade : COW bilateral trade, 1870-2014, top-quartile dyads per year, log(1+x)
#
# For each network two files are written to EMP_DATA:
#   <net>_series.rds : list(years, graph_layers)   native per-year igraph objects
#                      on the union node set (vertex names = COW ccodes)
#   <net>_union.rds  : list(layers, active, present, names)
#                      layers  = per-year N x N adjacency on the union node set
#                      active  = per-year logical mask, TRUE if the state has
#                                at least one tie that year (degree > 0)
#                      present = identical copy of `active` (see NOTE below)
#                      names   = union ccodes as character
#
# Consolidates: replication/extended/00_build_atop.R, 00_build_dca.R,
#               paper_scripts/build_igo.R, paper_scripts/build_trade.R
#
# Usage:  Rscript 05_build_networks.R            # builds all four
#         Rscript 05_build_networks.R <net>      # net in {atop, dca, igo, trade}
# Env:    DM_ROOT     project root (default getwd())
#         ATOP_YEARS  R expression for ATOP years (default "1816:2018")
#         IGO_THRESH  min shared-IGO count for an IGO tie (default 1)
# Input:  $DM_ROOT/data/DCAD-v1.0-dyadic.csv  (DCA only; Kinne DCAD v1.0)
# =============================================================================
set.seed(123)
suppressMessages({ library(peacesciencer); library(dplyr); library(igraph) })

ROOT     <- Sys.getenv("DM_ROOT", unset = getwd())
EMP_DATA <- file.path(ROOT, "output", "empirical_data")
dir.create(EMP_DATA, recursive = TRUE, showWarnings = FALSE)
stopifnot(dir.exists(EMP_DATA))

args <- commandArgs(trailingOnly = TRUE)
NETS <- if (length(args) >= 1) args[1] else c("atop", "dca", "igo", "trade")
stopifnot(all(NETS %in% c("atop", "dca", "igo", "trade")))

# NOTE: The fitter (06) and the scorers (07) consume the per-layer node mask
# under the name `present` (emp_fix.R, multinet_fix.R, 29/32), while the four
# original builders only ever saved it as `active`. No script in the sources
# constructs a separate `present` mask, so the two are the same object: the
# degree > 0 mask. It is saved under both names so every consumer works.

# ---------------------------------------------------------------------------
# ATOP: defense or offense pact, undirected, 1816-2018 (00_build_atop.R)
# ---------------------------------------------------------------------------
if ("atop" %in% NETS) {
  YEARS <- eval(parse(text = Sys.getenv("ATOP_YEARS", "1816:2018")))
  stopifnot(is.numeric(YEARS), length(YEARS) >= 2)
  dy  <- create_dyadyears(subset_years = YEARS, directed = FALSE)
  dy  <- add_atop_alliance(dy)
  stopifnot(all(c("atop_defense", "atop_offense") %in% names(dy)))
  def <- ifelse(is.na(dy$atop_defense), 0L, dy$atop_defense)
  off <- ifelse(is.na(dy$atop_offense), 0L, dy$atop_offense)
  dy$tie <- as.integer(def == 1L | off == 1L)
  E <- dy[dy$tie == 1L, c("ccode1", "ccode2", "year")]
  stopifnot(nrow(E) > 0)

  ucodes <- sort(unique(c(E$ccode1, E$ccode2))); N <- length(ucodes)
  idx <- setNames(seq_len(N), as.character(ucodes))
  graph_layers <- vector("list", length(YEARS))
  layers <- vector("list", length(YEARS)); active <- vector("list", length(YEARS))
  for (k in seq_along(YEARS)) {
    e <- E[E$year == YEARS[k], ]
    A <- matrix(0L, N, N, dimnames = list(as.character(ucodes), as.character(ucodes)))
    if (nrow(e) > 0) { ii <- idx[as.character(e$ccode1)]; jj <- idx[as.character(e$ccode2)]
      A[cbind(ii, jj)] <- 1L; A[cbind(jj, ii)] <- 1L }
    layers[[k]] <- A; active[[k]] <- rowSums(A) > 0
    g <- graph_from_adjacency_matrix(A, mode = "undirected", diag = FALSE)
    igraph::V(g)$name <- as.character(ucodes); graph_layers[[k]] <- g
  }
  stopifnot(all(sapply(layers, nrow) == N), all(sapply(active, length) == N))
  saveRDS(list(years = YEARS, graph_layers = graph_layers), file.path(EMP_DATA, "atop_series.rds"))
  saveRDS(list(layers = layers, active = active, present = active, names = as.character(ucodes)),
          file.path(EMP_DATA, "atop_union.rds"))
  stopifnot(file.exists(file.path(EMP_DATA, "atop_series.rds")),
            file.exists(file.path(EMP_DATA, "atop_union.rds")))
  cat("BUILT atop | years", length(YEARS), "| union", N,
      "| mean active", round(mean(sapply(active, sum)), 1),
      "| mean edges", round(mean(sapply(layers, function(m) sum(m) / 2)), 1), "\n")
  rm(dy, E, def, off, YEARS, ucodes, N, idx, graph_layers, layers, active); invisible(gc())
}

# ---------------------------------------------------------------------------
# DCA: Kinne DCAD v1.0, any DCA (narrow coding dcaAnyV1 == 1), undirected
# (00_build_dca.R). Years are those with at least one DCA tie.
# ---------------------------------------------------------------------------
if ("dca" %in% NETS) {
  dca_file <- file.path(ROOT, "data", "DCAD-v1.0-dyadic.csv")
  if (!file.exists(dca_file))
    stop("DCA input not found: ", dca_file,
         "\nPlace Kinne's DCAD v1.0 dyadic file (DCAD-v1.0-dyadic.csv) at $DM_ROOT/data/.",
         call. = FALSE)
  d <- read.csv(dca_file, stringsAsFactors = FALSE)
  stopifnot(all(c("ccode1", "ccode2", "year", "dcaAnyV1") %in% names(d)))
  tie <- ifelse(is.na(d$dcaAnyV1), 0L, d$dcaAnyV1) == 1L
  E <- d[tie, c("ccode1", "ccode2", "year")]
  stopifnot(nrow(E) > 0)
  YEARS <- sort(unique(E$year))
  stopifnot(length(YEARS) >= 2)

  ucodes <- sort(unique(c(E$ccode1, E$ccode2))); N <- length(ucodes)
  idx <- setNames(seq_len(N), as.character(ucodes))
  graph_layers <- vector("list", length(YEARS))
  layers <- vector("list", length(YEARS)); active <- vector("list", length(YEARS))
  for (k in seq_along(YEARS)) {
    e <- E[E$year == YEARS[k], ]
    A <- matrix(0L, N, N, dimnames = list(as.character(ucodes), as.character(ucodes)))
    if (nrow(e) > 0) { ii <- idx[as.character(e$ccode1)]; jj <- idx[as.character(e$ccode2)]
      A[cbind(ii, jj)] <- 1L; A[cbind(jj, ii)] <- 1L }
    layers[[k]] <- A; active[[k]] <- rowSums(A) > 0
    g <- graph_from_adjacency_matrix(A, mode = "undirected", diag = FALSE)
    igraph::V(g)$name <- as.character(ucodes); graph_layers[[k]] <- g
  }
  stopifnot(all(sapply(layers, nrow) == N), all(sapply(active, length) == N))
  saveRDS(list(years = YEARS, graph_layers = graph_layers), file.path(EMP_DATA, "dca_series.rds"))
  saveRDS(list(layers = layers, active = active, present = active, names = as.character(ucodes)),
          file.path(EMP_DATA, "dca_union.rds"))
  stopifnot(file.exists(file.path(EMP_DATA, "dca_series.rds")),
            file.exists(file.path(EMP_DATA, "dca_union.rds")))
  cat("BUILT dca | years", length(YEARS), sprintf("(%d-%d)", min(YEARS), max(YEARS)), "| union", N,
      "| mean active", round(mean(sapply(active, sum)), 1),
      "| mean edges", round(mean(sapply(layers, function(m) sum(m) / 2)), 1), "\n")
  rm(d, E, tie, YEARS, ucodes, N, idx, graph_layers, layers, active); invisible(gc())
}

# ---------------------------------------------------------------------------
# IGO: shared-IGO membership counts, 1816-2014, tie if count >= IGO_THRESH,
# weighted by the count (build_igo.R). Years are those with at least one tie.
# ---------------------------------------------------------------------------
if ("igo" %in% NETS) {
  THRESH <- as.integer(Sys.getenv("IGO_THRESH", "1"))
  stopifnot(!is.na(THRESH), THRESH >= 1L)
  dy <- create_dyadyears(subset_years = 1816:2014, directed = FALSE)
  dy <- add_igos(dy)
  stopifnot("dyadigos" %in% names(dy))
  dy$w <- ifelse(is.na(dy$dyadigos), 0, dy$dyadigos)
  E <- dy[dy$w >= THRESH, c("ccode1", "ccode2", "year", "w")]
  stopifnot(nrow(E) > 0)
  YEARS <- sort(unique(E$year))
  cat("IGO years with data:", length(YEARS), "range", paste(range(YEARS), collapse = "-"), "\n")
  stopifnot(length(YEARS) >= 2)

  ucodes <- sort(unique(c(E$ccode1, E$ccode2))); N <- length(ucodes)
  idx <- setNames(seq_len(N), as.character(ucodes))
  graph_layers <- vector("list", length(YEARS))
  layers <- vector("list", length(YEARS)); active <- vector("list", length(YEARS))
  for (k in seq_along(YEARS)) {
    e <- E[E$year == YEARS[k], ]
    A <- matrix(0, N, N, dimnames = list(as.character(ucodes), as.character(ucodes)))
    if (nrow(e) > 0) { ii <- idx[as.character(e$ccode1)]; jj <- idx[as.character(e$ccode2)]
      A[cbind(ii, jj)] <- e$w; A[cbind(jj, ii)] <- e$w }
    layers[[k]] <- A; active[[k]] <- rowSums(A > 0) > 0
    g <- graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
    igraph::V(g)$name <- as.character(ucodes); graph_layers[[k]] <- g
  }
  stopifnot(all(sapply(layers, nrow) == N), all(sapply(active, length) == N))
  saveRDS(list(years = YEARS, graph_layers = graph_layers), file.path(EMP_DATA, "igo_series.rds"))
  saveRDS(list(layers = layers, active = active, present = active, names = as.character(ucodes)),
          file.path(EMP_DATA, "igo_union.rds"))
  stopifnot(file.exists(file.path(EMP_DATA, "igo_series.rds")),
            file.exists(file.path(EMP_DATA, "igo_union.rds")))
  cat("BUILT igo | years", length(YEARS), "| union", N,
      "| mean active", round(mean(sapply(active, sum)), 1),
      "| mean edges", round(mean(sapply(layers, function(m) sum(m > 0) / 2)), 0),
      "| mean wt", round(mean(unlist(lapply(layers, function(m) m[m > 0]))), 1), "\n")
  rm(dy, E, YEARS, THRESH, ucodes, N, idx, graph_layers, layers, active); invisible(gc())
}

# ---------------------------------------------------------------------------
# TRADE: COW bilateral trade, 1870-2014; per year keep dyads with total flow at
# or above the 75th percentile of positive flows (years with < 4 positive dyads
# dropped); weight = log(total + 1) (build_trade.R).
# ---------------------------------------------------------------------------
if ("trade" %in% NETS) {
  dy <- create_dyadyears(subset_years = 1870:2014, directed = FALSE)
  dy <- add_cow_trade(dy)
  stopifnot(all(c("flow1", "flow2") %in% names(dy)))
  f1 <- ifelse(is.na(dy$flow1), 0, dy$flow1); f2 <- ifelse(is.na(dy$flow2), 0, dy$flow2)
  dy$tot <- f1 + f2
  YEARS0 <- sort(unique(dy$year)); keep <- rep(FALSE, nrow(dy))
  for (y in YEARS0) {
    i <- which(dy$year == y & dy$tot > 0); if (length(i) < 4) next
    thr <- as.numeric(quantile(dy$tot[i], 0.75)); keep[i[dy$tot[i] >= thr]] <- TRUE
  }
  E <- dy[keep, c("ccode1", "ccode2", "year")]; E$w <- log(dy$tot[keep] + 1)
  stopifnot(nrow(E) > 0, all(E$w > 0))
  YEARS <- sort(unique(E$year))
  stopifnot(length(YEARS) >= 2)

  ucodes <- sort(unique(c(E$ccode1, E$ccode2))); N <- length(ucodes)
  idx <- setNames(seq_len(N), as.character(ucodes))
  graph_layers <- vector("list", length(YEARS))
  layers <- vector("list", length(YEARS)); active <- vector("list", length(YEARS))
  for (k in seq_along(YEARS)) {
    e <- E[E$year == YEARS[k], ]
    A <- matrix(0, N, N, dimnames = list(as.character(ucodes), as.character(ucodes)))
    if (nrow(e) > 0) { ii <- idx[as.character(e$ccode1)]; jj <- idx[as.character(e$ccode2)]
      A[cbind(ii, jj)] <- e$w; A[cbind(jj, ii)] <- e$w }
    layers[[k]] <- A; active[[k]] <- rowSums(A > 0) > 0
    g <- graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
    igraph::V(g)$name <- as.character(ucodes); graph_layers[[k]] <- g
  }
  stopifnot(all(sapply(layers, nrow) == N), all(sapply(active, length) == N))
  saveRDS(list(years = YEARS, graph_layers = graph_layers), file.path(EMP_DATA, "trade_series.rds"))
  saveRDS(list(layers = layers, active = active, present = active, names = as.character(ucodes)),
          file.path(EMP_DATA, "trade_union.rds"))
  stopifnot(file.exists(file.path(EMP_DATA, "trade_series.rds")),
            file.exists(file.path(EMP_DATA, "trade_union.rds")))
  cat("BUILT trade | years", length(YEARS), "| range", paste(range(YEARS), collapse = "-"),
      "| union", N,
      "| mean active", round(mean(sapply(active, sum)), 1),
      "| mean edges", round(mean(sapply(layers, function(m) sum(m > 0) / 2)), 0), "\n")
  rm(dy, E, f1, f2, YEARS0, YEARS, keep, ucodes, N, idx, graph_layers, layers, active); invisible(gc())
}

cat("BUILD_DONE:", paste(NETS, collapse = " "), "->", EMP_DATA, "\n")
