# Shared fixture loaders for the 1.4.0 tests.
load_fixture_layers <- function(name, n_nodes) {
  ed <- utils::read.csv(testthat::test_path("fixtures", name))
  lapply(sort(unique(ed$layer)), function(t) {
    e <- ed[ed$layer == t, ]
    g <- igraph::make_empty_graph(n = n_nodes, directed = FALSE)
    igraph::V(g)$name <- as.character(seq_len(n_nodes))
    igraph::add_edges(g, as.vector(t(as.matrix(e[, c("from", "to")]))), weight = e$weight)
  })
}
sbm140_layers <- function() load_fixture_layers("multislice_sbm_edges.csv", 30)
chain140_layers <- function() load_fixture_layers("chain_panel_edges.csv", 60)
twostage_fitters <- function() list(
  jaccard = dynamicmultiplex::fit_multilayer_jaccard,
  overlap = dynamicmultiplex::fit_multilayer_overlap,
  weighted_jaccard = dynamicmultiplex::fit_multilayer_weighted_jaccard,
  weighted_overlap = dynamicmultiplex::fit_multilayer_weighted_overlap)
