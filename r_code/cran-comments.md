# dynamicmultiplex 1.3.2

Update of the CRAN release 1.3.1 (2026-09-24). No changes were requested by
CRAN; this is a maintainer release with bug fixes. NEWS.md lists every change.

Summary of changes:

* Default community algorithm is now `"leiden"` (was `"louvain"`) in all fit
  functions, `bootstrap_multilayer()` and `simulate_and_fit_multilayer()`,
  matching the Python package. `algorithm = "louvain"` reproduces the old
  default.
* Shared input policy: self-loops and zero-weight edges are dropped for both
  matrix and igraph input; layers with no vertices stop with an error.
* `fit_multilayer_hungarian()` matches nodes by vertex name and gains
  `allow_unequal_nodes`.
* Plotting fixes for graphs with named vertices.
* `fit_multilayer_identity_ties(directed = TRUE)` now keeps arc direction in
  the multislice stage (directed null model); undirected results unchanged.
* No new dependencies. Imports remain clue, igraph (>= 2.0.0), rlang.

## Test environments

* Ubuntu 22.04, R 4.3.3 (R CMD check --as-cran, examples with --run-donttest)
* local macOS, R 4.5.x
* win-builder (R-devel)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Downstream dependencies

There are currently no downstream dependencies for this package.
