# dynamicmultiplex 1.4.0

Update of the CRAN release 1.3.1 (2026-09-24). No changes were requested by
CRAN. NEWS.md lists every change; it includes the unreleased 1.3.2 changes.

Summary of changes:

* Breaking: `co_assignment_ci()` and `community_est()` are removed and
  `partition_stability()` returns only the stability score; the bundled
  calibration table in inst/extdata is removed.
* Breaking: fit functions default to `seed = 123` (was `NULL`).
* New: chain repair (`repair_chains`), past-only tracking
  (`tracking = "filtered"`), `layers_from_membership()`, and
  `allow_unequal_nodes` in `bootstrap_multilayer()`.
* Fixes: exact stay score in the undirected multislice mover; directed
  multislice slices keep arc direction; Leiden default; shared input policy.
* No new dependencies. Imports remain clue, igraph (>= 2.0.0), rlang.

## Test environments

* Ubuntu 22.04, R 4.3.3 (R CMD check --as-cran)
* local macOS (maintainer)
* win-builder (R-devel)

## R CMD check results

0 errors | 0 warnings | 0 notes
