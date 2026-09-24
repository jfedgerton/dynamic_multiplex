# dynamicmultiplex 1.3.1

Update of the CRAN release 1.1.0 (2026-08-07). No changes were requested
by CRAN; this is a maintainer release. NEWS.md lists every change since
1.1.0 (the intermediate 1.2.1 and 1.3.0 were internal and never submitted).

Summary of changes:

* New exported function `partition_stability()` (stability score and a
  calibrated accuracy floor for a tracked partition), with its calibration
  table shipped in `inst/extdata/` (31 rows, 2 KB).
* `bootstrap_multilayer()` returns two additional list elements consumed by
  `partition_stability()`; existing elements are unchanged.
* New `allow_unequal_nodes` argument (default `FALSE`, previous behaviour)
  on the `fit_multilayer_*()` functions.
* Bug fixes to `fit_multilayer_identity_ties()` (multislice null model) and
  to the two weighted similarity helpers; regression tests added.
* No new dependencies. Imports remain clue, igraph (>= 2.0.0), rlang.

## Test environments

* Ubuntu 22.04, R 4.3.3 (R CMD check --as-cran, examples with --run-donttest)
* local macOS, R 4.5.x
* win-builder (R-devel)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Downstream dependencies

There are currently no downstream dependencies for this package.
