# Archived: node co-membership confidence intervals (package 1.2.x, Sept 2026)

These scripts produced the coverage study that an earlier draft reported as a
95% node co-membership interval with a "reliability region". They are kept
for the record and are not part of the replication pipeline.

Why they were retired (Sept 20, 2026):

* The Wilson interval's coverage of the fresh-data co-assignment propensity
  was near nominal only for pairs whose bootstrap share is near 0 or 1 and
  fell to 0.03-0.15 for ambiguous pairs (share 0.3-0.8). The width gate
  selected the trivially decided pairs (99.96% of gated pairs at share < 0.1
  or > 0.9), which is why the pooled gated coverage looked like 0.96.
* A calibrated interval (conditional quantiles of the propensity given the
  share; `14_calibration_table.R`) was honest but uninformative mid-range
  (width 0.8-0.9).
* `05_alt_bootstrap.R` + `15_alternatives.R` tested two other bootstraps
  (degree-corrected, edge rewiring) and ten conditioning features under a
  pre-registered criterion (mid-range width <= 0.3 at coverage >= 0.90); the
  best achieved width 0.567. `diag_calibration_conditional.R` is the
  conditional-coverage diagnostic that surfaced the problem.
* `06_stability.R` + `16_stability.R` are the first version of the test that
  replaced the intervals (partition stability -> calibrated accuracy floor);
  the production version is `replication/sim/02_stability.R` and
  `replication/post/12_stability.R`.

Also note: the `identity` specification in `02_coverage_grid.R` and
`04_coverage_misspec.R` ran the pre-1.2.1 multislice implementation (plain
modularity on the stacked graph), so those rows are not meaningful.
