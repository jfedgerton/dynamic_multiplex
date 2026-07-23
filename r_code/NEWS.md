# dynamicmultiplex 1.1.0

Uncertainty-quantification overhaul, motivated by a large-scale simulation
study of empirical CI coverage (planted-partition multilayer networks;
n = 50-400 nodes, 3-10 communities, 5-15 layers, 5 switching rates,
3 density regimes, Louvain and Leiden, all similarity couplings).

## Breaking changes

* `community_ci()` is renamed `community_est()` and no longer returns a
  community-count confidence interval. The percentile `community_count_ci`
  (estimate/lower/upper) is replaced by a `community_count` data frame that
  reports, per layer, the observed-network `estimate` and its bootstrap
  `reproducibility` (the share of replicates whose community count equals
  the estimate), plus a plain-language `report`. Reason: the interval
  covered the truth at ~0.99 on well-specified planted-partition networks
  but collapsed to ~0.62 under strong community-size skew, and no
  observable diagnostic separated the reliable cases from the rest. The
  reproducibility summary is descriptive and makes no coverage claim; the
  package's validated interval is `co_assignment_ci()`. The raw per-replicate
  community counts are no longer exposed: `bootstrap_multilayer()` now returns
  `community_count_reproducibility` (the per-layer share of replicates
  matching the observed count) in place of the old `community_count_samples`
  vector.

* `community_ci()` no longer returns `modularity_ci`. The study showed its
  empirical coverage is never close to the nominal 95% at any network size
  (~0.40 at n = 50, 0.00 at n = 100, vacuously 1.00 at n = 200): community
  detection maximizes modularity, so the bootstrap interval concentrates
  around an optimized, upwardly biased estimate. An interval that is
  anti-conservative when detection is imperfect and uninformative when it
  is exact is not a confidence interval; it has been removed rather than
  documented around. Raw bootstrap draws remain available in
  `bootstrap_multilayer()$modularity_samples` for descriptive use.

## New features

* `bootstrap_multilayer()` now uses a parametric network bootstrap:
  within/between-community edge probabilities and edge-weight pools are
  estimated from the observed network using the point-estimate partition,
  and each replicate redraws the full edge set. The previous scheme
  (Bayesian bootstrap on edge weights with fixed topology) was REMOVED
  rather than kept as an option: because it never varies which edges
  exist, it understates the variability of fresh data, and simulation
  studies showed intervals built from it undercover badly (~45-48% at
  nominal 95% for pairwise co-assignment) with no recovery as networks
  grow. Under the project's design principle, a misleading option is
  worse than no option.

* `co_assignment_ci()`: Wilson score intervals for node-pair co-assignment
  propensities, computed from the bootstrap co-assignment probabilities.
  Label-invariant, so it avoids the label-switching problem that makes
  per-node membership intervals ill-defined.

## Improved warnings

* `community_est()` no longer issues the small-network community-count CI
  warning. With the count interval removed there is no coverage claim to
  qualify; the reproducibility summary is descriptive at any network size.
  The `co_assignment_ci()` guidance to interpret cautiously below 100 nodes
  remains.

* `directed = TRUE` with `algorithm = "leiden"` previously errored
  (igraph's `cluster_leiden()` supports undirected graphs only). Directed
  layers are now collapsed to weighted undirected graphs — matching the
  existing Louvain behavior — with a single warning per fit pointing users
  who need directed-aware detection to the Python package (leidenalg).

# dynamicmultiplex 1.0.0

Initial CRAN release.

* Multilayer community detection via Louvain or Leiden on evolving networks,
  with interlayer coupling built from Jaccard similarity, overlap coefficient,
  node-strength weighted variants, or direct node identity links
  (`fit_multilayer_jaccard()`, `fit_multilayer_overlap()`,
  `fit_multilayer_weighted_jaccard()`, `fit_multilayer_weighted_overlap()`,
  `fit_multilayer_identity_ties()`).
* User-specified layer connectivity via `layer_links`, enabling adjacent-only
  temporal coupling.
* Explicit objective selection (CPM or modularity) via the `objective`
  argument; negative edge weights are rejected for modularity-based methods.
* Bayesian bootstrap confidence intervals for community assignments
  (`bootstrap_multilayer()`, `community_ci()`).
* Simulation utilities for benchmarking (`simulate_and_fit_multilayer()`).
* Visualization helpers: layer series plots (`plot_multilayer_series()`),
  alluvial diagrams (`plot_multilayer_alluvial()`), and animated GIFs
  (`animate_multilayer_gif()`).
