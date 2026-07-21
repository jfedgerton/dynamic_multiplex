# dynamicmultiplex 1.1.0

Uncertainty-quantification overhaul, motivated by a large-scale simulation
study of empirical CI coverage (planted-partition multilayer networks;
n = 50-400 nodes, 3-10 communities, 5-15 layers, 5 switching rates,
3 density regimes, Louvain and Leiden, all similarity couplings).

## Breaking changes

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

* `co_assignment_ci()`: Wilson score intervals for node-pair co-assignment
  propensities, computed from the bootstrap co-assignment probabilities.
  Label-invariant, so it avoids the label-switching problem that makes
  per-node membership intervals ill-defined.

## Improved warnings

* `community_ci()` now warns at runtime when layers have fewer than 100
  nodes: at n = 50 the nominal 95% `community_count_ci` covered the true
  community count in only ~40% of simulations. Coverage recovers to at or
  above nominal for n >= 100. Documentation gained a prominent Warning
  section with the study numbers.

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
