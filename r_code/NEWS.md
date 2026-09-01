# dynamicmultiplex 1.1.1

## Bug fixes

* `weighted_overlap_similarity()` errored (`subscript out of bounds`) when a
  community contained a node with no recorded strength in that layer (e.g.,
  an isolate). Missing strengths now contribute 0, matching
  `weighted_jaccard_similarity()` and the Python implementation. This affected
  `fit_multilayer_weighted_overlap()` on layers with isolates.

# dynamicmultiplex 1.1.0

Uncertainty-quantification overhaul, motivated by a large-scale simulation
study of empirical CI coverage (planted-partition multilayer networks;
n = 50-400 nodes, 3-10 communities, 5-15 layers, 5 switching rates,
3 density regimes, Louvain and Leiden, all similarity couplings).

## New features

* **New method: `fit_multilayer_hungarian()`.** A two-stage
  snapshot-and-match tracker that detects communities independently per layer
  and then aligns their labels across consecutive layers with the Hungarian
  (optimal linear-sum assignment) algorithm on the community overlap matrix
  (`clue::solve_LSAP`; SciPy's `linear_sum_assignment` in the Python package).
  This is the standard "independent detection + optimal label matching"
  baseline, now exported as a first-class method with the same
  `multilayer_community_fit` interface (`layer_communities`,
  `meta_communities`, `extract_meta_membership()`). It complements the
  coupling-based fits: those jointly optimise on a coupled supra-graph, this
  one matches post hoc. `interlayer_ties` is `NULL` for this method.

* **Cross-layer community tracking (second-stage detection).** The
  `fit_multilayer_*()` functions now run a second community detection on the
  interlayer/community graph (the similarity ties plus self-loops) and return
  `meta_communities`: the tracked partition grouping per-layer communities into
  cross-layer meta-communities. This is what makes custom `layer_links` and the
  interlayer coupling actually affect the returned membership; previously the
  output equalled independent per-layer detection. Access it with the new
  exported `extract_meta_membership()`. The per-layer `layer_communities` are
  unchanged and still returned. For the Jaccard/overlap specifications the
  second stage runs on the community-similarity graph; for the **identity**
  specification it is **Mucha (2010) multislice modularity** — the layers are
  stacked into one supra-graph (intra-layer edges = each layer's own adjacency)
  with interlayer identity edges, and a single detection is run on the whole
  supra-graph. Custom `layer_links` set the coupling in both cases.

* **Multislice omega / resolution-limit controls.**
  `fit_multilayer_identity_ties()` gains an `omega` argument (Mucha's
  interlayer coupling strength) alongside `resolution_parameter`, both
  forwarded to the multislice supra-graph detection. `omega` scales the
  interlayer identity-edge weights: small values decouple toward independent
  per-layer detection, `omega = 1` recovers the planted partition on
  well-separated networks, and large values drive the known multislice
  degeneracy (temporal threads fragment or slices collapse). Together with
  `resolution_parameter` (Mucha's modularity resolution) these expose the
  omega-sensitivity and resolution-limit behavior of multislice modularity as
  tunable knobs. Mirrored in the Python package.

* **Uncertainty is now quantified on the meta-communities.**
  `bootstrap_multilayer()`, `co_assignment_ci()`, and `community_est()` operate
  on the cross-layer meta-communities rather than the independent per-layer
  partitions. Co-assignment now answers "do two nodes belong to the same
  *persistent* community," and the community count is the number of
  meta-communities per layer. (Coverage of these meta-community intervals must
  be re-validated in simulation; the earlier per-layer coverage numbers do not
  transfer.)

* **`seed` argument for reproducible detection.** All `fit_multilayer_*()`
  functions gain a `seed` argument. When supplied, the global RNG state is
  saved, the RNG is seeded for the duration of the call, and the previous
  state is restored on exit, so seeded detection never disturbs the caller's
  random number stream (e.g. `bootstrap_multilayer()` resampling). The default
  (`NULL`) preserves the previous behavior: detection inherits the caller's
  RNG stream, so existing scripts using `set.seed()` are unaffected. (In the
  Python package the default is `seed=123`, matching its previous internally
  fixed detection seed.)

* **Node-universe validation.** Layers must now share the same node set
  (compared on vertex names when present, on vertex indices otherwise);
  unequal universes error with instructions instead of silently misaligning
  interlayer ties. `fit_multilayer_identity_ties()` gains
  `allow_unequal_nodes = TRUE` to retain its documented support for nodes
  entering or exiting the system.

## Breaking changes

* Layers with different node sets are now rejected by all fit functions
  (previously assumed aligned without checking). For identity ties, opt back
  in with `allow_unequal_nodes = TRUE`.

* `animate_multilayer_gif()` no longer has a default `output_file`
  (previously `"multilayer_animation.gif"` in the working directory). An
  explicit path is now required, per CRAN policy that functions must not
  write to the user's filespace by default; use
  `tempfile(fileext = ".gif")` for a temporary location.

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
