# dynamicmultiplex 1.4.0

Includes every change listed under 1.3.2 below, which was not released
separately. The same changes are in the Python package `dynamic-multiplex`
1.4.0.

## Breaking changes

* The bootstrap is kept only for the stability score s used to compare two
  fits of the same network (for example DynMux against multislice).
  `co_assignment_ci()` and `community_est()` are removed, and
  `partition_stability()` now returns only `stability`,
  `stability_by_layer`, `stability_mc_se`, `metric` and `report`: the
  calibrated accuracy floor, node floors, decided-pair outputs, the
  `decided` and `calibration_table` arguments and the bundled calibration
  table are removed. Their coverage was not reliable.

* R fit functions now default to `seed = 123`, as in the Python package
  (previously `NULL`, the caller's random number stream). Pass `seed = NULL`
  for unseeded detection. Because bootstrap refits are now seeded, R values of
  s differ slightly from 1.3.x (for example 0.916 to 0.921 on a test network,
  within Monte Carlo error).

* `allow_unequal_nodes = TRUE` now requires vertex names on every layer even
  when layers have the same size, since nodes can enter and leave without the
  size changing. Previously same-size unnamed layers were matched by
  position. (Python: adjacency matrices are rejected with
  `allow_unequal_nodes=True`.)

## New features

* `repair_chains` and `repair_min_similarity` (two-stage fit functions,
  default off). Stage-2 modularity tends to cut a long-lived community into
  time segments (a resolution-limit effect on chains of communities). With
  `repair_chains = TRUE` a group that ends at layer s and one that starts at a
  linked layer t are joined when their communities are mutual best matches
  (Kawadia and Sreenivasan 2012) with similarity at least
  `repair_min_similarity` (default 0.2; Greene et al. 2010 use a matching
  threshold). Groups that coexist are never merged. Joins are logged in
  `chain_repairs`.

* `tracking = "filtered"` (two-stage fit functions; default `"smoothed"`, the
  previous behaviour). Stage 2 runs forward in time: the meta ids of layers
  1..t depend only on layers 1..t, so fitted labels carry no information from
  later layers. Works with Leiden and Louvain and with `repair_chains`.
  Needs a seed (`seed = NULL` is an error), because unseeded stage-1 fits of
  later layers would shift the random stream that stage 2 uses.

* `layers_from_membership()` builds one weighted layer per time slice from
  group memberships (wide 0/1 matrices or long node-group tables), with
  `weighting = "count"`, `"size"` or `"kumar"` (Kumar et al. hypergraph
  reweighting, as in hypernetx 2.4.3). There is no default weighting, each
  slice is built on its own, and `rescale = "mean"` puts slices on a common
  scale for the multislice coupling.

* `bootstrap_multilayer()` gains `allow_unequal_nodes` (nodes entering and
  leaving) and passes `tracking`, `repair_chains` and
  `repair_min_similarity` to the refits. Layers with different named node
  sets now stop with an error unless `allow_unequal_nodes = TRUE`.

## Python only (by design)

* The two packages deliberately diverge here, because the libraries these
  rely on exist in Python and not in R. `partition_membership_bipartite()`
  partitions membership slices by Barber's bipartite modularity (leidenalg's
  bipartite construction, Leiden only); `partition_membership_strict()` by
  strict (all-or-nothing) hypergraph modularity (hypernetx reference
  pipeline, optional extra `[hypergraph]`, which needs pandas < 3);
  `fit_multilayer_from_partitions()` tracks supplied per-slice partitions
  with the DynMux second stage or Hungarian matching. They are comparison
  specifications; in pilot simulations both did worse than the projections
  in `layers_from_membership()`.

## Bug fixes

* `fit_multilayer_identity_ties()` (undirected): the multislice local mover
  scored staying in a community the vertex has no tie into as 0, which could
  refuse an improving move (the same flaw was fixed for directed slices in
  1.3.2). The stay score is now exact. Results on the package's regression
  networks are unchanged.

* Python: `bootstrap_multilayer()` silently paired nodes by position when
  networkx layers had the same size but different nodes; it now stops unless
  `allow_unequal_nodes=True`, which matches nodes by label.

# dynamicmultiplex 1.3.2

Not released separately; included in 1.4.0.

## Changes in defaults

* All fit functions, `bootstrap_multilayer()` and
  `simulate_and_fit_multilayer()` now default to `algorithm = "leiden"`
  (previously `"louvain"`), matching the Python package. Pass
  `algorithm = "louvain"` to reproduce earlier default results.

* Input policy (shared with the Python package): within-layer self-loops and
  zero-weight edges are dropped for both adjacency-matrix and `igraph` input.
  Previously a loop-bearing `igraph` layer kept its loops while the
  equivalent matrix did not. Layers with no vertices now stop with a clear
  error.

## Bug fixes

* `fit_multilayer_identity_ties()` with `directed = TRUE`: the multislice stage
  built every slice as undirected, so each arc became an undirected edge and a
  reciprocated pair counted twice. A fully reciprocated directed network was
  therefore fitted as if `omega` were halved relative to its undirected twin,
  and one-way arcs were symmetrised. Each slice now keeps arc direction, with
  Leicht and Newman's directed null model k_out k_in / m_s; the interlayer
  identity ties stay undirected. With `directed = TRUE` every layer must be
  directed (error otherwise). Directed layers reaching the multislice stage
  with `directed = FALSE` are symmetrised as before, now with a warning.
  Results for undirected networks are unchanged (regression-tested). Same fix
  in the Python package.

* `fit_*()` with directed `igraph` layers and `directed = FALSE` stopped with
  an igraph error in per-layer detection. Such layers are now collapsed to
  weighted undirected graphs (reciprocal arc weights summed), with a warning.

* Python package: the multislice fit could merge two different node labels
  that print the same (for example `1` and `"1"`) when layers have unequal
  node sets; supra-vertices are now keyed by (layer, node).

* `fit_multilayer_hungarian()` matched nodes by position, so two layers with
  the same named vertices stored in different orders produced spurious label
  changes. Matching now uses vertex names (positions only for unnamed
  layers), and the function gains `allow_unequal_nodes` for layers with
  different named node sets.

* `plot_multilayer_series()` and `animate_multilayer_gif()` dropped every edge
  of graphs with named vertices (coordinates were indexed by name). Both now
  color nodes by the tracked `meta_communities` rather than the layer-local
  first-stage labels, interpolate the palette when there are more than eight
  communities (previously an error), and the series plot uses free panel
  scales.

* `plot_multilayer_alluvial()` converted actor names to layer-local
  positions, so it could connect different actors across layers. Flows now
  follow actors by name, use the tracked `meta_communities`, and are drawn
  with `ggalluvial::geom_flow()` (no warning when an actor changes community).

* Documentation: `objective` sets only the first-stage (per-layer) objective;
  the second stage always maximizes modularity (multislice modularity for
  `fit_multilayer_identity_ties()`). README corrections (package name,
  identity-tie description, five simulator fit types).

* `fit_multilayer_identity_ties()` (generalized Louvain multislice): the
  restart-selection score omitted a factor of 2 on the within-community edge
  weight, so the best of the `n_starts` restarts was ranked as if the
  resolution were 2 rather than `resolution_parameter`. Each restart already
  optimised the correct objective; only the choice among restarts changes.

# dynamicmultiplex 1.3.1

This release collects every change since 1.1.0 (CRAN, 2026-08-07). Versions
1.2.1 and 1.3.0 were internal and never submitted.

## New features

* **`partition_stability()`**: a reliability report for a tracked partition
  built from a `bootstrap_multilayer()` result. It returns (i) the stability
  score `s`, the mean NMI (or ARI) between each bootstrap replicate's
  meta-partition and the point estimate; (ii) a calibrated accuracy floor,
  the 5th percentile of accuracy against the planted partition among
  simulated fits with the same stability, so that accuracy exceeded the floor
  in 95 percent of held-out fits; (iii) node-level Jaccard stability with its
  own floor (reported only above stability 0.9); and (iv) decided node pairs
  (co-assignment share below/above user thresholds). The calibration table
  ships in `inst/extdata/stability_calibration_table.csv` and can be
  overridden with `calibration_table`.

* **`bootstrap_multilayer()`** now records `stability_samples` (replicate x
  layer NMI and ARI to the point estimate) and `node_jaccard_stability`
  (per-node mean Jaccard overlap with the point-estimate community), which
  `partition_stability()` summarises.

* **`allow_unequal_nodes`** argument for every `fit_multilayer_*()` function.
  Layers may now contain different node sets (nodes entering and exiting the
  system); nodes are matched across layers by vertex name and only nodes
  present in both layers of a linked pair contribute to the interlayer
  similarity. Default `FALSE` keeps the previous requirement of a shared node
  universe.

## Bug fixes

* **`fit_multilayer_identity_ties()`** now optimises Mucha et al. (2010)
  multislice modularity with a generalized Louvain (local moving plus
  aggregation), using each slice's own configuration null model
  `k_is k_js / 2 m_s` and null-free interlayer identity ties. Previously the
  function ran single-graph modularity on the stacked supra-graph, whose
  pooled null model is roughly `T` times too small and penalises cross-slice
  co-membership; at `omega <= 1` that returned each slice as one community.
  `algorithm` is kept for API compatibility but the optimiser is always the
  generalized Louvain. The optimiser uses 10 random restarts (3 in the
  internal 1.3.0); on small change-point series joint NMI varied by up to 0.2
  across runs at 3.

* **`weighted_jaccard_similarity()`** and **`weighted_overlap_similarity()`**
  used layer-wide node strengths, so a node's strength counted toward a
  community it was not a member of. Two disjoint communities could score 1.0
  whenever their members had non-zero strength in both layers. A node's
  strength now counts toward community `a` only when the node is a member of
  `a` (and likewise for `b`).

* **`weighted_overlap_similarity()`** errored (`subscript out of bounds`)
  when a community contained a node with no recorded strength in that layer
  (e.g., an isolate). Missing strengths now contribute 0, matching
  `weighted_jaccard_similarity()` and the Python implementation. This
  affected `fit_multilayer_weighted_overlap()` on layers with isolates.

## Changes

* **`co_assignment_ci()`** is now documented as a descriptive Monte Carlo
  interval for the bootstrap co-assignment share, not a calibrated confidence
  interval. The simulation study found pair-level intervals are not
  calibratable for ambiguous pairs; use `partition_stability()` for a
  calibrated reliability measure of the tracked partition. Function signature
  and output are unchanged.

* Documentation regenerated with roxygen2 7.3.3; `partition_stability` help
  page added; cross-references synced.

* Regression tests added for the weighted similarity fixes, the multislice
  null model, and `partition_stability()`.

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
