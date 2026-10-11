"""Directed within-slice edges in the multislice (Mucha) fit.

Covers ``detect_multislice_communities`` / ``fit_multilayer_identity_ties``
only. Seed 123 throughout. The known-answer networks are hard coded.

* Reciprocal equivalence: a directed network in which every arc is
  reciprocated must give the same fit as its undirected twin at the same
  omega. Before 1.3.2 the slice graphs were always built undirected, so the
  two arcs became parallel edges of double weight, which is the same as
  fitting the undirected network at omega / 2.
* Directed null model: with ``directed=True`` each slice uses Leicht and
  Newman's k_out * k_in / m null model. A 2-slice, 4-node directed network
  is solved by brute force over all 4,140 partitions of the 8 supra-nodes.
* Regression: undirected inputs, and directed layers fitted with
  ``directed=False``, reproduce the 1.3.1 memberships exactly.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import networkx as nx
import numpy as np
import pytest

pytest.importorskip("leidenalg")
pytest.importorskip("igraph")

from dynamic_multiplex.fit_multilayer_identity_ties import fit_multilayer_identity_ties
from dynamic_multiplex.multilayer_utils import detect_multislice_communities

FIXTURES = Path(__file__).parent / "fixtures"
SEED = 123


# ---------------------------------------------------------------------------
# Hard-coded networks
# ---------------------------------------------------------------------------
def _block_edges(assign: list[int]) -> list[tuple[int, int]]:
    n = len(assign)
    edges = [(i, j) for i in range(n) for j in range(i + 1, n) if assign[i] == assign[j]]
    return edges + [(0, 11)]


def planted12(directed: bool) -> list[nx.Graph | nx.DiGraph]:
    """3 slices, two 6-node cliques plus one bridge; nodes 5 and 6 swap blocks in slice 3."""
    base = [0] * 6 + [1] * 6
    s3 = base.copy()
    s3[5], s3[6] = 1, 0
    layers = []
    for assign in (base, base, s3):
        edges = _block_edges(assign)
        if directed:
            g = nx.DiGraph()
            g.add_nodes_from(range(12))
            g.add_edges_from(edges, weight=1.0)
            g.add_edges_from([(v, u) for u, v in edges], weight=1.0)
        else:
            g = nx.Graph()
            g.add_nodes_from(range(12))
            g.add_edges_from(edges, weight=1.0)
        layers.append(g)
    return layers


# 2 slices x 4 nodes, one-way arcs. Under the directed null model the unique
# optimum at omega = 0.25, gamma = 1 splits {0, 3} from {1, 2} in both slices
# (Q = 46/15); symmetrising the arcs makes one all-node community optimal.
DIRECTED_ARCS = [
    [(1, 0), (1, 2), (1, 3), (2, 0), (2, 1), (3, 0)],
    [(0, 1), (0, 2), (0, 3), (2, 1), (3, 2)],
]
DIRECTED_OMEGA = 0.25
DIRECTED_BEST_Q = 46.0 / 15.0
DIRECTED_BEST_PARTITION = [1, 2, 2, 1, 1, 2, 2, 1]  # slice 1 nodes 1..4, then slice 2


def directed4() -> list[nx.DiGraph]:
    layers = []
    for arcs in DIRECTED_ARCS:
        g = nx.DiGraph()
        g.add_nodes_from(range(4))
        g.add_edges_from(arcs, weight=1.0)
        layers.append(g)
    return layers


def sbm_layers(drop_node_from_layer2: bool = False) -> list[nx.Graph]:
    """30 nodes, 4 layers, 3 blocks, weighted; fixture written once with seed 123."""
    layers = {t: nx.Graph() for t in range(1, 5)}
    for g in layers.values():
        g.add_nodes_from(range(1, 31))
    with open(FIXTURES / "multislice_sbm_edges.csv", newline="") as f:
        for row in csv.DictReader(f):
            layers[int(row["layer"])].add_edge(
                int(row["from"]), int(row["to"]), weight=float(row["weight"])
            )
    out = [layers[t] for t in range(1, 5)]
    if drop_node_from_layer2:
        out[1].remove_node(30)
    return out


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def canon(meta: list) -> list[int]:
    """Relabel a stacked membership by order of first appearance."""
    seen: dict = {}
    out = []
    for x in np.concatenate([np.asarray(m) for m in meta]):
        seen.setdefault(int(x), len(seen) + 1)
        out.append(seen[int(x)])
    return out


def all_partitions(n: int):
    """All set partitions of n items as restricted growth strings."""
    a = [0] * n

    def rec(i: int, mx: int):
        if i == n:
            yield tuple(a)
            return
        for v in range(mx + 2):
            a[i] = v
            yield from rec(i + 1, max(mx, v))

    yield from rec(1, 0)


def directed_multislice_q(memb, arcs_by_slice, n, omega, gamma=1.0) -> float:
    """Unnormalised multislice modularity with directed slices (ordered pairs)."""
    q = 0.0
    for s, arcs in enumerate(arcs_by_slice):
        a = np.zeros((n, n))
        for u, v in arcs:
            a[u, v] += 1.0
        k_out, k_in, m = a.sum(1), a.sum(0), a.sum()
        for i in range(n):
            for j in range(n):
                if memb[s * n + i] == memb[s * n + j]:
                    q += a[i, j] - gamma * k_out[i] * k_in[j] / m
    for i in range(n):
        if memb[i] == memb[n + i]:
            q += 2.0 * omega
    return q


# ---------------------------------------------------------------------------
# Regression cases (fixtures recorded from 1.3.1 / pkg-1.3.2 before the change)
# ---------------------------------------------------------------------------
def regression_cases() -> dict:
    sbm = sbm_layers()
    return {
        "planted12_omega1": dict(layers=planted12(False), omega=1.0),
        "planted12_omega6": dict(layers=planted12(False), omega=6.0),
        "sbm_omega1": dict(layers=sbm, omega=1.0),
        "sbm_omega05_gamma12": dict(layers=sbm, omega=0.5, resolution_parameter=1.2),
        "sbm_all_pairs": dict(
            layers=sbm,
            omega=1.0,
            layer_links=[
                {"from": a, "to": b, "weight": 1.0 / (b - a)}
                for a in range(1, 5)
                for b in range(a + 1, 5)
            ],
        ),
        "sbm_unequal_nodes": dict(
            layers=sbm_layers(drop_node_from_layer2=True), omega=1.0, allow_unequal_nodes=True
        ),
        # directed layers fitted with directed=False keep the symmetrised
        # (A + A^T) behaviour of 1.3.1 exactly
        "planted12_digraph_undirected_fit": dict(layers=planted12(True), omega=6.0),
    }


def run_case(case: dict) -> list[list[int]]:
    kw = {k: v for k, v in case.items() if k != "layers"}
    fit = fit_multilayer_identity_ties(case["layers"], seed=SEED, **kw)
    return [[int(x) for x in m] for m in fit["meta_communities"]]


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
def test_leidenalg_quality_is_ordered_pair_and_direction_aware():
    """Basis of the fix: leidenalg sums qualities over ordered pairs.

    A directed graph with every arc reciprocated has the same RBConfiguration
    and CPM(0) quality as the undirected graph, so no omega rescaling is
    needed when slices are built as directed graphs.
    """
    import igraph as ig
    import leidenalg

    edges = [(0, 1), (1, 2), (0, 2), (3, 4), (4, 5), (3, 5), (2, 3)]
    gu = ig.Graph(n=6, edges=edges, directed=False)
    gd = ig.Graph(n=6, edges=edges + [(v, u) for u, v in edges], directed=True)
    gu.es["weight"] = [1.0] * gu.ecount()
    gd.es["weight"] = [1.0] * gd.ecount()
    memb = [0, 0, 0, 1, 1, 1]
    qu = leidenalg.RBConfigurationVertexPartition(gu, initial_membership=memb, weights="weight")
    qd = leidenalg.RBConfigurationVertexPartition(gd, initial_membership=memb, weights="weight")
    assert qu.quality() == pytest.approx(5.0)  # hand-computed, ordered pairs
    assert qd.quality() == pytest.approx(qu.quality())
    # one-way arcs: k_out * k_in / m null model
    g1 = ig.Graph(n=6, edges=edges, directed=True)
    g1.es["weight"] = [1.0] * g1.ecount()
    q1 = leidenalg.RBConfigurationVertexPartition(g1, initial_membership=memb, weights="weight")
    assert q1.quality() == pytest.approx(18.0 / 7.0)


@pytest.mark.parametrize("omega", [0.5, 1.0, 3.0, 5.0, 6.0, 8.0, 10.0])
def test_reciprocal_directed_equals_undirected(omega):
    und = fit_multilayer_identity_ties(planted12(False), omega=omega, seed=SEED)
    dig = fit_multilayer_identity_ties(planted12(True), omega=omega, directed=True, seed=SEED)
    assert canon(dig["meta_communities"]) == canon(und["meta_communities"])


def test_reciprocal_directed_no_longer_halves_omega():
    """At omega = 6 the undirected fit keeps nodes 6 and 7 with their slice-1
    blocks in slice 3; the old symmetrised directed fit behaved like omega = 3
    and moved them."""
    und6 = canon(fit_multilayer_identity_ties(planted12(False), omega=6.0, seed=SEED)["meta_communities"])
    und3 = canon(fit_multilayer_identity_ties(planted12(False), omega=3.0, seed=SEED)["meta_communities"])
    assert und6 != und3  # the network separates the two omegas
    dig6 = canon(
        fit_multilayer_identity_ties(planted12(True), omega=6.0, directed=True, seed=SEED)[
            "meta_communities"
        ]
    )
    assert dig6 == und6


def test_directed_null_model_matches_brute_force():
    parts = list(all_partitions(8))
    assert len(parts) == 4140
    qs = np.array([directed_multislice_q(p, DIRECTED_ARCS, 4, DIRECTED_OMEGA) for p in parts])
    best = qs.max()
    assert best == pytest.approx(DIRECTED_BEST_Q)
    assert int(np.isclose(qs, best).sum()) == 1  # unique optimum
    fit = fit_multilayer_identity_ties(directed4(), omega=DIRECTED_OMEGA, directed=True, seed=SEED)
    got = canon(fit["meta_communities"])
    assert got == DIRECTED_BEST_PARTITION
    q_got = directed_multislice_q([g - 1 for g in got], DIRECTED_ARCS, 4, DIRECTED_OMEGA)
    assert q_got == pytest.approx(DIRECTED_BEST_Q)


def test_symmetrised_fit_of_directed4_differs():
    """With directed=False the same arcs are symmetrised (documented
    behaviour) and the optimum is a single community: direction matters."""
    with pytest.warns(UserWarning, match="symmetrised"):
        fit = fit_multilayer_identity_ties(directed4(), omega=DIRECTED_OMEGA, directed=False, seed=SEED)
    assert canon(fit["meta_communities"]) == [1] * 8


def test_mixed_directedness_raises():
    layers = planted12(True)
    layers[1] = nx.Graph(layers[1])
    with pytest.raises(ValueError, match="requires every layer to be a directed"):
        detect_multislice_communities(
            graph_layers=layers,
            interlayer_ties=None,
            omega=1.0,
            directed=True,
        )


def test_undirected_layers_with_directed_true_raise():
    with pytest.raises(ValueError, match="requires every layer to be a directed"):
        fit_multilayer_identity_ties(planted12(False), omega=1.0, directed=True, seed=SEED)


def test_directed_slice_without_edges_runs():
    layers = directed4()
    layers.append(nx.DiGraph())
    layers[2].add_nodes_from(range(4))  # slice 3 has nodes but no arcs (m_s = 0)
    fit = fit_multilayer_identity_ties(layers, omega=DIRECTED_OMEGA, directed=True, seed=SEED)
    meta = fit["meta_communities"]
    assert len(meta) == 3
    assert all(np.all(np.isfinite(np.asarray(m, dtype=float))) for m in meta)


def test_directed_fit_is_seed_reproducible():
    a = fit_multilayer_identity_ties(planted12(True), omega=6.0, directed=True, seed=SEED)
    b = fit_multilayer_identity_ties(planted12(True), omega=6.0, directed=True, seed=SEED)
    assert [list(m) for m in a["meta_communities"]] == [list(m) for m in b["meta_communities"]]


def test_fit_object_fields_unchanged():
    fit = fit_multilayer_identity_ties(directed4(), omega=DIRECTED_OMEGA, directed=True, seed=SEED)
    for key in ("algorithm", "layer_communities", "layer_links", "interlayer_ties", "directed"):
        assert key in fit
    assert fit["directed"] is True


@pytest.mark.filterwarnings("ignore:Directed layers were passed with directed=False")
def test_regression_undirected_and_symmetrised_paths_unchanged():
    expected = json.loads((FIXTURES / "multislice_regression_py.json").read_text())
    cases = regression_cases()
    assert set(expected) == set(cases)
    for name, case in cases.items():
        assert run_case(case) == expected[name], name
