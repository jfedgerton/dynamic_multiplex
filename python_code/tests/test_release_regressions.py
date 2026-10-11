"""Regression guards added for the 1.3.2 release (pre-release review).

* Default ``fit_multilayer_jaccard`` output must equal 1.3.1 exactly. The
  fixture was recorded with the 1.3.1 code (commit 8e9c150) on a network
  (numpy default_rng(3)) where an extra Leiden iteration in the
  community-graph stage changes the partition.
* Multislice supra-vertices must keep node labels of different types apart:
  1 and "1" are different networkx nodes.
"""
from __future__ import annotations

import json
from pathlib import Path

import networkx as nx
import numpy as np
import pytest

pytest.importorskip("leidenalg")
pytest.importorskip("igraph")

from dynamic_multiplex.fit_multilayer_identity_ties import fit_multilayer_identity_ties
from dynamic_multiplex.fit_multilayer_jaccard import fit_multilayer_jaccard

FIXTURES = Path(__file__).parent / "fixtures"


def _seed3_layers() -> list[nx.Graph]:
    rng = np.random.default_rng(3)
    layers = []
    for _ in range(8):
        a = rng.integers(1, 13, 100)
        g = nx.Graph()
        g.add_nodes_from(range(100))
        g.add_edges_from((i, j) for i in range(100) for j in range(i + 1, 100) if a[i] == a[j])
        layers.append(g)
    return layers


def test_default_jaccard_matches_1_3_1_exactly():
    expected = json.loads((FIXTURES / "jaccard_default_regression_py.json").read_text())
    fit = fit_multilayer_jaccard(_seed3_layers())
    got = [[int(x) for x in m] for m in fit["meta_communities"]]
    assert got == expected
    assert len({x for m in got for x in m}) == 37


def test_multislice_keeps_int_and_str_labels_apart():
    layers = [nx.empty_graph([1, "1", "x"]), nx.empty_graph([1, "1", "x", "y"])]
    fit = fit_multilayer_identity_ties(layers, allow_unequal_nodes=True, omega=0)
    meta = [[int(x) for x in m] for m in fit["meta_communities"]]
    # no arcs and no coupling: every supra-vertex is its own community
    assert len({x for m in meta for x in m}) == 7
