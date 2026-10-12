"""Shared fixture loaders for the 1.4.0 tests."""
from pathlib import Path
import networkx as nx
import pandas as pd

FIX = Path(__file__).parent / "fixtures"


def load_layers(name, n_nodes):
    d = pd.read_csv(FIX / name)
    layers = []
    for t in sorted(d["layer"].unique()):
        g = nx.Graph()
        g.add_nodes_from(range(1, n_nodes + 1))
        e = d[d["layer"] == t]
        g.add_weighted_edges_from(zip(e["from"], e["to"], e["weight"]))
        layers.append(g)
    return layers


def sbm_layers():
    return load_layers("multislice_sbm_edges.csv", 30)


def chain_layers():
    return load_layers("chain_panel_edges.csv", 60)


def nmi(a, b):
    """Normalized mutual information, arithmetic normalization (sklearn's
    default), written out so the tests need no scikit-learn."""
    import numpy as np
    _, ia = np.unique(np.asarray(a), return_inverse=True)
    _, ib = np.unique(np.asarray(b), return_inverse=True)
    n = len(ia)
    c = np.zeros((ia.max() + 1, ib.max() + 1))
    np.add.at(c, (ia, ib), 1)
    p = c / n
    pa, pb = p.sum(1), p.sum(0)
    nz = p > 0
    mi = float((p[nz] * np.log(p[nz] / np.outer(pa, pb)[nz])).sum())
    ha = float(-(pa[pa > 0] * np.log(pa[pa > 0])).sum())
    hb = float(-(pb[pb > 0] * np.log(pb[pb > 0])).sum())
    if ha == 0 and hb == 0:
        return 1.0
    return mi / ((ha + hb) / 2)
