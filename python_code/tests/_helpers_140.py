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
