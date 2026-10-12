"""1.4.0: bootstrapping is kept only for the stability score s (manuscript decision tree).
Interval and floor outputs are removed. s must reproduce 1.3.2 exactly. Seed 123."""
import importlib.resources
import importlib.util
import inspect
import json

import networkx as nx
import numpy as np
import pytest

import dynamic_multiplex
from _helpers_140 import FIX, sbm_layers
from dynamic_multiplex import bootstrap_multilayer, partition_stability

pytestmark = pytest.mark.skipif(
    importlib.util.find_spec("leidenalg") is None or importlib.util.find_spec("community") is None,
    reason="needs leidenalg and python-louvain",
)


def test_interval_and_floor_functions_are_removed():
    for name in ("co_assignment_ci", "community_est"):
        assert not hasattr(dynamic_multiplex, name), name
        assert name not in dynamic_multiplex.__all__, name
    import dynamic_multiplex.bootstrap_multilayer as mod
    for name in ("co_assignment_ci", "community_est", "_load_stability_table", "_floor_lookup"):
        assert not hasattr(mod, name), name
    assert "calibration_table" not in inspect.signature(partition_stability).parameters
    data = importlib.resources.files("dynamic_multiplex") / "data" / "stability_calibration_table.csv"
    assert not data.is_file()


def test_partition_stability_returns_s_only():
    mats = [nx.to_numpy_array(g, nodelist=range(1, 31)) for g in sbm_layers()]
    p = partition_stability(bootstrap_multilayer(mats, n_boot=3, seed=123))
    assert set(p) == {"stability", "stability_by_layer", "stability_mc_se", "metric", "report"}
    assert "floor" not in p["report"] and "calibrat" not in p["report"]


@pytest.mark.parametrize("key", ["jaccard|leiden", "jaccard|louvain", "identity|leiden"])
def test_stability_score_matches_132(key):
    rec = json.load(open(FIX / "other_regression_py.json"))
    assert rec["recorded_from"] == "b55c66c"
    ft, alg = key.split("|")
    mats = [nx.to_numpy_array(g, nodelist=range(1, 31)) for g in sbm_layers()]
    p = partition_stability(bootstrap_multilayer(mats, fit_type=ft, algorithm=alg, n_boot=8, seed=123))
    exp = rec["stability"][key]
    assert p["stability"] == exp["stability"]
    assert list(map(float, p["stability_by_layer"])) == exp["stability_by_layer"]
    assert p["stability_mc_se"] == exp["stability_mc_se"]
