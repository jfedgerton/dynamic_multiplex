"""Python counterpart of test_package_regressions.R (dynamic_multiplex >= 1.2.1).
Usage: python replication/tests/test_package_regressions.py
"""
import numpy as np
from sklearn.metrics import normalized_mutual_info_score as nmi
from dynamic_multiplex import (extract_meta_membership, fit_multilayer_identity_ties,
                               fit_multilayer_weighted_jaccard)
from dynamic_multiplex.multilayer_utils import weighted_jaccard_similarity

w = {i: 5.0 for i in range(1, 11)}
assert weighted_jaccard_similarity([1, 2, 3, 4, 5], [6, 7, 8, 9, 10], w, w) == 0.0
assert abs(weighted_jaccard_similarity([1, 2, 3, 4, 5], [4, 5, 6, 7, 8], w, w) - 0.25) < 1e-12
print("weighted Jaccard: ok")

rng = np.random.default_rng(123); n, K = 100, 4
mem = rng.integers(1, K + 1, n)
def build(m, p_in, p_out):
    P = np.where(m[:, None] == m[None, :], p_in, p_out); np.fill_diagonal(P, 0)
    A = np.triu((rng.random((n, n)) < P).astype(float), 1); return A + A.T
L = [build(mem, 0.3, 0.04) for _ in range(5)]
fit = fit_multilayer_identity_ties(L, algorithm="leiden", omega=1.0, seed=123)
m = extract_meta_membership(fit)
assert all(len(set(x)) == K for x in m), [len(set(x)) for x in m]
assert len(set(np.concatenate(m))) == K
assert all(nmi(x, mem) > 0.95 for x in m)
fit0 = fit_multilayer_identity_ties(L, algorithm="leiden", omega=0.0, seed=123)
assert len(set(np.concatenate(extract_meta_membership(fit0)))) == K * 5
print("multislice identity ties: ok")
fw = fit_multilayer_weighted_jaccard(L, algorithm="leiden", seed=123)
mw = extract_meta_membership(fw)
assert all(len(set(x)) == K for x in mw)
assert np.mean([nmi(x, mem) for x in mw]) > 0.85
print("weighted Jaccard fitter: ok\nALL PACKAGE REGRESSION TESTS PASSED")
