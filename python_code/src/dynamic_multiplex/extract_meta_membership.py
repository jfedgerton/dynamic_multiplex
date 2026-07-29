from __future__ import annotations


def extract_meta_membership(fit: dict):
    """Extract cross-layer meta-community membership.

    Returns the tracked (cross-layer) community assignment produced by the
    second-stage detection: for each layer, an array giving every node's
    meta-community. Unlike the per-layer ``layer_communities`` (detected
    independently, ignoring the coupling), the meta-communities are the
    partition that reflects the interlayer ties and any custom
    ``layer_links``. This is the membership that ``bootstrap_multilayer``
    validates.

    Parameters
    ----------
    fit : dict
        A fit object from one of the ``fit_multilayer_*`` functions.

    Returns
    -------
    list[np.ndarray]
        One array per layer giving each node's meta-community assignment
        (node order).

    Raises
    ------
    ValueError
        If ``fit`` has no ``meta_communities`` (refit with a current
        ``fit_multilayer_*`` function).
    """
    meta = fit.get("meta_communities")
    if meta is None:
        raise ValueError(
            "`fit` has no meta_communities; refit with a current "
            "fit_multilayer_* function."
        )
    return meta
