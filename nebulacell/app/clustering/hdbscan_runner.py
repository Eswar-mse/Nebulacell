"""
Runs HDBSCAN on spatial coordinates to identify tissue domains,
independent of expression-based (e.g. Leiden) clustering.
"""

import numpy as np
import hdbscan
import anndata as ad


def suggest_params(n_cells: int) -> tuple[int, int]:
    """
    Scale default HDBSCAN params to dataset size.

    Visium's hex-spot grid has natural density gaps (fixed spot spacing,
    gaps between rows) that a small min_cluster_size can pick up on even
    at ~2-3k cells. Imaging-based assays (MERFISH, Xenium) are dense,
    roughly uniform continuous scans with no such gaps — a fixed small
    min_cluster_size just finds "the whole tissue" as one blob at scale.
    There's no single constant that's right for both; this scales the
    default with cell count as a starting point, not a guaranteed fix —
    the UI lets you override and re-run cheaply from there.
    """
    min_cluster_size = int(np.clip(n_cells * 0.002, 10, 300))
    min_samples = max(3, min_cluster_size // 4)
    return min_cluster_size, min_samples


def run_spatial_clustering(
    adata: ad.AnnData,
    min_cluster_size: int | None = None,
    min_samples: int | None = None,
) -> ad.AnnData:
    """
    Cluster cells by spatial coordinates only (adata.obsm['spatial']).
    Writes results to adata.obs['spatial_domain'] (cluster id, -1 = noise)
    and adata.obs['spatial_domain_prob'] (membership confidence).

    If min_cluster_size / min_samples aren't given, they're auto-scaled
    to dataset size via suggest_params().

    Returns the same AnnData object (mutated in place) for convenience.
    """
    if "spatial" not in adata.obsm:
        raise ValueError("AnnData object has no 'spatial' key in .obsm")

    coords = adata.obsm["spatial"]

    if min_cluster_size is None or min_samples is None:
        auto_size, auto_samples = suggest_params(len(coords))
        min_cluster_size = min_cluster_size or auto_size
        min_samples = min_samples or auto_samples

    clusterer = hdbscan.HDBSCAN(
        min_cluster_size=min_cluster_size,
        min_samples=min_samples,
    )
    labels = clusterer.fit_predict(coords)

    adata.obs["spatial_domain"] = labels.astype(str)
    adata.obs["spatial_domain"] = adata.obs["spatial_domain"].astype("category")
    adata.obs["spatial_domain_prob"] = clusterer.probabilities_

    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise = int(np.sum(labels == -1))

    print(
        f"[hdbscan] min_cluster_size={min_cluster_size}, min_samples={min_samples} -> "
        f"{n_clusters} spatial domains, {n_noise} noise points out of {len(labels)} cells"
    )

    return adata
