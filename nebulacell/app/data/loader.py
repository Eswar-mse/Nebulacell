"""
AnnData loading utilities.

- Small datasets: load fully into memory.
- Large datasets (matrix size above MAX_ELEMENTS): load in backed mode ('r'),
  which keeps the expression matrix on disk and only reads slices on demand.

The backed-mode decision is based on actual matrix size (cells * genes),
not cell count alone. A dataset with many cells but a small targeted gene
panel (e.g. MERFISH, ~150 genes) is tiny in memory even at 70k+ cells, and
forcing it through the backed-mode path breaks downstream code that expects
a dense/in-memory array.
"""

import scanpy as sc
import squidpy as sq
import anndata as ad

MAX_ELEMENTS = 50_000_000  # cells * genes; the real backed-mode trigger


def load_demo_dataset():
    """Load Squidpy's built-in Visium fluorescence dataset for dev/testing."""
    adata = sq.datasets.visium_fluo_adata()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    return adata


def load_h5ad(path: str) -> ad.AnnData:
    """
    Load a real .h5ad file. Peeks at shape/metadata cheaply first to decide
    in-memory vs backed mode.
    """
    peek = sc.read_h5ad(path, backed="r")

    if "spatial" not in peek.obsm:
        peek.file.close()
        raise ValueError(
            f"'{path}' has no 'spatial' coordinates in .obsm — "
            "this doesn't look like a spatial transcriptomics dataset."
        )

    n_elements = peek.n_obs * peek.n_vars

    if n_elements <= MAX_ELEMENTS:
        peek.file.close()
        adata = sc.read_h5ad(path)  # full load — small enough regardless of cell count
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        return adata

    # Genuinely large dataset (many cells AND many genes): stay backed.
    # Caller is responsible for pulling only the slices/subsets it needs.
    return peek


def get_spatial_coords(adata: ad.AnnData):
    """Return the (x, y) spatial coordinates as a plain array."""
    if "spatial" not in adata.obsm:
        raise ValueError("AnnData object has no 'spatial' key in .obsm")
    return adata.obsm["spatial"]
