"""
Helpers for live hover tooltips: top expressed genes per cell,
and expression lookup for a user-selected gene.
"""

import numpy as np
import anndata as ad


def compute_top_genes(adata: ad.AnnData, n_top: int = 3) -> list[str]:
    """
    For each cell, return a formatted string of its top N expressed genes.
    Vectorized with numpy (no per-cell Python loop) — both faster and
    avoids iterating sparse matrix rows directly, which raises
    "sparse array length is ambiguous" on some backed/sparse types.
    """
    X = adata.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X)

    gene_names = adata.var_names.to_numpy()
    n_cells, n_genes = X.shape
    n_top_eff = min(n_top, n_genes)

    top_idx = np.argpartition(X, -n_top_eff, axis=1)[:, -n_top_eff:]
    row_idx = np.arange(n_cells)[:, None]
    order = np.argsort(-X[row_idx, top_idx], axis=1)
    sorted_top_idx = np.take_along_axis(top_idx, order, axis=1)

    gene_matrix = gene_names[sorted_top_idx]
    return [", ".join(row) for row in gene_matrix]


def get_gene_expression(adata: ad.AnnData, gene: str) -> np.ndarray:
    """Return the expression vector (one value per cell) for a given gene name."""
    if gene not in adata.var_names:
        raise ValueError(f"Gene '{gene}' not found in dataset")
    col = adata[:, gene].X
    if hasattr(col, "toarray"):
        col = col.toarray()
    return np.asarray(col).flatten()


def get_gene_options(adata: ad.AnnData, limit: int = 500) -> list[str]:
    """
    Return a manageable list of gene names for a dropdown.
    Prefers highly_variable genes if available, else first `limit` genes.
    """
    if "highly_variable" in adata.var:
        genes = adata.var_names[adata.var["highly_variable"]].tolist()
        if genes:
            return sorted(genes)[:limit]
    return sorted(adata.var_names.tolist())[:limit]
