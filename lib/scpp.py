# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.3'
#       jupytext_version: 0.8.3
# ---

"""
module containing helper functions for proprocessing
scingle cell datasets
"""
import scanpy as sc
import numpy as np
import warnings

def norm_log(adata):
    """
    Normalize and log-transform adata object.

    Store the raw (input) adata in `adata.uns["raw_counts"]`
    and per-cell normalized, log-transformed data
    in `adata.raw` so that it can be used for
    visualization.

    Args:
        adata: the adata object, containing raw, unfiltered
            non-log-transformed (UMI) counts.

    Returns:
        Nothing, but modifies adata inplace.

    Note:
        Stores a flag in the adata object that it has
        already been normalized. If this flag is detected
        this function doesn't do anything but raises a
        Warning instead.
    """
    if "norm_log" in adata.uns and adata.uns["norm_log"]:
        warnings.warn("This adata object has already been processed"
                      " using the `norm_log` function. Will do"
                      " nothing this time. ")
    else:
        adata.uns["norm_log"] = True
        adata.uns["raw_counts"] = adata.X.copy()

        sc.pp.filter_genes(adata, min_cells=1)
        sc.pp.normalize_per_cell(adata, counts_per_cell_after=1000)
        sc.pp.log1p(adata)

        adata.raw = adata.copy()

