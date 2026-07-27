"""
pseudobulk_utils.py — donor-level pseudobulk aggregation. Counts only, no statistics.
=====================================================================================
The biology is done here, NOT on any embedding: raw integer UMI counts are summed to
donor x tissue x frozen-label pseudobulk and strata below the cell floor are dropped.

This module deliberately runs **no** differential expression. DE is the R side of the
seam — `02_analysis/scripts/03b_pseudobulk_de.R` (edgeR/limma-voom, paired model
`~ donor + tissue` where both arms share >= 2 donors) — which reads the CSVs written by
`03a_pseudobulk_export.py` and emits the sign-preserving ranked lists keyed by HGNC
symbol and ordered on the moderated `t`.

Used by 03a_pseudobulk_export.py.
"""
from __future__ import annotations

from typing import List, Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sp


def aggregate_pseudobulk(adata, groupby: List[str], layer: str = "counts"
                         ) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Sum raw counts over `groupby` strata.

    Returns (counts, coldata):
      counts   — samples (strata) x genes, integer DataFrame (index = stratum id).
      coldata  — one row per stratum with the `groupby` columns + `n_cells`.
    Stratum id = the groupby values joined by '|'.
    """
    obs = adata.obs
    X = adata.layers[layer]
    X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(np.asarray(X))

    keys = obs[groupby].astype(str).agg("|".join, axis=1)
    uniq = keys.unique()
    rows = []
    meta = []
    for k in uniq:
        idx = np.where((keys == k).to_numpy())[0]
        rows.append(np.asarray(X[idx].sum(axis=0)).ravel())
        vals = k.split("|")
        rec = dict(zip(groupby, vals))
        rec["n_cells"] = int(len(idx))
        rec["stratum"] = k
        meta.append(rec)

    counts = pd.DataFrame(np.vstack(rows), index=uniq, columns=adata.var_names.astype(str))
    counts = counts.round().astype(int)
    coldata = pd.DataFrame(meta).set_index("stratum")
    coldata = coldata.loc[counts.index]
    return counts, coldata


def filter_strata(counts: pd.DataFrame, coldata: pd.DataFrame, min_cells: int
                  ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Drop strata below `min_cells`. Returns (counts, coldata, dropped_table)."""
    keep = coldata["n_cells"] >= min_cells
    dropped = coldata.loc[~keep, ["n_cells"]].copy()
    dropped["reason"] = f"n_cells < {min_cells}"
    return counts.loc[keep], coldata.loc[keep], dropped.reset_index()
