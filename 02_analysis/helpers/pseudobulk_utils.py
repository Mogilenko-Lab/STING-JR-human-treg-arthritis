"""
pseudobulk_utils.py — donor-level pseudobulk aggregation + integer-count DE.
===========================================================================
The biology is done here, NOT on any embedding. Raw integer UMI
counts are summed to donor x tissue x frozen-label pseudobulk, then SF-vs-PB DE
is run per population with pydeseq2 (paired model `~ donor + tissue`).

Emits sign-preserving ranked lists (by the Wald stat) for pre-ranked fgsea.
Used by 03_pseudobulk_de.py.
"""
from __future__ import annotations

from typing import Dict, List, Tuple

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


def run_de_sf_vs_pb(counts: pd.DataFrame, coldata: pd.DataFrame,
                    tissue_key: str, donor_key: str,
                    num: str, den: str,
                    min_cells_per_gene: int = 10,
                    gene_symbols: Dict[str, str] | None = None) -> pd.DataFrame:
    """Paired SF-vs-PB DE for ONE population's strata via pydeseq2.

    Model is chosen from the strata present: `~ donor + tissue` when both tissue
    arms share >=2 donors (paired); falls back to `~ tissue` otherwise. A donor
    seen in only one tissue (e.g. p3) is absorbed by its donor coefficient and
    contributes nothing to the tissue contrast (effectively unpaired for it).

    Returns a DE table (index = ensembl_id) with columns:
      gene_symbol, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, model.
    """
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    cd = coldata.copy()
    cd = cd[cd[tissue_key].isin([num, den])]
    c = counts.loc[cd.index]

    # Gene filter: keep genes detected in >= min_cells_per_gene pseudobulk samples.
    detected = (c > 0).sum(axis=0)
    c = c.loc[:, detected >= min(min_cells_per_gene, c.shape[0])]

    # Decide model: is the donor factor estimable alongside tissue?
    tab = cd.groupby(tissue_key)[donor_key].nunique()
    donors_per_arm = tab.reindex([num, den]).fillna(0)
    shared = set(cd.loc[cd[tissue_key] == num, donor_key]) & set(cd.loc[cd[tissue_key] == den, donor_key])
    paired = (donors_per_arm.min() >= 2) and (len(shared) >= 2)
    design = "~ donor + tissue" if paired else "~ tissue"
    # Normalise the formula factor names to the actual columns.
    design = design.replace("donor", donor_key).replace("tissue", tissue_key)

    cd[tissue_key] = pd.Categorical(cd[tissue_key], categories=[den, num])  # ref = denominator
    cd[donor_key] = cd[donor_key].astype("category")

    dds = DeseqDataSet(counts=c, metadata=cd, design=design, quiet=True)
    dds.deseq2()
    stats = DeseqStats(dds, contrast=[tissue_key, num, den], quiet=True)
    stats.summary()
    res = stats.results_df.copy()
    res.index.name = "ensembl_id"
    if gene_symbols is not None:
        res["gene_symbol"] = [gene_symbols.get(g, g) for g in res.index]
    res["model"] = design
    res["n_paired_donors"] = len(shared)
    return res


def ranked_by_stat(de: pd.DataFrame, symbol_col: str = "gene_symbol",
                   stat_col: str = "stat") -> pd.DataFrame:
    """Sign-preserving ranked list keyed by human symbol (max |stat| on collision).

    Returns a 2-col DataFrame (symbol, stat) sorted descending — the .rnk for fgsea.
    """
    d = de.dropna(subset=[stat_col, symbol_col]).copy()
    d[symbol_col] = d[symbol_col].astype(str)
    d = d[d[symbol_col] != "nan"]
    d["abs_stat"] = d[stat_col].abs()
    d = d.sort_values("abs_stat", ascending=False).drop_duplicates(symbol_col)
    out = d[[symbol_col, stat_col]].rename(columns={symbol_col: "symbol", stat_col: "stat"})
    return out.sort_values("stat", ascending=False).reset_index(drop=True)
