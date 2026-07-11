#!/usr/bin/env python
"""
01_qc_filter.py — COMPUTE (no plotting).
========================================
MAD-based adaptive QC computed PER GSM (so small samples aren't clipped by pooled
thresholds), species-aware MT/RP/HB annotation on the HGNC gene_symbol (var_names
are Ensembl), scrublet doublet flagging per GSM, then a FIRST unsupervised
embedding (normalize->log1p->HVG->PCA->neighbors->UMAP->leiden) stored as
`X_umap_unsupervised` for usability review ONLY (biology is never read off it).

Outputs:
  objects/01_qc.h5ad          (filtered; lognorm X + counts layer + X_umap_unsupervised)
  03_results/01_qc/tables/qc_metrics_per_cell.csv
  03_results/01_qc/tables/qc_thresholds_per_gsm.csv
  03_results/01_qc/tables/cells_kept_dropped.csv
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "02_analysis"))
os.chdir(ROOT)

from config import PATHS, PARAMS, SPECIES_DB, CONFIG  # noqa: E402
from helpers.anndata_utils import annotate_gene_classes  # noqa: E402

STAGE = "01_qc"
QC_METRICS = ["log1p_total_counts", "log1p_n_genes_by_counts"]


def mad_outlier(series: pd.Series, n_mads: float, upper_only: bool = False) -> pd.Series:
    """Boolean outlier mask: value beyond median +/- n_mads*MAD (median abs deviation)."""
    med = series.median()
    mad = (series - med).abs().median()
    if mad == 0:
        mad = 1e-9
    hi = series > med + n_mads * mad
    lo = series < med - n_mads * mad
    return hi if upper_only else (hi | lo)


def main() -> None:
    adata = sc.read_h5ad(PATHS.object("00_build"))
    adata.X = adata.layers["counts"].copy()  # QC on raw counts
    annotate_gene_classes(adata, SPECIES_DB)

    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], percent_top=None,
                               log1p=True, inplace=True)

    n_mads = float(PARAMS.qc_n_mads)
    pct_mt_max = float(PARAMS.qc_pct_mt_max)
    min_genes = int(PARAMS.qc_min_genes)
    # Mito-specific MAD multiplier: None -> ceiling-only (retain stressed high-mito cells).
    _n_mads_mt = PARAMS.get("qc_n_mads_mt", None)
    n_mads_mt = None if _n_mads_mt is None else float(_n_mads_mt)
    mt_policy = "ceiling_only" if n_mads_mt is None else f"mad{n_mads_mt:g}+ceiling"
    # GSMs the human hard-excludes at BREAKPOINT 01 (decisions.qc.drop_gsms).
    drop_gsms = set((CONFIG.get("decisions", {}).get("qc", {}) or {}).get("drop_gsms", []) or [])

    # --- per-GSM MAD thresholds ---
    thr_rows = []
    outlier = pd.Series(False, index=adata.obs_names)
    for gsm, idx in adata.obs.groupby("gsm", observed=True).groups.items():
        sub = adata.obs.loc[idx]
        gsm_out = pd.Series(False, index=idx)
        rec = {"gsm": gsm, "n_cells_in": len(idx)}
        for m in QC_METRICS:
            o = mad_outlier(sub[m], n_mads)
            gsm_out |= o
            med, mad = sub[m].median(), (sub[m] - sub[m].median()).abs().median()
            rec[f"{m}_lo"] = float(med - n_mads * mad)
            rec[f"{m}_hi"] = float(med + n_mads * mad)
        # %mt gate: per-GSM MAD-upper only if n_mads_mt set, else ceiling alone.
        if n_mads_mt is not None:
            o_mt = mad_outlier(sub["pct_counts_mt"], n_mads_mt, upper_only=True)
            mt_med = sub["pct_counts_mt"].median()
            mt_mad = (sub["pct_counts_mt"] - mt_med).abs().median()
            rec["pct_mt_hi"] = float(min(mt_med + n_mads_mt * mt_mad, pct_mt_max))
        else:
            o_mt = pd.Series(False, index=idx)
            rec["pct_mt_hi"] = pct_mt_max
        o_mt = o_mt | (sub["pct_counts_mt"] > pct_mt_max)
        gsm_out |= o_mt
        rec["pct_mt_ceiling"] = pct_mt_max
        rec["mt_policy"] = mt_policy
        rec["excluded_gsm"] = gsm in drop_gsms
        rec["n_flagged_outlier"] = int(gsm_out.sum())
        thr_rows.append(rec)
        outlier.loc[idx] = gsm_out.values

    adata.obs["mad_outlier"] = outlier
    adata.obs["low_genes"] = adata.obs["n_genes_by_counts"] < min_genes
    adata.obs["excluded_gsm"] = adata.obs["gsm"].isin(drop_gsms)
    if drop_gsms:
        print(f"[01_qc] hard-excluding {len(drop_gsms)} GSM(s) (decisions.qc.drop_gsms): "
              f"{sorted(drop_gsms)} = {int(adata.obs['excluded_gsm'].sum())} cells; mt_policy={mt_policy}")

    # --- doublets (scrublet per GSM, on raw counts) ---
    # NB: scanpy's sc.pp.scrublet(batch_key=...) is all-or-nothing — it runs scrublet per
    # batch in a list comprehension, so a SINGLE GSM with fewer cells than n_prin_comps
    # (default 30) makes its internal PCA raise ("n_components=30 must be between 1 and
    # min(n_samples,n_features)"), which kills doublet detection for EVERY batch (the whole
    # call throws → doublet_score all-NaN dataset-wide). We run it per-GSM instead so one
    # tiny sample can't zero out the rest; tiny batches are skipped (left NaN/False).
    edr = float(PARAMS.scrublet_expected_doublet_rate)
    n_prin_comps = 30
    adata.obs["doublet_score"] = np.nan
    adata.obs["predicted_doublet"] = False
    for gsm, idx in adata.obs.groupby("gsm", observed=True).groups.items():
        sub = adata[idx].copy()
        if sub.n_obs <= n_prin_comps + 1:  # too few cells for the scrublet PCA
            print(f"[01_qc] scrublet: skipping GSM {gsm} (n={sub.n_obs} <= n_prin_comps={n_prin_comps})")
            continue
        try:
            sc.pp.scrublet(sub, expected_doublet_rate=edr, n_prin_comps=n_prin_comps)
            adata.obs.loc[idx, "doublet_score"] = sub.obs["doublet_score"].to_numpy()
            adata.obs.loc[idx, "predicted_doublet"] = sub.obs["predicted_doublet"].to_numpy()
        except Exception as exc:  # per-batch guard; a bad batch no longer sinks the rest
            print(f"[01_qc] scrublet warning for GSM {gsm}: {exc}; leaving NaN/False for this batch")
    n_scored = int(adata.obs["doublet_score"].notna().sum())
    print(f"[01_qc] scrublet scored {n_scored}/{adata.n_obs} cells; "
          f"flagged {int(adata.obs['predicted_doublet'].sum())} doublets")
    adata.obs["predicted_doublet"] = adata.obs["predicted_doublet"].astype(bool)

    # --- pass/fail ---
    fail = (adata.obs["mad_outlier"] | adata.obs["low_genes"]
            | adata.obs["predicted_doublet"] | adata.obs["excluded_gsm"])
    adata.obs["pass_qc"] = ~fail

    # --- tables ---
    tdir = PATHS.tables(STAGE)
    per_cell = adata.obs[[
        "gsm", "donor", "tissue", "population", "population_short",
        "total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo",
        "doublet_score", "predicted_doublet", "mad_outlier", "low_genes", "excluded_gsm", "pass_qc",
    ]].copy()
    per_cell.to_csv(tdir / "qc_metrics_per_cell.csv")
    pd.DataFrame(thr_rows).to_csv(tdir / "qc_thresholds_per_gsm.csv", index=False)

    kept_dropped = (
        adata.obs.groupby(["gsm", "donor", "tissue", "population"], observed=True)["pass_qc"]
        .agg(n_cells="size", n_kept="sum").reset_index()
    )
    kept_dropped["n_dropped"] = kept_dropped["n_cells"] - kept_dropped["n_kept"]
    kept_dropped["frac_kept"] = kept_dropped["n_kept"] / kept_dropped["n_cells"]
    kept_dropped.sort_values(["population", "tissue", "donor"]).to_csv(
        tdir / "cells_kept_dropped.csv", index=False)
    print(f"[01_qc] kept {int(adata.obs['pass_qc'].sum())}/{adata.n_obs} cells "
          f"({adata.obs['pass_qc'].mean():.1%})")

    # --- subset to passing cells + gene filter ---
    adata = adata[adata.obs["pass_qc"]].copy()
    sc.pp.filter_genes(adata, min_cells=int(PARAMS.qc_min_cells_per_gene))

    # --- first UNSUPERVISED embedding (usability only) ---
    adata.layers["counts"] = adata.layers["counts"]  # keep raw
    adata.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata  # log-norm snapshot for scoring convenience
    sc.pp.highly_variable_genes(adata, n_top_genes=int(PARAMS.hvg_n_top))
    emb = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(emb, max_value=10)
    sc.tl.pca(emb, n_comps=int(PARAMS.n_pcs))
    sc.pp.neighbors(emb, n_pcs=int(PARAMS.n_pcs))
    sc.tl.umap(emb)
    sc.tl.leiden(emb, resolution=float(PARAMS.leiden_resolution), flavor="igraph",
                 n_iterations=2, directed=False)
    adata.obsm["X_umap_unsupervised"] = emb.obsm["X_umap"]
    adata.obs["leiden_unsupervised"] = emb.obs["leiden"].values

    out = PATHS.object(STAGE)
    adata.write_h5ad(out)
    print(f"[01_qc] wrote checkpoint {out}: {adata.n_obs} cells x {adata.n_vars} genes")


if __name__ == "__main__":
    main()
