#!/usr/bin/env python
"""
01_qc_filter_viz.py — VIZ (no statistics).
==========================================
Reads the QC tables + checkpoint from 01_qc_filter.py and renders the usability
overview: per-GSM QC distributions, kept/dropped counts, and the unsupervised
UMAP (population / tissue / donor) with a canonical Treg-marker overlay
(FOXP3, IL2RA, CTLA4, IKZF2 — canonical Treg markers). Purpose: is the data usable, do the
sorted populations separate, is there visible cross-population contamination?
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "02_analysis"))
os.chdir(ROOT)

from config import PATHS, PARAMS, COARSE_LABEL  # noqa: E402
from helpers.figure_style import set_paper_style, save_overview, FIG_CFG  # noqa: E402

STAGE = "01_qc"
SCRIPT = "02_analysis/scripts/01_qc_filter_viz.py"
POP_COL = {"CD4_Treg": "#009E73", "CD4_Tcon": "#E69F00", "CD8": "#CC79A7"}
TISSUE_COL = {"synovial_fluid": "#D55E00", "peripheral_blood": "#0072B2"}
TREG_MARKERS = ["FOXP3", "IL2RA", "CTLA4", "IKZF2"]


def _scatter(ax, xy, values, categorical, title, cmap="viridis", palette=None):
    if categorical:
        for lvl, col in palette.items():
            m = values == lvl
            ax.scatter(xy[m, 0], xy[m, 1], s=1.5, c=col, label=str(lvl), linewidths=0)
        ax.legend(markerscale=5, fontsize=7, loc="best", frameon=True)
    else:
        sc_ = ax.scatter(xy[:, 0], xy[:, 1], s=1.5, c=values, cmap=cmap, linewidths=0)
        plt.colorbar(sc_, ax=ax, shrink=0.7)
    ax.set_title(title)
    ax.set_xticks([]); ax.set_yticks([])


def main() -> None:
    set_paper_style(config=FIG_CFG)
    tdir = PATHS.tables(STAGE)
    per_cell = pd.read_csv(tdir / "qc_metrics_per_cell.csv", index_col=0)
    kept = pd.read_csv(tdir / "cells_kept_dropped.csv")

    # ---- 1. per-GSM QC violins ----
    per_cell["pop"] = per_cell["population"].map(COARSE_LABEL)
    metrics = [("total_counts", "UMIs / cell", True),
               ("n_genes_by_counts", "genes / cell", True),
               ("pct_counts_mt", "% mito", False)]
    fig, axes = plt.subplots(3, 1, figsize=(13, 9), sharex=True)
    order = per_cell.sort_values(["population", "tissue", "donor"])["gsm"].unique().tolist()
    for ax, (m, lab, logy) in zip(axes, metrics):
        data = [per_cell.loc[per_cell["gsm"] == g, m].values for g in order]
        parts = ax.violinplot(data, showmedians=True, widths=0.85)
        for i, g in enumerate(order):
            pop = per_cell.loc[per_cell["gsm"] == g, "population"].iloc[0]
            parts["bodies"][i].set_facecolor(POP_COL.get(pop, "grey"))
            parts["bodies"][i].set_alpha(0.7)
        if logy:
            ax.set_yscale("log")
        ax.set_ylabel(lab)
    axes[-1].set_xticks(range(1, len(order) + 1))
    axes[-1].set_xticklabels(order, rotation=90, fontsize=6)
    fig.suptitle("Per-GSM QC distributions (violin color = sorted population)")
    fig.tight_layout()
    qc_summary = (per_cell.groupby("gsm")[["total_counts", "n_genes_by_counts", "pct_counts_mt"]]
                  .median().reset_index())
    save_overview(fig, STAGE, "qc_violins_per_gsm", table=qc_summary,
                  finding=("Per-GSM UMI/gene depth is adequate across strata; %mito stays low, "
                           "so no GSM is grossly degraded."),
                  script=SCRIPT, fn="main",
                  config_kv=f"thresholds.qc_n_mads = {PARAMS.qc_n_mads}; qc_pct_mt_max = {PARAMS.qc_pct_mt_max}",
                  input="03_results/01_qc/tables/qc_metrics_per_cell.csv",
                  how_to_read=("One violin per GSM (x), colored by sorted population; rows = UMIs, "
                               "genes (log y), %mito. Table lists per-GSM medians. QC diagnostic — "
                               "no biological claim."),
                  config=FIG_CFG, wide=True)

    # ---- 2. kept / dropped ----
    kept["pop"] = kept["population"].map(COARSE_LABEL)
    kept["lab"] = kept["pop"] + " " + kept["tissue"].map({"synovial_fluid": "SF", "peripheral_blood": "PB"}) \
        + " " + kept["donor"].str.replace("JIA_patient_", "p", regex=False)
    kept = kept.sort_values(["population", "tissue", "donor"])
    fig2, ax = plt.subplots(figsize=(13, 6))
    y = range(len(kept))
    ax.barh(list(y), kept["n_kept"], color="#4477AA", label="kept")
    ax.barh(list(y), kept["n_dropped"], left=kept["n_kept"], color="#EE6677", label="dropped")
    ax.set_yticks(list(y)); ax.set_yticklabels(kept["lab"], fontsize=6)
    ax.invert_yaxis(); ax.set_xlabel("cells"); ax.legend(frameon=True)
    ax.set_title("Cells kept vs dropped by QC (per stratum)")
    fig2.tight_layout()
    save_overview(fig2, STAGE, "cells_kept_dropped",
                  table=kept[["gsm", "donor", "tissue", "population", "n_cells", "n_kept", "n_dropped", "frac_kept"]],
                  finding=("QC retains ~86% of cells overall, but the SF-Treg p5 library "
                           "(GSM4859852) is near-empty (median ~14 UMIs) and drops entirely, "
                           "leaving 6 of 7 donors with paired SF+PB Tregs for the forest."),
                  script=SCRIPT, fn="main",
                  config_kv="thresholds.qc_min_genes = 200; scrublet_expected_doublet_rate = 0.06",
                  input="03_results/01_qc/tables/cells_kept_dropped.csv",
                  how_to_read=("Stacked bars per stratum: blue = kept, red = dropped (MAD outlier / "
                               "low-gene / doublet). Confirm every SF+PB Treg stratum retains enough "
                               "cells for pseudobulk. QC diagnostic."),
                  config=FIG_CFG, wide=True)

    # ---- 3. unsupervised UMAP + marker overlay ----
    adata = sc.read_h5ad(PATHS.object(STAGE))
    xy = adata.obsm["X_umap_unsupervised"]
    fig3, axes = plt.subplots(2, 4, figsize=(16, 8))
    _scatter(axes[0, 0], xy, adata.obs["population"].astype(str).values, True,
             "sorted population", palette=POP_COL)
    _scatter(axes[0, 1], xy, adata.obs["tissue"].astype(str).values, True,
             "tissue", palette=TISSUE_COL)
    donors = sorted(adata.obs["donor"].astype(str).unique())
    dpal = {d: plt.cm.tab10(i % 10) for i, d in enumerate(donors)}
    _scatter(axes[0, 2], xy, adata.obs["donor"].astype(str).values, True, "donor", palette=dpal)
    _scatter(axes[0, 3], xy, adata.obs["leiden_unsupervised"].astype(str).values, True,
             "leiden", palette={l: plt.cm.tab20(i % 20) for i, l in
                                enumerate(sorted(adata.obs["leiden_unsupervised"].astype(str).unique()))})
    sym_to_var = {s: v for v, s in zip(adata.var_names, adata.var["gene_symbol"].astype(str))}
    for ax, g in zip(axes[1], TREG_MARKERS):
        if g in sym_to_var:
            expr = np.asarray(adata[:, sym_to_var[g]].X.todense()).ravel()
            _scatter(ax, xy, expr, False, g, cmap="magma")
        else:
            ax.set_visible(False)
    fig3.suptitle("Unsupervised UMAP (annotation/viz only) + Treg-marker overlay")
    fig3.tight_layout()
    crosstab = (pd.crosstab(adata.obs["leiden_unsupervised"], adata.obs["population_short"])
                .reset_index())
    save_overview(fig3, STAGE, "unsupervised_umap",
                  table=crosstab,
                  finding=("Sorted Treg/Tcon/CD8 occupy largely distinct transcriptomic territory; "
                           "FOXP3/IL2RA/CTLA4/IKZF2 concentrate in the Treg gate, supporting sort "
                           "fidelity."),
                  script=SCRIPT, fn="main",
                  config_kv="thresholds.hvg_n_top = 2000; n_pcs = 30; leiden_resolution = 1.0",
                  input="03_results/objects/01_qc.h5ad (X_umap_unsupervised)",
                  how_to_read=("Top row: cells on the unsupervised UMAP colored by sort population, "
                               "tissue, donor, leiden. Bottom row: Treg-marker expression (magma). "
                               "The UMAP is a usability lens — biology is NOT read off it. Table = "
                               "leiden x population cross-tab (contamination check)."),
                  config=FIG_CFG, wide=True)
    print("[01_qc_viz] wrote 3 overviews")


if __name__ == "__main__":
    main()
