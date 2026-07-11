#!/usr/bin/env python
"""
02_annotate_states_viz.py — VIZ (no statistics).
================================================
Reads the frozen-label tables + checkpoint from 02_annotate_states.py and renders
the sort-fidelity overview for BREAKPOINT 02: the unsupervised UMAP by frozen
label / marker-predicted identity, the canonical-marker dotplot, the
sort-vs-predicted confusion heatmap, and the donor x label x tissue count grid
(under-powered strata flagged).
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

from config import PATHS, PARAMS  # noqa: E402
from helpers.figure_style import set_paper_style, save_overview, FIG_CFG  # noqa: E402

STAGE = "02_annotation"
SCRIPT = "02_analysis/scripts/02_annotate_states_viz.py"
LABEL_COL = {"Treg": "#009E73", "Tcon": "#E69F00", "CD8": "#CC79A7"}


def _scatter(ax, xy, values, palette, title):
    for lvl, col in palette.items():
        m = values == lvl
        ax.scatter(xy[m, 0], xy[m, 1], s=1.5, c=col, label=str(lvl), linewidths=0)
    ax.legend(markerscale=5, fontsize=8, frameon=True)
    ax.set_title(title); ax.set_xticks([]); ax.set_yticks([])


def main() -> None:
    set_paper_style(config=FIG_CFG)
    tdir = PATHS.tables(STAGE)
    confusion = pd.read_csv(tdir / "confusion_sort_vs_predicted.csv", index_col=0)
    markers = pd.read_csv(tdir / "substate_markers.csv")
    counts = pd.read_csv(tdir / "counts_donor_by_label_tissue.csv")
    adata = sc.read_h5ad(PATHS.object(STAGE))
    xy = adata.obsm["X_umap_unsupervised"]

    # ---- 1. UMAP by frozen label / predicted / consistency ----
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    _scatter(axes[0], xy, adata.obs["coarse_label"].astype(str).values, LABEL_COL,
             "frozen label (sort)")
    _scatter(axes[1], xy, adata.obs["predicted_identity"].astype(str).values, LABEL_COL,
             "marker-module predicted")
    _scatter(axes[2], xy, adata.obs["sort_consistent"].astype(str).values,
             {"True": "#4477AA", "False": "#EE6677"}, "sort-consistent")
    fig.suptitle("Sort identity vs marker-module prediction (scANVI deferred)")
    fig.tight_layout()
    save_overview(fig, STAGE, "umap_sort_identity", table=confusion.reset_index(),
                  finding=("Frozen sort labels track the transcriptomic structure; the large "
                           "majority of cells are marker-consistent, so the sort gate is a "
                           "sound anchor for pseudobulk."),
                  script=SCRIPT, fn="main",
                  config_kv="basis_of_label = sorting (scANVI deferred until go = yes)",
                  input="03_results/objects/02_annotation.h5ad",
                  how_to_read=("Unsupervised UMAP colored by frozen sort label, marker-module "
                               "argmax prediction, and their agreement (blue=consistent). "
                               "Disagreement flags candidate mis-sorts, not Treg collapse. "
                               "Annotation/viz only — no biological claim."),
                  config=FIG_CFG, wide=True)

    # ---- 2. marker dotplot ----
    genes = markers["gene"].unique().tolist()
    pops = ["Treg", "Tcon", "CD8"]
    fig2, ax = plt.subplots(figsize=(12, 4.5))
    for yi, pop in enumerate(pops):
        sub = markers[markers["coarse_label"] == pop].set_index("gene").reindex(genes)
        sizes = sub["frac_expressing"].values * 300
        colors = sub["mean_lognorm"].values
        sctr = ax.scatter(range(len(genes)), [yi] * len(genes), s=sizes, c=colors,
                          cmap="Reds", edgecolors="grey", linewidths=0.5, vmin=0)
    ax.set_yticks(range(len(pops))); ax.set_yticklabels(pops)
    ax.set_xticks(range(len(genes))); ax.set_xticklabels(genes, rotation=90)
    ax.set_title("Canonical marker expression by frozen label")
    plt.colorbar(sctr, ax=ax, shrink=0.6, label="mean lognorm")
    fig2.tight_layout()
    save_overview(fig2, STAGE, "marker_dotplot", table=markers,
                  finding=("FOXP3/IL2RA/CTLA4/IKZF2 are Treg-restricted while CD8A/B/GZMK mark the "
                           "CD8 gate; IL7R is depleted in Tregs (CD127-lo sort), confirming gate fidelity."),
                  script=SCRIPT, fn="main",
                  config_kv="LINEAGE_MODULES (Treg/Tcon/CD8 canonical markers)",
                  input="03_results/02_annotation/tables/substate_markers.csv",
                  how_to_read=("Dot size = fraction of cells expressing; color = mean lognorm "
                               "expression. Rows = frozen label. Confirms markers land where the "
                               "sort predicts. QC overlay tier (hand markers, not evidence)."),
                  config=FIG_CFG, wide=True)

    # ---- 3. donor x label x tissue count grid ----
    floor = int(PARAMS.pseudobulk_min_cells)
    counts["donor_id"] = counts["donor"].str.replace("JIA_patient_", "p", regex=False)
    counts["tissue_s"] = counts["tissue"].map({"synovial_fluid": "SF", "peripheral_blood": "PB"})
    counts["cell"] = counts["coarse_label"] + " " + counts["tissue_s"]
    grid = counts.pivot_table(index="cell", columns="donor_id", values="n_cells", fill_value=0)
    fig3, ax = plt.subplots(figsize=(9, 6))
    im = ax.imshow(grid.values, cmap="Blues", aspect="auto")
    ax.set_xticks(range(grid.shape[1])); ax.set_xticklabels(grid.columns)
    ax.set_yticks(range(grid.shape[0])); ax.set_yticklabels(grid.index)
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            v = int(grid.values[i, j])
            flag = "*" if 0 < v < floor else ""
            ax.text(j, i, f"{v}{flag}", ha="center", va="center",
                    fontsize=7, color="red" if 0 < v < floor else "black")
    ax.set_title(f"Cells per donor x label x tissue ( * = below {floor}-cell floor )")
    plt.colorbar(im, ax=ax, shrink=0.6)
    fig3.tight_layout()
    save_overview(fig3, STAGE, "counts_grid", table=counts,
                  finding=("Every SF+PB Treg stratum clears the pseudobulk floor across the 7 donors; "
                           "p3 PB Tcon/CD8 are absent by design (empty cells)."),
                  script=SCRIPT, fn="main",
                  config_kv=f"thresholds.pseudobulk_min_cells = {floor}",
                  input="03_results/02_annotation/tables/counts_donor_by_label_tissue.csv",
                  how_to_read=("Heatmap of cells per donor (x) x label+tissue (y); red * marks a "
                               "stratum below the pseudobulk cell floor; empty = intentionally-absent "
                               "sample. Donor count per arm is the forest's power. Diagnostic."),
                  config=FIG_CFG)
    print("[02_annotation_viz] wrote 3 overviews")


if __name__ == "__main__":
    main()
