#!/usr/bin/env python
"""
02_annotate_states.py — COMPUTE (no plotting).
==============================================
Sort-identity check + FROZEN coarse labels for a SORTED dataset. On the critical
path, scVI->scANVI sub-state resolution is DEFERRED
until the go/no-go returns positive: a sorted dataset needs only the sort gate +
canonical markers + the unsupervised embedding to freeze Treg/Tcon/CD8.

What this does in place of scANVI:
  - freeze `coarse_label` from the sort `population` (basis = "sorting");
  - score canonical lineage modules (Treg / Tcon / CD8) per cell;
  - assign a marker-module `predicted_identity` (argmax of z-scored modules) and
    build the sort-vs-predicted confusion table — the BREAKPOINT-02 sort-fidelity
    check (a lightweight stand-in for the deferred scANVI confusion);
  - flag grossly inconsistent cells (`sort_consistent = False`).

Outputs:
  objects/02_annotation.h5ad   (frozen labels; both counts + lognorm retained)
  03_results/02_annotation/tables/confusion_sort_vs_predicted.csv
  03_results/02_annotation/tables/substate_markers.csv
  03_results/02_annotation/tables/counts_donor_by_label_tissue.csv
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

from config import PATHS, COARSE_LABEL  # noqa: E402

STAGE = "02_annotation"

LINEAGE_MODULES = {
    "Treg": ["FOXP3", "IL2RA", "CTLA4", "IKZF2", "TIGIT"],
    "Tcon": ["IL7R", "CD40LG", "ANK3", "LEF1"],
    "CD8": ["CD8A", "CD8B", "NKG7", "GZMK", "GZMA"],
}
# markers shown in the BREAKPOINT-02 dotplot
PANEL_MARKERS = ["FOXP3", "IL2RA", "CTLA4", "IKZF2", "TIGIT", "IL7R", "CD40LG",
                 "CD4", "CD8A", "CD8B", "GZMK", "NKG7"]


def main() -> None:
    adata = sc.read_h5ad(PATHS.object("01_qc"))

    # --- freeze coarse labels from the sort gate (sorting is the anchor) ---
    adata.obs["coarse_label"] = adata.obs["population"].map(COARSE_LABEL).astype("category")
    adata.obs["basis_of_label"] = "sorting"

    # --- canonical lineage module scores (z-scored, argmax = predicted identity) ---
    sym_to_var = {s: v for v, s in zip(adata.var_names, adata.var["gene_symbol"].astype(str))}
    module_cols = []
    for lineage, genes in LINEAGE_MODULES.items():
        present = [sym_to_var[g] for g in genes if g in sym_to_var]
        col = f"module_{lineage}"
        sc.tl.score_genes(adata, gene_list=present, score_name=col, use_raw=False, random_state=0)
        module_cols.append(col)
    Z = adata.obs[module_cols].apply(lambda s: (s - s.mean()) / (s.std() + 1e-9))
    pred = Z.idxmax(axis=1).str.replace("module_", "", regex=False)
    adata.obs["predicted_identity"] = pd.Categorical(pred, categories=["Treg", "Tcon", "CD8"])
    adata.obs["sort_consistent"] = (adata.obs["predicted_identity"].astype(str)
                                    == adata.obs["coarse_label"].astype(str))

    # --- tables ---
    tdir = PATHS.tables(STAGE)
    order = ["Treg", "Tcon", "CD8"]
    confusion = pd.crosstab(adata.obs["coarse_label"], adata.obs["predicted_identity"],
                            rownames=["sort"], colnames=["predicted"])
    # Force identical row+column order so the diagonal is the sort=predicted agreement.
    confusion = confusion.reindex(index=order, columns=order, fill_value=0)
    confusion.to_csv(tdir / "confusion_sort_vs_predicted.csv")

    # marker mean expression (lognorm) per sorted population -> dotplot source
    rows = []
    for g in PANEL_MARKERS:
        if g not in sym_to_var:
            continue
        expr = np.asarray(adata[:, sym_to_var[g]].X.todense()).ravel()
        for pop in ["Treg", "Tcon", "CD8"]:
            m = (adata.obs["coarse_label"].astype(str) == pop).to_numpy()
            rows.append({"gene": g, "coarse_label": pop,
                         "mean_lognorm": float(expr[m].mean()),
                         "frac_expressing": float((expr[m] > 0).mean())})
    pd.DataFrame(rows).to_csv(tdir / "substate_markers.csv", index=False)

    counts = (adata.obs.groupby(["donor", "coarse_label", "tissue"], observed=True)
              .size().rename("n_cells").reset_index()
              .sort_values(["coarse_label", "tissue", "donor"]))
    counts.to_csv(tdir / "counts_donor_by_label_tissue.csv", index=False)

    print("[02_annotation] sort-vs-predicted confusion:\n", confusion.to_string())
    print(f"[02_annotation] sort_consistent: {adata.obs['sort_consistent'].mean():.1%} of cells")

    out = PATHS.object(STAGE)
    adata.write_h5ad(out)
    print(f"[02_annotation] wrote checkpoint {out}: {adata.n_obs} cells (labels frozen)")


if __name__ == "__main__":
    main()
