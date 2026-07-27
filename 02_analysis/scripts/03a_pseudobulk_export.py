#!/usr/bin/env python
"""
03a_pseudobulk_export.py — COMPUTE (no plotting).
=================================================
Aggregate RAW integer UMI counts to donor x tissue x
frozen-label pseudobulk. Emits CSVs for downstream R (limma-voom) DE.

Outputs:
  03_results/03_pseudobulk/tables/pseudobulk_counts.csv
  03_results/03_pseudobulk/tables/pseudobulk_coldata.csv
  03_results/03_pseudobulk/tables/gene_symbols.csv
  03_results/03_pseudobulk/tables/strata_dropped.csv

The counts matrix is keyed by Ensembl id, but every downstream consumer of the ranked
lists (the mouse-projection signatures, the Hallmark/HSR references) matches on HGNC
symbol. The map travels with the counts so the R seam can rename without reopening the
AnnData.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pandas as pd
import scanpy as sc

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "02_analysis"))
os.chdir(ROOT)

from config import PATHS, PARAMS, TISSUE_KEY, DONOR_KEY  # noqa: E402
from helpers.anndata_utils import counts_are_integer  # noqa: E402
from helpers.pseudobulk_utils import aggregate_pseudobulk, filter_strata  # noqa: E402

STAGE = "03_pseudobulk"

def main() -> None:
    adata = sc.read_h5ad(PATHS.object("02_annotation"))
    assert counts_are_integer(adata, "counts"), "counts layer not integer"

    groupby = [DONOR_KEY, TISSUE_KEY, "coarse_label"]
    counts, coldata = aggregate_pseudobulk(adata, groupby, layer="counts")
    counts, coldata, dropped = filter_strata(counts, coldata, int(PARAMS.pseudobulk_min_cells))
    print(f"[03a_pseudobulk_export] {counts.shape[0]} strata x {counts.shape[1]} genes; "
          f"dropped {len(dropped)} below floor")

    tdir = PATHS.tables(STAGE)
    counts.to_csv(tdir / "pseudobulk_counts.csv")
    coldata.to_csv(tdir / "pseudobulk_coldata.csv")

    # Ensembl -> HGNC map for the R seam. Restricted to the exported columns and in their
    # order, so a mismatch between counts and map is impossible by construction.
    gene_map = pd.DataFrame({
        "ensembl_id": adata.var_names.astype(str),
        "gene_symbol": adata.var["gene_symbol"].astype(str),
    }).set_index("ensembl_id").loc[counts.columns]
    gene_map.to_csv(tdir / "gene_symbols.csv")
    n_unmapped = int((gene_map["gene_symbol"].isin(["nan", "None", ""])).sum())
    print(f"[03a_pseudobulk_export] gene_map: {len(gene_map)} ids, "
          f"{gene_map['gene_symbol'].nunique()} distinct symbols, {n_unmapped} unmapped")

    if len(dropped):
        dropped.to_csv(tdir / "strata_dropped.csv", index=False)

    print("[03a_pseudobulk_export] done")


if __name__ == "__main__":
    main()
