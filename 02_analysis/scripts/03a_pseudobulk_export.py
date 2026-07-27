#!/usr/bin/env python
"""
03a_pseudobulk_export.py — COMPUTE (no plotting).
=================================================
Aggregate RAW integer UMI counts to donor x tissue x
frozen-label pseudobulk. Emits CSVs for downstream R (limma-voom) DE.

Outputs:
  03_results/03_pseudobulk/tables/pseudobulk_counts.csv
  03_results/03_pseudobulk/tables/pseudobulk_coldata.csv
  03_results/03_pseudobulk/tables/strata_dropped.csv
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

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
    if len(dropped):
        dropped.to_csv(tdir / "strata_dropped.csv", index=False)
        
    print("[03a_pseudobulk_export] done")


if __name__ == "__main__":
    main()
