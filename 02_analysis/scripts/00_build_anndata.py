#!/usr/bin/env python
"""
00_build_anndata.py — COMPUTE (no plotting).
============================================
Ingest the 40 per-GSM CellRanger filtered 10x H5 files (GSE160097), assemble a
gene-union AnnData with donor/tissue/population provenance from samples.csv,
preserve raw integer UMI counts in `layers['counts']`, and flag the two
intentionally-absent samples (p3 PB Tcon/CD8) + the low-input SF-Treg p5.

TCR ingest is DEFERRED (a later critical-path step) — no scirpy join here.

Outputs:
  objects/00_build.h5ad
  03_results/00_build/tables/cells_per_gsm.csv
  03_results/00_build/tables/genes_union_summary.csv

Figures live in 00_build_anndata_viz.py.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "02_analysis"))
os.chdir(ROOT)

from config import PATHS, SPECIES_DB, COARSE_LABEL  # noqa: E402
from helpers.anndata_utils import build_pooled, counts_are_integer  # noqa: E402

STAGE = "00_build"


def main() -> None:
    samples = pd.read_csv(PATHS.samples_csv)
    print(f"[00_build] {len(samples)} GSMs in samples.csv; reading H5 from {PATHS.raw}")

    adata = build_pooled(samples, PATHS.raw)
    adata.obs["population_short"] = adata.obs["population"].map(COARSE_LABEL).astype("category")
    for c in ["gsm", "title", "donor", "condition", "tissue", "population"]:
        adata.obs[c] = adata.obs[c].astype("category")

    assert counts_are_integer(adata, "counts"), "counts layer is not integer-valued"
    print(f"[00_build] pooled: {adata.n_obs} cells x {adata.n_vars} genes (union)")

    # --- tables ---
    tdir = PATHS.tables(STAGE)
    cells_per_gsm = (
        adata.obs.groupby(["gsm", "donor", "tissue", "population"], observed=True)
        .size().rename("n_cells").reset_index()
        .sort_values(["population", "tissue", "donor"])
    )
    cells_per_gsm.to_csv(tdir / "cells_per_gsm.csv", index=False)

    # Which genes are detected (union coverage): per-gene number of GSMs w/ >0 counts is
    # expensive; summarise union size + non-empty gene count instead.
    genes_detected = int(((adata.layers["counts"] > 0).sum(axis=0) > 0).sum())
    genes_union_summary = pd.DataFrame([{
        "n_genes_union": int(adata.n_vars),
        "n_genes_detected_any_cell": genes_detected,
        "n_gsm": int(adata.obs["gsm"].nunique()),
        "n_cells_total": int(adata.n_obs),
        "species_db": SPECIES_DB,
    }])
    genes_union_summary.to_csv(tdir / "genes_union_summary.csv", index=False)

    # The feature union by NAME, not only by count. Downstream a gene set member absent
    # from the count matrix has two opposite causes — the symbol vintage differs, or the
    # gene is genuinely not expressed in sorted T cells — and only this list separates
    # them. EGFR, EPCAM and INHBA are all absent from the post-QC matrix and all present
    # in the CellRanger reference: they are a detection fact. MB21D1 and TMEM173 are
    # present under names no current gene set uses: that is a vocabulary fact. A count
    # cannot tell those apart, so the symbols themselves are persisted.
    pd.DataFrame({
        "ensembl_id": adata.var_names.astype(str),
        "gene_symbol": adata.var["gene_symbol"].astype(str),
    }).to_csv(tdir / "reference_feature_symbols.csv", index=False)

    # Design-completeness table: expected 42 (3 pops x 2 tissues x 7 donors) minus 2 absent.
    design = (
        adata.obs.groupby(["population", "tissue"], observed=True)["donor"].nunique()
        .rename("n_donors").reset_index()
    )
    design.to_csv(tdir / "design_completeness.csv", index=False)
    print("[00_build] design completeness:\n", design.to_string(index=False))

    # --- checkpoint ---
    out = PATHS.object(STAGE)
    adata.write_h5ad(out)
    print(f"[00_build] wrote checkpoint {out}")


if __name__ == "__main__":
    main()
