#!/usr/bin/env python
"""
03_pseudobulk_de.py — COMPUTE (no plotting).
============================================
The biology: aggregate RAW integer UMI counts to donor x tissue x
frozen-label pseudobulk, then run paired SF-vs-PB DE per population (Treg primary;
Tcon + CD8 are the Treg-specificity control, on-path). Emits sign-preserving
ranked lists (signed Wald stat) per population for pre-ranked fgsea (stage 05).

Outputs:
  objects/03_pseudobulk.h5ad-free — matrices written as CSV under tables/:
  03_results/03_pseudobulk/tables/pseudobulk_counts.csv
  03_results/03_pseudobulk/tables/pseudobulk_coldata.csv
  03_results/03_pseudobulk/tables/de_SFvsPB_{treg,tcon,cd8}.csv
  03_results/03_pseudobulk/tables/ranked_{treg,tcon,cd8}.tsv
  03_results/03_pseudobulk/tables/de_summary.csv
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

from config import PATHS, PARAMS, TISSUE_KEY, DONOR_KEY, TISSUE_NUM, TISSUE_DEN  # noqa: E402
from helpers.anndata_utils import counts_are_integer  # noqa: E402
from helpers.pseudobulk_utils import (aggregate_pseudobulk, filter_strata,  # noqa: E402
                                      run_de_sf_vs_pb, ranked_by_stat)

STAGE = "03_pseudobulk"
POP_TAG = {"Treg": "treg", "Tcon": "tcon", "CD8": "cd8"}


def main() -> None:
    adata = sc.read_h5ad(PATHS.object("02_annotation"))
    assert counts_are_integer(adata, "counts"), "counts layer not integer"

    groupby = [DONOR_KEY, TISSUE_KEY, "coarse_label"]
    counts, coldata = aggregate_pseudobulk(adata, groupby, layer="counts")
    counts, coldata, dropped = filter_strata(counts, coldata, int(PARAMS.pseudobulk_min_cells))
    print(f"[03_pseudobulk] {counts.shape[0]} strata x {counts.shape[1]} genes; "
          f"dropped {len(dropped)} below floor")

    tdir = PATHS.tables(STAGE)
    counts.to_csv(tdir / "pseudobulk_counts.csv")
    coldata.to_csv(tdir / "pseudobulk_coldata.csv")
    if len(dropped):
        dropped.to_csv(tdir / "strata_dropped.csv", index=False)

    gene_symbols = dict(zip(adata.var_names.astype(str), adata.var["gene_symbol"].astype(str)))

    summary = []
    for pop, tag in POP_TAG.items():
        strata = coldata.index[coldata["coarse_label"] == pop]
        cd = coldata.loc[strata]
        n_sf = int((cd[TISSUE_KEY] == TISSUE_NUM).sum())
        n_pb = int((cd[TISSUE_KEY] == TISSUE_DEN).sum())
        if min(n_sf, n_pb) < int(PARAMS.pseudobulk_min_donors):
            print(f"[03_pseudobulk] {pop}: underpowered (SF={n_sf}, PB={n_pb}) — skipping DE")
            summary.append({"population": pop, "n_sf": n_sf, "n_pb": n_pb,
                            "model": "skipped", "n_sig_de": 0})
            continue
        de = run_de_sf_vs_pb(
            counts.loc[strata], cd, tissue_key=TISSUE_KEY, donor_key=DONOR_KEY,
            num=TISSUE_NUM, den=TISSUE_DEN, gene_symbols=gene_symbols)
        de = de.sort_values("pvalue")
        de.to_csv(tdir / f"de_SFvsPB_{tag}.csv")

        ranked = ranked_by_stat(de)
        ranked.to_csv(tdir / f"ranked_{tag}.tsv", sep="\t", header=False, index=False)

        n_sig = int(((de["padj"] < float(PARAMS.de_fdr)) &
                     (de["log2FoldChange"].abs() >= float(PARAMS.de_logfc))).sum())
        model = de["model"].iloc[0]
        print(f"[03_pseudobulk] {pop}: model {model}; {n_sig} sig DE "
              f"(padj<{PARAMS.de_fdr}, |lfc|>={PARAMS.de_logfc}); ranked {len(ranked)} genes")
        summary.append({"population": pop, "n_sf": n_sf, "n_pb": n_pb, "model": model,
                        "n_sig_de": n_sig, "n_ranked": len(ranked)})

    pd.DataFrame(summary).to_csv(tdir / "de_summary.csv", index=False)
    print("[03_pseudobulk] done")


if __name__ == "__main__":
    main()
