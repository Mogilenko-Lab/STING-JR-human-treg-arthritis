#!/usr/bin/env python
"""
export_explorers.py — build compact interactive-explorer tables (utility, not a stage).
=======================================================================================
Materializes small parquet tables under `03_results/interactive/` that the jscatter
breakpoint notebooks (`02_analysis/notebooks/*_explore/`) load instantly — so the live
interactive kernel never has to open a multi-GB checkpoint. Read-only projection of
already-computed checkpoints; recomputes NO biology.

Re-run whenever the QC / annotation / scoring checkpoints change:
    python 02_analysis/scripts/export_explorers.py

Outputs (03_results/interactive/):
  01_qc_explore.parquet          x,y + QC obs + Treg/CD8 markers   (from 01_qc.h5ad)
  02_annotation_explore.parquet  x,y + frozen labels + markers     (from 02_annotation.h5ad)
  05_gonogo_explore.parquet      x,y + labels + per-cell WT_heat    (02_annotation.h5ad + scores)
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

from config import PATHS  # noqa: E402
from helpers.geneset_utils import derive_etreg_from_xlsx, load_signature  # noqa: E402

ETREG_XLSX = ROOT / "00_data/GSE161426_eTreg-signature/raw/GSE161426_Gene_expression_table_log2.xlsx"

OBSM = "X_umap_unsupervised"
MARKERS = ["FOXP3", "IL2RA", "CTLA4", "IKZF2", "TIGIT", "IL7R", "CD40LG",
           "CD4", "CD8A", "CD8B", "GZMK", "NKG7"]

# QC-overlay module scores (Tier-3 hand markers — eyeball only, never evidence). Used to
# discriminate a stressed-Treg pocket (HSP/proteostasis up, identity retained) from a
# low-quality/dying tail (apoptosis up, FOXP3 lost). Human HGNC symbols.
MODULES = {
    "score_HSP": ["HSPA1A", "HSPA1B", "HSPA6", "HSPH1", "HSPB1", "DNAJB1", "DNAJA1",
                  "BAG3", "HSPD1", "HSPE1", "SERPINH1", "AHSA1", "FKBP4"],
    "score_apoptosis": ["BAX", "BAK1", "BBC3", "PMAIP1", "BCL2L11", "CASP3", "CASP7",
                        "CASP8", "CASP9", "CYCS", "APAF1", "FAS", "TNFRSF10B", "DIABLO", "XAF1"],
}


def _coords(adata) -> pd.DataFrame:
    xy = np.asarray(adata.obsm[OBSM])
    return pd.DataFrame({"x": xy[:, 0], "y": xy[:, 1]}, index=adata.obs_names.astype(str))


def _add_obs(df, adata, cols):
    for c in cols:
        if c not in adata.obs:
            continue
        col = adata.obs[c]
        num = (pd.api.types.is_numeric_dtype(col)
               and not isinstance(col.dtype, pd.CategoricalDtype)
               and not pd.api.types.is_bool_dtype(col))  # bool -> string (jscatter needs str/category)
        df[c] = col.to_numpy() if num else col.astype(str).to_numpy()


def _sym_map(adata):
    m = {}
    for vn, sym in zip(adata.var_names.astype(str), adata.var["gene_symbol"].astype(str)):
        if sym and sym != "nan" and sym not in m:
            m[sym] = vn
    return m


def _add_genes(df, adata, genes):
    sym_to_var = _sym_map(adata)
    for g in genes:
        vn = sym_to_var.get(g)
        if vn is None:
            continue
        x = adata[:, vn].X
        df[g] = np.asarray(x.todense()).ravel() if hasattr(x, "todense") else np.asarray(x).ravel()


def _add_modules(df, adata, modules):
    """Score each gene set via scanpy score_genes (on lognorm X) and add to df."""
    sym_to_var = _sym_map(adata)
    for name, symbols in modules.items():
        present = [sym_to_var[s] for s in symbols if s in sym_to_var]
        if len(present) < 3:
            continue
        sc.tl.score_genes(adata, gene_list=present, score_name=name, use_raw=False, random_state=0)
        df[name] = adata.obs[name].to_numpy()


def _add_zscore_module(df, adata, name, symbols):
    """Safe small-set score: mean of per-gene z-scored lognorm expression (no random
    background, so it's stable for tiny gene sets like the 7-gene Interaction set)."""
    sym_to_var = _sym_map(adata)
    present = [sym_to_var[s] for s in symbols if s in sym_to_var]
    if not present:
        return
    X = adata[:, present].X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    Z = (X - X.mean(axis=0)) / (X.std(axis=0) + 1e-9)
    df[name] = Z.mean(axis=1)


def _slim(adata):
    """Free the heavy bits we don't need for the explorer (raw snapshot + counts layer);
    keep X (lognorm), var, obs, obsm. Halves peak memory so scoring doesn't OOM."""
    adata.raw = None
    for k in list(adata.layers.keys()):
        del adata.layers[k]
    return adata


def _persist_etreg(sig) -> None:
    d = ROOT / "00_data" / "references" / "etreg_GSE161426"
    d.mkdir(parents=True, exist_ok=True)
    (d / "eTreg_up.txt").write_text("\n".join(sig["up"]) + "\n")
    (d / "eTreg_down.txt").write_text("\n".join(sig["down"]) + "\n")
    sig["ranked"].to_csv(d / "eTreg_ranked.tsv", sep="\t", index=False, header=False)
    (d / "PROVENANCE.md").write_text(
        "# eTreg signature — derived from GSE161426 (EDA overlay)\n\n"
        "Contrast: SF Treg vs PB Treg on the log2 matrix "
        "(GSE161426_Gene_expression_table_log2.xlsx). Mean-difference + Welch t "
        f"(SF n={sig['n_sf']} vs PB n={sig['n_pb']}); up/down gated |log2FC|>=1 & p<0.05, "
        "capped 200. Source: Mijnheer/Lutter 2021, Nat Commun 12:2710.\n\n"
        "**Tier: derived-for-EDA.** A rigorous version (limma/DESeq2 on GEO counts, or the "
        "paper Supp core-signature list) is the post-go/no-go deliverable. Regenerated by "
        "`02_analysis/scripts/export_explorers.py`.\n")


def main() -> None:
    out = PATHS.interactive_dir()

    # eTreg signature (GSE161426 SF-vs-PB Treg) — derived once, persisted, scored per cell.
    etreg_mod = {}
    if ETREG_XLSX.exists():
        try:
            sig = derive_etreg_from_xlsx(ETREG_XLSX)
            etreg_mod = {"score_eTreg": sig["up"]}
            _persist_etreg(sig)
            print(f"[export] eTreg: {len(sig['up'])} up / {len(sig['down'])} down "
                  f"(SF n={sig['n_sf']} vs PB n={sig['n_pb']})")
        except Exception as exc:  # openpyxl / parse issues shouldn't break the explorers
            print(f"[export] eTreg derivation skipped: {exc}")
    else:
        print(f"[export] eTreg xlsx not found at {ETREG_XLSX} — skipping eTreg score")

    # Comparator mouse contrasts from the frozen contract: KO_heat (cGAS-KO heat, normal-size
    # set -> score_genes) and Interaction (cGAS-DEPENDENT heat = 7 IFN/STING genes; too small
    # for a stable module score, so we add a safe z-composite + each gene individually).
    contract = PATHS.signature_contract
    ko = load_signature(contract, "KO_heat")
    inter = load_signature(contract, "Interaction")
    inter_genes = list(inter["up"])   # IFI16 IFIT1B IRF7 MX1 RTP4 TRIM5 XAF1
    ko_mod = {"score_KO_heat_up": ko["up"], "score_KO_heat_down": ko["down"]}
    print(f"[export] comparators: KO_heat {len(ko['up'])}up/{len(ko['down'])}down; "
          f"Interaction {len(inter_genes)} genes {inter_genes}")

    # Per-cell WT_heat scores (from stage 05), joined by barcode where present.
    scores_csv = PATHS.tables("05_scoring") / "per_cell_scores.csv"
    scores = None
    if scores_csv.exists():
        scores = pd.read_csv(scores_csv, index_col=0)
        scores.index = scores.index.astype(str)

    def _join_scores(df):
        if scores is None:
            return
        # The per-cell WT_heat lens is now AUCell (stage 05, rank-based [0,1]); expose it to
        # the QC/07 maps under the historical overlay names so the map pathway stays coherent.
        auc = {"WT_heat_up": "WT_heat_up_AUCell", "WT_heat_down": "WT_heat_down_AUCell"}
        for name, col in auc.items():
            if col in scores.columns:
                df[name] = scores[col].reindex(df.index).to_numpy()
        if {"WT_heat_up", "WT_heat_down"} <= set(df.columns):
            df["WT_heat_updown"] = df["WT_heat_up"] - df["WT_heat_down"]

    # --- 01_qc explorer (QC + markers + mouse-stress score in one table) ---
    a1 = _slim(sc.read_h5ad(PATHS.object("01_qc")))
    df1 = _coords(a1)
    _add_obs(df1, a1, ["population", "population_short", "tissue", "donor",
                       "leiden_unsupervised", "pct_counts_mt", "total_counts",
                       "n_genes_by_counts", "doublet_score"])
    _add_genes(df1, a1, MARKERS + inter_genes)          # markers + the 7 Interaction genes one-by-one
    _add_modules(df1, a1, {**MODULES, **etreg_mod, **ko_mod})
    _add_zscore_module(df1, a1, "score_Interaction_z", inter_genes)   # safe 7-gene composite
    _join_scores(df1)                                    # WT_heat_up/down/updown (from stage 05)
    if {"score_KO_heat_up", "score_KO_heat_down"} <= set(df1.columns):
        df1["score_KO_heat_updown"] = df1["score_KO_heat_up"] - df1["score_KO_heat_down"]
    df1.to_parquet(out / "01_qc_explore.parquet")
    print(f"[export] 01_qc_explore.parquet {df1.shape}")
    del a1

    # --- 02_annotation + 05_gonogo explorers (share the annotation checkpoint) ---
    a2 = _slim(sc.read_h5ad(PATHS.object("02_annotation")))
    df2 = _coords(a2)
    _add_obs(df2, a2, ["coarse_label", "predicted_identity", "sort_consistent",
                       "tissue", "donor", "population"])
    _add_genes(df2, a2, MARKERS)
    _add_modules(df2, a2, etreg_mod)
    df2.to_parquet(out / "02_annotation_explore.parquet")
    print(f"[export] 02_annotation_explore.parquet {df2.shape}")

    # go/no-go: same coords + labels + per-cell WT_heat + eTreg scores + FOXP3/%mt cross-checks
    keep5 = ["x", "y", "coarse_label", "tissue", "donor", "FOXP3", "CTLA4", "IL2RA"]
    if "score_eTreg" in df2.columns:
        keep5.append("score_eTreg")
    df5 = df2[keep5].copy()
    _join_scores(df5)
    df5.to_parquet(out / "05_gonogo_explore.parquet")
    print(f"[export] 05_gonogo_explore.parquet {df5.shape}")

    print(f"[export] wrote 3 explorer tables to {out}")


if __name__ == "__main__":
    main()
