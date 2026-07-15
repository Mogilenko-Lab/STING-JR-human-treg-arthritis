#!/usr/bin/env python
"""
10_hsr_lens.py — COMPUTE ONLY. Curated HSR second lens.
=======================================================
Purpose: carry the mouse-anchor curated, activation-free human HSR lens into the
JIA SF/PB T-cell compartment alongside empirical WT_heat_up, then ask whether the
SF-vs-PB separation survives a cleaner proteostasis lens and whether the same SF
cells are WT_heat_up-high and HSR-high.

Inputs:
  - 03_results/objects/02_annotation.h5ad
  - 00_data/references/temp_hsr_lens/{HSR_core,HSR_sensitivity}.txt
  - 03_results/03_pseudobulk/tables/ranked_{treg,tcon,cd8}.tsv
  - 03_results/05_scoring/tables/per_cell_scores.csv
  - ../mouse_anchor/03_results/human_projection/signatures/WT_heat/WT_heat_up.txt

Outputs:
  - 03_results/10_hsr_lens/tables/per_cell_hsr_scores.csv
  - 03_results/10_hsr_lens/tables/hsr_gsea_{treg,tcon,cd8}.csv
  - 03_results/10_hsr_lens/tables/hsr_lens_nes.csv
  - 03_results/10_hsr_lens/tables/hsr_colocalization.csv
  - 03_results/10_hsr_lens/tables/hsr_wtheatup_overlap.csv

Tier note: ANNOTATION / secondary tier, firewalled from the confirmatory WT_heat
claim spine. This script appends HSR NES rows only to master_gsea_pseudobulk and
does not write any row to effect_sizes_treg_arthritis.csv.

Honest ceiling: even the clean HSR core is proteotoxic-stress-general, not
fever-specific; only the mouse 37/39 contrast can measure thermal-ness. In JIA
we carry the lens and read it correlatively.

Sign convention: NES > 0 means the HSR set is enriched toward genes up in
synovial fluid vs paired peripheral blood.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "02_analysis"))
os.chdir(ROOT)

from config import DONOR_KEY, PARAMS, PATHS, TISSUE_KEY, TISSUE_NUM  # noqa: E402
from helpers.figure_style import FIG_CFG, append_master_table, round_numeric_cols  # noqa: E402
from helpers.geneset_utils import (  # noqa: E402
    _symbol_to_varname,
    load_signature,
    score_cells_aucell_ucell,
)

STAGE = "10_hsr_lens"
PRIMARY = "WT_heat"
DATASET = "GSE160097"
POP_TAG = {"Treg": "treg", "Tcon": "tcon", "CD8": "cd8"}
FGSEA_R = "02_analysis/helpers/fgsea_prerank.R"
HSR_DIR = ROOT / "00_data" / "references" / "temp_hsr_lens"
HSR_TERMS = ("HSR_core", "HSR_sensitivity")


def read_gene_list(path: Path) -> list[str]:
    return [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]


def write_gene_list(path: Path, genes: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(list(genes)) + "\n")


def load_hsr_sets() -> dict[str, list[str]]:
    missing = [term for term in HSR_TERMS if not (HSR_DIR / f"{term}.txt").exists()]
    if missing:
        raise FileNotFoundError(
            f"frozen HSR set(s) absent: {missing} under {HSR_DIR}. "
            "Regenerate with: Rscript 02_analysis/scripts/freeze_hsr_lens.R"
        )
    return {term: read_gene_list(HSR_DIR / f"{term}.txt") for term in HSR_TERMS}


def prepare_hsr_signature_dir(tables_dir: Path, hsr: dict[str, list[str]]) -> Path:
    sig_dir = tables_dir / "_signatures_hsr"
    if sig_dir.exists():
        shutil.rmtree(sig_dir)
    sig_dir.mkdir(parents=True, exist_ok=True)
    for term in HSR_TERMS:
        write_gene_list(sig_dir / f"{term}.txt", sorted(set(hsr[term])))
    return sig_dir


def run_fgsea(ranked_path: Path, out_csv: Path, contrast: str, sig_dir: Path) -> pd.DataFrame:
    cmd = [
        "Rscript",
        FGSEA_R,
        str(ranked_path),
        str(out_csv),
        contrast,
        str(PARAMS.gsea_min_size),
        str(PARAMS.gsea_max_size),
        str(PARAMS.gsea_seed),
        str(PARAMS.gsea_nperm),
        f"HSR_core={sig_dir / 'HSR_core.txt'}",
        f"HSR_sensitivity={sig_dir / 'HSR_sensitivity.txt'}",
    ]
    subprocess.run(cmd, check=True)
    return pd.read_csv(out_csv)


def drop_legacy_hsr_master_rows() -> None:
    """Remove this stage's pre-database-stamped HSR rows before re-appending.

    The fgsea helper writes `database=mouse_projection`; this stage owns only rows
    with stage=10_hsr_lens and pathway_id in HSR_TERMS. Removing that narrow slice
    lets append_master_table re-add rows under database=HSR_{tag}.
    """
    master = PATHS.master_file("master_gsea_pseudobulk.csv")
    if not master.exists():
        return
    df = pd.read_csv(master)
    if not {"stage", "pathway_id"}.issubset(df.columns):
        return
    owned = df["stage"].astype(str).eq(STAGE) & df["pathway_id"].astype(str).isin(HSR_TERMS)
    if owned.any():
        df.loc[~owned].to_csv(master, index=False)
        print(f"[10_hsr_lens] removed {int(owned.sum())} prior stage-10 HSR row(s) from master_gsea_pseudobulk.csv")


def leading_edge_from_row(row: pd.Series) -> str:
    for col in ("leading_edge", "core_enrichment"):
        if col in row and pd.notna(row[col]):
            return str(row[col]).replace("/", ";")
    return ""


def per_cell_hsr_scores(adata, hsr: dict[str, list[str]], tables_dir: Path) -> pd.DataFrame:
    sym_to_var = _symbol_to_varname(adata, "gene_symbol")
    for term, genes in hsr.items():
        print(f"[10_hsr_lens] {term} coverage: {sum(g in sym_to_var for g in genes)}/{len(genes)} genes present")

    scores = score_cells_aucell_ucell(
        adata,
        hsr,
        layer=None,
        symbol_col="gene_symbol",
        n_cores=int(PARAMS.get("percell_score_ncores", 4)),
    )
    score_cols = list(scores.columns)
    per_cell = adata.obs[[DONOR_KEY, TISSUE_KEY, "coarse_label"]].copy()
    for col in score_cols:
        per_cell[col] = scores[col].reindex(adata.obs_names.astype(str)).to_numpy()
    per_cell.to_csv(tables_dir / "per_cell_hsr_scores.csv")
    return per_cell


def hsr_nes(tables_dir: Path, sig_dir: Path) -> pd.DataFrame:
    drop_legacy_hsr_master_rows()
    ranked_dir = PATHS.tables("03_pseudobulk")
    rows = []
    for pop, tag in POP_TAG.items():
        ranked_path = ranked_dir / f"ranked_{tag}.tsv"
        if not ranked_path.exists():
            print(f"[10_hsr_lens] {pop}: no ranked list at {ranked_path} — skipping")
            continue
        gsea = run_fgsea(
            ranked_path,
            tables_dir / f"hsr_gsea_{tag}.csv",
            f"SF_vs_PB_{pop}_HSR",
            sig_dir,
        )
        gsea = gsea.copy()
        gsea["database"] = f"HSR_{tag}"
        gsea["evidence_tier"] = "secondary_annotation"
        append_master_table(
            gsea,
            database=f"HSR_{tag}",
            stage=STAGE,
            name="master_gsea_pseudobulk",
            config=FIG_CFG,
        )
        for _, r in gsea.iterrows():
            rows.append(
                {
                    "population": pop,
                    "signature": r["pathway_id"],
                    "contrast": "SF_vs_PB",
                    "nes": r["nes"],
                    "pvalue": r["pvalue"],
                    "padj": r["padj"],
                    "set_size": r["set_size"],
                    "leading_edge": leading_edge_from_row(r),
                    "evidence_tier": "secondary_annotation",
                }
            )
    out = pd.DataFrame(rows)
    out = out[
        [
            "population",
            "signature",
            "contrast",
            "nes",
            "pvalue",
            "padj",
            "set_size",
            "leading_edge",
            "evidence_tier",
        ]
    ]
    out = round_numeric_cols(out)
    out.to_csv(tables_dir / "hsr_lens_nes.csv", index=False)
    return out


def corr_pair(x: pd.Series, y: pd.Series, method: str) -> tuple[float, float]:
    ok = x.notna() & y.notna()
    if ok.sum() < 3:
        return np.nan, np.nan
    if method == "spearman":
        res = stats.spearmanr(x[ok], y[ok])
    elif method == "pearson":
        res = stats.pearsonr(x[ok], y[ok])
    else:
        raise ValueError(method)
    return float(res.statistic), float(res.pvalue)


def load_joined_cell_scores(hsr_per_cell: pd.DataFrame) -> pd.DataFrame:
    hsr = hsr_per_cell.copy()
    hsr.index = hsr.index.astype(str)
    if hsr.index.duplicated().any():
        raise ValueError("HSR per-cell score barcodes are not unique")

    heat = pd.read_csv(PATHS.tables("05_scoring") / "per_cell_scores.csv", index_col=0)
    heat.index = heat.index.astype(str)
    if heat.index.duplicated().any():
        raise ValueError("stage-05 per-cell score barcodes are not unique")
    if "WT_heat_up_AUCell" not in heat.columns:
        raise ValueError("stage-05 per-cell score table lacks WT_heat_up_AUCell")

    common = hsr.index.intersection(heat.index)
    if len(common) != len(hsr.index) or len(common) != len(heat.index):
        raise ValueError(
            f"cell-score join would drop cells: HSR={len(hsr.index)}, WT_heat={len(heat.index)}, overlap={len(common)}"
        )

    joined = hsr.join(heat[["WT_heat_up_AUCell"]], how="inner", validate="one_to_one")
    if len(joined) != len(hsr):
        raise ValueError("cell-score join changed row count")
    return joined


def hsr_colocalization(tables_dir: Path, hsr_per_cell: pd.DataFrame) -> pd.DataFrame:
    cells = load_joined_cell_scores(hsr_per_cell)
    sf = cells[cells[TISSUE_KEY].astype(str) == TISSUE_NUM].copy()
    rows = []
    for pop in POP_TAG:
        d = sf[sf["coarse_label"].astype(str) == pop]
        for hsr_col, hsr_term in (
            ("HSR_core_AUCell", "HSR_core"),
            ("HSR_sensitivity_AUCell", "HSR_sensitivity"),
        ):
            for method in ("spearman", "pearson"):
                r, p = corr_pair(d["WT_heat_up_AUCell"], d[hsr_col], method)
                rows.append(
                    {
                        "population": pop,
                        "hsr_term": hsr_term,
                        "level": "cell",
                        "method": method,
                        "r": r,
                        "pvalue": p,
                        "n": int(d[["WT_heat_up_AUCell", hsr_col]].dropna().shape[0]),
                        "evidence_tier": "secondary_percell",
                    }
                )

            dm = (
                d.groupby(DONOR_KEY, observed=True)[["WT_heat_up_AUCell", hsr_col]]
                .mean()
                .reset_index()
            )
            for method in ("spearman", "pearson"):
                r, p = corr_pair(dm["WT_heat_up_AUCell"], dm[hsr_col], method)
                rows.append(
                    {
                        "population": pop,
                        "hsr_term": hsr_term,
                        "level": "donor_sf_mean",
                        "method": method,
                        "r": r,
                        "pvalue": p,
                        "n": int(dm[["WT_heat_up_AUCell", hsr_col]].dropna().shape[0]),
                        "evidence_tier": "secondary_percell",
                    }
                )
    out = round_numeric_cols(pd.DataFrame(rows))
    out.to_csv(tables_dir / "hsr_colocalization.csv", index=False)
    return out


def hsr_wtheatup_overlap(tables_dir: Path, hsr: dict[str, list[str]]) -> pd.DataFrame:
    sig = load_signature(PATHS.signature_contract, PRIMARY)
    wt_up = set(sig["up"])
    rows = []
    for term, genes in hsr.items():
        hsr_set = set(genes)
        inter = sorted(wt_up & hsr_set)
        union = wt_up | hsr_set
        rows.append(
            {
                "set_a": "WT_heat_up",
                "set_b": term,
                "n_a": len(wt_up),
                "n_b": len(hsr_set),
                "n_intersect": len(inter),
                "jaccard": len(inter) / len(union) if union else np.nan,
                "genes_intersect": ";".join(inter),
            }
        )
    out = round_numeric_cols(pd.DataFrame(rows))
    out.to_csv(tables_dir / "hsr_wtheatup_overlap.csv", index=False)
    return out


def main() -> None:
    tables_dir = PATHS.tables(STAGE)
    adata = sc.read_h5ad(PATHS.object("02_annotation"))
    print(f"[10_hsr_lens] annotation object: {adata.n_obs} cells x {adata.n_vars} genes")

    hsr = load_hsr_sets()
    print("[10_hsr_lens] HSR sets: " + ", ".join(f"{k}={len(v)}" for k, v in hsr.items()))

    hsr_per_cell = per_cell_hsr_scores(adata, hsr, tables_dir)
    sig_dir = prepare_hsr_signature_dir(tables_dir, hsr)
    nes = hsr_nes(tables_dir, sig_dir)
    coloc = hsr_colocalization(tables_dir, hsr_per_cell)
    overlap = hsr_wtheatup_overlap(tables_dir, hsr)

    treg_core = nes[(nes["population"] == "Treg") & (nes["signature"] == "HSR_core")]
    treg_spear = coloc[
        (coloc["population"] == "Treg")
        & (coloc["hsr_term"] == "HSR_core")
        & (coloc["level"] == "cell")
        & (coloc["method"] == "spearman")
    ]
    if len(treg_core):
        print(
            "[10_hsr_lens] headline HSR_core Treg SF-vs-PB NES: "
            f"{float(treg_core['nes'].iloc[0]):.3g}, padj={float(treg_core['padj'].iloc[0]):.3g}"
        )
    if len(treg_spear):
        print(
            "[10_hsr_lens] headline within-SF Treg WT_heat_up vs HSR_core Spearman: "
            f"r={float(treg_spear['r'].iloc[0]):.3g}, n={int(treg_spear['n'].iloc[0])}"
        )
    print("[10_hsr_lens] WT_heat_up overlap:\n", overlap.to_string(index=False))
    print("[10_hsr_lens] done — annotation tier only; effect_sizes_treg_arthritis.csv untouched")


if __name__ == "__main__":
    main()
