#!/usr/bin/env python
"""
09_heat_hypoxia.py -- COMPUTE ONLY.

Ask whether the JIA SF-vs-PB mouse WT_heat signal survives a simple hypoxia
gene purge, then keep two per-cell reads as secondary annotation:

1. Re-run the stage-05 fgsea helper on full WT_heat and on WT_heat with
   HALLMARK_HYPOXIA-overlapping genes removed.
2. Correlate WT_heat_up_AUCell with HALLMARK_HYPOXIA_AUCell within SF cells.
3. Classify the full WT_heat_up leading edge by heat-proteostatic genes,
   hypoxia-overlap genes, and other genes.

No plotting. Outputs are tables under 03_results/09_heat_hypoxia/tables.
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
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "02_analysis"))
os.chdir(ROOT)

from config import DONOR_KEY, PARAMS, PATHS, TISSUE_DEN, TISSUE_KEY, TISSUE_NUM  # noqa: E402
from helpers.figure_style import FIG_CFG, append_master_table  # noqa: E402
from helpers.geneset_utils import load_signature  # noqa: E402

STAGE = "09_heat_hypoxia"
PRIMARY = "WT_heat"
DATASET = "GSE160097"
POP_TAG = {"Treg": "treg", "Tcon": "tcon", "CD8": "cd8"}
FGSEA_R = "02_analysis/helpers/fgsea_prerank.R"
HYPOXIA_PATH = Path("00_data/references/msigdb_hallmark/HALLMARK_HYPOXIA.txt")

# Leading-edge gene taxonomy: a frozen, provenance-stamped classification of the
# WT_heat_up leading-edge genes into biological programs (heat_shock_proteostasis,
# hypoxia_HIF, immediate_early_stress, effector_activation, other), produced by an
# external large-context model (agy/Gemini) via delegate-cli and committed as a
# reference so the composition tally is reproducible. See its PROVENANCE.md.
TAXONOMY_PATH = Path("00_data/references/heat_leadingedge_taxonomy/leadingedge_gene_taxonomy.csv")
LE_CATEGORIES = [
    "heat_shock_proteostasis",
    "hypoxia_HIF",
    "immediate_early_stress",
    "effector_activation",
    "other",
]


def read_gene_list(path: Path) -> list[str]:
    return [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]


def write_gene_list(path: Path, genes: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(list(genes)) + "\n")


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
        f"{PRIMARY}_up={sig_dir / f'{PRIMARY}_up.txt'}",
        f"{PRIMARY}_down={sig_dir / f'{PRIMARY}_down.txt'}",
    ]
    subprocess.run(cmd, check=True)
    return pd.read_csv(out_csv)


def row_for(gsea: pd.DataFrame, pathway_id: str) -> pd.Series:
    hit = gsea[gsea["pathway_id"] == pathway_id]
    if len(hit) != 1:
        raise ValueError(f"expected one {pathway_id} row, found {len(hit)}")
    return hit.iloc[0]


def prepare_signature_dirs(tables_dir: Path, sig: dict[str, object], hypoxia: set[str]) -> tuple[Path, Path, dict[str, list[str]]]:
    full_dir = tables_dir / "_signatures_full"
    purged_dir = tables_dir / "_signatures_purged"
    for d in (full_dir, purged_dir):
        if d.exists():
            shutil.rmtree(d)
        d.mkdir(parents=True, exist_ok=True)

    up = list(sig["up"])
    down = list(sig["down"])
    removed = {
        "up": sorted(set(up) & hypoxia),
        "down": sorted(set(down) & hypoxia),
    }
    purged = {
        "up": [g for g in up if g not in hypoxia],
        "down": [g for g in down if g not in hypoxia],
    }
    write_gene_list(full_dir / f"{PRIMARY}_up.txt", up)
    write_gene_list(full_dir / f"{PRIMARY}_down.txt", down)
    write_gene_list(purged_dir / f"{PRIMARY}_up.txt", purged["up"])
    write_gene_list(purged_dir / f"{PRIMARY}_down.txt", purged["down"])
    return full_dir, purged_dir, {
        "full_up": up,
        "full_down": down,
        "purged_up": purged["up"],
        "purged_down": purged["down"],
        "removed_up": removed["up"],
        "removed_down": removed["down"],
    }


def gene_purge_nes(tables_dir: Path, sig: dict[str, object], hypoxia: set[str]) -> tuple[pd.DataFrame, dict[str, Path]]:
    full_dir, purged_dir, genes = prepare_signature_dirs(tables_dir, sig, hypoxia)
    ranked_dir = PATHS.tables("03_pseudobulk")
    runsum_paths: dict[str, Path] = {}
    rows = []
    for pop, tag in POP_TAG.items():
        ranked_path = ranked_dir / f"ranked_{tag}.tsv"
        if not ranked_path.exists():
            print(f"[09_heat_hypoxia] {pop}: no ranked list at {ranked_path} — skipping")
            continue
        full = run_fgsea(
            ranked_path,
            tables_dir / f"gsea_full_{tag}.csv",
            f"SF_vs_PB_{pop}_full",
            full_dir,
        )
        purged = run_fgsea(
            ranked_path,
            tables_dir / f"gsea_purged_{tag}.csv",
            f"SF_vs_PB_{pop}_hypoxia_purged",
            purged_dir,
        )
        full_up = row_for(full, f"{PRIMARY}_up")
        purged_up = row_for(purged, f"{PRIMARY}_up")
        rows.append(
            {
                "population": pop,
                "signature": f"{PRIMARY}_up",
                "contrast": "SF_vs_PB",
                "NES_full": full_up["nes"],
                "pvalue_full": full_up["pvalue"],
                "padj_full": full_up["padj"],
                "NES_purged": purged_up["nes"],
                "pvalue_purged": purged_up["pvalue"],
                "padj_purged": purged_up["padj"],
                "delta_NES_purged_minus_full": purged_up["nes"] - full_up["nes"],
                "n_genes_full": len(genes["full_up"]),
                "n_genes_purged": len(genes["purged_up"]),
                "n_genes_removed": len(genes["removed_up"]),
                "genes_removed": ";".join(genes["removed_up"]),
                "set_size_full_in_ranked": full_up["set_size"],
                "set_size_purged_in_ranked": purged_up["set_size"],
                "evidence_tier": "primary_pseudobulk",
            }
        )
        runsum_paths[pop] = tables_dir / f"runsum_interactive_gsea_full_{tag}_{PRIMARY}_up.csv"
    out = pd.DataFrame(rows)
    out.to_csv(tables_dir / "gene_purge_nes_comparison.csv", index=False)
    return out, runsum_paths


def load_cell_scores() -> pd.DataFrame:
    harvest = pd.read_parquet(PATHS.interactive_dir() / "08_harvest_readout.parquet")
    if "WT_heat_up_AUCell" in harvest.columns:
        cells = harvest.copy()
    else:
        scores = pd.read_csv(PATHS.tables("05_scoring") / "per_cell_scores.csv")
        scores = scores.rename(columns={scores.columns[0]: "barcode"})
        need = ["barcode", DONOR_KEY, TISSUE_KEY, "coarse_label", "WT_heat_up_AUCell"]
        missing = [c for c in need if c not in scores.columns]
        if missing:
            raise ValueError(f"stage-05 per-cell score table missing columns: {missing}")
        if scores["barcode"].duplicated().any():
            raise ValueError("stage-05 per-cell score barcodes are not unique")
        if harvest["barcode"].duplicated().any():
            raise ValueError("harvest readout barcodes are not unique")

        cells = harvest.merge(
            scores[need],
            on="barcode",
            how="left",
            suffixes=("", "_score"),
            validate="one_to_one",
        )
        if cells["WT_heat_up_AUCell"].isna().any():
            raise ValueError("WT_heat_up_AUCell join left missing values")
        for col in (DONOR_KEY, TISSUE_KEY, "coarse_label"):
            other = f"{col}_score"
            if other in cells.columns:
                bad = cells[col].astype(str) != cells[other].astype(str)
                if bad.any():
                    raise ValueError(f"{col} disagrees after WT_heat score join for {int(bad.sum())} cells")
                cells = cells.drop(columns=[other])
    return cells


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


def heat_hypoxia_colocalization(tables_dir: Path) -> pd.DataFrame:
    cells = load_cell_scores()
    sf = cells[cells[TISSUE_KEY] == TISSUE_NUM].copy()
    rows = []
    for pop in POP_TAG:
        d = sf[sf["coarse_label"].astype(str) == pop]
        for method in ("spearman", "pearson"):
            r, p = corr_pair(d["WT_heat_up_AUCell"], d["HALLMARK_HYPOXIA_AUCell"], method)
            rows.append(
                {
                    "population": pop,
                    "level": "cell",
                    "method": method,
                    "r": r,
                    "pvalue": p,
                    "n": int(d[["WT_heat_up_AUCell", "HALLMARK_HYPOXIA_AUCell"]].dropna().shape[0]),
                    "evidence_tier": "secondary_percell",
                }
            )

        dm = (
            d.groupby(DONOR_KEY, observed=True)[["WT_heat_up_AUCell", "HALLMARK_HYPOXIA_AUCell"]]
            .mean()
            .reset_index()
        )
        for method in ("spearman", "pearson"):
            r, p = corr_pair(dm["WT_heat_up_AUCell"], dm["HALLMARK_HYPOXIA_AUCell"], method)
            rows.append(
                {
                    "population": pop,
                    "level": "donor_sf_mean",
                    "method": method,
                    "r": r,
                    "pvalue": p,
                    "n": int(dm[["WT_heat_up_AUCell", "HALLMARK_HYPOXIA_AUCell"]].dropna().shape[0]),
                    "evidence_tier": "secondary_percell",
                }
            )
    out = pd.DataFrame(rows)
    out.to_csv(tables_dir / "heat_hypoxia_colocalization.csv", index=False)
    return out


def leadingedge_composition(tables_dir: Path, runsum_paths: dict[str, Path]) -> pd.DataFrame:
    """Tally each population's WT_heat_up leading edge by biological program, using the
    frozen agy/Gemini gene taxonomy (heat_shock_proteostasis / hypoxia_HIF /
    immediate_early_stress / effector_activation / other). Secondary annotation tier —
    it interprets WHAT the enriching genes are, never a confirmatory statistic."""
    tax = pd.read_csv(TAXONOMY_PATH)
    cat_of = dict(zip(tax["gene"].astype(str), tax["category"].astype(str)))
    tax_source = str(tax["source"].iloc[0]) if "source" in tax.columns and len(tax) else "unknown"
    rows = []
    for pop, path in runsum_paths.items():
        runsum = pd.read_csv(path)
        le_col = runsum["leading_edge"]
        # Robust to R emitting real bool vs "TRUE"/"FALSE" strings (a bare astype(bool)
        # on strings would silently mark the WHOLE ranked list as leading edge).
        if le_col.dtype == bool:
            le_mask = le_col
        else:
            le_mask = le_col.astype(str).str.strip().str.upper().isin({"TRUE", "T", "1"})
        genes = (runsum.loc[le_mask, "gene"].astype(str).drop_duplicates().tolist())
        n = len(genes)
        by_cat = {c: [g for g in genes if cat_of.get(g) == c] for c in LE_CATEGORIES}
        unclassified = [g for g in genes if g not in cat_of]
        rec = {"population": pop, "signature": f"{PRIMARY}_up", "n_leading_edge": n}
        for c in LE_CATEGORIES:
            rec[f"n_{c}"] = len(by_cat[c])
            rec[f"frac_{c}"] = len(by_cat[c]) / n if n else np.nan
        rec["n_unclassified"] = len(unclassified)
        for c in LE_CATEGORIES:
            rec[f"genes_{c}"] = ";".join(by_cat[c])
        rec["genes_unclassified"] = ";".join(unclassified)
        rec["taxonomy_source"] = tax_source
        rec["evidence_tier"] = "secondary_exploratory"
        rows.append(rec)
    out = pd.DataFrame(rows)
    out.to_csv(tables_dir / "leadingedge_composition.csv", index=False)
    total_unclassified = int(out["n_unclassified"].sum())
    if total_unclassified:
        print(f"[09_heat_hypoxia] WARNING: {total_unclassified} leading-edge gene(s) "
              "absent from the taxonomy (counted as unclassified) — regenerate the taxonomy.")
    return out


def update_effect_sizes(gene_purge: pd.DataFrame) -> pd.DataFrame:
    eff_path = PATHS.master_file("effect_sizes_treg_arthritis.csv")
    existing = pd.read_csv(eff_path) if eff_path.exists() else pd.DataFrame()
    primary_existing = pd.DataFrame()
    if not existing.empty:
        primary_existing = existing[
            (existing["signature"] == f"{PRIMARY}_up")
            & (existing["contrast"] == "SF_vs_PB")
            & (existing["effect_metric"] == "pseudobulk_nes")
            & (existing["evidence_tier"] == "primary_pseudobulk")
        ].copy()

    rows = []
    for _, r in gene_purge.iterrows():
        nes = r["NES_purged"]
        ref = primary_existing[primary_existing["cell_state"] == r["population"]]
        n_donors = ref["n_donors"].iloc[0] if len(ref) else np.nan
        n_cells = ref["n_cells"].iloc[0] if len(ref) else np.nan
        rows.append(
            {
                "dataset": DATASET,
                "signature": f"{PRIMARY}_up_purged_hypoxia",
                "cell_state": r["population"],
                "contrast": "SF_vs_PB_heat_purged",
                "effect_metric": "pseudobulk_nes",
                "evidence_tier": "primary_pseudobulk",
                "estimate": nes,
                "se": np.nan,
                "ci_low": np.nan,
                "ci_high": np.nan,
                "direction": "up" if pd.notna(nes) and nes > 0 else "down" if pd.notna(nes) else "ns",
                "pvalue": r["pvalue_purged"],
                "padj": r["padj_purged"],
                "n_donors": n_donors,
                "n_cells": n_cells,
            }
        )
    new = pd.DataFrame(rows)
    if existing.empty:
        combined = new
    else:
        key_cols = ["dataset", "signature", "cell_state", "contrast", "effect_metric", "evidence_tier"]
        old_key = existing[key_cols].astype(str).agg("\t".join, axis=1)
        new_key = new[key_cols].astype(str).agg("\t".join, axis=1)
        combined = pd.concat([existing.loc[~old_key.isin(set(new_key))], new], ignore_index=True)
        combined = combined[existing.columns.tolist()]
    combined.to_csv(eff_path, index=False)

    combined_for_master = combined.copy()
    combined_for_master["stage"] = np.where(
        combined_for_master["contrast"].astype(str).eq("SF_vs_PB_heat_purged"),
        STAGE,
        "05_scoring",
    )
    append_master_table(
        combined_for_master,
        database=DATASET,
        stage=STAGE,
        name="master_effect_sizes",
        config=FIG_CFG,
    )
    return new


def main() -> None:
    tables_dir = PATHS.tables(STAGE)
    hypoxia = set(read_gene_list(HYPOXIA_PATH))
    sig = load_signature(PATHS.signature_contract, PRIMARY)

    gene_purge, runsum_paths = gene_purge_nes(tables_dir, sig, hypoxia)
    colocalization = heat_hypoxia_colocalization(tables_dir)
    leadingedge = leadingedge_composition(tables_dir, runsum_paths)
    eff = update_effect_sizes(gene_purge)

    print("[09_heat_hypoxia] gene purge NES:")
    print(gene_purge[["population", "NES_full", "padj_full", "NES_purged", "padj_purged", "n_genes_removed"]].to_string(index=False))
    print("[09_heat_hypoxia] SF cell-level Spearman:")
    print(colocalization[(colocalization["level"] == "cell") & (colocalization["method"] == "spearman")][["population", "r", "n"]].to_string(index=False))
    print("[09_heat_hypoxia] leading-edge composition:")
    print(leadingedge[["population", "n_leading_edge", "frac_effector_activation",
                       "frac_immediate_early_stress", "frac_hypoxia_HIF",
                       "frac_heat_shock_proteostasis"]].to_string(index=False))
    print("[09_heat_hypoxia] effect-size rows:")
    print(eff[["cell_state", "contrast", "estimate", "padj"]].to_string(index=False))
    print("[09_heat_hypoxia] done")


if __name__ == "__main__":
    main()
