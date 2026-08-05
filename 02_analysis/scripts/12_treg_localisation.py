#!/usr/bin/env python
"""
12_treg_localisation.py -- COMPUTE ONLY (no plotting).
=============================================================================
Per-cell AUCell scoring of five gene signatures on the sorted JIA dataset:
1. WT_heat_up (202 genes nominal, 177 in dataset)
2. Interaction_up (7 genes nominal, 7 in dataset)
3. Interaction_fdrOnly_up (19 genes nominal, 19 in dataset)
4. HALLMARK_HYPOXIA (200 genes nominal, 178 in dataset)
5. WT_heat_up_purged_hypoxia (184 genes nominal, 161 in dataset)

Subsets the 99,915-cell substrate by sort gate coarse_label == 'Treg' (FACS-sorted,
not score-selected) and calculates per score x tissue summary statistics and
power-band classifications.

Note: Interaction_down.txt and Interaction_fdrOnly_down.txt are 0-byte files
in mouse_anchor/03_results/human_projection/signatures/Interaction/. Those arms are
structurally absent at nominal zero (a property of the mouse contrast, not of JIA).

Outputs:
  03_results/12_treg_localisation/tables/treg_per_cell_scores.csv
  03_results/12_treg_localisation/tables/treg_localisation_summary.csv
  03_results/12_treg_localisation/tables/_overview/treg_localisation.csv
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

COMPARTMENT_ROOT = Path(__file__).resolve().parents[2]
REPO_ROOT = COMPARTMENT_ROOT.parent

sys.path.insert(0, str(COMPARTMENT_ROOT))
sys.path.insert(0, str(COMPARTMENT_ROOT / "02_analysis"))
os.chdir(COMPARTMENT_ROOT)

from config import (CONFIG, PATHS, PARAMS, DONOR_KEY, TISSUE_KEY, TISSUE_NUM,  # noqa: E402
                    TISSUE_DEN, COARSE_LABEL)
from helpers.geneset_utils import (load_alias_map, resolve_symbols,  # noqa: E402
                                   score_cells_aucell_ucell, _symbol_to_varname)

STAGE = "12_treg_localisation"
DATASET = "GSE160097"

# Signature definition file paths
SIG_PATHS = {
    "WT_heat_up": REPO_ROOT / "mouse_anchor/03_results/human_projection/signatures/WT_heat/WT_heat_up.txt",
    "Interaction_up": REPO_ROOT / "mouse_anchor/03_results/human_projection/signatures/Interaction/Interaction_up.txt",
    "Interaction_fdrOnly_up": REPO_ROOT / "mouse_anchor/03_results/human_projection/signatures/Interaction/Interaction_fdrOnly_up.txt",
    "HALLMARK_HYPOXIA": COMPARTMENT_ROOT / "00_data/references/msigdb_hallmark/HALLMARK_HYPOXIA.txt",
    "WT_heat_up_purged_hypoxia": COMPARTMENT_ROOT / "03_results/09_heat_hypoxia/tables/_signatures_purged/WT_heat_up.txt",
}

# 0-byte down arms to verify and log
DOWN_ARM_PATHS = {
    "Interaction_down": REPO_ROOT / "mouse_anchor/03_results/human_projection/signatures/Interaction/Interaction_down.txt",
    "Interaction_fdrOnly_down": REPO_ROOT / "mouse_anchor/03_results/human_projection/signatures/Interaction/Interaction_fdrOnly_down.txt",
}


def read_gene_list(path: Path) -> list[str]:
    if not path.exists():
        raise FileNotFoundError(f"Signature file absent at {path}")
    return [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]


def assign_power_band(n_nominal: int) -> str:
    if n_nominal >= 15:
        return "testable"
    elif n_nominal >= 5:
        return "underpowered_reported"
    else:
        return "untestable"


def main() -> None:
    # 1. Verify input signature files and 0-byte down arms
    gene_sets = {}
    nominal_sizes = {}
    for sig_name, p in SIG_PATHS.items():
        genes = read_gene_list(p)
        gene_sets[sig_name] = genes
        nominal_sizes[sig_name] = len(genes)
        print(f"[12_treg_localisation] {sig_name}: {len(genes)} nominal genes from {p.relative_to(REPO_ROOT)}")

    for down_name, p in DOWN_ARM_PATHS.items():
        if p.exists():
            size = len(read_gene_list(p))
            print(f"[12_treg_localisation] {down_name}: {size} bytes / {size} genes (structurally absent at nominal zero)")
        else:
            print(f"[12_treg_localisation] WARNING: {down_name} file missing at {p}")

    # 2. Load AnnData checkpoint
    h5ad_path = PATHS.object("02_annotation")
    print(f"[12_treg_localisation] Loading AnnData checkpoint from {h5ad_path}...")
    adata = sc.read_h5ad(h5ad_path)
    print(f"[12_treg_localisation] AnnData loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    # The sets ship current HGNC symbols and this object carries the hg19-era vintage, so
    # they are resolved into it before the effective size is counted. Nominal size is left
    # alone above: it is a fact about the source file, not about this matrix.
    sym_to_var = _symbol_to_varname(adata, "gene_symbol")
    alias_map = load_alias_map(CONFIG["symbol_alias"]["map_path"])
    effective_sizes = {}
    for sig_name, genes in list(gene_sets.items()):
        genes_r, applied = resolve_symbols(genes, alias_map, set(sym_to_var))
        gene_sets[sig_name] = genes_r
        in_data = [g for g in genes_r if g in sym_to_var]
        effective_sizes[sig_name] = len(in_data)
        print(f"[12_treg_localisation] {sig_name}: {len(in_data)}/{len(genes)} genes present in AnnData"
              + (f" (+{len(applied)} via alias: {' '.join(f'{a}->{b}' for a, b in applied)})"
                 if applied else ""))

    # 3. Score all 5 gene sets using established AUCell method (one parallel Rscript execution)
    print(f"[12_treg_localisation] Running score_cells_aucell_ucell for {len(gene_sets)} signatures...")
    n_cores = int(PARAMS.get("percell_score_ncores", 4))
    scores_df = score_cells_aucell_ucell(
        adata,
        gene_sets,
        layer=None,
        symbol_col="gene_symbol",
        n_cores=n_cores,
    )

    # 4. Assemble complete per-cell table
    per_cell = pd.DataFrame(
        {
            "barcode": adata.obs_names,
            DONOR_KEY: adata.obs[DONOR_KEY].values,
            TISSUE_KEY: adata.obs[TISSUE_KEY].values,
            "coarse_label": adata.obs["coarse_label"].values,
        },
        index=adata.obs_names,
    )
    for col in scores_df.columns:
        per_cell[col] = scores_df[col].values

    # Join embedding coordinates from hook_factor_substrate.parquet if available
    embed_path = PATHS.tables("07_embedding") / "hook_factor_substrate.parquet"
    if embed_path.exists():
        emb = pd.read_parquet(embed_path)
        if "x" in emb.columns and "y" in emb.columns:
            per_cell["x"] = emb["x"].values
            per_cell["y"] = emb["y"].values

    # 5. Verify reproduction of heat_hypoxia_colocalization.csv (Treg SF cell-level correlation)
    sf_tregs = per_cell[(per_cell[TISSUE_KEY] == TISSUE_NUM) & (per_cell["coarse_label"] == "Treg")]
    r_sp, p_sp = stats.spearmanr(sf_tregs["WT_heat_up_AUCell"], sf_tregs["HALLMARK_HYPOXIA_AUCell"])
    r_pe, p_pe = stats.pearsonr(sf_tregs["WT_heat_up_AUCell"], sf_tregs["HALLMARK_HYPOXIA_AUCell"])
    print(f"[12_treg_localisation] SF Treg per-cell colocalization reproduction check:")
    print(f"  Spearman r = {r_sp:.6f} (target in 09_heat_hypoxia: ~0.196072), N = {len(sf_tregs)}")
    print(f"  Pearson r  = {r_pe:.6f} (target in 09_heat_hypoxia: ~0.220789), N = {len(sf_tregs)}")

    # 6. Save full per-cell scores table
    out_dir = PATHS.tables(STAGE)
    out_dir.mkdir(parents=True, exist_ok=True)
    per_cell_out = out_dir / "treg_per_cell_scores.csv"
    per_cell.to_csv(per_cell_out, index=False)
    print(f"[12_treg_localisation] Saved per-cell scores ({len(per_cell)} rows, {per_cell.shape[1]} cols) to {per_cell_out}")

    # 7. Subset to Treg cells alone for summary statistics
    tregs = per_cell[per_cell["coarse_label"] == "Treg"].copy()
    print(f"[12_treg_localisation] Subsetted to Treg sort gate: {len(tregs)} cells across {tregs[DONOR_KEY].nunique()} donors")

    summary_rows = []
    for sig_name in SIG_PATHS:
        col = f"{sig_name}_AUCell"
        n_nom = nominal_sizes[sig_name]
        n_eff = effective_sizes[sig_name]
        band = assign_power_band(n_nom)
        is_interpretable = (n_nom >= 15)

        for tissue in [TISSUE_NUM, TISSUE_DEN]:
            sub = tregs[tregs[TISSUE_KEY] == tissue]
            vals = sub[col].dropna().values
            n_c = len(vals)
            n_d = sub[DONOR_KEY].nunique()

            mean_val = float(np.mean(vals)) if n_c > 0 else np.nan
            std_val = float(np.std(vals, ddof=1)) if n_c > 1 else np.nan
            median_val = float(np.median(vals)) if n_c > 0 else np.nan
            q25_val = float(np.percentile(vals, 25)) if n_c > 0 else np.nan
            q75_val = float(np.percentile(vals, 75)) if n_c > 0 else np.nan
            iqr_val = q75_val - q25_val if n_c > 0 else np.nan

            summary_rows.append({
                "signature": sig_name,
                "tissue": tissue,
                "coarse_label": "Treg",
                "n_cells": n_c,
                "n_donors": n_d,
                "mean_auc": mean_val,
                "std_auc": std_val,
                "median_auc": median_val,
                "q25_auc": q25_val,
                "q75_auc": q75_val,
                "iqr_auc": iqr_val,
                "set_size_nominal": n_nom,
                "set_size_effective": n_eff,
                "power_band": band,
                "interpretable_at_percell": is_interpretable,
                "evidence_tier": "secondary_percell",
            })

    summary_df = pd.DataFrame(summary_rows)
    summary_path = out_dir / "treg_localisation_summary.csv"
    summary_df.to_csv(summary_path, index=False)
    print(f"[12_treg_localisation] Saved summary statistics to {summary_path}")

    # Mirror to _overview directory for save_overview() consumption
    overview_dir = out_dir / "_overview"
    overview_dir.mkdir(parents=True, exist_ok=True)
    overview_csv = overview_dir / "treg_localisation.csv"
    summary_df.to_csv(overview_csv, index=False)
    print(f"[12_treg_localisation] Mirrored summary table to {overview_csv}")
    print("[12_treg_localisation] COMPUTE DONE.")


if __name__ == "__main__":
    main()
