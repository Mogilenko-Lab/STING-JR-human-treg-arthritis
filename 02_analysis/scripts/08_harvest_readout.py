#!/usr/bin/env python
"""
08_harvest_readout.py — COMPUTE (no plotting). Exploratory READOUT harvest.
===========================================================================
Adds two MSigDB Hallmark per-cell readouts — HALLMARK_HYPOXIA and
HALLMARK_UNFOLDED_PROTEIN_RESPONSE — to the frozen JIA SF/PB T-cell annotation,
aligned column-for-column with the sting compartment's Phase-1 readout so the two
are directly comparable. These are READOUTS, not causal claims: a hypoxia score is
"consistent with a low-O2 / metabolically-stressed state", NEVER a HIF-causality
statement. This whole stage is SECONDARY / annotation tier — per-cell, rank-based,
composition-robust — and is NEVER pooled with the confirmatory donor-pseudobulk NES
spine (stage 05) or the OR-gated POI harvest (stage 07). It feeds the reactive
review and the cross-species harvest question only.

Engine: `score_cells_aucell_ucell` (helpers/geneset_utils.py → percell_score.R), the
same AUCell+UCell rank-based scorer stage 05 and the sting compartment use. AUCell is
the canonical readout; UCell rides alongside as a cross-check. The already-derived
per-cell readouts of THIS compartment — score_HSP, score_eTreg, WT_heat_updown — are
NOT re-derived here; they are carried in verbatim from the frozen explorer parquets so
the review has one tidy table.

Outputs:
  03_results/interactive/08_harvest_readout.parquet          (per-cell; gitignored/regenerable)
  03_results/08_harvest_readout/tables/harvest_readout_summary.csv
  03_results/08_harvest_readout/README.md                    (caption; written by hand, not here)

Run in-container from the compartment root (or anywhere, via `from config import`):
    python 02_analysis/scripts/08_harvest_readout.py
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

from config import (PATHS, PARAMS, TISSUE_KEY, DONOR_KEY, TISSUE_NUM,  # noqa: E402
                    TISSUE_DEN)
from helpers.geneset_utils import score_cells_aucell_ucell, _symbol_to_varname  # noqa: E402

STAGE = "08_harvest_readout"
COARSE = "coarse_label"

# Frozen Hallmark sets (freeze_hallmark_sets.R). Scored per-cell → <set>_{AUCell,UCell}.
HALLMARK_DIR = ROOT / "00_data" / "references" / "msigdb_hallmark"
HALLMARK_SETS = ["HALLMARK_HYPOXIA", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"]

# Readout name → source column in the assembled per-cell table. AUCell is canonical
# for the two new Hallmark sets (aligned with the sting compartment).
READOUTS = {
    "WT_heat_updown": "WT_heat_updown",                       # mouse 39C anchor (annotation only)
    "score_eTreg": "score_eTreg",                             # effector-Treg (GSE161426)
    "score_HSP": "score_HSP",                                 # heat-shock / proteostasis
    "HALLMARK_HYPOXIA": "HALLMARK_HYPOXIA_AUCell",
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE": "HALLMARK_UNFOLDED_PROTEIN_RESPONSE_AUCell",
}


def _read_gene_list(path: Path) -> list[str]:
    with open(path) as fh:
        return [ln.strip() for ln in fh if ln.strip()]


def load_hallmark_sets() -> dict[str, list[str]]:
    missing = [s for s in HALLMARK_SETS if not (HALLMARK_DIR / f"{s}.txt").exists()]
    if missing:
        raise FileNotFoundError(
            f"frozen Hallmark set(s) absent: {missing} under {HALLMARK_DIR}. "
            "Regenerate with: Rscript 02_analysis/scripts/freeze_hallmark_sets.R")
    return {s: _read_gene_list(HALLMARK_DIR / f"{s}.txt") for s in HALLMARK_SETS}


def carry_frozen_readouts(index: pd.Index) -> pd.DataFrame:
    """Join already-derived per-cell readouts (score_HSP, score_eTreg, WT_heat_updown)
    from the frozen explorer parquets — reindexed to the annotation object's barcodes.
    These are carried verbatim, NEVER re-derived (guardrail)."""
    idir = PATHS.interactive_dir()
    out = pd.DataFrame(index=index)

    gonogo = pd.read_parquet(idir / "05_gonogo_explore.parquet")
    gonogo.index = gonogo.index.astype(str)
    for col in ("WT_heat_updown", "score_eTreg"):
        out[col] = gonogo[col].reindex(index).to_numpy() if col in gonogo else np.nan

    qc = pd.read_parquet(idir / "01_qc_explore.parquet")
    qc.index = qc.index.astype(str)
    out["score_HSP"] = qc["score_HSP"].reindex(index).to_numpy() if "score_HSP" in qc else np.nan
    return out


def build_parquet(adata, scored: pd.DataFrame, carried: pd.DataFrame) -> pd.DataFrame:
    """One row per cell: barcode + UMAP + coarse_label/tissue/donor + %mt + the two
    new Hallmark AUCell/UCell scores + the three carried readouts."""
    idx = adata.obs_names.astype(str)
    df = pd.DataFrame(index=idx)
    df["barcode"] = idx
    um = np.asarray(adata.obsm["X_umap_unsupervised"])
    df["x"] = um[:, 0]
    df["y"] = um[:, 1]
    df[COARSE] = adata.obs[COARSE].astype(str).to_numpy()
    df[TISSUE_KEY] = adata.obs[TISSUE_KEY].astype(str).to_numpy()
    df[DONOR_KEY] = adata.obs[DONOR_KEY].astype(str).to_numpy()
    df["pct_counts_mt"] = adata.obs["pct_counts_mt"].to_numpy()
    for col in scored.columns:                       # HALLMARK_*_{AUCell,UCell}
        df[col] = scored[col].reindex(idx).to_numpy()
    for col in carried.columns:                      # WT_heat_updown, score_eTreg, score_HSP
        df[col] = carried[col].reindex(idx).to_numpy()
    return df.reset_index(drop=True)


def summarise(df: pd.DataFrame) -> pd.DataFrame:
    """Per (coarse_label × tissue × readout): mean/spread + mean %mt + n, plus a compact
    verdict — fraction above the readout's global P90 (high-pocket signal), and the
    per-cell SF-vs-PB shift (mean difference + standardized d), attached per cell-state.
    Long/tidy, mirroring the sting Phase-1 summary. Descriptive, correlative only."""
    # Global P90 per readout (across all cells) → the high-pocket threshold.
    p90 = {name: np.nanpercentile(df[col].to_numpy(dtype=float), 90)
           for name, col in READOUTS.items()}

    rows = []
    for label, gl in df.groupby(COARSE, observed=True):
        # per cell-state SF-vs-PB per-cell shift (mean diff + Cohen's d), per readout.
        sf_all = gl[gl[TISSUE_KEY] == TISSUE_NUM]
        pb_all = gl[gl[TISSUE_KEY] == TISSUE_DEN]
        shift = {}
        for name, col in READOUTS.items():
            sf = sf_all[col].to_numpy(dtype=float)
            pb = pb_all[col].to_numpy(dtype=float)
            sf = sf[np.isfinite(sf)]
            pb = pb[np.isfinite(pb)]
            mean_shift = (float(np.mean(sf)) - float(np.mean(pb))
                          if len(sf) and len(pb) else np.nan)
            d = np.nan
            if len(sf) >= 2 and len(pb) >= 2:
                n1, n2 = len(sf), len(pb)
                pooled = np.sqrt(((n1 - 1) * sf.std(ddof=1) ** 2
                                  + (n2 - 1) * pb.std(ddof=1) ** 2) / (n1 + n2 - 2))
                d = (sf.mean() - pb.mean()) / (pooled + 1e-12)
            shift[name] = (mean_shift, d)

        for tissue, gt in gl.groupby(TISSUE_KEY, observed=True):
            for name, col in READOUTS.items():
                v = gt[col].to_numpy(dtype=float)
                v = v[np.isfinite(v)]
                rows.append(dict(
                    coarse_label=label, tissue=tissue, readout=name,
                    n_cells=int(len(v)),
                    mean=float(np.mean(v)) if len(v) else np.nan,
                    sd=float(np.std(v, ddof=1)) if len(v) > 1 else np.nan,
                    median=float(np.median(v)) if len(v) else np.nan,
                    frac_above_p90=(float(np.mean(v > p90[name])) if len(v) else np.nan),
                    mean_pct_counts_mt=float(np.nanmean(gt["pct_counts_mt"].to_numpy(dtype=float)))
                    if len(gt) else np.nan,
                    sf_minus_pb_mean=shift[name][0],
                    sf_minus_pb_smd=shift[name][1],
                    global_p90=float(p90[name])))
    out = pd.DataFrame(rows)
    order = {"Treg": 0, "Tcon": 1, "CD8": 2}
    out = out.sort_values(
        by=["coarse_label", "readout", "tissue"],
        key=lambda s: s.map(order) if s.name == "coarse_label" else s
    ).reset_index(drop=True)
    return out


def main() -> None:
    n_cores = int(PARAMS.get("percell_score_ncores", 4))
    hallmark = load_hallmark_sets()
    print(f"[08_harvest] Hallmark sets: "
          + ", ".join(f"{k}={len(v)}" for k, v in hallmark.items()))

    adata = sc.read_h5ad(PATHS.object("02_annotation"))
    print(f"[08_harvest] annotation object: {adata.n_obs} cells x {adata.n_vars} genes")

    # per-cell coverage (how many set genes are present in the object)
    sym_to_var = _symbol_to_varname(adata, "gene_symbol")
    for k, v in hallmark.items():
        print(f"[08_harvest] {k} coverage: "
              f"{sum(g in sym_to_var for g in v)}/{len(v)} genes present")

    # --- score the two Hallmark sets per cell (AUCell + UCell on lognorm X) ---
    scored = score_cells_aucell_ucell(adata, hallmark, layer=None,
                                      symbol_col="gene_symbol", n_cores=n_cores)
    auc_cols = [c for c in scored.columns if c.endswith("_AUCell")]
    for c in auc_cols:
        vals = scored[c].to_numpy(dtype=float)
        vals = vals[np.isfinite(vals)]
        assert vals.min() >= 0.0 and vals.max() <= 1.0, f"{c} out of [0,1]: {vals.min()}..{vals.max()}"
        print(f"[08_harvest] {c} in [{vals.min():.4f}, {vals.max():.4f}] (n={len(vals)}) OK")

    # --- carry already-derived readouts from the frozen explorers (never re-derived) ---
    carried = carry_frozen_readouts(adata.obs_names.astype(str))
    for col in carried.columns:
        print(f"[08_harvest] carried {col}: NaN frac {carried[col].isna().mean():.4f}")

    # --- per-cell parquet (gitignored/regenerable) ---
    df = build_parquet(adata, scored, carried)
    pq_path = PATHS.interactive_dir() / "08_harvest_readout.parquet"
    df.to_parquet(pq_path, index=False)
    print(f"[08_harvest] wrote {pq_path.name} ({df.shape[0]} cells, {df.shape[1]} cols)")

    # --- summary (per coarse_label × tissue × readout) ---
    summary = summarise(df)
    sdir = PATHS.tables(STAGE)
    summary.to_csv(sdir / "harvest_readout_summary.csv", index=False)
    print(f"[08_harvest] wrote harvest_readout_summary.csv ({summary.shape[0]} rows)")
    print("\n[08_harvest] summary (mean, frac>P90, SF-PB shift):\n",
          summary[["coarse_label", "tissue", "readout", "n_cells", "mean",
                   "frac_above_p90", "sf_minus_pb_mean", "sf_minus_pb_smd"]]
          .to_string(index=False))
    print("\n[08_harvest] DONE. Exploratory / annotation tier — never pooled with the "
          "pseudobulk NES spine.")


if __name__ == "__main__":
    main()
