#!/usr/bin/env python
"""
07_embedding.py — COMPUTE (no plotting). Treg-compartment embedding atlas.
==========================================================================
Builds the per-cell HOOK-FACTOR substrate that the `07_embedding_viz.py`
figures render, implementing the bounded-OR harvest design.

WHAT THIS IS (and is NOT): this is a VISUALISATION substrate for the multi-hook
harvest-design preview — "show the various signatures we are
drafting for and where they land on their own and where they land under OR
condition". It is NOT the statistical readout and NOT an anchor-score-gated
selection. Per the umbrella guardrails, embeddings are annotation/visualisation
ONLY; the mouse `WT_heat` anchor score is carried here as an ANNOTATION, never a
selection predicate (the anchor-orthogonal selection safeguard). No cells are
lassoed/subset to a file here — harvest selection design is deferred.

Substrate: the frozen explorer parquets (per-cell; NOT the multi-GB h5ad):
  03_results/interactive/01_qc_explore.parquet      (richest: coords + scores + markers)
  03_results/interactive/02_annotation_explore.parquet  (frozen coarse_label)

HOOK FACTORS (booleans, kept SEPARATE so downstream contrasts stay clean):
  hook_lineage       — sorted Treg identity (FACS gate; the strongest orthogonal
                       hook, fully independent of the anchor score).
  hook_effector      — effector-state high: score_eTreg >= dataset-internal top
                       decile (P90). A DEFINED MINORITY, admits effector cells
                       beyond the sort gate; never the whole set.
  hook_mthi_viable   — the mt-hi-but-VIABLE Treg pocket (cluster-6-like), mirroring
                       01_qc_mthi_characterize.py rule B + an explicit viability
                       gate: Treg AND pct_counts_mt >= within-Treg P97.5 AND
                       score_eTreg >= within-Treg median (excludes the eTreg-LOW
                       cluster-16 junk) AND n_genes_by_counts >= 2x the 200-gene
                       floor (real cells, not debris).
Annotation factors (NEVER selection predicates; safeguard S2):
  anno_heat_hi/_lo   — WT_heat_up top/bottom decile (matched hi/lo baseline).
  anno_stingspecific — placeholder: the SAVI STING-specific axis is not joined in
                       this compartment yet (see NOTE below); column is left NA.
OR-union: hook_or_union = hook_lineage | hook_effector | hook_mthi_viable.
  Matched-lo baselines (safeguard S3): baseline_heat_lo (anno_heat_lo & ~union),
  baseline_effector_lo (score_eTreg <= P10 & ~union). The union MUST stay a
  bounded, contrastable minority — reported as a fraction of all cells.

Outputs (03_results/07_embedding/tables/):
  hook_factor_definitions.csv   — one row per factor: kind, definition, threshold, n, frac
  or_union_membership.csv       — union membership categories: n, frac_all, frac_of_union
  hook_per_lineage_summary.csv  — per sort-lineage: factor counts + score medians
  signatures_per_lineage.csv    — per lineage: median/mean of each candidate signature (fig source)
  markers_per_lineage.csv       — per lineage: median expr + frac expressing per POI marker (fig source)
  hook_factor_substrate.parquet — per-cell substrate (coords, labels, booleans, scores, markers) for viz
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "02_analysis"))
os.chdir(ROOT)

from config import PATHS, PARAMS  # noqa: E402

STAGE = "07_embedding"
SCRIPT = "02_analysis/scripts/07_embedding.py"

QC_PARQUET = PATHS.interactive / "01_qc_explore.parquet"
ANNO_PARQUET = PATHS.interactive / "02_annotation_explore.parquet"

# --- threshold constants (dataset-internal; reported in hook_factor_definitions) ---
EFFECTOR_Q = 0.90          # effector-HIGH = top decile of score_eTreg across all cells
MT_PCTILE = 97.5           # within-Treg %mt percentile that separates the mt-hi clusters {6,16}
HEAT_HI_Q = 0.90           # anno_heat_hi = top decile of WT_heat_up (annotation only)
HEAT_LO_Q = 0.10           # anno_heat_lo = bottom decile of WT_heat_up (matched lo baseline)
EFFECTOR_LO_Q = 0.10       # effector-lo matched baseline

# Candidate signatures we are "drafting for" (continuous annotations shown on the UMAP).
SIGNATURE_COLS = ["WT_heat_up", "score_eTreg", "score_HSP"]
# POI lineage markers (curated; orthogonal to the anchor score).
MARKER_GENES = ["FOXP3", "IL2RA", "CTLA4", "IKZF2", "CD8A", "IL7R"]


def main() -> None:
    tdir = PATHS.tables(STAGE)

    # ---- load frozen explorer parquets (NOT the h5ad) --------------------- #
    qc = pd.read_parquet(QC_PARQUET)
    anno = pd.read_parquet(ANNO_PARQUET)
    # coarse_label is the frozen annotation label; index is aligned across explorers.
    df = qc.copy()
    df["coarse_label"] = anno.reindex(df.index)["coarse_label"]
    n_all = len(df)

    treg = df["population_short"] == "Treg"

    # ---- selection HOOKS (anchor-orthogonal) ----------------------------- #
    # hook_lineage — sorted Treg identity (strongest orthogonal hook).
    df["hook_lineage"] = treg.to_numpy()

    # hook_effector — effector-state high (top decile of score_eTreg, all cells).
    eff_thr = float(df["score_eTreg"].quantile(EFFECTOR_Q))
    df["hook_effector"] = (df["score_eTreg"] >= eff_thr).to_numpy()

    # hook_mthi_viable — mt-hi-but-viable Treg pocket (mirror mthi rule B + viability gate).
    tg = df[treg]
    mt_thr = float(np.percentile(tg["pct_counts_mt"], MT_PCTILE))          # ~10% (splits clusters 6/16)
    etreg_med_treg = float(tg["score_eTreg"].median())                     # excludes eTreg-LOW cluster 16
    viability_floor = 2 * int(PARAMS.qc_min_genes)                         # n_genes "well above" the 200 floor
    df["hook_mthi_viable"] = (
        treg
        & (df["pct_counts_mt"] >= mt_thr)
        & (df["score_eTreg"] >= etreg_med_treg)
        & (df["n_genes_by_counts"] >= viability_floor)
    ).to_numpy()

    # ---- ANNOTATION factors (NEVER selection predicates; safeguard S2) ---- #
    heat_hi_thr = float(df["WT_heat_up"].quantile(HEAT_HI_Q))
    heat_lo_thr = float(df["WT_heat_up"].quantile(HEAT_LO_Q))
    df["anno_heat_hi"] = (df["WT_heat_up"] >= heat_hi_thr).to_numpy()
    df["anno_heat_lo"] = (df["WT_heat_up"] <= heat_lo_thr).to_numpy()
    # anno_stingspecific — placeholder: SAVI STING-specific axis not joined in this compartment yet.
    df["anno_stingspecific"] = pd.Series(pd.NA, index=df.index, dtype="object")

    # ---- OR-union + matched-lo baselines (safeguard S3) ------------------ #
    df["hook_or_union"] = (
        df["hook_lineage"] | df["hook_effector"] | df["hook_mthi_viable"]
    ).to_numpy()
    eff_lo_thr = float(df["score_eTreg"].quantile(EFFECTOR_LO_Q))
    df["effector_lo"] = (df["score_eTreg"] <= eff_lo_thr).to_numpy()
    df["baseline_heat_lo"] = (df["anno_heat_lo"] & ~df["hook_or_union"]).to_numpy()
    df["baseline_effector_lo"] = (df["effector_lo"] & ~df["hook_or_union"]).to_numpy()

    # ---- membership category (which hook(s) a cell satisfies; for the union view) ---- #
    def _membership(r) -> str:
        if not r["hook_or_union"]:
            return "baseline (not in union)"
        if r["hook_mthi_viable"]:
            return "mt-hi viable pocket"       # M is a subset of L; highlight distinctly
        if r["hook_lineage"] and r["hook_effector"]:
            return "lineage + effector"
        if r["hook_lineage"]:
            return "lineage only"
        return "effector only"                 # effector-high non-Treg
    df["hook_membership"] = df.apply(_membership, axis=1)

    # matched-lo baseline map (for the baseline panel of the union view).
    def _baseline(r) -> str:
        if r["hook_or_union"]:
            return "in union"
        hl, el = bool(r["anno_heat_lo"]), bool(r["effector_lo"])
        if hl and el:
            return "heat-lo & effector-lo"
        if hl:
            return "heat-lo baseline"
        if el:
            return "effector-lo baseline"
        return "other (not in union)"
    df["baseline_map"] = df.apply(_baseline, axis=1)

    # ---- TABLE: hook-factor definitions (kind, definition, threshold, n, frac) ---- #
    def _n(col: str) -> int:
        return int(df[col].sum())
    defs = pd.DataFrame([
        {"factor": "hook_lineage", "kind": "selection_hook",
         "definition": "sorted Treg identity (FACS gate; population_short == 'Treg')",
         "threshold": "population_short=='Treg'", "n_cells": _n("hook_lineage")},
        {"factor": "hook_effector", "kind": "selection_hook",
         "definition": f"effector-high: score_eTreg >= dataset P{int(EFFECTOR_Q*100)} (top decile)",
         "threshold": f"score_eTreg>={eff_thr:.4f}", "n_cells": _n("hook_effector")},
        {"factor": "hook_mthi_viable", "kind": "selection_hook",
         "definition": (f"mt-hi VIABLE Treg pocket (cluster-6-like): Treg AND pct_counts_mt>=within-Treg "
                        f"P{MT_PCTILE} AND score_eTreg>=within-Treg median (excludes eTreg-LOW cluster 16) "
                        f"AND n_genes_by_counts>=2x qc_min_genes floor (viability gate)"),
         "threshold": f"mt>={mt_thr:.2f} & eTreg>={etreg_med_treg:.4f} & n_genes>={viability_floor}",
         "n_cells": _n("hook_mthi_viable")},
        {"factor": "hook_or_union", "kind": "or_union",
         "definition": "hook_lineage OR hook_effector OR hook_mthi_viable (anchor score NEVER a disjunct)",
         "threshold": "L|E|M", "n_cells": _n("hook_or_union")},
        {"factor": "anno_heat_hi", "kind": "annotation_only",
         "definition": f"WT_heat_up >= dataset P{int(HEAT_HI_Q*100)} (anchor annotation; NOT a selection gate)",
         "threshold": f"WT_heat_up>={heat_hi_thr:.4f}", "n_cells": _n("anno_heat_hi")},
        {"factor": "anno_heat_lo", "kind": "annotation_baseline",
         "definition": f"WT_heat_up <= dataset P{int(HEAT_LO_Q*100)} (matched heat-lo baseline)",
         "threshold": f"WT_heat_up<={heat_lo_thr:.4f}", "n_cells": _n("anno_heat_lo")},
        {"factor": "anno_stingspecific", "kind": "annotation_placeholder",
         "definition": "SAVI STING-specific axis NOT joined in this compartment yet (placeholder, NA)",
         "threshold": "n/a", "n_cells": 0},
        {"factor": "baseline_heat_lo", "kind": "matched_baseline",
         "definition": "anno_heat_lo AND NOT in OR-union (matched lo baseline for factorial contrasts)",
         "threshold": "anno_heat_lo & ~union", "n_cells": _n("baseline_heat_lo")},
        {"factor": "baseline_effector_lo", "kind": "matched_baseline",
         "definition": f"score_eTreg <= P{int(EFFECTOR_LO_Q*100)} AND NOT in OR-union (matched effector-lo baseline)",
         "threshold": f"score_eTreg<={eff_lo_thr:.4f} & ~union", "n_cells": _n("baseline_effector_lo")},
    ])
    defs["frac_all_cells"] = defs["n_cells"] / n_all
    defs["n_all_cells"] = n_all
    defs["note"] = "VISUALISATION substrate; anchor score is annotation-only, never a selection predicate."
    defs.to_csv(tdir / "hook_factor_definitions.csv", index=False)

    # ---- TABLE: OR-union membership breakdown ----------------------------- #
    n_union = _n("hook_or_union")
    mem = (df["hook_membership"].value_counts().rename_axis("membership_category")
           .reset_index(name="n_cells"))
    mem["frac_all_cells"] = mem["n_cells"] / n_all
    mem["frac_of_union"] = np.where(
        mem["membership_category"] == "baseline (not in union)", np.nan,
        mem["n_cells"] / n_union)
    mem["n_all_cells"] = n_all
    mem["n_union"] = n_union
    mem["union_frac_all_cells"] = n_union / n_all
    mem = mem.sort_values("n_cells", ascending=False).reset_index(drop=True)
    mem.to_csv(tdir / "or_union_membership.csv", index=False)

    # ---- TABLE: per-lineage hook + score summary -------------------------- #
    rows = []
    for lab, g in df.groupby("coarse_label", observed=True):
        rec = {"lineage": lab, "n_cells": len(g)}
        for col in ["hook_lineage", "hook_effector", "hook_mthi_viable",
                    "hook_or_union", "anno_heat_hi", "anno_heat_lo"]:
            rec[f"n_{col}"] = int(g[col].sum())
            rec[f"frac_{col}"] = float(g[col].mean())
        for s in SIGNATURE_COLS + ["pct_counts_mt", "n_genes_by_counts"]:
            rec[f"median_{s}"] = float(g[s].median())
        rows.append(rec)
    perlin = pd.DataFrame(rows).sort_values("n_cells", ascending=False)
    perlin.to_csv(tdir / "hook_per_lineage_summary.csv", index=False)

    # ---- TABLE: signatures per lineage (fig source for umap_signatures) --- #
    rows = []
    for lab, g in df.groupby("coarse_label", observed=True):
        rec = {"lineage": lab, "n_cells": len(g)}
        for s in SIGNATURE_COLS + ["pct_counts_mt"]:
            rec[f"median_{s}"] = float(g[s].median())
            rec[f"mean_{s}"] = float(g[s].mean())
        rows.append(rec)
    pd.DataFrame(rows).sort_values("n_cells", ascending=False).to_csv(
        tdir / "signatures_per_lineage.csv", index=False)

    # ---- TABLE: markers per lineage (fig source for umap_markers) --------- #
    rows = []
    for lab, g in df.groupby("coarse_label", observed=True):
        for gene in MARKER_GENES:
            rows.append({
                "lineage": lab, "gene": gene, "n_cells": len(g),
                "median_expr": float(g[gene].median()),
                "mean_expr": float(g[gene].mean()),
                "frac_expressing": float((g[gene] > 0).mean()),
            })
    pd.DataFrame(rows).to_csv(tdir / "markers_per_lineage.csv", index=False)

    # ---- per-cell substrate parquet (viz input; NOT a figure source table) ---- #
    keep = (["x", "y", "population_short", "coarse_label", "leiden_unsupervised",
             "tissue", "donor", "pct_counts_mt", "n_genes_by_counts",
             "hook_lineage", "hook_effector", "hook_mthi_viable", "hook_or_union",
             "anno_heat_hi", "anno_heat_lo", "anno_stingspecific",
             "baseline_heat_lo", "baseline_effector_lo",
             "hook_membership", "baseline_map"]
            + SIGNATURE_COLS + MARKER_GENES)
    df[keep].to_parquet(tdir / "hook_factor_substrate.parquet")

    # ---- console summary -------------------------------------------------- #
    frac = n_union / n_all
    print(f"[07_embedding] cells={n_all}; hook_lineage={_n('hook_lineage')} "
          f"({_n('hook_lineage')/n_all*100:.1f}%), hook_effector={_n('hook_effector')} "
          f"({_n('hook_effector')/n_all*100:.1f}%), hook_mthi_viable={_n('hook_mthi_viable')}")
    print(f"[07_embedding] OR-UNION = {n_union} cells = {frac*100:.1f}% of all cells "
          f"({'BOUNDED — a contrastable minority' if frac < 0.6 else 'WARNING: approaching whole dataset — TIGHTEN'})")
    print(f"[07_embedding] baselines: heat_lo={_n('baseline_heat_lo')}, "
          f"effector_lo={_n('baseline_effector_lo')}")
    print("[07_embedding] wrote 5 tables + substrate parquet to", tdir)


if __name__ == "__main__":
    main()
