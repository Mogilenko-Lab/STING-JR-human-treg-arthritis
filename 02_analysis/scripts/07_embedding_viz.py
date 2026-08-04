#!/usr/bin/env python
"""
07_embedding_viz.py — VIZ (no statistics). Treg-compartment embedding atlas.
============================================================================
Renders the multi-hook HARVEST-DESIGN preview — "show the
various signatures we are drafting for and where they land on their own and
where they land under OR condition" — from the substrate written by
`07_embedding.py`. These are SUPPLEMENTARY / VISUALISATION figures: embeddings
are annotation/visualisation ONLY, never the statistical readout (umbrella
guardrail). The mouse `WT_heat` anchor score appears as an ANNOTATION overlay,
NOT as a selection gate. No cells are lassoed/subset here.

Figures (each PDF+PNG + same-stem source table + README caption, via save_overview).
They render in the canonical drafted-subset order — annotation overview →
quad-marker context → signatures we draft for → what we end up with:
  umap_annotation_treg  — annotation-overview UMAP coloured by the frozen coarse_label
                          (Treg/Tcon/CD8), with a tissue (SF/PB) companion panel.
                          Establishes what is on the map first.
  umap_quadmarkers_treg — 2x2 context patch: pct_counts_mt + three canonical POI genes
                          for JIA sorted Tregs (FOXP3 / CTLA4 / IKZF2 — core Treg
                          identity + suppressive/activation axis), continuous colormaps.
  umap_signatures_treg  — each candidate signature landing INDIVIDUALLY on the UMAP
                          (WT_heat_up anchor annotation, score_eTreg effector,
                          score_HSP stress) + coarse_label reference. What we draft for.
  umap_drafted_treg      — the FINAL drafted-subset view: cells the multi-hook OR-union
                          rules would harvest, highlighted vs background. What we end up with.
  umap_or_union_treg    — supporting hook-composition view: cells coloured by which
                          hook-factor(s) they satisfy, with matched heat-lo / effector-lo
                          baselines (the "where they land together" panel).
  umap_markers_treg     — POI lineage markers (FOXP3/IL2RA/CTLA4/IKZF2/CD8A/IL7R)
                          + coarse_label categorical reference (retained artifact).
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "02_analysis"))
os.chdir(ROOT)

from config import PATHS, POPULATION_COLORS, TISSUE_COLORS  # noqa: E402
from helpers.figure_style import set_paper_style, save_overview, FIG_CFG  # noqa: E402

STAGE = "07_embedding"
SCRIPT = "02_analysis/scripts/07_embedding_viz.py"

# Okabe-Ito categorical palette (colourblind-safe), read from analysis_config.yaml.
# Populations come from the one shared palette at colors.populations.
LINEAGE_COL = POPULATION_COLORS
TISSUE_COL = TISSUE_COLORS
# Three canonical POI/disease genes for the quad-marker context patch (JIA sorted
# Tregs): FOXP3 (master Treg-lineage TF), CTLA4 (suppressive effector molecule),
# IKZF2/Helios (stable/tTreg identity + activation). Core Treg identity/activation
# axis — the most-motivated trio for this compartment; documented in the caption.
QUAD_GENES = ["FOXP3", "CTLA4", "IKZF2"]
MEMBERSHIP_COL = {
    "lineage only": "#009E73",
    "lineage + effector": "#0072B2",
    "mt-hi viable pocket": "#D55E00",
    "effector only": "#E69F00",
    "baseline (not in union)": "#D9D9D9",
}
BASELINE_COL = {
    "in union": "#B0B0B0",
    "heat-lo baseline": "#56B4E9",
    "effector-lo baseline": "#CC79A7",
    "heat-lo & effector-lo": "#0072B2",
    "other (not in union)": "#ECECEC",
}
# Drafted-subset view: the concrete "these are the cells our OR-union rules would
# harvest" highlight (warm vermillion) vs light-grey non-drafted background.
DRAFTED_COL = {
    "not drafted (background)": "#D9D9D9",
    "drafted (rules would harvest)": "#D55E00",
}
PT = 2.0            # scatter point size (dense embedding)
CMAP = "viridis"    # perceptually-uniform continuous colormap


def _scatter_cont(ax, d, col, title):
    """Continuous overlay on the UMAP, robustly clipped to 2-98th percentile."""
    v = d[col].to_numpy(dtype=float)
    lo, hi = np.nanpercentile(v, [2, 98])
    order = np.argsort(v)  # draw high-value cells on top
    sc = ax.scatter(d["x"].to_numpy()[order], d["y"].to_numpy()[order],
                    c=v[order], s=PT, cmap=CMAP, vmin=lo, vmax=hi,
                    linewidths=0, rasterized=True)
    ax.set_title(title)
    ax.set_xticks([]); ax.set_yticks([])
    cb = ax.figure.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)
    cb.ax.tick_params(labelsize=8)


def _scatter_cat(ax, d, col, palette, title, order=None, legend=True):
    """Categorical overlay; grey/baseline levels drawn first so signal sits on top."""
    cats = order or list(palette.keys())
    present = [c for c in cats if c in set(d[col].unique())]
    for c in present:
        m = (d[col] == c).to_numpy()
        ax.scatter(d["x"].to_numpy()[m], d["y"].to_numpy()[m], s=PT,
                   color=palette[c], linewidths=0, rasterized=True, label=c)
    ax.set_title(title)
    ax.set_xticks([]); ax.set_yticks([])
    if legend:
        handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor=palette[c],
                          markeredgecolor="none", markersize=7, label=c) for c in present]
        ax.legend(handles=handles, frameon=True, fontsize=8, loc="best",
                  markerscale=1.2, handletextpad=0.3)


def main() -> None:
    set_paper_style(config=FIG_CFG)
    tdir = PATHS.tables(STAGE)
    d = pd.read_parquet(tdir / "hook_factor_substrate.parquet")

    # ================================================================== #
    # 1. umap_annotation_treg — annotation-overview (what's on the map)    #
    #    coarse_label (Treg/Tcon/CD8) + a clean tissue SF/PB companion.     #
    # ================================================================== #
    fig, axes = plt.subplots(1, 2, figsize=(15, 7))
    _scatter_cat(axes[0], d, "coarse_label", LINEAGE_COL,
                 "Frozen cell-type annotation (coarse_label)",
                 order=["Treg", "Tcon", "CD8"])
    _scatter_cat(axes[1], d, "tissue", TISSUE_COL,
                 "Tissue of origin (SF / PB companion)",
                 order=["synovial_fluid", "peripheral_blood"])
    fig.suptitle("Annotation overview — the frozen sorted T-cell map (Treg / Tcon / CD8)",
                 fontsize=16)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    # Source table: per-annotation cell composition by tissue (describes the figure).
    anno_tbl = (d.groupby(["coarse_label", "tissue"], observed=True).size()
                .rename("n_cells").reset_index())
    anno_tbl["frac_of_total"] = anno_tbl["n_cells"] / len(d)
    save_overview(
        fig, STAGE, "umap_annotation_treg",
        table=anno_tbl,
        finding=("Annotation-overview UMAP: every cell coloured by its frozen sort-lineage "
                 "annotation (coarse_label = Treg / Tcon / CD8), with a tissue (synovial-fluid / "
                 "peripheral-blood) companion. Establishes WHAT IS ON THE MAP before any "
                 "signature/hook overlay is read against it. VISUALISATION only — the frozen "
                 "annotation, not a statistical readout."),
        script=SCRIPT, fn="main",
        config_kv=("colour=coarse_label (Okabe-Ito Treg/Tcon/CD8); companion=tissue (SF/PB); "
                   "categorical scatter, baseline levels drawn first"),
        input="03_results/07_embedding/tables/hook_factor_substrate.parquet",
        how_to_read=("Left: the frozen sorted-lineage annotation the whole atlas is built on "
                     "(Treg green, Tcon orange, CD8 pink). Right: the same cells coloured by tissue "
                     "of origin (synovial fluid vs peripheral blood). Read as the reference layout "
                     "every downstream signature and hook overlay sits on top of. Source table = "
                     "per-annotation x tissue cell counts. Correlative, annotation only."),
        width=15, height=7, config=FIG_CFG)

    # ================================================================== #
    # 2. umap_quadmarkers_treg — 2x2 context patch (POI/disease markers)   #
    #    pct_counts_mt + FOXP3 / CTLA4 / IKZF2 (core Treg identity axis).   #
    # ================================================================== #
    fig, axes = plt.subplots(2, 2, figsize=(13, 11))
    _scatter_cont(axes[0, 0], d, "pct_counts_mt", "Mitochondrial fraction (pct_counts_mt)")
    _scatter_cont(axes[0, 1], d, QUAD_GENES[0], f"{QUAD_GENES[0]} (Treg-lineage TF)")
    _scatter_cont(axes[1, 0], d, QUAD_GENES[1], f"{QUAD_GENES[1]} (suppressive effector)")
    _scatter_cont(axes[1, 1], d, QUAD_GENES[2], f"{QUAD_GENES[2]} / Helios (stable-Treg identity)")
    fig.suptitle("Context markers on the sorted T-cell UMAP — %mt + core Treg identity/activation genes"
                 " (FOXP3 / CTLA4 / IKZF2)", fontsize=15)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    # Source table: assembled from the compute stage's precomputed per-lineage CSVs.
    sig_pl = pd.read_csv(tdir / "signatures_per_lineage.csv")
    mk_pl = pd.read_csv(tdir / "markers_per_lineage.csv")
    mt_rows = (sig_pl[["lineage", "median_pct_counts_mt"]]
               .rename(columns={"median_pct_counts_mt": "median_value"}))
    mt_rows.insert(1, "channel", "pct_counts_mt")
    mt_rows["frac_expressing"] = np.nan
    gene_rows = (mk_pl[mk_pl["gene"].isin(QUAD_GENES)]
                 [["lineage", "gene", "median_expr", "frac_expressing"]]
                 .rename(columns={"gene": "channel", "median_expr": "median_value"}))
    quad_tbl = pd.concat([mt_rows, gene_rows], ignore_index=True)
    save_overview(
        fig, STAGE, "umap_quadmarkers_treg",
        table=quad_tbl,
        finding=("Quadruple context patch: mitochondrial fraction plus the three canonical "
                 "person-of-interest / disease genes for JIA sorted Tregs — FOXP3 (master "
                 "Treg-lineage transcription factor), CTLA4 (suppressive effector molecule) and "
                 "IKZF2/Helios (stable/thymic-Treg identity marker). Chosen as the core Treg "
                 "identity + suppressive-activation axis, the most-motivated trio for this "
                 "compartment; %mt anchors QC/viability context. VISUALISATION only."),
        script=SCRIPT, fn="main",
        config_kv=(f"channels: pct_counts_mt + {', '.join(QUAD_GENES)}; continuous viridis, clipped "
                   "2-98th pct; genes = core Treg identity/activation (FOXP3 TF, CTLA4 effector, "
                   "IKZF2/Helios stable-Treg)"),
        input="03_results/07_embedding/tables/hook_factor_substrate.parquet",
        how_to_read=("Four continuous panels (viridis, clipped 2-98th pct, high on top): top-left is "
                     "mitochondrial fraction (QC/viability context); the other three are the canonical "
                     "Treg genes — FOXP3 and CTLA4 concentrate on the Treg gate, IKZF2/Helios marks the "
                     "stable/thymic-Treg fraction. Read as the marker context the drafted subset is "
                     "designed against. Source table = per-lineage median values (+ frac expressing for "
                     "genes) from the compute stage. Correlative, annotation only."),
        width=13, height=11, config=FIG_CFG)

    # ================================================================== #
    # 3. umap_signatures_treg — each candidate signature landing alone     #
    # ================================================================== #
    fig, axes = plt.subplots(2, 2, figsize=(13, 11))
    _scatter_cont(axes[0, 0], d, "WT_heat_up",
                  "Mouse WT_heat_up (anchor annotation)")
    _scatter_cont(axes[0, 1], d, "score_eTreg", "Effector-Treg score (score_eTreg)")
    _scatter_cont(axes[1, 0], d, "score_HSP", "Heat-shock / stress score (score_HSP)")
    _scatter_cat(axes[1, 1], d, "coarse_label", LINEAGE_COL,
                 "Sort lineage (reference)", order=["Treg", "Tcon", "CD8"])
    fig.suptitle("Candidate harvest signatures on the sorted T-cell UMAP — each landing individually",
                 fontsize=16)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    save_overview(
        fig, STAGE, "umap_signatures_treg",
        table=pd.read_csv(tdir / "signatures_per_lineage.csv"),
        finding=("Previews the multi-hook harvest design: where each candidate signature we are "
                 "drafting for (mouse WT_heat anchor, effector-Treg, heat-shock/stress) lands ON ITS "
                 "OWN across the sorted Treg/Tcon/CD8 embedding, next to the frozen sort-lineage "
                 "reference. VISUALISATION only — not the statistical evidence."),
        script=SCRIPT, fn="main",
        config_kv="signatures: WT_heat_up (anchor, annotation-only), score_eTreg, score_HSP; cmap=viridis",
        input="03_results/07_embedding/tables/hook_factor_substrate.parquet",
        how_to_read=("Three continuous panels colour every cell by one candidate signature (viridis, "
                     "clipped 2-98th pct; high on top); the fourth shows the frozen sort lineage. Read "
                     "as WHERE each signature concentrates on the embedding, alone. The WT_heat anchor "
                     "is shown as an ANNOTATION, never a selection gate; the source table gives per-"
                     "lineage medians. Correlative preview of the harvest design, not an effect size."),
        width=13, height=11, config=FIG_CFG)

    # ================================================================== #
    # 4. umap_or_union_treg — the OR-union money view (+ matched baselines) #
    # ================================================================== #
    fig, axes = plt.subplots(1, 2, figsize=(15, 7))
    _scatter_cat(axes[0], d, "hook_membership", MEMBERSHIP_COL,
                 "OR-union membership: which hook(s) admit a cell",
                 order=["baseline (not in union)", "lineage only", "lineage + effector",
                        "effector only", "mt-hi viable pocket"])
    _scatter_cat(axes[1], d, "baseline_map", BASELINE_COL,
                 "Matched lo baselines (heat-lo / effector-lo)",
                 order=["other (not in union)", "in union", "heat-lo baseline",
                        "effector-lo baseline", "heat-lo & effector-lo"])
    union_frac = float(d["hook_or_union"].mean())
    fig.suptitle(f"OR-union harvest preview — bounded to {union_frac*100:.0f}% of cells "
                 "(lineage OR effector OR mt-hi viable); anchor score NOT a disjunct",
                 fontsize=15)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    save_overview(
        fig, STAGE, "umap_or_union_treg",
        table=pd.read_csv(tdir / "or_union_membership.csv"),
        finding=("The 'where they land together' panel: cells coloured by which anchor-orthogonal "
                 "hook-factor(s) admit them under the bounded OR-union (lineage / effector / mt-hi "
                 "viable / multiple), with the matched heat-lo and effector-lo baseline regions shown. "
                 "Previews the harvest design's factorial contrastability — VISUALISATION, not "
                 "selection-to-file and not a statistical readout."),
        script=SCRIPT, fn="main",
        config_kv=("OR-union = hook_lineage | hook_effector | hook_mthi_viable; anchor score NEVER a "
                   f"disjunct; union = {union_frac*100:.1f}% of all cells (bounded minority)"),
        input="03_results/07_embedding/tables/hook_factor_substrate.parquet",
        how_to_read=("Left: each cell coloured by the hook(s) it satisfies (grey = not in union). Right: "
                     "the matched-lo baselines (heat-lo, effector-lo) that give every factorial contrast a "
                     "defined negative arm. The union is a BOUNDED minority of cells (fraction in the "
                     "title + source table), so an OR sweep over these hooks does not take in most of "
                     "the dataset. No cells are lassoed/subset; harvest selection is deferred. "
                     "Correlative."),
        width=15, height=7, config=FIG_CFG)

    # ================================================================== #
    # 5. umap_markers_treg — POI lineage markers + lineage reference        #
    # ================================================================== #
    markers = ["FOXP3", "IL2RA", "CTLA4", "IKZF2", "CD8A", "IL7R"]
    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    axf = axes.ravel()
    for ax, gene in zip(axf, markers):
        _scatter_cont(ax, d, gene, gene)
    _scatter_cat(axf[6], d, "coarse_label", LINEAGE_COL,
                 "Sort lineage (reference)", order=["Treg", "Tcon", "CD8"])
    axf[7].axis("off")
    fig.suptitle("POI lineage markers on the sorted T-cell UMAP (curated, anchor-orthogonal hooks)",
                 fontsize=16)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    save_overview(
        fig, STAGE, "umap_markers_treg",
        table=pd.read_csv(tdir / "markers_per_lineage.csv"),
        finding=("Curated POI lineage markers (FOXP3/IL2RA/CTLA4/IKZF2 Treg identity; CD8A; IL7R) that "
                 "define the anchor-orthogonal lineage/marker hooks, next to the frozen sort-lineage "
                 "reference — confirming the hooks track real lineage biology, not the anchor score."),
        script=SCRIPT, fn="main",
        config_kv="markers: FOXP3, IL2RA, CTLA4, IKZF2, CD8A, IL7R; cmap=viridis (log-norm expression)",
        input="03_results/07_embedding/tables/hook_factor_substrate.parquet",
        how_to_read=("Each panel colours cells by one marker's log-normalised expression (viridis, "
                     "clipped 2-98th pct); last panel is the frozen sort lineage. Read as: the Treg-"
                     "identity markers concentrate on the Treg gate, CD8A on CD8, IL7R on conventional "
                     "T — the curated hooks are lineage-faithful and independent of the anchor score. "
                     "Source table = per-lineage median expression + fraction expressing. Correlative."),
        width=20, height=10, config=FIG_CFG)

    # ================================================================== #
    # 6. umap_drafted_treg — the EXPLICIT drafted-subset view              #
    #    Primary panel: drafted (hook_or_union) highlight vs grey backdrop; #
    #    secondary panel: drafted cells faceted by which hook drafted them. #
    # ================================================================== #
    drafted = d["hook_or_union"].to_numpy(dtype=bool)
    n_drafted = int(drafted.sum())
    n_total = int(len(d))
    draft_frac = n_drafted / n_total
    d = d.assign(
        drafted_view=np.where(drafted, "drafted (rules would harvest)",
                              "not drafted (background)"),
        # which hook drafted a cell (drafted only); background greyed out.
        drafted_by_hook=np.where(drafted, d["hook_membership"].astype(str),
                                 "not drafted (background)"),
    )
    hook_facet_col = {**MEMBERSHIP_COL,
                      "not drafted (background)": "#D9D9D9"}

    fig, axes = plt.subplots(1, 2, figsize=(15, 7))
    _scatter_cat(axes[0], d, "drafted_view", DRAFTED_COL,
                 "Drafted subset: cells the OR-union rules would harvest",
                 order=["not drafted (background)", "drafted (rules would harvest)"])
    _scatter_cat(axes[1], d, "drafted_by_hook", hook_facet_col,
                 "Drafted cells, by which hook drafted them",
                 order=["not drafted (background)", "lineage only", "lineage + effector",
                        "effector only", "mt-hi viable pocket"])
    fig.suptitle(f"Drafted-subset preview — {n_drafted:,} / {n_total:,} cells "
                 f"({draft_frac*100:.0f}%) the multi-hook OR-union rules would draft "
                 "(design preview, not a committed cohort)", fontsize=15)
    fig.tight_layout(rect=(0, 0, 1, 0.97))

    # Source table: the drafted-vs-background split + per-hook drafted counts.
    draft_rows = pd.DataFrame({
        "group": ["drafted (hook_or_union)", "not drafted (background)",
                  "drafted: hook_lineage", "drafted: hook_effector",
                  "drafted: hook_mthi_viable"],
        "n_cells": [n_drafted, n_total - n_drafted,
                    int(d["hook_lineage"].sum()), int(d["hook_effector"].sum()),
                    int(d["hook_mthi_viable"].sum())],
        "frac_of_total": [draft_frac, 1 - draft_frac,
                          float(d["hook_lineage"].mean()), float(d["hook_effector"].mean()),
                          float(d["hook_mthi_viable"].mean())],
    })
    save_overview(
        fig, STAGE, "umap_drafted_treg",
        table=draft_rows,
        finding=("The explicit 'these are the cells our rules would harvest' view: the Treg-"
                 "compartment UMAP with the drafted subset (hook_or_union == TRUE, "
                 f"{n_drafted:,} cells / {draft_frac*100:.0f}%) highlighted against the light-grey "
                 "non-drafted background, plus a facet of the drafted cells by which hook drafted "
                 "them. The subset the multi-hook rules would draft (design preview / visualization, "
                 "not a committed cohort or statistical evidence); anchor heat score is annotation, "
                 "never a selection gate."),
        script=SCRIPT, fn="main",
        config_kv=(f"drafted = hook_or_union (lineage|effector|mt-hi viable); {n_drafted} / {n_total} "
                   f"cells ({draft_frac*100:.1f}%); highlight={DRAFTED_COL['drafted (rules would harvest)']}, "
                   "background=#D9D9D9"),
        input="03_results/07_embedding/tables/hook_factor_substrate.parquet",
        how_to_read=("Left (primary): every cell coloured by whether the OR-union harvest rules would "
                     "draft it (vermillion = drafted, light grey = background) — the concrete subset "
                     "preview. Right: the drafted cells split by which hook admitted them (lineage / "
                     "effector / mt-hi viable / multiple); non-drafted cells stay grey. This is a "
                     "VISUALISATION of the harvest design, not a committed cohort and not a statistical "
                     "readout; no cells are lassoed/subset to file and the anchor heat score is never a "
                     "selection gate. Source table = drafted-vs-background split + per-hook counts. "
                     "Correlative."),
        width=15, height=7, config=FIG_CFG)

    print("[07_embedding_viz] wrote 6 overviews (annotation, quadmarkers, signatures, "
          "OR-union, markers, drafted)")


if __name__ == "__main__":
    main()
