#!/usr/bin/env python
"""
12_treg_localisation_viz.py -- VIZ ONLY (no computing).
=============================================================================
Draws the Treg-only niche localisation figure showing per-cell AUCell score
distributions across synovial fluid versus peripheral blood for five signatures:
1. WT_heat_up (mouse 39 °C-derived up arm)
2. Interaction_up (cGAS-dependent up arm, 7-gene gate)
3. Interaction_fdrOnly_up (cGAS-dependent up arm, 18-gene gate)
4. HALLMARK_HYPOXIA (MSigDB Hallmark hypoxia program)
5. WT_heat_up_purged_hypoxia (mouse 39 °C-derived up arm purged of hypoxia genes)

Tier discipline:
- Headline is a question / localisation statement (never an answer).
- Declared on face: HYPOTHESIS-GENERATING TIER · per-cell AUCell score distributions.
- Cites donor-level pseudobulk differential expression as the confirmatory tier.
- Explains why Treg subsetting is legitimate (GSE160097 is FACS-sorted; coarse_label == 'Treg'
  is a sort gate, not a score-derived selection).
- Set sizes and power bands declared on face. Plain statement that the 7-gene gate
  is dominated by detection noise and not interpretable at per-cell resolution.

Outputs:
  03_results/12_treg_localisation/figures/_overview/treg_localisation.pdf
  03_results/12_treg_localisation/figures/_overview/treg_localisation.png
  03_results/12_treg_localisation/tables/_overview/treg_localisation.csv
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

COMPARTMENT_ROOT = Path(__file__).resolve().parents[2]
REPO_ROOT = COMPARTMENT_ROOT.parent

sys.path.insert(0, str(COMPARTMENT_ROOT))
sys.path.insert(0, str(COMPARTMENT_ROOT / "02_analysis"))
os.chdir(COMPARTMENT_ROOT)

from config import (CONFIG, PATHS, PARAMS, TISSUE_KEY, DONOR_KEY, TISSUE_NUM,  # noqa: E402
                    TISSUE_DEN)
from helpers.figure_style import FIG_CFG, save_overview, set_paper_style  # noqa: E402

STAGE = "12_treg_localisation"
SCRIPT = "02_analysis/scripts/12_treg_localisation_viz.py"
STEM = "treg_localisation"

SIG_ORDER = [
    "WT_heat_up",
    "Interaction_up",
    "Interaction_fdrOnly_up",
    "HALLMARK_HYPOXIA",
    "WT_heat_up_purged_hypoxia",
]

SIG_META = {
    "WT_heat_up": {
        "title": "WT_heat_up",
        "description": "mouse 39 °C-derived up arm",
        "pseudobulk_cite": "Primary pseudobulk: NES +2.18, FDR 3.2e-4",
    },
    "Interaction_up": {
        "title": "Interaction_up (7-gene gate)",
        "description": "cGAS-dependent up arm",
        "pseudobulk_cite": "7-gene gate: underpowered & dominated by detection noise",
    },
    "Interaction_fdrOnly_up": {
        "title": "Interaction_fdrOnly_up (18-gene gate)",
        "description": "cGAS-dependent up arm",
        "pseudobulk_cite": "Primary pseudobulk: FDR-gated cGAS-dependent arm",
    },
    "HALLMARK_HYPOXIA": {
        "title": "HALLMARK_HYPOXIA",
        "description": "MSigDB Hallmark hypoxia program",
        "pseudobulk_cite": "Primary pseudobulk: NES +2.05, FDR < 0.001 (08_harvest_readout)",
    },
    "WT_heat_up_purged_hypoxia": {
        "title": "WT_heat_up (purged of hypoxia)",
        "description": "mouse 39 °C-derived up arm minus hypoxia genes",
        "pseudobulk_cite": "Primary pseudobulk: NES +1.95, FDR 1.2e-3 (09_heat_hypoxia)",
    },
}


def build_figure(per_cell: pd.DataFrame, summary: pd.DataFrame) -> plt.Figure:
    """Draw a 5-panel small-multiple figure comparing per-cell AUCell score distributions across niches."""
    palette = FIG_CFG["colors"]
    color_sf = palette["okabe_ito"]["vermillion"]
    color_pb = palette["okabe_ito"]["sky_blue"]
    text_color = palette["okabe_ito"]["black"]

    tregs = per_cell[per_cell["coarse_label"] == "Treg"].copy()

    fig, axes = plt.subplots(1, 5, figsize=(21, 6.0), sharey=False)
    fig.subplots_adjust(wspace=0.35, top=0.74, bottom=0.15, left=0.06, right=0.98)

    for i, sig in enumerate(SIG_ORDER):
        ax = axes[i]
        col = f"{sig}_AUCell"

        sig_sum = summary[summary["signature"] == sig]
        sf_sum = sig_sum[sig_sum["tissue"] == TISSUE_NUM].iloc[0]
        pb_sum = sig_sum[sig_sum["tissue"] == TISSUE_DEN].iloc[0]

        n_nom = int(sf_sum["set_size_nominal"])
        n_eff = int(sf_sum["set_size_effective"])
        band = str(sf_sum["power_band"])
        meta = SIG_META[sig]

        sub = tregs[[col, TISSUE_KEY, DONOR_KEY]].dropna()
        sub_sf = sub[sub[TISSUE_KEY] == TISSUE_NUM][col]
        sub_pb = sub[sub[TISSUE_KEY] == TISSUE_DEN][col]

        # Draw boxplots / violins
        data_to_plot = [sub_pb, sub_sf]
        positions = [0, 1]
        bp = ax.boxplot(
            data_to_plot,
            positions=positions,
            widths=0.55,
            patch_artist=True,
            showmeans=True,
            meanline=True,
            showfliers=False,
            boxprops=dict(linewidth=1.2),
            medianprops=dict(color="black", linewidth=1.8),
            meanprops=dict(color="darkred", linewidth=1.5, linestyle="--"),
            whiskerprops=dict(linewidth=1.2),
            capprops=dict(linewidth=1.2),
        )

        bp["boxes"][0].set_facecolor(color_pb)
        bp["boxes"][0].set_alpha(0.7)
        bp["boxes"][1].set_facecolor(color_sf)
        bp["boxes"][1].set_alpha(0.7)

        ax.set_xticks(positions)
        ax.set_xticklabels([f"PB\n(n={len(sub_pb):,})", f"SF\n(n={len(sub_sf):,})"], fontsize=10)

        # Title & Subtitle for each panel
        title_wrapped = textwrap.fill(meta["title"], width=26)
        ax.set_title(
            f"{title_wrapped}\n"
            f"{n_nom} genes ({n_eff} in data)\n{band}",
            fontsize=11,
            fontweight="bold",
            pad=12,
        )

        # Annotations
        delta_median = sf_sum["median_auc"] - pb_sum["median_auc"]
        ax.text(
            0.5,
            0.96,
            f"Δ median: {delta_median:+.4f}\nSF med: {sf_sum['median_auc']:.4f}\nPB med: {pb_sum['median_auc']:.4f}",
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=8.5,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8, edgecolor="gray"),
        )

        # Plain statement for 7-gene gate on face
        if sig == "Interaction_up":
            ax.text(
                0.5,
                0.62,
                "Caveat: 7-gene gate\nis underpowered &\ndominated by detection\nnoise at cell level;\nnot interpretable.",
                transform=ax.transAxes,
                ha="center",
                va="top",
                fontsize=8.0,
                color="darkred",
                fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="#FFEEEE", alpha=0.9, edgecolor="red"),
            )

        ax.set_ylabel("AUCell Score", fontsize=10)
        ax.grid(axis="y", linewidth=0.5, alpha=0.4)
        ax.set_axisbelow(True)

        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)

    # Super titles and tier declaration
    fig.suptitle(
        "Where does each signature localise across the synovial fluid vs peripheral blood split in JIA Tregs?",
        fontsize=14,
        fontweight="bold",
        y=0.98,
    )

    fig.text(
        0.5,
        0.915,
        "HYPOTHESIS-GENERATING TIER · per-cell AUCell score distributions across tissue niches (7 JIA donors; 6 paired)",
        ha="center",
        va="bottom",
        fontsize=11,
        color=palette["diverging"]["down"],
        fontweight="bold",
    )

    fig.text(
        0.5,
        0.875,
        "Confirmatory claims are carried by donor-level pseudobulk differential expression. "
        "GSE160097 is FACS-sorted: coarse_label == 'Treg' is a sort gate (non-circular).",
        ha="center",
        va="bottom",
        fontsize=9.5,
        color=text_color,
        style="italic",
    )

    fig.text(
        0.5,
        0.02,
        "Dashed red line = mean; solid black line = median. Power bands: testable (≥15 genes), underpowered_reported (5–14 genes). "
        "Interaction_down arms (0 genes) are structurally absent at nominal zero.",
        ha="center",
        va="top",
        fontsize=9.0,
        color=text_color,
    )

    return fig


def main() -> None:
    set_paper_style(config=FIG_CFG)

    scores_path = PATHS.tables(STAGE) / "treg_per_cell_scores.csv"
    summary_path = PATHS.tables(STAGE) / "treg_localisation_summary.csv"

    if not scores_path.exists() or not summary_path.exists():
        raise FileNotFoundError(f"Compute artifacts missing; run {STAGE} compute script first.")

    per_cell = pd.read_csv(scores_path)
    summary = pd.read_csv(summary_path)

    fig = build_figure(per_cell, summary)

    # Save dual PDF + PNG via save_overview()
    save_overview(
        fig,
        STAGE,
        STEM,
        table=summary,
        finding=(
            "Across sorted JIA Tregs, per-cell AUCell scores for WT_heat_up, HALLMARK_HYPOXIA, and "
            "WT_heat_up_purged_hypoxia are consistently higher in synovial fluid than in peripheral blood; "
            "the cGAS-dependent Interaction_fdrOnly_up (18-gene gate) shows modest SF elevation, whereas the "
            "7-gene Interaction_up gate is underpowered and dominated by detection noise at per-cell resolution."
        ),
        script=SCRIPT,
        fn="build_figure",
        config_kv="gsea_fdr = 0.05, percell_score_ncores = 4",
        input="03_results/12_treg_localisation/tables/treg_per_cell_scores.csv",
        how_to_read=(
            "Box plots show per-cell AUCell score distributions for sorted CD4_Treg cells in peripheral blood (PB, blue) "
            "versus synovial fluid (SF, vermillion). Solid lines indicate medians; dashed red lines indicate means. "
            "Set sizes (nominal and in-dataset effective) and power-band classifications are declared for each panel. "
            "This panel is hypothesis-generating tier; primary statistical claims are carried by donor-level pseudobulk DE."
        ),
        config=FIG_CFG,
        wide=True,
        height=6.2,
    )

    print("[12_treg_localisation_viz] VIZ DONE.")


if __name__ == "__main__":
    main()
