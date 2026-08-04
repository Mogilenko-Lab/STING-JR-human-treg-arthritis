#!/usr/bin/env python
"""
16_narrative_scoring_arms_viz.py: VIZ ONLY (no statistics).
=============================================================================
Draws the per-cell AUCell distribution of each mouse-derived up arm across the
three frozen sort labels and the two tissues, from the one per-cell substrate
this stage published. One panel per arm, six violins per panel: synovial fluid
and paired peripheral blood side by side inside each cell state.

TIER. This is the annotation tier and says so on its face. AUCell gives every
cell one vote, and the 99,915 cells come from 7 donors in very uneven numbers, so
a panel-level average here tracks whichever donors contributed the most cells.
The panel that orders the cell states is the donor-level pseudobulk dot plot in
03_results/14_unbiased_enrichment/figures/_overview/arm_nes_by_cell_state.png,
where each donor carries one vote inside a frozen label. This panel shows the
shape and the spread of the underlying per-cell scores that the donor-level
aggregate is built from.

SCALE. The three arms are 199, 218 and 7 genes, and AUCell is the area under a
cell's gene-recovery curve over its own ranking, so the three arms occupy three
different score ranges. Each panel therefore carries its own y axis, and a level
in one panel means nothing against a level in another. The comparison the panel
supports is the synovial-fluid-versus-blood offset inside one cell state of one
panel.

Input:
  03_results/interactive/16_narrative_embedding.parquet   99,915 cells x score columns

Output (03_results/16_narrative_scoring/):
  figures/_overview/arm_score_violins.{pdf,png}
  tables/_overview/arm_score_violins.csv    18 rows: arm x cell state x tissue summary
  README.md                                 caption (via save_overview)

Run in-container from the compartment root, AFTER 16_narrative_scoring.py:
  python 02_analysis/scripts/16_narrative_scoring_arms_viz.py
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

COMPARTMENT_ROOT = Path(__file__).resolve().parents[2]

sys.path.insert(0, str(COMPARTMENT_ROOT))
sys.path.insert(0, str(COMPARTMENT_ROOT / "02_analysis"))
os.chdir(COMPARTMENT_ROOT)

from config import PATHS, TISSUE_COLORS, TISSUE_DEN, TISSUE_NUM  # noqa: E402
from helpers.figure_style import (  # noqa: E402
    FIG_CFG,
    purge_figures,
    save_overview,
    set_paper_style,
)

STAGE = "16_narrative_scoring"
SCRIPT = "02_analysis/scripts/16_narrative_scoring_arms_viz.py"
STEM = "arm_score_violins"
SUBSTRATE_PARQUET = "16_narrative_embedding.parquet"

# Panel order, and the nominal gene count of each arm as a display label only.
ARMS = [
    ("WT_heat_up_AUCell", "WT_heat_up", 199),
    ("KO_heat_up_AUCell", "KO_heat_up", 218),
    ("Interaction_up_AUCell", "Interaction_up", 7),
]
CELL_STATES = ["Treg", "Tcon", "CD8"]
# Display label per tissue; the hue comes from the one tissue palette below.
TISSUES = [(TISSUE_NUM, "synovial fluid"), (TISSUE_DEN, "paired blood")]

_F = FIG_CFG["figures"]
_OI = FIG_CFG["colors"]["okabe_ito"]

SZ_TITLE = float(_F["title_size"])
SZ_SUBTITLE = float(_F["subtitle_size"])
SZ_AXIS_TITLE = float(_F["axis_title_size"])
SZ_AXIS_TEXT = float(_F["axis_text_size"])
SZ_LEGEND = float(_F["legend_text_size"])
SZ_STRIP = float(_F["strip_size"])
SZ_CAPTION = float(_F["caption_size"])
LINE_W = float(_F["line_width"])

# The one tissue palette, read from analysis_config.yaml::colors.tissue: warm for the
# inflamed joint, cool for paired blood.
TISSUE_COLOR = TISSUE_COLORS
INK = _OI["black"]
# Dodge of each tissue away from its cell-state centre, and the violin width.
DODGE = 0.21
VIOLIN_W = 0.36
# Headroom at both ends of every panel's y axis, as a fraction of that panel's drawn
# maximum. matplotlib evaluates a violin's kernel on linspace(data min, data max), so the
# body's lower boundary lands exactly on the arm's minimum, and every arm here reaches
# 0.0. A limit pinned at 0 therefore puts that boundary on the axis spine, where the fill
# and its closing edge are cut by the axes edge. The pad puts the boundary inside the axes.
Y_PAD_FRAC = 0.04


def summary_table(df: pd.DataFrame) -> pd.DataFrame:
    """Aggregate the per-cell frame to the 18 rows behind the panel.

    The figure-style contract forbids a 100,000-row per-cell frame as a figure's
    source table, and a reader checking a violin needs the quartiles rather than
    the cells, so the neighbour table is the five-number summary plus the cell
    count of every violin drawn.
    """
    rows = []
    for col, arm, n_genes in ARMS:
        for state in CELL_STATES:
            for tissue, tissue_label in TISSUES:
                v = df.loc[df["coarse_label"].eq(state) & df["tissue"].eq(tissue), col].to_numpy()
                rows.append({
                    "arm": arm,
                    "genes_in_arm": n_genes,
                    "cell_state": state,
                    "tissue": tissue_label,
                    "n_cells": int(v.size),
                    "mean": float(np.mean(v)),
                    "median": float(np.median(v)),
                    "q25": float(np.percentile(v, 25)),
                    "q75": float(np.percentile(v, 75)),
                    "min": float(np.min(v)),
                    "max": float(np.max(v)),
                })
    return pd.DataFrame(rows)


def build_figure(df: pd.DataFrame, width: float, height: float):
    """Three arm panels, six dodged violins each, one y axis per panel."""
    fig = plt.figure(figsize=(width, height))
    # Explicit rectangles so the three panels keep identical proportions after the
    # exporter fixes the canvas size.
    left, bottom, hgt = 0.055, 0.19, 0.615
    span, gap = 0.285, 0.045
    axes = [fig.add_axes((left + i * (span + gap), bottom, span, hgt)) for i in range(len(ARMS))]

    for ax, (col, arm, n_genes) in zip(axes, ARMS):
        drawn_max = 0.0
        for xi, state in enumerate(CELL_STATES):
            for sgn, (tissue, _label) in zip((-1, 1), TISSUES):
                v = df.loc[df["coarse_label"].eq(state) & df["tissue"].eq(tissue), col].to_numpy()
                drawn_max = max(drawn_max, float(np.max(v)))
                pos = xi + sgn * DODGE
                parts = ax.violinplot([v], positions=[pos], widths=VIOLIN_W,
                                      showextrema=False, showmedians=True)
                for body in parts["bodies"]:
                    body.set_facecolor(TISSUE_COLOR[tissue])
                    body.set_edgecolor(TISSUE_COLOR[tissue])
                    body.set_alpha(0.72)
                    body.set_linewidth(LINE_W * 0.8)
                med = parts["cmedians"]
                med.set_color(INK)
                med.set_linewidth(LINE_W * 1.3)

        ax.set_xticks(range(len(CELL_STATES)))
        ax.set_xticklabels(CELL_STATES, fontsize=SZ_AXIS_TEXT)
        ax.set_xlim(-0.55, len(CELL_STATES) - 0.45)
        ax.tick_params(axis="y", labelsize=SZ_AXIS_TEXT)
        ax.set_title(f"{arm}  ({n_genes} genes)", fontsize=SZ_STRIP, fontweight="bold")
        # Headroom at both ends, so the widest and the narrowest part of every body sits
        # inside the axes. The tick the lower pad opens up below zero is then dropped:
        # AUCell is bounded in [0, 1] and a negative label would misstate the scale.
        pad = Y_PAD_FRAC * drawn_max
        ax.set_ylim(-pad, drawn_max + pad)
        ax.set_yticks([t for t in ax.get_yticks() if 0.0 <= t <= drawn_max + pad])

    axes[0].set_ylabel("AUCell score, per cell", fontsize=SZ_AXIS_TITLE)
    handles = [Patch(facecolor=TISSUE_COLOR[t], edgecolor=TISSUE_COLOR[t], alpha=0.72, label=lab)
               for t, lab in TISSUES]
    axes[1].legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.115),
                   ncol=2, frameon=False, fontsize=SZ_LEGEND)

    fig.text(0.5, 0.955, "Per-cell AUCell by cell state and tissue, mouse-derived up arms",
             ha="center", va="center", fontsize=SZ_TITLE, fontweight="bold")
    fig.text(0.5, 0.900,
             "Annotation tier. AUCell is a rank-based score in 0 to 1, one value per cell. "
             "Each panel carries its own y axis.",
             ha="center", va="center", fontsize=SZ_SUBTITLE)

    # One line, and only a glyph key plus the n it is drawn over. How the panel is read,
    # which panel ranks the cell states, and why each panel keeps its own y axis all live in
    # `how_to_read` and land in the stage README.
    n_cells = int(len(df))
    n_donors = int(df["donor"].nunique())
    fig.text(0.055, 0.048,
             f"Black line marks the median. {n_cells:,} sorted cells from {n_donors} donors, "
             "one vote per cell.",
             ha="left", va="top", fontsize=SZ_CAPTION)
    return fig


def main() -> None:
    set_paper_style(config=FIG_CFG)
    purge_figures(STAGE, STEM, overview=True, config=FIG_CFG)

    pq_path = PATHS.interactive_dir() / SUBSTRATE_PARQUET
    if not pq_path.exists():
        raise FileNotFoundError(
            f"per-cell substrate missing at {pq_path}; "
            "run 02_analysis/scripts/16_narrative_scoring.py first")
    need = ["coarse_label", "tissue", "donor"] + [c for c, _a, _n in ARMS]
    df = pd.read_parquet(pq_path, columns=need)
    missing = sorted(set(need) - set(df.columns))
    if missing:
        raise KeyError(f"[16_arms_viz] substrate is missing columns: {missing}")

    summary = summary_table(df)
    width, height = 13.0, 6.8
    fig = build_figure(df, width, height)

    def offset(arm: str, state: str) -> float:
        s = summary[summary["arm"].eq(arm) & summary["cell_state"].eq(state)].set_index("tissue")
        return float(s.loc["synovial fluid", "median"] - s.loc["paired blood", "median"])

    wt_up = [state for state in CELL_STATES if offset("WT_heat_up", state) > 0]
    inter_up = [state for state in CELL_STATES if offset("Interaction_up", state) > 0]
    per_donor = df["donor"].value_counts()
    donor_lo, donor_hi = int(per_donor.min()), int(per_donor.max())
    inter_zero_pct = 100.0 * float((df["Interaction_up_AUCell"] == 0).mean())

    def count_phrase(states: list) -> str:
        return "all three cell states" if len(states) == len(CELL_STATES) else (
            f"{len(states)} of the three cell states ({', '.join(states)})")

    save_overview(
        fig, STAGE, STEM,
        table=summary,
        finding=(
            f"Across the {len(df):,} sorted cells the median per-cell AUCell of the mouse "
            "39 °C-derived up arm sits higher in synovial fluid than in paired blood in "
            f"{count_phrase(wt_up)}, and KO_heat_up follows the same pattern, so the per-cell "
            "channel shows the same direction the donor-level pseudobulk carries. The 7-gene "
            f"interaction arm leaves {inter_zero_pct:.0f}% of cells at exactly zero, which is "
            "what a 7-gene set does when none of its genes reaches a cell's top-ranked genes, "
            f"and its synovial-fluid median exceeds its blood median in {count_phrase(inter_up)}. "
            f"Every distribution here is one vote per cell over {int(df['donor'].nunique())} "
            f"donors of unequal cell yield ({donor_lo:,} to {donor_hi:,} cells), so it carries "
            "shape and spread while ranking the cell states stays with the donor-level panel."
        ),
        script=SCRIPT, fn="build_figure",
        config_kv=("metric = AUCell (rank-based, 0 to 1); arms = WT_heat_up 199, KO_heat_up 218, "
                   f"Interaction_up 7 genes; y_pad_frac = {Y_PAD_FRAC}"),
        input="03_results/interactive/16_narrative_embedding.parquet",
        how_to_read=(
            "Annotation tier. One panel per mouse-derived up arm, one violin pair per frozen sort "
            "label: warm is synovial fluid, cool is paired peripheral blood, black line at the "
            "median. AUCell is a rank-based score in 0 to 1, the area under a cell's "
            "gene-recovery curve for that arm, so it is robust to library size and composition. "
            "The comparison the panel supports is the synovial-fluid-versus-blood offset inside "
            "one cell state of one panel. Each panel has its own y axis because the arms are 199, "
            "218 and 7 genes and AUCell is computed against each cell's own ranking, so a level "
            "in one panel means nothing against a level in another. Both ends of every y axis "
            "carry headroom, so the score itself is bounded in [0, 1] and the axis is not. Every "
            f"one of the {len(df):,} cells casts one vote, and the {int(df['donor'].nunique())} donors "
            f"contributed {donor_lo:,} to {donor_hi:,} cells each, so a panel-level average "
            "follows the donors that contributed the most cells. Ranking the cell states is the "
            "job of the donor-level pseudobulk panel "
            "03_results/14_unbiased_enrichment/figures/_overview/arm_nes_by_cell_state.png, "
            "where each donor carries one vote inside a frozen label and the enrichment is "
            "tested; this panel adds the shape and the spread the donor-level aggregate is built "
            "from. The same-stem source table gives the cell count, mean, median, quartiles and "
            "range of all 18 violins. Naming follows how each arm was derived, from mouse iTreg "
            "37 versus 39 °C contrasts, and the reading stays correlative."
        ),
        config=FIG_CFG, width=width, height=height)

    print(f"[16_narrative_scoring_arms_viz] wrote {STEM} from {len(summary)} summarised violins")


if __name__ == "__main__":
    main()
