#!/usr/bin/env python
"""
16_narrative_scoring_viz.py — VIZ ONLY (no computing).
=============================================================================
Draws the provenance seam check for the narrative per-cell substrate: for each
score that already existed in the published per-cell readout, the newly computed
AUCell value against the published one, cell by cell.

This is a provenance check, not a result. It carries no biology and makes no
claim about JIA. Its only job is to let a reader see, before trusting any
colouring built on `interactive/16_narrative_embedding.parquet`, whether that
substrate sits on the same scale as the published readout. So the panel is
deliberately plain: grey points, one identity line, the correlations on the face.

Two panels land on the identity line at r = 1.000000 — the AUCell scorer is
bit-for-bit reproducible. The third does not, and that is the point of drawing it:
`published_WT_heat_up` is a stale mean-centred scanpy `score_genes` module score,
not AUCell, so the two axes are different metrics and no identity line applies.
The same mouse arm reproduces at r = 1.000000 against stage 05's canonical AUCell
column, stated in the figure-level note beneath the panels (it sat inside the panel
until it was found to cover the point cloud and the y tick labels it explains). The
mouse up arm must therefore be coloured with `WT_heat_up_AUCell`, never with
`published_WT_heat_up`.

Outputs:
  03_results/16_narrative_scoring/figures/_overview/narrative_score_seam_check.pdf
  03_results/16_narrative_scoring/figures/_overview/narrative_score_seam_check.png
  03_results/16_narrative_scoring/tables/_overview/narrative_score_seam_check.csv
  03_results/16_narrative_scoring/README.md                (caption, via save_overview)

Run in-container from the compartment root, AFTER 16_narrative_scoring.py:
    python 02_analysis/scripts/16_narrative_scoring_viz.py
"""
from __future__ import annotations

import os
import sys
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

COMPARTMENT_ROOT = Path(__file__).resolve().parents[2]

sys.path.insert(0, str(COMPARTMENT_ROOT))
sys.path.insert(0, str(COMPARTMENT_ROOT / "02_analysis"))
os.chdir(COMPARTMENT_ROOT)

from config import PATHS  # noqa: E402
from helpers.figure_style import FIG_CFG, save_overview, set_paper_style  # noqa: E402

STAGE = "16_narrative_scoring"
SCRIPT = "02_analysis/scripts/16_narrative_scoring_viz.py"
STEM = "narrative_score_seam_check"
SUBSTRATE_PARQUET = "16_narrative_embedding.parquet"

# The three sets the brief asks to see, in panel order. Each names the new AUCell column
# and the published column it is checked against.
PANELS = [
    ("HALLMARK_HYPOXIA", "HALLMARK_HYPOXIA_AUCell", "published_HALLMARK_HYPOXIA_AUCell"),
    ("HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE_AUCell",
     "published_HALLMARK_UNFOLDED_PROTEIN_RESPONSE_AUCell"),
    ("WT_heat_up", "WT_heat_up_AUCell", "published_WT_heat_up"),
]


def _seam_row(seam: pd.DataFrame, new_col: str, ref_col: str) -> pd.Series:
    hit = seam[(seam["new_column"] == new_col) & (seam["reference_column"] == ref_col)]
    if hit.empty:
        raise ValueError(f"seam-check table has no row for {new_col} vs {ref_col}")
    return hit.iloc[0]


def build_figure(df: pd.DataFrame, seam: pd.DataFrame) -> plt.Figure:
    """One scatter per checked set: newly computed AUCell (y) against the published column
    (x), every cell drawn. Same-metric panels get an identity line; the cross-metric panel
    does not, because an identity line between two different metrics would assert a
    comparability that is not there."""
    palette = FIG_CFG["colors"]
    ink = palette["okabe_ito"]["black"]
    ok_hue = palette["okabe_ito"]["bluish_green"]
    warn_hue = palette["okabe_ito"]["vermillion"]
    ident_hue = palette["diverging"]["down"]

    # Type sizes come from the project figure contract, never from a local guess — those are
    # FLOORS, so overlaps are resolved with canvas and spacing, never by shrinking text.
    fz = FIG_CFG["figures"]
    SZ_TITLE = fz["title_size"] + 1          # 17 — figure question
    SZ_PANEL = fz["base_size"]               # 14 — panel titles
    SZ_AXIS = fz["axis_title_size"]          # 13 — axis titles
    SZ_ANNOT = fz["legend_text_size"]        # 11 — in-panel statistics + legend
    SZ_CAPTION = fz["caption_size"] + 1      # 10 — figure-level notes

    # Canvas matches the wide preset save_figure() will enforce, so the layout drawn here
    # is the layout that ships (a mismatch silently re-squeezes every panel).
    fig, axes = plt.subplots(1, len(PANELS), figsize=(13, 8.2))
    fig.subplots_adjust(wspace=0.38, top=0.66, bottom=0.28, left=0.085, right=0.985)

    cross_note = ""
    for ax, (set_name, new_col, ref_col) in zip(np.atleast_1d(axes), PANELS):
        row = _seam_row(seam, new_col, ref_col)
        same_metric = str(row["comparison_kind"]) == "same_metric"
        hue = ok_hue if same_metric else warn_hue

        x = df[ref_col].to_numpy(dtype=float)
        y = df[new_col].to_numpy(dtype=float)
        ok = np.isfinite(x) & np.isfinite(y)

        # Rasterize the point cloud: ~100k markers must not become 100k vector paths.
        ax.scatter(x[ok], y[ok], s=1.6, c="#4D4D4D", alpha=0.10, linewidths=0,
                   rasterized=True)

        if same_metric:
            lo = float(min(x[ok].min(), y[ok].min()))
            hi = float(max(x[ok].max(), y[ok].max()))
            ax.plot([lo, hi], [lo, hi], color=ident_hue, linewidth=1.4, linestyle="--",
                    zorder=3, label="identity (y = x)")
            ax.legend(loc="lower right", frameon=True, fontsize=SZ_ANNOT)

        # Set name on line 1, the comparability verdict on line 2. Line 2 stays to two words:
        # WHICH metrics are being compared is already stated on both axis titles, and spelling
        # the pair out here wrapped the title onto an orphan third line.
        ax.set_title(
            textwrap.fill(set_name.replace("_", " "), width=28) + "\n"
            + ("same metric" if same_metric else "DIFFERENT METRICS"),
            fontsize=SZ_PANEL, fontweight="bold", pad=10,
            color=ink if same_metric else warn_hue)
        # Short, semantically complete axis titles. The exact column names live in the
        # same-stem seam-check table, and the panel title already names the set, so spelling
        # them out on both axes only makes neighbouring panels overprint each other.
        ax.set_xlabel(f"published score  ({row['reference_metric']})\n"
                      f"{Path(str(row['reference_source'])).name}", fontsize=SZ_AXIS)
        ax.set_ylabel("newly computed score  (AUCell)\nthis stage", fontsize=SZ_AXIS)

        verdict = ("reproduces" if bool(row["passes_floor"])
                   else "FAILS the r >= %.2f floor" % float(row["r_floor"]))
        ax.text(0.03, 0.97,
                f"Pearson r = {float(row['pearson_r']):.6f}\n"
                f"Spearman r = {float(row['spearman_r']):.6f}\n"
                f"n = {int(row['n_shared_cells']):,} cells\n{verdict}",
                transform=ax.transAxes, ha="left", va="top", fontsize=SZ_ANNOT,
                fontweight="bold", color=hue,
                bbox=dict(boxstyle="round,pad=0.35", facecolor="white", alpha=0.88,
                          edgecolor=hue))

        if not same_metric:
            # No extra in-panel flag: the vermillion title, the vermillion "FAILS the floor"
            # verdict and the vermillion footer note already carry it three times over, and a
            # fourth label only collides with the statistics box it is repeating.
            anchor = seam[(seam["new_column"] == new_col)
                          & (seam["comparison_kind"] == "same_metric")]
            anchor_txt = ""
            if not anchor.empty:
                a = anchor.iloc[0]
                anchor_txt = (f" Against the canonical AUCell column for the same arm "
                              f"({a['reference_source']}) it reproduces at Pearson "
                              f"r = {float(a['pearson_r']):.6f}.")
            cross_note = (
                f"Note on {set_name}: `{ref_col}` is a mean-centred scanpy score_genes module "
                "score, NOT AUCell — its values run negative, whereas AUCell is bounded in "
                f"[0, 1]. The r = {float(row['pearson_r']):.3f} is therefore a metric "
                "difference, not a scoring drift." + anchor_txt)

        ax.grid(linewidth=0.5, alpha=0.35)
        ax.set_axisbelow(True)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)

    fig.suptitle("Does the new per-cell substrate sit on the same scale as the published "
                 "per-cell readout?", fontsize=SZ_TITLE, fontweight="bold", y=0.985)
    fig.text(0.5, 0.905,
             "PROVENANCE CHECK · no biology, no claim — one point per cell "
             "(99,915 JIA SF/PB T cells), AUCell recomputed against the published columns",
             ha="center", va="bottom", fontsize=SZ_ANNOT, color=ident_hue, fontweight="bold")
    fig.text(0.5, 0.865,
             "Both genuine-AUCell comparisons land on the identity line at r = 1.000000. "
             "The third is a metric mismatch in the published readout, not a drift here.",
             ha="center", va="bottom", fontsize=SZ_CAPTION, color=ink, style="italic")
    # The cross-metric explanation lives HERE, at figure level, not inside panel 3 — in the
    # panel it covered the point cloud and the y tick labels it was explaining.
    if cross_note:
        fig.text(0.5, 0.135, textwrap.fill(cross_note, width=140),
                 ha="center", va="top", fontsize=SZ_CAPTION, color=warn_hue,
                 fontweight="bold")
    fig.text(0.5, 0.032,
             "Colour the mouse 39 °C-derived up arm with `WT_heat_up_AUCell`, never with "
             "`published_WT_heat_up`. The stale column is carried in the substrate only so "
             "this discrepancy stays visible.",
             ha="center", va="top", fontsize=SZ_CAPTION, color=ink)
    return fig


def main() -> None:
    set_paper_style(config=FIG_CFG)

    pq_path = PATHS.interactive_dir() / SUBSTRATE_PARQUET
    seam_path = PATHS.tables(STAGE) / "_overview" / f"{STEM}.csv"
    for p in (pq_path, seam_path):
        if not p.exists():
            raise FileNotFoundError(
                f"compute artifact missing at {p}; run 02_analysis/scripts/16_narrative_scoring.py first")

    df = pd.read_parquet(pq_path)
    seam = pd.read_csv(seam_path)
    # save_overview() rewrites this same table, so normalise the cell count on the way in.
    # It still lands as 99915.0 in the CSV (the figure-style table writer floats every numeric
    # column), but the written value no longer depends on what a previous render left on disk.
    seam["n_shared_cells"] = seam["n_shared_cells"].astype("int64")

    fig = build_figure(df, seam)

    save_overview(
        fig, STAGE, STEM,
        table=seam,
        finding=(
            "Re-deriving the two genuine AUCell columns of the published per-cell readout "
            "reproduces them exactly (Pearson and Spearman r = 1.000000 over all 99,915 "
            "cells), so this substrate sits on the published AUCell scale; the third "
            "comparison reaches only r = 0.755 because the published `WT_heat_up` column is "
            "a stale mean-centred scanpy score_genes module score rather than AUCell, and "
            "the same mouse up arm reproduces at r = 1.000000 against the canonical AUCell "
            "column for that arm."
        ),
        script=SCRIPT,
        fn="build_figure",
        config_kv="percell_score_ncores = 8",
        input="03_results/interactive/16_narrative_embedding.parquet",
        how_to_read=(
            "One grey point per cell: the published score on x, the score newly computed in "
            "this stage on y. Green statistics and a dashed blue identity line mark a "
            "same-metric comparison (AUCell vs AUCell), where landing on the line means the "
            "scorer reproduced the published column; vermillion statistics and a vermillion "
            "title mark a comparison between two DIFFERENT metrics, where no identity line "
            "applies and the scatter is expected to spread. Correlations are Pearson and "
            "Spearman over the shared barcodes, with the r >= 0.98 pass floor stated as a "
            "verdict in each box. This is a provenance panel and carries no biological claim; "
            "confirmatory claims are made by donor-level pseudobulk differential expression."
        ),
        config=FIG_CFG,
        wide=True,
        height=8.2,
    )

    print(f"[{STAGE}_viz] VIZ DONE.")


if __name__ == "__main__":
    main()
