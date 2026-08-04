#!/usr/bin/env python
"""
17_treg_reembedding_viz.py: VIZ ONLY (no statistics).
=============================================================================
Static, print-ready UMAP of the Treg-only re-embedding: the 27,175 cells of the
frozen `coarse_label == "Treg"` sort gate laid out on their own map, drawn from
the committed per-cell substrate `03_results/interactive/17_treg_reembedding.parquet`.

  umap_treg_reembedding: tissue of origin, the mouse WT 39 °C-derived up arm,
                          and the curated Hallmark hypoxia lens, over ONE frame
                          so every panel holds the identical cells at identical
                          coordinates.

WHICH COORDINATES. `x`/`y` in this substrate are the HARMONY-CORRECTED Treg-only
coordinates and `x_uncorrected`/`y_uncorrected` are the raw pair, which is the
opposite of the convention in `16_narrative_embedding.parquet`. The corrected
pair is what this figure draws, and the figure says so under its panel row.
Justification lives in the caption: at k = 30 the same-donor neighbour fraction
is 0.661 on the uncorrected map and 0.201 after Harmony over donor, against
0.146 expected from the donor proportions.

COLUMN DISCIPLINE. Score panels colour on the `*_AUCell` columns only. The
stale `published_*` columns carried alongside them are never drawn.

TIER. An embedding places cells against their labels, and Harmony reshapes the
space it corrects. Claims in this compartment rest on donor-level pseudobulk
differential expression within the frozen cell states.

Run in-container from the compartment root, AFTER 17_treg_reembedding.py:
    python 02_analysis/scripts/17_treg_reembedding_viz.py
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
from matplotlib.lines import Line2D  # noqa: E402

COMPARTMENT_ROOT = Path(__file__).resolve().parents[2]

sys.path.insert(0, str(COMPARTMENT_ROOT))
sys.path.insert(0, str(COMPARTMENT_ROOT / "02_analysis"))
os.chdir(COMPARTMENT_ROOT)

from config import PATHS  # noqa: E402
from helpers.figure_style import (  # noqa: E402
    FIG_CFG, save_overview, set_paper_style, write_caption,
)

STAGE = "17_treg_reembedding"
SOURCE_STAGE = "16_narrative_scoring"
SCRIPT = "02_analysis/scripts/17_treg_reembedding_viz.py"
SUBSTRATE = "17_treg_reembedding.parquet"
SUBSTRATE_REL = f"03_results/interactive/{SUBSTRATE}"
SUMMARY_REL = "03_results/16_narrative_scoring/tables/narrative_score_summary.csv"
MIXING_REL = "03_results/17_treg_reembedding/tables/treg_reembedding_mixing.csv"

# --- drawing parameters (same contract as 16_narrative_embedding_viz.py) -----
SAMPLE_N = None         # 27,175 Treg cells is a drawable density, so every cell is drawn
SAMPLE_SEED = 0
PT = 2.4                # scatter point size at this density
CLIP = (2, 98)
CMAP = "viridis"        # the SAME sequential scale every AUCell panel in this bundle uses
BOX_PAD = 0.03

# Panel geometry is placed by hand in figure coordinates rather than by tight_layout, so
# every panel box is EXACTLY square and the slack left in each cell is a known amount of
# room for that panel's colourbar. Same canvas as the full-object figures of
# 02_analysis/scripts/16_narrative_embedding_viz.py, so the bundle reads as one set.
FIG_W, FIG_H = 15.5, 6.0
ROW_BOTTOM, ROW_TOP = 0.185, 0.855   # the band the square panels occupy
ROW_LEFT, ROW_RIGHT = 0.030, 0.995
CB_W_IN, CB_PAD_IN = 0.14, 0.11      # colourbar bar width and its gap from the panel
FOOTER_Y, FOOTER_STEP = 0.035, 0.040  # standing lines under the row

_OI = (FIG_CFG.get("colors", {}) or {}).get("okabe_ito", {}) or {}
TISSUE_COL = {"synovial_fluid": _OI["vermillion"], "peripheral_blood": _OI["blue"]}
TISSUE_LABEL = {"synovial_fluid": "synovial fluid", "peripheral_blood": "paired blood"}

COORD_LINE = ("Coordinates: Treg-only UMAP with Harmony over donor applied to the PCA "
              "(columns x / y of the substrate).")
TIER_LINE = ("Annotation tier. Claims in this compartment rest on donor-level pseudobulk "
             "differential expression within the frozen cell states.")

# Panel titles carry the set identifier and its size and nothing else. Both identifiers
# already name the set by how it was derived, and the derivation is spelled out in the
# README caption rather than on the panel.
PANELS = [
    ("WT_heat_up_AUCell", "WT_heat_up (199 genes)"),
    ("HALLMARK_HYPOXIA_AUCell", "HALLMARK_HYPOXIA (200 genes)"),
]
PANEL_SETS = ["WT_heat_up", "HALLMARK_HYPOXIA"]


# =============================================================================
# Shared panel primitives, same behaviour as the embedding-atlas primitives in
# 02_analysis/scripts/07_embedding_viz.py (2nd-to-98th percentile clip, argsort
# so high-scoring cells draw last, rasterized data layer), extended with a
# shared square bounding box and a colourbar slot on every panel.
# =============================================================================
def _fs(key: str) -> float:
    """One font size, read from the config `figures:` block (never a literal)."""
    return float((FIG_CFG.get("figures", {}) or {})[key])


def sample_frame(df: pd.DataFrame, n: int | None, seed: int = SAMPLE_SEED) -> pd.DataFrame:
    """Draw the ONE frame every panel shares, then shuffle its row order."""
    d = df if (n is None or n >= len(df)) else df.sample(n=n, random_state=seed)
    return d.sample(frac=1.0, random_state=seed).reset_index(drop=True)


def square_box(d: pd.DataFrame, pad: float = BOX_PAD):
    """One square bounding box covering the frame, so panels are directly comparable."""
    x = d["x"].to_numpy(dtype=float)
    y = d["y"].to_numpy(dtype=float)
    cx, cy = (x.min() + x.max()) / 2.0, (y.min() + y.max()) / 2.0
    half = max(x.max() - x.min(), y.max() - y.min()) / 2.0 * (1.0 + pad)
    return (cx - half, cx + half), (cy - half, cy + half)


def row_axes(fig, n: int = 3):
    """Place `n` exactly-square panel boxes across the reserved band, left to right."""
    band = ROW_TOP - ROW_BOTTOM
    side_frac = band * FIG_H / FIG_W          # square: width in inches == height in inches
    cell = (ROW_RIGHT - ROW_LEFT) / n
    return [fig.add_axes([ROW_LEFT + i * cell, ROW_BOTTOM, side_frac, band])
            for i in range(n)]


def _frame_panel(ax, xlim, ylim, title: str):
    """Tick-free UMAP panel on the shared square bounding box."""
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(title, fontsize=_fs("base_size"), pad=10)


def _colorbar(fig, ax, sc):
    """A colourbar in the slack beside its panel, labelled by the unit it carries."""
    band = ROW_TOP - ROW_BOTTOM
    pos = ax.get_position()
    cax = fig.add_axes([pos.x1 + CB_PAD_IN / FIG_W, ROW_BOTTOM + 0.08 * band,
                        CB_W_IN / FIG_W, 0.84 * band])
    cb = fig.colorbar(sc, cax=cax)
    cb.set_label("AUCell", fontsize=_fs("axis_title_size"))
    cb.ax.tick_params(labelsize=_fs("axis_text_size"))
    cb.outline.set_visible(False)


def scatter_cont(ax, d: pd.DataFrame, col: str, title: str, xlim, ylim):
    """Continuous AUCell panel: robust clip, high values drawn last, colourbar reads AUCell."""
    _frame_panel(ax, xlim, ylim, title)
    v = d[col].to_numpy(dtype=float)
    lo, hi = np.nanpercentile(v, list(CLIP))
    order = np.argsort(v)
    sc = ax.scatter(d["x"].to_numpy()[order], d["y"].to_numpy()[order], c=v[order],
                    s=PT, cmap=CMAP, vmin=lo, vmax=hi, linewidths=0, rasterized=True)
    _colorbar(ax.figure, ax, sc)


def scatter_cat(ax, d: pd.DataFrame, col: str, palette: dict, title: str,
                xlim, ylim, order: list, labels: dict | None = None, ncol: int = 2):
    """Categorical panel: one shuffled scatter call, legend in the reserved band below."""
    _frame_panel(ax, xlim, ylim, title)
    values = d[col].astype(str).to_numpy()
    colours = np.array([palette[v] for v in values])
    ax.scatter(d["x"].to_numpy(), d["y"].to_numpy(), c=colours, s=PT,
               linewidths=0, rasterized=True)
    present = [c for c in order if c in set(values)]
    show = labels or {}
    handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor=palette[c],
                      markeredgecolor="none", markersize=9, label=show.get(c, c))
               for c in present]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.02),
              ncol=ncol, frameon=False, fontsize=_fs("legend_text_size"),
              handletextpad=0.25, columnspacing=1.2, borderpad=0.0)


def _dress(fig, suptitle: str, subtitle: str, footer_lines: list):
    """Title block above the row, standing lines below it."""
    fig.suptitle(suptitle, fontsize=_fs("title_size"), y=0.993)
    fig.text(0.5, 0.945, subtitle, ha="center", va="top", fontsize=_fs("subtitle_size"))
    y = FOOTER_Y + FOOTER_STEP * (len(footer_lines) - 1)
    for line in footer_lines:
        fig.text(0.5, y, line, ha="center", va="center",
                 fontsize=_fs("axis_title_size"), color=_OI["black"])
        y -= FOOTER_STEP


# =============================================================================
# Source table. The committed per (set x cell state x tissue) summary, restricted
# to the Treg gate and the two sets drawn.
# =============================================================================
def score_table(summary: pd.DataFrame, sets: list, label: str = "Treg") -> pd.DataFrame:
    tbl = summary[(summary["set_name"].isin(sets)) & (summary["coarse_label"] == label)].copy()
    tbl["set_name"] = pd.Categorical(tbl["set_name"], categories=sets, ordered=True)
    return tbl.sort_values(["set_name", "tissue"]).reset_index(drop=True)


def _mean_of(summary: pd.DataFrame, set_name: str, tissue: str, label: str = "Treg") -> float:
    hit = summary[(summary["set_name"] == set_name) & (summary["coarse_label"] == label)
                  & (summary["tissue"] == tissue)]
    if hit.empty:
        raise ValueError(f"no summary row for {set_name} / {label} / {tissue}")
    return float(hit["mean"].iloc[0])


def _mixing(mix: pd.DataFrame, embedding: str, key: str) -> float:
    hit = mix[(mix["embedding"] == embedding) & (mix["grouping_key"] == key)
              & (mix["group"] == "_all_")]
    if hit.empty:
        raise ValueError(f"no mixing row for {embedding} / {key} / _all_")
    return float(hit["observed_same_frac"].iloc[0])


def _expected(mix: pd.DataFrame, embedding: str, key: str) -> float:
    hit = mix[(mix["embedding"] == embedding) & (mix["grouping_key"] == key)
              & (mix["group"] == "_all_")]
    return float(hit["expected_same_frac"].iloc[0])


# =============================================================================
def figure_treg(drawn: pd.DataFrame, xlim, ylim) -> plt.Figure:
    """Tissue, the mouse WT 39 °C-derived up arm, and the Hallmark hypoxia lens."""
    fig = plt.figure(figsize=(FIG_W, FIG_H))
    axes = row_axes(fig, 3)
    scatter_cat(axes[0], drawn, "tissue", TISSUE_COL, "Tissue of origin", xlim, ylim,
                order=["synovial_fluid", "peripheral_blood"], labels=TISSUE_LABEL, ncol=2)
    for ax, (col, title) in zip(axes[1:], PANELS):
        scatter_cont(ax, drawn, col, title, xlim, ylim)
    _dress(fig,
           "Treg-only map, Harmony corrected over donor",
           f"{len(drawn):,} sorted Treg cells, one frame shared by all three panels. Score "
           "colour is per-cell AUCell, each set clipped to its own 2nd and 98th percentile.",
           [COORD_LINE, TIER_LINE])
    return fig


def main() -> None:
    set_paper_style(config=FIG_CFG)

    treg = pd.read_parquet(PATHS.interactive / SUBSTRATE)
    summary = pd.read_csv(PATHS.tables(SOURCE_STAGE) / "narrative_score_summary.csv")
    mix = pd.read_csv(PATHS.tables(STAGE) / "treg_reembedding_mixing.csv")
    print(f"[17_treg_reembedding_viz] substrate {treg.shape[0]:,} Treg cells x {treg.shape[1]} cols, "
          f"{treg['donor'].nunique()} donors, cell states {sorted(treg['coarse_label'].unique())}")

    donor_raw = _mixing(mix, "treg_only", "donor")
    donor_harmony = _mixing(mix, "treg_only_harmony", "donor")
    donor_chance = _expected(mix, "treg_only_harmony", "donor")
    tissue_harmony = _mixing(mix, "treg_only_harmony", "tissue")
    tissue_chance = _expected(mix, "treg_only_harmony", "tissue")
    print(f"[17_treg_reembedding_viz] same-donor neighbours at k=30: uncorrected {donor_raw:.3f}, "
          f"Harmony {donor_harmony:.3f}, chance {donor_chance:.3f}")

    drawn = sample_frame(treg, SAMPLE_N)
    xlim, ylim = square_box(drawn)

    fig = figure_treg(drawn, xlim, ylim)
    save_overview(
        fig, STAGE, "umap_treg_reembedding",
        table=score_table(summary, PANEL_SETS),
        finding=("On the Treg-only map, drawn on the Harmony-corrected coordinates, the synovial-"
                 f"fluid and paired-blood cells still occupy distinct territory ({tissue_harmony:.3f} "
                 f"same-tissue neighbours at k = 30 against {tissue_chance:.3f} expected) after the "
                 f"same-donor neighbour fraction has dropped from {donor_raw:.3f} to "
                 f"{donor_harmony:.3f} against {donor_chance:.3f} expected, and both the mouse WT "
                 "39 \u00b0C-derived up arm and the curated Hallmark hypoxia lens colour the synovial-"
                 "fluid territory brighter (per-cell AUCell mean "
                 f"{_mean_of(summary, 'WT_heat_up', 'peripheral_blood'):.4f} to "
                 f"{_mean_of(summary, 'WT_heat_up', 'synovial_fluid'):.4f} for WT_heat_up and "
                 f"{_mean_of(summary, 'HALLMARK_HYPOXIA', 'peripheral_blood'):.4f} to "
                 f"{_mean_of(summary, 'HALLMARK_HYPOXIA', 'synovial_fluid'):.4f} for the hypoxia "
                 "lens)."),
        script=SCRIPT, fn="figure_treg",
        config_kv=(f"coordinates = x / y (Harmony over donor), all {len(drawn):,} cells drawn, "
                   f"point_size = {PT}, cmap = {CMAP}, clip_percentiles = {list(CLIP)}, "
                   "figures.dpi = 300, figures.rasterized_dpi = 600, "
                   "columns = WT_heat_up_AUCell, HALLMARK_HYPOXIA_AUCell"),
        input=f"{SUBSTRATE_REL}, {SUMMARY_REL}, {MIXING_REL}",
        how_to_read=(
            "Three panels over ONE frame of the same 27,175 sorted Treg cells at the same "
            "coordinates, sharing one square bounding box. Left is tissue of origin, synovial "
            "fluid in vermillion and paired blood in blue, drawn in shuffled order. Middle and "
            "right colour every cell by per-cell AUCell of one gene set, on the scale the "
            "full-object figures use, clipped to the 2nd and 98th percentile with the "
            "highest-scoring cells drawn last. Panel titles carry the set identifier and its "
            "size: WT_heat_up is the up arm of the mouse WT iTreg 39 versus 37 °C contrast in "
            "human projection, HALLMARK_HYPOXIA the curated MSigDB Hallmark program. The two sets "
            "are unrelated and their ranges differ, so each keeps its own colour scale. "
            "The coordinates are the Harmony-corrected pair, because at k = 30 the same-donor "
            f"neighbour fraction is {donor_raw:.3f} on the uncorrected Treg-only map and "
            f"{donor_harmony:.3f} after Harmony over donor against {donor_chance:.3f} expected, "
            f"while same-tissue neighbours hold at {tissue_harmony:.3f} against "
            f"{tissue_chance:.3f} expected. Harmony reshapes the space it corrects, so this map "
            "is annotation, and cells are pooled across donors, making a tissue difference read "
            "off the colouring pseudoreplicated and descriptive. Claims in this compartment rest "
            "on donor-level pseudobulk differential expression within the frozen cell states."),
        width=FIG_W, height=FIG_H, config=FIG_CFG)

    write_caption(
        STAGE, "tables/_overview/umap_treg_reembedding.csv",
        finding=("Per-cell AUCell summaries of the two sets drawn on the Treg-only map, "
                 "WT_heat_up and HALLMARK_HYPOXIA, restricted to the Treg gate, one row per "
                 "tissue, so the colouring can be read as numbers."),
        script=SCRIPT, fn="score_table",
        config_kv="rows = 2 sets x Treg x 2 tissues, metric = AUCell",
        input=SUMMARY_REL,
        how_to_read=("A restriction of the narrative scoring summary table to the Treg gate and "
                     "the two sets this figure draws. One row per (`set_name` x `tissue`) with "
                     "the mean, median and standard deviation of the per-cell AUCell score and "
                     "the cell and donor counts behind it. AUCell is bounded in [0, 1] and its "
                     "scale depends on set size, so values compare across tissue within a "
                     "`set_name`. Cells are pooled across donors, so the unit of replication is "
                     "the cell and the tissue difference here is pseudoreplicated. "
                     "`evidence_tier` reads `secondary_percell` throughout."),
        config=FIG_CFG)

    print("[17_treg_reembedding_viz] wrote 1 overview (umap_treg_reembedding) + 1 table caption")


if __name__ == "__main__":
    main()
