#!/usr/bin/env python
"""
16_narrative_embedding_viz.py: VIZ ONLY (no statistics).
=============================================================================
Static, print-ready UMAP grids of the frozen 99,915-cell sorted JIA T-cell
annotation (GSE160097), drawn from the committed per-cell substrate
`03_results/interactive/16_narrative_embedding.parquet`. Three figures, one
sampled frame each, so every panel in a figure holds the identical cells at
identical coordinates and a point sits in the same place across the row:

  umap_full_reference: the reference layout, tissue of origin, the frozen sort
                        gate (Treg / Tcon / CD8), and donor.
  umap_full_arms:      the three mouse-derived, human-projected UP arms,
                        WT_heat_up (199), KO_heat_up (218), Interaction_up (7).
  umap_full_programs:  three curated, anchor-independent lenses,
                        HALLMARK_HYPOXIA, the 21 published IFN-independent
                        STING-activation genes, the 200-gene generic type-I
                        interferon axis.

TIER. An embedding places cells against their labels. Claims in this
compartment rest on donor-level pseudobulk differential expression within the
frozen cell states, so every figure carries that standing line under its panel
row and the fuller statement in its README caption.

COLUMN DISCIPLINE. Score panels colour on the `*_AUCell` columns only. The
substrate also carries stale `published_*` columns, one of which
(`published_WT_heat_up`) is a mean-centred scanpy `score_genes` module score on
a different scale, and those are never drawn.

Run in-container from the compartment root, AFTER 16_narrative_scoring.py:
    python 02_analysis/scripts/16_narrative_embedding_viz.py
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

STAGE = "16_narrative_scoring"
SCRIPT = "02_analysis/scripts/16_narrative_embedding_viz.py"
SUBSTRATE = "16_narrative_embedding.parquet"
SUBSTRATE_REL = f"03_results/interactive/{SUBSTRATE}"
SUMMARY_REL = "03_results/16_narrative_scoring/tables/narrative_score_summary.csv"

# --- drawing parameters (shared with 17_treg_reembedding_viz.py) -------------
SAMPLE_N = 60_000       # cells drawn per figure; sampled ONCE and reused by every panel
SAMPLE_SEED = 0         # fixed so the frame is identical across runs and across figures
PT = 1.6                # scatter point size at this density
CLIP = (2, 98)          # robust percentile clip for every continuous panel
CMAP = "viridis"        # ONE perceptually-uniform sequential scale for every AUCell panel
BOX_PAD = 0.03          # fraction of the square bounding box added as margin

# Panel geometry is placed by hand in figure coordinates rather than by tight_layout,
# so every panel box is EXACTLY square and the slack left over in each cell is a known
# amount of room for that panel's colourbar.
FIG_W, FIG_H = 15.5, 6.0   # one canvas for every three-panel row in this bundle
ROW_BOTTOM, ROW_TOP = 0.165, 0.855   # the band the square panels occupy
ROW_LEFT, ROW_RIGHT = 0.030, 0.995
CB_W_IN, CB_PAD_IN = 0.14, 0.11      # colourbar bar width and its gap from the panel
FOOTER_Y, FOOTER_STEP = 0.035, 0.040  # standing lines under the row

# Categorical palettes, read from the config Okabe-Ito block by NAME (no raw hex here).
_OI = (FIG_CFG.get("colors", {}) or {}).get("okabe_ito", {}) or {}
TISSUE_COL = {"synovial_fluid": _OI["vermillion"], "peripheral_blood": _OI["blue"]}
STATE_COL = {"Treg": _OI["bluish_green"], "Tcon": _OI["orange"], "CD8": _OI["reddish_purple"]}
# Seven donors, seven Okabe-Ito hues. Orange is left out rather than yellow, because
# orange and vermillion are the one confusable pair in this palette at scatter-point size.
DONOR_HUES = ["vermillion", "sky_blue", "bluish_green", "yellow",
              "blue", "reddish_purple", "black"]

TISSUE_LABEL = {"synovial_fluid": "synovial fluid", "peripheral_blood": "paired blood"}

TIER_LINE = ("Annotation tier. Claims in this compartment rest on donor-level pseudobulk "
             "differential expression within these frozen cell states.")

# Panel titles carry the set identifier and its size and nothing else. Every one of
# these identifiers already names the set by how it was derived, and the derivation is
# spelled out in the README caption rather than on the panel.
ARM_PANELS = [
    ("WT_heat_up_AUCell", "WT_heat_up (199 genes)"),
    ("KO_heat_up_AUCell", "KO_heat_up (218 genes)"),
    ("Interaction_up_AUCell", "Interaction_up (7 genes)"),
]
PROGRAM_PANELS = [
    ("HALLMARK_HYPOXIA_AUCell", "HALLMARK_HYPOXIA (200 genes)"),
    ("sting_specific_published_AUCell", "sting_specific_published (21 genes)"),
    ("ifn_generic_axis_AUCell", "ifn_generic_axis (200 genes)"),
]
# WT_heat_up and KO_heat_up share 182 of their genes and are meant to be read against
# each other, so they go on ONE colour scale. Interaction_up spans a different range by
# an order of magnitude and keeps its own bar.
ARM_SHARED_SCALE = ["WT_heat_up_AUCell", "KO_heat_up_AUCell"]
ARM_SETS = ["WT_heat_up", "KO_heat_up", "Interaction_up"]
PROGRAM_SETS = ["HALLMARK_HYPOXIA", "sting_specific_published", "ifn_generic_axis"]


# =============================================================================
# Shared panel primitives.
# Same behaviour as the embedding-atlas primitives already used in this
# compartment (02_analysis/scripts/07_embedding_viz.py): a 2nd-to-98th
# percentile clip, an argsort so high-scoring cells draw last, and a rasterized
# data layer. Extended here with a shared square bounding box and a colourbar
# slot on every panel so the panels of a row stay the same size.
# =============================================================================
def _fs(key: str) -> float:
    """One font size, read from the config `figures:` block (never a literal)."""
    return float((FIG_CFG.get("figures", {}) or {})[key])


def sample_frame(df: pd.DataFrame, n: int | None, seed: int = SAMPLE_SEED) -> pd.DataFrame:
    """Draw the ONE frame a figure's panels all share, then shuffle its row order.

    Shuffling matters for the categorical panels: every point is drawn in a single
    scatter call in shuffled order, so no group is systematically painted over another.
    """
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


def shared_clip(d: pd.DataFrame, cols: list) -> "tuple[float, float]":
    """One 2nd-to-98th percentile clip over the pooled values of several columns.

    Used where two sets are meant to be read against each other: putting them on one
    scale is what makes the comparison possible.
    """
    pooled = np.concatenate([d[c].to_numpy(dtype=float) for c in cols])
    lo, hi = np.nanpercentile(pooled, list(CLIP))
    return float(lo), float(hi)


def scatter_cont(ax, d: pd.DataFrame, col: str, title: str, xlim, ylim, vlim=None):
    """Continuous AUCell panel: robust clip, high values drawn last, colourbar reads AUCell.

    `vlim` overrides the per-panel clip with a scale shared across panels.
    """
    _frame_panel(ax, xlim, ylim, title)
    v = d[col].to_numpy(dtype=float)
    lo, hi = vlim if vlim is not None else np.nanpercentile(v, list(CLIP))
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


def _row_figure(n_panels: int = 3):
    """One row of square panels on the shared canvas, with bands reserved top and bottom."""
    fig = plt.figure(figsize=(FIG_W, FIG_H))
    return fig, row_axes(fig, n_panels)


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
# Source tables. Aggregates only: a per-cell frame is not a readable source CSV.
# =============================================================================
def reference_table(full: pd.DataFrame, drawn: pd.DataFrame) -> pd.DataFrame:
    """Per (cell state x tissue x donor) cell counts, full object and as drawn."""
    keys = ["coarse_label", "tissue", "donor"]
    tbl = full.groupby(keys, observed=True).size().rename("n_cells").reset_index()
    got = drawn.groupby(keys, observed=True).size().rename("n_cells_drawn").reset_index()
    tbl = tbl.merge(got, on=keys, how="left")
    tbl["n_cells_drawn"] = tbl["n_cells_drawn"].fillna(0).astype(int)
    tbl["frac_of_total"] = tbl["n_cells"] / len(full)
    return tbl.sort_values(keys).reset_index(drop=True)


def score_table(summary: pd.DataFrame, sets: list, labels: list | None = None) -> pd.DataFrame:
    """The committed per (set x cell state x tissue) summary, restricted to the panels drawn."""
    tbl = summary[summary["set_name"].isin(sets)].copy()
    if labels is not None:
        tbl = tbl[tbl["coarse_label"].isin(labels)]
    tbl["set_name"] = pd.Categorical(tbl["set_name"], categories=sets, ordered=True)
    return tbl.sort_values(["set_name", "coarse_label", "tissue"]).reset_index(drop=True)


def _mean_of(summary: pd.DataFrame, set_name: str, label: str, tissue: str) -> float:
    hit = summary[(summary["set_name"] == set_name) & (summary["coarse_label"] == label)
                  & (summary["tissue"] == tissue)]
    if hit.empty:
        raise ValueError(f"no summary row for {set_name} / {label} / {tissue}")
    return float(hit["mean"].iloc[0])


def _shift(summary: pd.DataFrame, set_name: str, label: str) -> str:
    pb = _mean_of(summary, set_name, label, "peripheral_blood")
    sf = _mean_of(summary, set_name, label, "synovial_fluid")
    return f"{label} {pb:.4f} to {sf:.4f}"


# =============================================================================
# Figures
# =============================================================================
def figure_reference(full: pd.DataFrame, drawn: pd.DataFrame, xlim, ylim) -> plt.Figure:
    """The reference layout: tissue, frozen sort gate, donor, over one shared frame."""
    donors = sorted(drawn["donor"].unique())
    donor_col = {d: _OI[h] for d, h in zip(donors, DONOR_HUES)}
    donor_label = {d: d.replace("_", " ") for d in donors}

    fig, axes = _row_figure()
    scatter_cat(axes[0], drawn, "tissue", TISSUE_COL, "Tissue of origin", xlim, ylim,
                order=["synovial_fluid", "peripheral_blood"], labels=TISSUE_LABEL, ncol=2)
    scatter_cat(axes[1], drawn, "coarse_label", STATE_COL, "Frozen cell state",
                xlim, ylim, order=["Treg", "Tcon", "CD8"], ncol=3)
    scatter_cat(axes[2], drawn, "donor", donor_col, "Donor", xlim, ylim,
                order=donors, labels=donor_label, ncol=4)
    _dress(fig,
           "Sorted JIA T cells: tissue, cell state and donor",
           f"UMAP of the frozen annotation. {len(drawn):,} of {len(full):,} cells drawn, "
           "one sampled frame shared by all three panels.",
           [TIER_LINE])
    return fig


def figure_arms(drawn: pd.DataFrame, xlim, ylim, shared) -> plt.Figure:
    """The three mouse-derived human-projected up arms on the shared frame.

    WT_heat_up and KO_heat_up are drawn on ONE colour scale so the two can be compared
    directly. Interaction_up spans a range an order of magnitude wider and keeps its own.
    """
    fig, axes = _row_figure()
    for ax, (col, title) in zip(axes, ARM_PANELS):
        scatter_cont(ax, drawn, col, title, xlim, ylim,
                     vlim=shared if col in ARM_SHARED_SCALE else None)
    _dress(fig,
           "Mouse 39 °C-derived up arms on the JIA T-cell map",
           "Colour is per-cell AUCell. WT_heat_up and KO_heat_up share one scale "
           f"({shared[0]:.3f} to {shared[1]:.3f}). Interaction_up is scaled on its own.",
           [TIER_LINE])
    return fig


def figure_programs(drawn: pd.DataFrame, xlim, ylim) -> plt.Figure:
    """Three curated, anchor-independent lenses on the same frame as the arms."""
    fig, axes = _row_figure()
    for ax, (col, title) in zip(axes, PROGRAM_PANELS):
        scatter_cont(ax, drawn, col, title, xlim, ylim)
    _dress(fig,
           "Curated program lenses on the JIA T-cell map",
           "Colour is per-cell AUCell, clipped to the 2nd and 98th percentile within a panel.",
           [TIER_LINE])
    return fig


# =============================================================================
def main() -> None:
    set_paper_style(config=FIG_CFG)

    full = pd.read_parquet(PATHS.interactive / SUBSTRATE)
    summary = pd.read_csv(PATHS.tables(STAGE) / "narrative_score_summary.csv")
    print(f"[16_narrative_embedding_viz] substrate {full.shape[0]:,} cells x {full.shape[1]} cols, "
          f"{full['donor'].nunique()} donors, {full['coarse_label'].nunique()} frozen cell states")

    drawn = sample_frame(full, SAMPLE_N)
    xlim, ylim = square_box(drawn)
    print(f"[16_narrative_embedding_viz] one shared frame of {len(drawn):,} cells, "
          f"square box x {xlim[0]:.2f}..{xlim[1]:.2f}")

    arm_shared = shared_clip(drawn, ARM_SHARED_SCALE)
    print(f"[16_narrative_embedding_viz] WT_heat_up + KO_heat_up shared colour scale "
          f"{arm_shared[0]:.4f} to {arm_shared[1]:.4f}")

    cfg_common = (f"figures.dpi = 300, figures.rasterized_dpi = 600, sample_n = {SAMPLE_N}, "
                  f"sample_seed = {SAMPLE_SEED}, point_size = {PT}")

    # --- 1. reference layout -------------------------------------------------
    fig = figure_reference(full, drawn, xlim, ylim)
    save_overview(
        fig, STAGE, "umap_full_reference",
        table=reference_table(full, drawn),
        finding=("One sampled frame of the frozen 99,915-cell sorted JIA T-cell map coloured "
                 "three ways, by tissue of origin, by the frozen FACS sort gate and by donor, "
                 "establishing the layout that every score colouring of this substrate is read "
                 "against: all seven JIA donors contribute cells to both the synovial-fluid and "
                 "the paired-blood side of the map, and the sort gates occupy largely distinct "
                 "territory within each tissue."),
        script=SCRIPT, fn="figure_reference",
        config_kv=(f"{cfg_common}, colours = colors.okabe_ito (tissue vermillion/blue, "
                   "state green/orange/pink, donor 7 hues)"),
        input=SUBSTRATE_REL,
        how_to_read=(
            "Three panels over ONE sampled frame of the same cells at the same coordinates, so a "
            "cell sits in the same place in all three. Left is tissue of origin, synovial fluid in "
            "vermillion and paired blood in blue. Middle is the frozen FACS sort gate the "
            "compartment is built on, Treg, Tcon and CD8. Right is donor, one hue per JIA "
            "patient. Points are drawn in "
            "shuffled order so overlapping groups paint evenly, the axes are UMAP coordinates "
            "without units, and all three panels share one square bounding box, so the row is "
            "comparable panel to panel and comparable to the score figures beside it. "
            f"{SAMPLE_N:,} of the 99,915 cells are drawn with a fixed seed, and the source table "
            "gives the full per cell state, tissue and donor counts next to the counts drawn. This "
            "is annotation. Claims in this compartment rest on donor-level pseudobulk differential "
            "expression within these frozen cell states, meta-analysed as effect sizes with "
            "confidence intervals."),
        width=FIG_W, height=FIG_H, config=FIG_CFG)

    # --- 2. mouse-derived up arms -------------------------------------------
    fig = figure_arms(drawn, xlim, ylim, arm_shared)
    save_overview(
        fig, STAGE, "umap_full_arms",
        table=score_table(summary, ARM_SETS),
        finding=("All three mouse 39 °C-derived up arms colour the synovial-fluid side of the map "
                 "brighter than the paired-blood side in every frozen sort gate, and the per-cell "
                 "AUCell means behind that colouring put the WT_heat_up shift at "
                 f"{_shift(summary, 'WT_heat_up', 'Treg')} in Treg, "
                 f"{_shift(summary, 'WT_heat_up', 'Tcon')} in Tcon and "
                 f"{_shift(summary, 'WT_heat_up', 'CD8')} in CD8, so the colouring tracks tissue "
                 "across the whole sorted compartment and reads as gate-shared."),
        script=SCRIPT, fn="figure_arms",
        config_kv=(f"{cfg_common}, cmap = {CMAP}, clip_percentiles = {list(CLIP)}, "
                   "columns = WT_heat_up_AUCell, KO_heat_up_AUCell, Interaction_up_AUCell, "
                   f"shared_scale = WT_heat_up + KO_heat_up at {arm_shared[0]:.4f}-{arm_shared[1]:.4f}"),
        input=f"{SUBSTRATE_REL}, {SUMMARY_REL}",
        how_to_read=(
            "Three panels over the SAME sampled frame and square bounding box as the reference "
            "figure. Panel titles carry the set identifier and its size. WT_heat_up is the up arm "
            "of the mouse WT iTreg 39 versus 37 °C contrast, 199 human symbols. KO_heat_up is the "
            "same contrast in cGAS-knockout iTregs, 218 symbols. Interaction_up is the mouse "
            "genotype by temperature up arm, 7 symbols, small enough that one gene moves the "
            f"score, so read it for location and treat its spread as noise. WT_heat_up and "
            f"KO_heat_up share 182 genes and are drawn on ONE colour scale, {arm_shared[0]:.4f} to "
            f"{arm_shared[1]:.4f}, so the two panels can be compared pixel for pixel. "
            "Interaction_up spans a range an order of magnitude wider, so a common scale would "
            "flatten both, and it carries its own bar. Colour is per-cell AUCell with the "
            "highest-scoring cells drawn last, and the limits are the 2nd and 98th percentile of "
            "the values on that scale. AUCell is bounded in [0, 1] and scales with set size, so "
            "the source table carries mean, median, standard deviation and cell and donor counts "
            "for any comparison the colour cannot make. Cells are pooled across donors, making a "
            "tissue difference here pseudoreplicated and descriptive. Claims in this compartment "
            "rest on donor-level pseudobulk differential expression within the frozen cell "
            "states."),
        width=FIG_W, height=FIG_H, config=FIG_CFG)

    # --- 3. curated program lenses ------------------------------------------
    fig = figure_programs(drawn, xlim, ylim)
    save_overview(
        fig, STAGE, "umap_full_programs",
        table=score_table(summary, PROGRAM_SETS),
        finding=("The curated hypoxia lens and the generic type-I interferon axis colour the "
                 "synovial-fluid side brighter in all three sort gates, while the 21 published "
                 "IFN-independent STING-activation genes sit far lower in the Treg gate than in "
                 "Tcon or CD8 in both tissues (per-cell AUCell mean "
                 f"{_mean_of(summary, 'sting_specific_published', 'Treg', 'peripheral_blood'):.4f} "
                 "in Treg blood against "
                 f"{_mean_of(summary, 'sting_specific_published', 'Tcon', 'peripheral_blood'):.4f} "
                 "in Tcon and "
                 f"{_mean_of(summary, 'sting_specific_published', 'CD8', 'peripheral_blood'):.4f} "
                 "in CD8), so that panel reports a sort-gate difference alongside a tissue one."),
        script=SCRIPT, fn="figure_programs",
        config_kv=(f"{cfg_common}, cmap = {CMAP}, clip_percentiles = {list(CLIP)}, "
                   "columns = HALLMARK_HYPOXIA_AUCell, sting_specific_published_AUCell, "
                   "ifn_generic_axis_AUCell"),
        input=f"{SUBSTRATE_REL}, {SUMMARY_REL}",
        how_to_read=(
            "Three panels over the SAME sampled frame and bounding box as the mouse-arm figure, "
            "on the same sequential colormap. Panel titles carry the set identifier and its size. "
            "Each set here is curated, versioned and derived without reference to the mouse "
            "anchor: "
            "HALLMARK_HYPOXIA from MSigDB Hallmark, sting_specific_published the 21 published "
            "IFN-independent STING-activation genes, and ifn_generic_axis a 200-gene generic "
            "type-I interferon axis of which 116 symbols match this object, the thinnest "
            "intersection in the panel. Colour is per-cell AUCell, clipped to the 2nd and 98th "
            "percentile within each panel, highest-scoring cells drawn last. These three sets "
            "are unrelated to each other and their ranges differ, so each keeps its own scale and "
            "brightness compares tissues within a panel while the source table carries the "
            "cross-panel numbers. Hypoxia and temperature are both imposed by the "
            "inflamed joint and stay entangled in cross-sectional human data, so this hypoxia "
            "panel is one lens on that niche. Claims in this compartment rest on donor-level "
            "pseudobulk differential expression within the frozen cell states."),
        width=FIG_W, height=FIG_H, config=FIG_CFG)

    # --- 4. captions for the three same-stem source tables -------------------
    write_caption(
        STAGE, "tables/_overview/umap_full_reference.csv",
        finding=("The 39 populated (cell state x tissue x donor) strata behind the reference map: "
                 "all seven donors appear in both tissues, and the three strata that are absent "
                 "are single sort gates in one donor arm rather than a missing donor."),
        script=SCRIPT, fn="reference_table",
        config_kv=f"sample_n = {SAMPLE_N}, sample_seed = {SAMPLE_SEED}",
        input=SUBSTRATE_REL,
        how_to_read=("One row per (`coarse_label` x `tissue` x `donor`) stratum present in the "
                     "frozen annotation. `n_cells` is the full-object count and `n_cells_drawn` "
                     "the count in the sampled frame the figure draws, so the two say how faithful "
                     "the drawn frame is to the object. `frac_of_total` is `n_cells` over 99,915. "
                     "Absent combinations carry no row. Annotation tier, no test and no effect "
                     "size."),
        config=FIG_CFG)
    for stem, sets in (("umap_full_arms", ARM_SETS), ("umap_full_programs", PROGRAM_SETS)):
        write_caption(
            STAGE, f"tables/_overview/{stem}.csv",
            finding=(f"Per-cell AUCell summaries of the {len(sets)} sets drawn in "
                     f"`figures/_overview/{stem}.png` ({', '.join(sets)}), one row per cell state "
                     "and tissue, so the colouring can be read as numbers."),
            script=SCRIPT, fn="score_table",
            config_kv="rows = 3 sets x 3 frozen cell states x 2 tissues, metric = AUCell",
            input=SUMMARY_REL,
            how_to_read=("A restriction of the stage summary table to the sets this figure draws. "
                         "One row per (`set_name` x `coarse_label` x `tissue`) with the mean, "
                         "median and standard deviation of the per-cell AUCell score and the cell "
                         "and donor counts behind it. AUCell is bounded in [0, 1] and its scale "
                         "depends on set size, so values compare across tissue within a "
                         "`set_name` and the source of a cross-set comparison is the gene lists, "
                         "not these means. Cells are pooled across donors, so the unit of "
                         "replication is the cell and every tissue difference here is "
                         "pseudoreplicated. `evidence_tier` reads `secondary_percell` "
                         "throughout."),
            config=FIG_CFG)

    print("[16_narrative_embedding_viz] wrote 3 overviews "
          "(umap_full_reference, umap_full_arms, umap_full_programs) + 3 table captions")


if __name__ == "__main__":
    main()
