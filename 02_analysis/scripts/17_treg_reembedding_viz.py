#!/usr/bin/env python
"""
17_treg_reembedding_viz.py: VIZ ONLY (no statistics).
=============================================================================
Six-panel UMAP strips of the Treg-only re-embedding: the 27,175 cells of the frozen
`coarse_label == "Treg"` sort gate on their own map, drawn from
`03_results/interactive/17_treg_reembedding.parquet`.

  umap_treg_reembedding  tissue, four Treg identity genes, donor
  umap_treg_arms         counterpart of `16_narrative_scoring/umap_full_arms`
  umap_treg_programs     counterpart of `16_narrative_scoring/umap_full_programs`
  umap_treg_signatures   counterpart of `07_embedding/umap_signatures_treg`
  umap_treg_patchwork    the reference strip over the arm strip, for layout

GEOMETRY comes from `02_analysis/helpers/umap_grid.py`, shared with
`16_narrative_embedding_viz.py`, so a strip here stacks against its full-object twin with
panels at identical size and a column reads top to bottom.

ONE SCALE ACROSS EACH PAIR. Every counterpart panel takes its full-object twin's colour
limits, derived from the full-object substrate on the frame stage 16 draws, so a washed-out
panel states something about the Treg gate's range. The arm pair's shared limits are asserted
against the value stage 16 records, holding the two in step.

ONE BAR PER SCORE STRIP, as in stage 16. Each panel is rescaled across the twin's clip onto
[0, 1], so the strip carries one gradient. The clip is the twin's either way, so the picture and
the pairing hold; LEVEL moves to the caption's config line. The signature strip is the exception,
being module scores against a stage-07 twin drawn in its own units.

THE MARKERS ALSO SHARE ONE BAR, but in REAL units: four genes, one quantity, comparable ranges,
pooled on the full-object clip. Expression and AUCell stay on separate keys.

WHICH COORDINATES. `x`/`y` are the HARMONY-CORRECTED Treg-only coordinates and
`x_uncorrected`/`y_uncorrected` the raw pair, inverting the convention in
`16_narrative_embedding.parquet`. These strips draw the corrected pair and say so under the
panel row: at k = 30 the same-donor neighbour fraction is 0.661 raw and 0.201 after Harmony
over donor, against 0.146 expected from the donor proportions.

TIER. An embedding places cells against their labels, and Harmony reshapes the space it
corrects. Claims rest on donor-level pseudobulk differential expression within the frozen
cell states.

COLUMNS. Score panels colour on `*_AUCell`, sized from the stage-16 manifest. Marker genes
and any score column absent from this substrate are joined on barcode from their source file,
with the join asserted complete. The stale `published_*` columns ride along undrawn.

Run in-container from the compartment root, after 17_treg_reembedding.py and
16_narrative_scoring.py:
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

COMPARTMENT_ROOT = Path(__file__).resolve().parents[2]

sys.path.insert(0, str(COMPARTMENT_ROOT))
sys.path.insert(0, str(COMPARTMENT_ROOT / "02_analysis"))
os.chdir(COMPARTMENT_ROOT)

from config import (  # noqa: E402
    ARM_STRIP_SETS as CONFIG_ARM_STRIP_SETS,
    PATHS,
    PROGRAM_STRIP_SETS,
    TISSUE_COLORS,
)
from helpers.figure_style import (  # noqa: E402
    FIG_CFG, save_overview, set_paper_style, write_caption,
)
from helpers import umap_grid as ug  # noqa: E402

STAGE = "17_treg_reembedding"
SOURCE_STAGE = "16_narrative_scoring"
SCRIPT = "02_analysis/scripts/17_treg_reembedding_viz.py"
SUBSTRATE = "17_treg_reembedding.parquet"
SUBSTRATE_REL = f"03_results/interactive/{SUBSTRATE}"
SUMMARY_REL = "03_results/16_narrative_scoring/tables/narrative_score_summary.csv"
MANIFEST_REL = "03_results/16_narrative_scoring/tables/narrative_scoring_manifest.csv"
MIXING_REL = "03_results/17_treg_reembedding/tables/treg_reembedding_mixing.csv"
MARKER_SUBSTRATE_REL = "03_results/07_embedding/tables/hook_factor_substrate.parquet"

N_COL = 6
SAMPLE_N = None         # 27,175 Treg cells is a drawable density, so every cell is drawn
SAMPLE_SEED = 0

_OI = (FIG_CFG.get("colors", {}) or {}).get("okabe_ito", {}) or {}
TISSUE_COL = TISSUE_COLORS
TISSUE_LABEL = {"synovial_fluid": "synovial fluid", "peripheral_blood": "paired blood"}
DONOR_HUES = ["vermillion", "sky_blue", "bluish_green", "yellow",
              "blue", "reddish_purple", "black"]

COORD_LINE = ("Coordinates: Treg-only UMAP with Harmony over donor applied to the PCA "
              "(columns x / y of the substrate).")
COUNTERPART_LINE = ("Colour limits are the full-object figure's, so the pair is comparable "
                    "panel for panel.")
TIER_LINE = ("Annotation tier. Claims in this compartment rest on donor-level pseudobulk "
             "differential expression within the frozen cell states.")

# --- the reference strip: tissue and donor as bookends, four markers between ---
# The sort-gate panel of the full-object twin is one category here, so IKZF2 takes the slot.
MARKER_PANELS = [
    ("FOXP3", "FOXP3\nTreg lineage TF"),
    ("IL2RA", "IL2RA\nCD25"),
    ("CTLA4", "CTLA4\nsuppressive effector"),
    ("IKZF2", "IKZF2\nHelios"),
]
MARKER_CHANNELS = [c for c, _ in MARKER_PANELS]
REFERENCE_BANDS = [
    ("tissue", 0, 0, True),
    ("core Treg identity genes", 1, 4, True),
    ("donor of origin", 5, 5, False),
]

# --- the two counterpart strips: the same sets, in the same order, as the full object ---
# From `analysis_config.yaml::percell_map_panels`, so both strips draw one list.
ARM_STRIP_SETS = list(CONFIG_ARM_STRIP_SETS)
ARM_SETS = ARM_STRIP_SETS[:3]
ARM_LENS_SETS = ARM_STRIP_SETS[3:]
ARM_SHARED_SCALE = ["WT_heat_up_AUCell", "KO_heat_up_AUCell"]
ARM_BANDS = [
    ("mouse 39 °C-derived up arms — anchor-dependent", 0, 2, True),
    ("curated, anchor-independent", 3, 4, True),
    ("derived here from GSE161426", 5, 5, False),   # stage 16's wording, so the pair matches
]

PROGRAM_SETS = list(PROGRAM_STRIP_SETS)
PROGRAM_BANDS = [
    ("cGAS-STING and type-I interferon", 0, 2, True),
    ("inflammation and activation", 3, 5, False),
]

# --- the full-object twin each counterpart is scaled against -------------------
# The frame is stage 16's, cell for cell: same substrate, size and seed, so the limits derived
# here are the limits its colourbars carry.
FULL_SUBSTRATE = "16_narrative_embedding.parquet"
FULL_SUBSTRATE_REL = f"03_results/interactive/{FULL_SUBSTRATE}"
FULL_SAMPLE_N = 60_000
FULL_SAMPLE_SEED = 0
# The shared arm scale stage 16 records for umap_full_arms, to 4 decimals. Reproducing it is
# the check that this script reads the same frame; a drift means the pair has stopped sharing
# a scale.
FULL_ARM_SCALE_COMMITTED = (0.0023, 0.0328)

# The three continuous channels of 07_embedding/umap_signatures_treg: scanpy `score_genes`
# module scores, mean-centred and on a different scale from AUCell, living in the 07_embedding
# substrate upstream of this stage.
SIGNATURE_PANELS = [
    ("WT_heat_up", "WT_heat_up\nmouse anchor annotation"),
    ("score_eTreg", "score_eTreg\neffector-Treg"),
    ("score_HSP", "score_HSP\nheat-shock / stress"),
]
SIGNATURE_CHANNELS = [c for c, _ in SIGNATURE_PANELS]


# =============================================================================
# Inputs
# =============================================================================
def join_on_barcode(frame: pd.DataFrame, source: pd.DataFrame, cols: list,
                    source_rel: str) -> pd.DataFrame:
    """Add `cols` to `frame`, matched on barcode, halting on an incomplete join.

    The values come from the source file so a counterpart colours by the same numbers as its
    twin. A reindexed column would put one cell's value at another cell's coordinates.
    """
    take = source[cols].copy()
    take.index = source.index.astype(str)
    out = frame.join(take, on="barcode")
    missing = int(out[cols].isna().any(axis=1).sum())
    if missing:
        raise AssertionError(
            f"{missing} of {len(out)} cells carry no {source_rel} value for {cols}; the two "
            "substrates disagree on barcodes, so the panels cannot be drawn.")
    return out


def with_all_channels(frame: pd.DataFrame, full: pd.DataFrame,
                      markers: pd.DataFrame) -> pd.DataFrame:
    """Every channel these strips draw, on one frame: scores from stage 16, genes from 07."""
    needed = [f"{s}_AUCell" for s in ARM_STRIP_SETS + PROGRAM_SETS]
    absent = [c for c in needed if c not in frame.columns]
    if absent:
        frame = join_on_barcode(frame, full.set_index("barcode"), absent, FULL_SUBSTRATE_REL)
        print(f"[17_treg_reembedding_viz] joined {len(absent)} score column(s) from "
              f"{FULL_SUBSTRATE_REL}: {', '.join(absent)}")
    return join_on_barcode(frame, markers, MARKER_CHANNELS + SIGNATURE_CHANNELS,
                           MARKER_SUBSTRATE_REL)


def set_sizes(manifest: pd.DataFrame) -> dict:
    """{set_name: genes scored in the object}, so a pair's titles agree."""
    return dict(zip(manifest["set_name"], manifest["n_genes_found_in_object"].astype(int)))


def score_panels(sets: list, sizes: dict) -> list:
    """(column, title) per scored set, titles of uniform height across the strip."""
    return list(zip((f"{s}_AUCell" for s in sets), ug.set_titles(sets, sizes)))


def full_object_frame(markers: pd.DataFrame) -> pd.DataFrame:
    """The frame the full-object figures draw, reproduced cell for cell, markers joined.

    A clip taken from a different frame is a different clip, and the pair would stop sharing
    a scale.
    """
    full = pd.read_parquet(PATHS.interactive / FULL_SUBSTRATE)
    d = full.sample(n=FULL_SAMPLE_N, random_state=FULL_SAMPLE_SEED)
    d = d.sample(frac=1.0, random_state=FULL_SAMPLE_SEED).reset_index(drop=True)
    return join_on_barcode(d, markers, MARKER_CHANNELS, MARKER_SUBSTRATE_REL)


def arm_clips(frame: pd.DataFrame, panels: list) -> dict:
    """Per-column limits for the arm strip, WT and KO pooled onto one scale.

    The pooled pair is checked against the number stage 16 records. That assertion is the
    guarantee the pair shares a scale: a changed substrate or sample size would otherwise
    produce a second scale while the two figures still looked like a matched set.
    """
    shared = ug.clip_of(frame, ARM_SHARED_SCALE)
    want = FULL_ARM_SCALE_COMMITTED
    if (round(shared[0], 4), round(shared[1], 4)) != want:
        raise AssertionError(
            f"the WT_heat_up + KO_heat_up shared scale reproduces as "
            f"{shared[0]:.4f}-{shared[1]:.4f}, while this script records "
            f"{want[0]:.4f}-{want[1]:.4f} for umap_full_arms. The counterpart would stop "
            "sharing a colour scale with its full-object twin; reconcile the frame and "
            "update FULL_ARM_SCALE_COMMITTED before drawing.")
    return {col: (shared if col in ARM_SHARED_SCALE else ug.clip_of(frame, [col]))
            for col, _ in panels}


# =============================================================================
# Source tables
# =============================================================================
def score_table(summary: pd.DataFrame, sets: list, label: str = "Treg") -> pd.DataFrame:
    """The committed per (set x tissue) summary for one gate, restricted to the panels drawn."""
    tbl = summary[(summary["set_name"].isin(sets)) & (summary["coarse_label"] == label)].copy()
    tbl["set_name"] = pd.Categorical(tbl["set_name"], categories=sets, ordered=True)
    return tbl.sort_values(["set_name", "tissue"]).reset_index(drop=True)


def channel_table(drawn: pd.DataFrame, channels: list, metric: str) -> pd.DataFrame:
    """Per (channel x tissue) descriptive summary within the Treg gate.

    Mirrors the shape of `narrative_score_summary.csv` and carries `metric`, because module
    scores and AUCell share no scale.
    """
    rows = []
    for channel in channels:
        for tissue, grp in drawn.groupby("tissue", observed=True):
            v = grp[channel].to_numpy(dtype=float)
            rows.append({
                "set_name": channel, "coarse_label": "Treg", "tissue": tissue,
                "n_cells": len(grp), "n_donors": int(grp["donor"].nunique()),
                "mean": float(np.nanmean(v)), "median": float(np.nanmedian(v)),
                "sd": float(np.nanstd(v, ddof=1)),
                "metric": metric, "evidence_tier": "secondary_percell",
            })
    out = pd.DataFrame(rows)
    out["set_name"] = pd.Categorical(out["set_name"], categories=channels, ordered=True)
    return out.sort_values(["set_name", "tissue"]).reset_index(drop=True)


def _mean_of(summary: pd.DataFrame, set_name: str, tissue: str, label: str = "Treg") -> float:
    hit = summary[(summary["set_name"] == set_name) & (summary["coarse_label"] == label)
                  & (summary["tissue"] == tissue)]
    if hit.empty:
        raise ValueError(f"no summary row for {set_name} / {label} / {tissue}")
    return float(hit["mean"].iloc[0])


def _median_of(summary: pd.DataFrame, set_name: str, tissue: str, label: str = "Treg") -> float:
    hit = summary[(summary["set_name"] == set_name) & (summary["coarse_label"] == label)
                  & (summary["tissue"] == tissue)]
    if hit.empty:
        raise ValueError(f"no summary row for {set_name} / {label} / {tissue}")
    return float(hit["median"].iloc[0])


def _shift(summary: pd.DataFrame, set_name: str) -> str:
    pb = _mean_of(summary, set_name, "peripheral_blood")
    sf = _mean_of(summary, set_name, "synovial_fluid")
    return f"{set_name} {pb:.4f} to {sf:.4f}"


def _tissue_shift(tbl: pd.DataFrame, name: str) -> str:
    def m(tissue: str) -> float:
        hit = tbl[(tbl["set_name"] == name) & (tbl["tissue"] == tissue)]
        return float(hit["mean"].iloc[0])
    return f"{name} {m('peripheral_blood'):.4f} to {m('synovial_fluid'):.4f}"


def _n_zero_median(summary: pd.DataFrame, set_name: str) -> tuple[int, int]:
    """How many (cell state x tissue) rows of `set_name` carry a median of exactly zero."""
    rows = summary[summary["set_name"] == set_name]
    return int((rows["median"] == 0).sum()), int(len(rows))


def _mixing(mix: pd.DataFrame, embedding: str, key: str) -> tuple[float, float]:
    hit = mix[(mix["embedding"] == embedding) & (mix["grouping_key"] == key)
              & (mix["group"] == "_all_")]
    if hit.empty:
        raise ValueError(f"no mixing row for {embedding} / {key} / _all_")
    return (float(hit["observed_same_frac"].iloc[0]),
            float(hit["expected_same_frac"].iloc[0]))


# =============================================================================
# Strips
# =============================================================================
def _new_strip(n_row: int = 1, n_col: int = N_COL):
    fig = plt.figure(figsize=ug.canvas(n_col, n_row))
    return fig, ug.grid(fig, n_col, n_row)


def draw_reference_row(axes: list, drawn: pd.DataFrame, xlim, ylim, marker_clip: tuple):
    """Tissue, four Treg identity genes, donor — over one shared frame.

    The markers are one quantity, so they pool onto one clip — the full-object twin's, as every
    limit here is — and read on one bar in REAL units, where level compares between the genes.
    Stage 16 pools its three the same way.
    """
    donors = sorted(drawn["donor"].unique())
    donor_col = {d: _OI[h] for d, h in zip(donors, DONOR_HUES)}
    donor_label = {d: d.replace("JIA_patient_", "JIA ") for d in donors}

    ug.scatter_cat(axes[0], drawn, "tissue", TISSUE_COL, "Tissue of origin\nsampling site",
                   xlim, ylim, FIG_CFG, order=["synovial_fluid", "peripheral_blood"],
                   labels=TISSUE_LABEL, ncol=2)
    sc = None
    for ax, (col, title) in zip(axes[1:5], MARKER_PANELS):
        sc = ug.scatter_cont(ax, drawn, col, title, xlim, ylim, FIG_CFG,
                             vlim=marker_clip, unit=ug.UNIT_EXPR, colorbar=False)
    ug.row_colorbar(axes[4].figure, axes[4], sc, ug.UNIT_EXPR, FIG_CFG)
    ug.scatter_cat(axes[5], drawn, "donor", donor_col, "Donor\n7 JIA patients", xlim, ylim,
                   FIG_CFG, order=donors, labels=donor_label, ncol=4)


def draw_score_row(axes: list, drawn: pd.DataFrame, xlim, ylim, panels: list, clips: dict,
                   unit: str = ug.UNIT_AUCELL, rescale: bool = True):
    """One row of continuous panels, each clipped to the limits it was handed.

    `rescale` maps each panel across those limits onto [0, 1] and gives the row ONE bar, as
    stage 16 does. The limits stay the full-object twin's, so the pairing survives.

    The signature strip passes `rescale=False`: module scores against a stage-07 counterpart
    drawn in its own units.
    """
    sc = None
    for ax, (col, title) in zip(axes, panels):
        sc = ug.scatter_cont(ax, drawn, col, title, xlim, ylim, FIG_CFG,
                             vlim=clips[col], unit=ug.UNIT_RESCALED if rescale else unit,
                             rescale=rescale, colorbar=not rescale)
    if rescale:
        ug.row_colorbar(axes[-1].figure, axes[-1], sc, ug.UNIT_RESCALED, FIG_CFG)


def clip_note(panels: list, clips: dict) -> str:
    """The limits each panel was rescaled over. The bar reads 0 to 1, so the caption needs
    these to be interpretable."""
    return "; ".join(f"{col} {clips[col][0]:.4f}-{clips[col][1]:.4f}" for col, _ in panels)


# =============================================================================
def main() -> None:
    set_paper_style(config=FIG_CFG)

    treg = pd.read_parquet(PATHS.interactive / SUBSTRATE)
    full = pd.read_parquet(PATHS.interactive / FULL_SUBSTRATE)
    markers = pd.read_parquet(COMPARTMENT_ROOT / MARKER_SUBSTRATE_REL)
    summary = pd.read_csv(PATHS.tables(SOURCE_STAGE) / "narrative_score_summary.csv")
    manifest = pd.read_csv(PATHS.tables(SOURCE_STAGE) / "narrative_scoring_manifest.csv")
    mix = pd.read_csv(PATHS.tables(STAGE) / "treg_reembedding_mixing.csv")
    sizes = set_sizes(manifest)
    print(f"[17_treg_reembedding_viz] substrate {treg.shape[0]:,} Treg cells x {treg.shape[1]} "
          f"cols, {treg['donor'].nunique()} donors, "
          f"cell states {sorted(treg['coarse_label'].unique())}")

    donor_raw, _ = _mixing(mix, "treg_only", "donor")
    donor_harmony, donor_chance = _mixing(mix, "treg_only_harmony", "donor")
    tissue_harmony, tissue_chance = _mixing(mix, "treg_only_harmony", "tissue")
    print(f"[17_treg_reembedding_viz] same-donor neighbours at k=30: uncorrected {donor_raw:.3f}, "
          f"Harmony {donor_harmony:.3f}, chance {donor_chance:.3f}")

    treg = with_all_channels(treg, full, markers)
    drawn = ug.sample_frame(treg, SAMPLE_N, SAMPLE_SEED)
    xlim, ylim = ug.data_box(drawn)

    # Every colour limit comes from the full-object frame, so each panel matches its twin.
    twin = full_object_frame(markers)
    arm_panels = score_panels(ARM_STRIP_SETS, sizes)
    a_clips = arm_clips(twin, arm_panels)
    program_panels = score_panels(PROGRAM_SETS, sizes)
    p_clips = {col: ug.clip_of(twin, [col]) for col, _ in program_panels}
    marker_clip = ug.clip_of(twin, MARKER_CHANNELS)   # one clip over all four, one bar
    s_clips = {c: ug.clip_of(markers, [c]) for c in SIGNATURE_CHANNELS}
    print(f"[17_treg_reembedding_viz] full-object frame {len(twin):,} cells; arm shared scale "
          f"{a_clips['WT_heat_up_AUCell'][0]:.4f} to {a_clips['WT_heat_up_AUCell'][1]:.4f} "
          "reproduces the recorded value")

    fw, fh = ug.canvas(N_COL)
    cfg_common = (f"coordinates = x / y (Harmony over donor), all {len(drawn):,} cells drawn, "
                  f"panel = {ug.PANEL_W_IN:.2f} x {ug.PANEL_H_IN:.2f} in x {N_COL} columns, "
                  f"point_size = {ug.point_size(len(drawn))}, cmap = {ug.CMAP}, "
                  "figures.dpi = 300, figures.rasterized_dpi = 600, colour limits from "
                  f"{FULL_SUBSTRATE_REL} at sample_n = {FULL_SAMPLE_N}, "
                  f"sample_seed = {FULL_SAMPLE_SEED}, clip_percentiles = {list(ug.CLIP)}")
    pair_line = (
        "Every colour limit is the full-object figure's, from the same frame and seed, so the "
        "pair compares brightness for brightness and a washed-out counterpart panel states "
        "something real about the Treg gate's range. The pair shares its cells and its colour "
        "scale; the coordinates differ, because these are a re-embedding of the Treg cells "
        "alone. Cells are pooled across donors, so a tissue difference read off the colouring "
        "is pseudoreplicated, and Harmony reshapes the space it corrects, so this is a map. "
        "Claims in this compartment rest on donor-level pseudobulk differential expression "
        "within the frozen cell states.")
    pair_line_short = (
        "Limits, coordinates and tier follow the counterpart contract at the top of this page: "
        "limits are the full-object figure's, the coordinates are this map's own, cells pool "
        "across donors so a tissue difference is pseudoreplicated, and claims rest on "
        "donor-level pseudobulk differential expression within the frozen cell states.")

    # --- 1. reference strip --------------------------------------------------
    fig, rows = _new_strip()
    draw_reference_row(rows[0], drawn, xlim, ylim, marker_clip)
    ug.group_bands(fig, rows[0], REFERENCE_BANDS, FIG_CFG, colour=_OI["black"])
    ug.dress(fig, "Treg-only map, Harmony corrected over donor",
             f"{len(drawn):,} sorted Treg cells, one frame shared by all six panels. Marker "
             "colour is log-normalised expression on one bar, pooled over the four genes at the "
             "full-object figure's limits.",
             [COORD_LINE, TIER_LINE], FIG_CFG, colour=_OI["black"])
    ref_tbl = channel_table(drawn, MARKER_CHANNELS, "log_normalised_expression")
    save_overview(
        fig, STAGE, "umap_treg_reembedding",
        table=ref_tbl,
        finding=("On the Treg-only map, drawn on the Harmony-corrected coordinates, the synovial-"
                 f"fluid and paired-blood cells still occupy distinct territory — "
                 f"{tissue_harmony:.3f} same-tissue neighbours at k = 30 against "
                 f"{tissue_chance:.3f} expected — after the same-donor neighbour fraction has "
                 f"fallen from {donor_raw:.3f} to {donor_harmony:.3f} against "
                 f"{donor_chance:.3f} expected. The four Treg identity genes hold across the map "
                 f"rather than marking a corner of it ({_tissue_shift(ref_tbl, 'FOXP3')}, "
                 f"{_tissue_shift(ref_tbl, 'IKZF2')} in mean log-normalised expression), so the "
                 "layout separates tissue while the gate stays uniform."),
        script=SCRIPT, fn="draw_reference_row",
        config_kv=(f"{cfg_common}, columns = tissue, {', '.join(MARKER_CHANNELS)}, donor, "
                   f"marker_scale = one pooled full-object clip at {marker_clip[0]:.4f}-"
                   f"{marker_clip[1]:.4f} over all {len(MARKER_CHANNELS)} genes, one bar"),
        input=f"{SUBSTRATE_REL}, {MARKER_SUBSTRATE_REL}, {MIXING_REL}",
        how_to_read=(
            f"Six panels over one frame of the same {len(drawn):,} sorted Treg cells, "
            f"sharing one bounding box. Tissue and donor are the bookends, both drawn in "
            f"shuffled order so overlapping groups paint evenly. Between them, four Treg "
            f"identity genes in log-normalised expression: FOXP3 the lineage transcription "
            f"factor, IL2RA the CD25 chain, CTLA4 the suppressive effector, IKZF2 the Helios "
            f"subset marker. The full-object twin carries a Treg/Tcon/CD8 sort-gate panel in "
            f"the second slot, which is a single category here, so IKZF2 takes it. All four "
            f"pool onto ONE clip and ONE bar in real units, {marker_clip[0]:.2f} to "
            f"{marker_clip[1]:.2f}, so brightness compares between the genes as well as "
            f"within one; that clip is the full-object figure's, so a panel brighter here is "
            f"brighter, and the values are 07_embedding's own, joined on barcode. "
            f"Coordinates are the "
            f"Harmony-corrected pair: correction over donor takes same-donor neighbours from "
            f"{donor_raw:.3f} to {donor_harmony:.3f} at k = 30 against {donor_chance:.3f} "
            f"expected while same-tissue neighbours hold at {tissue_harmony:.3f} against "
            f"{tissue_chance:.3f}, so it acted on donor and left the tissue separation "
            f"standing. Harmony reshapes the space it corrects, so this map is annotation. "
            f"Cells pool across donors, so a tissue difference is pseudoreplicated. Claims "
            f"rest on donor-level pseudobulk differential expression within the frozen cell "
            f"states."),
        width=fw, height=fh, config=FIG_CFG)
    plt.close(fig)

    # --- 2. counterpart of the arm strip -------------------------------------
    fig, rows = _new_strip()
    draw_score_row(rows[0], drawn, xlim, ylim, arm_panels, a_clips)
    ug.group_bands(fig, rows[0], ARM_BANDS, FIG_CFG, colour=_OI["black"])
    ug.dress(fig, "Mouse 39 °C-derived up arms and their lenses, on the Treg-only map",
             f"{len(drawn):,} sorted Treg cells, one frame shared by all six panels. Colour is "
             f"per-cell AUCell rescaled to each panel's full-object clip, so all six read on one "
             f"0-to-1 bar; WT_heat_up and KO_heat_up are rescaled together "
             f"({a_clips['WT_heat_up_AUCell'][0]:.3f} to "
             f"{a_clips['WT_heat_up_AUCell'][1]:.3f}).",
             [COORD_LINE, COUNTERPART_LINE, TIER_LINE], FIG_CFG, colour=_OI["black"])
    save_overview(
        fig, STAGE, "umap_treg_arms",
        table=score_table(summary, ARM_STRIP_SETS),
        finding=("Read on the Treg gate's own map and on the full object's colour scale, all "
                 "three mouse 39 °C-derived up arms still colour the synovial-fluid territory "
                 "brighter than the paired-blood territory, so the tissue contrast the "
                 "full-object row shows survives viewing Treg apart from Tcon and CD8 — and the "
                 "three anchor-independent lenses beside them do the same. Per-cell AUCell means "
                 f"run {_shift(summary, 'WT_heat_up')}, {_shift(summary, 'KO_heat_up')} and "
                 f"{_shift(summary, 'Interaction_up')} for the arms, against "
                 f"{_shift(summary, 'HALLMARK_HYPOXIA')}, {_shift(summary, 'HSR_core')} and "
                 f"{_shift(summary, 'eTreg_up')} for the lenses."),
        script=SCRIPT, fn="draw_score_row",
        config_kv=(f"{cfg_common}, columns = " + ", ".join(c for c, _ in arm_panels)
                   + f", shared_scale = WT_heat_up + KO_heat_up at "
                     f"{a_clips['WT_heat_up_AUCell'][0]:.4f}-"
                     f"{a_clips['WT_heat_up_AUCell'][1]:.4f}"
                   + ", colour = rescaled to panel clip onto [0, 1], one bar for the row; "
                     "AUCell limits rescaled over: " + clip_note(arm_panels, a_clips)),
        input=f"{SUBSTRATE_REL}, {FULL_SUBSTRATE_REL}, {SUMMARY_REL}, {MANIFEST_REL}",
        how_to_read=(
            f"The Treg-only counterpart of "
            f"`16_narrative_scoring/figures/_overview/umap_full_arms.png`: same six sets, "
            f"same order, same geometry and colormap, on the {len(drawn):,} sorted Treg "
            f"cells alone. The vertical rules carry the provenance. The left three are "
            f"anchor-dependent: WT_heat_up, the mouse WT iTreg 39-versus-37 °C up arm in "
            f"human projection; KO_heat_up, the same contrast in cGAS-knockout iTregs; "
            f"Interaction_up, the genotype-by-temperature arm at {sizes['Interaction_up']} "
            f"genes, small enough that one gene moves the score, so read it for location and "
            f"treat its spread as noise. The middle two are curated and versioned, "
            f"HALLMARK_HYPOXIA and the activation-free HSR_core proteostasis lens. The right "
            f"one is eTreg_up, this compartment's own GSE161426 effector-Treg contrast, "
            f"derived for exploration and ruled off for that reason. ONE BAR SERVES THE "
            f"ROW, and it reads 0 to 1, not AUCell: each panel is clipped to its "
            f"full-object twin's limits as before and then rescaled across that clip, so "
            f"the picture is unchanged and the twin comparison holds panel for panel. "
            f"Brightness compares tissue WITHIN a panel and says nothing about level "
            f"BETWEEN panels; the limits are in the config line below. WT_heat_up and "
            f"KO_heat_up share 182 genes and one clip, so those two do compare pixel for "
            f"pixel. Titles give the symbols scored, and the source table carries "
            f"the means, spreads and counts the colour cannot show. {pair_line_short}"),
        width=fw, height=fh, config=FIG_CFG)
    plt.close(fig)

    # --- 3. counterpart of the curated-lens strip ----------------------------
    fig, rows = _new_strip()
    draw_score_row(rows[0], drawn, xlim, ylim, program_panels, p_clips)
    ug.group_bands(fig, rows[0], PROGRAM_BANDS, FIG_CFG, colour=_OI["black"])
    ug.dress(fig, "Curated program lenses on the Treg-only map",
             f"{len(drawn):,} sorted Treg cells, one frame shared by all six panels. Each panel "
             "carries its full-object twin's colour limits, rescaled onto one shared 0-to-1 bar; "
             "AUCell values are in the source table.",
             [COORD_LINE, COUNTERPART_LINE, TIER_LINE], FIG_CFG, colour=_OI["black"])
    n_zero, n_rows = _n_zero_median(summary, "sting_specific_published")
    save_overview(
        fig, STAGE, "umap_treg_programs",
        table=score_table(summary, PROGRAM_SETS),
        finding=("All six curated lenses colour synovial-fluid territory brighter than paired "
                 "blood on the Treg gate's own map, so that structure appears without the other "
                 f"two sort gates: per-cell AUCell means {_shift(summary, 'ifn_generic_axis')}, "
                 f"{_shift(summary, 'HALLMARK_INTERFERON_ALPHA_RESPONSE')}, "
                 f"{_shift(summary, 'HALLMARK_TNFA_SIGNALING_VIA_NFKB')} and "
                 f"{_shift(summary, 'sting_specific_published')}. The shared clip carries the "
                 "level, which the rescaled bar does not: the published STING panel is scored "
                 "over the whole Treg map against a full-object limit it barely reaches, Treg "
                 "blood mean "
                 f"{_mean_of(summary, 'sting_specific_published', 'peripheral_blood'):.4f} "
                 "against Tcon "
                 f"{_mean_of(summary, 'sting_specific_published', 'peripheral_blood', 'Tcon'):.4f} "
                 "and CD8 "
                 f"{_mean_of(summary, 'sting_specific_published', 'peripheral_blood', 'CD8'):.4f}. "
                 "Its median is exactly "
                 f"{_median_of(summary, 'sting_specific_published', 'peripheral_blood'):.3f}, the "
                 f"only {n_zero} of {n_rows} rows at zero, so at least half of Treg blood cells "
                 "score zero on that set and the Treg synovial-versus-blood difference on it "
                 "rests partly on a zero-inflated blood baseline."),
        script=SCRIPT, fn="draw_score_row",
        config_kv=(f"{cfg_common}, columns = " + ", ".join(c for c, _ in program_panels)
                   + "; per-panel full-object limits, no pooling across sets"
                   + ", colour = rescaled to panel clip onto [0, 1], one bar for the row; "
                     "AUCell limits rescaled over: " + clip_note(program_panels, p_clips)),
        input=f"{SUBSTRATE_REL}, {FULL_SUBSTRATE_REL}, {SUMMARY_REL}, {MANIFEST_REL}",
        how_to_read=(
            f"The Treg-only counterpart of "
            f"`16_narrative_scoring/figures/_overview/umap_full_programs.png`: same six "
            f"sets, same order, same geometry and colormap, on the {len(drawn):,} sorted "
            f"Treg cells alone. Every set here is curated, versioned and derived "
            f"independently of the mouse anchor, so a colouring stands apart from the "
            f"anchor. The rule splits two families, the cGAS-STING and type-I interferon "
            f"sets on the left and the inflammation and activation programs on the right. "
            f"Reading them together is the point: a synovial-high colouring shared by both "
            f"families is generic inflammation, and only a pattern the left family carries "
            f"and the right one lacks would be specific to STING or interferon. The six "
            f"differ in range, so each panel keeps its own set's limits; what is shared is "
            f"the object those limits come from, which makes each panel comparable to its "
            f"own twin. ONE BAR SERVES THE ROW, and it reads 0 to 1, not AUCell: each panel "
            f"is rescaled across the twin's clip it was already drawn on, so the picture is "
            f"unchanged, the twin comparison holds, and brightness compares tissue WITHIN a "
            f"panel while saying nothing about level BETWEEN panels. The limits are in the "
            f"config line below and the values in the source table. The published STING set "
            f"is 21 genes, "
            f"{sizes['sting_specific_published']} of them scored here, and its IFN-β "
            f"validation in the positive-control compartment is underpowered at three "
            f"donors, so a bright or dim panel there is consistent with STING pathway "
            f"activity and never proof of it. Hypoxia and temperature are both imposed by "
            f"the inflamed joint and stay entangled in cross-sectional human data, so these "
            f"lenses carry no HIF claim. {pair_line_short}"),
        width=fw, height=fh, config=FIG_CFG)
    plt.close(fig)

    # --- 4. counterpart of the 07_embedding candidate signatures -------------
    sw, sh = ug.canvas(len(SIGNATURE_PANELS))
    fig, rows = _new_strip(n_col=len(SIGNATURE_PANELS))
    draw_score_row(rows[0], drawn, xlim, ylim, SIGNATURE_PANELS, s_clips,
                   unit=ug.UNIT_MODULE, rescale=False)
    ug.dress(fig, "Candidate harvest signatures on the Treg-only map",
             f"{len(drawn):,} sorted Treg cells, one frame shared by all three panels. Colour is "
             "the scanpy module score, on the full-object figure's limits.",
             [COORD_LINE, COUNTERPART_LINE, TIER_LINE], FIG_CFG, colour=_OI["black"])
    sig_tbl = channel_table(drawn, SIGNATURE_CHANNELS, "scanpy_score_genes_module_score")
    save_overview(
        fig, STAGE, "umap_treg_signatures",
        table=sig_tbl,
        finding=("The three candidate harvest signatures on the Treg gate's own map and on the "
                 "full-object figure's limits. Within the Treg gate the effector-Treg and "
                 "heat-shock module scores still separate synovial fluid from paired blood "
                 f"({_tissue_shift(sig_tbl, 'score_eTreg')}, "
                 f"{_tissue_shift(sig_tbl, 'score_HSP')}), alongside the mouse anchor "
                 f"annotation ({_tissue_shift(sig_tbl, 'WT_heat_up')}), so the structure the "
                 "full-object figure shows for these three channels appears without viewing "
                 "the Treg gate against Tcon and CD8."),
        script=SCRIPT, fn="draw_score_row",
        config_kv=(f"{cfg_common}, columns = " + ", ".join(SIGNATURE_CHANNELS)
                   + f"; joined on barcode from {MARKER_SUBSTRATE_REL}; limits from all "
                     f"{len(markers):,} cells of that substrate, which is the frame "
                     "07_embedding_viz.py draws; metric = scanpy score_genes module score"),
        input=f"{SUBSTRATE_REL}, {MARKER_SUBSTRATE_REL}",
        how_to_read=(
            f"The Treg-only counterpart of "
            f"`07_embedding/figures/_overview/umap_signatures_treg.png`, which draws these "
            f"same three channels across all three sort gates. Two differences matter first. "
            f"The unit: these are scanpy `score_genes` module scores, mean-centred against a "
            f"sampled background and signed, so they share a scale neither with the AUCell "
            f"panels here nor with each other, and a value near zero means at background — "
            f"the bounded AUCell readings of the same two programs are the WT_heat_up and "
            f"eTreg_up panels of `umap_treg_arms`. The panel count: the full-object twin "
            f"carries a fourth panel, a Treg/Tcon/CD8 sort-gate reference, which is a single "
            f"category here, and the tissue reference sits in `umap_treg_reembedding`. "
            f"WT_heat_up here is the mouse WT 39-versus-37 °C up arm carried as annotation "
            f"only, never a selection predicate, and the harvest design it was previewed for "
            f"is frozen as implemented. Limits come from all {len(markers):,} cells of the "
            f"07_embedding substrate, the frame that figure draws. {pair_line_short}"),
        width=sw, height=sh, config=FIG_CFG)
    plt.close(fig)

    # --- 5. the two strips on one canvas, for layout -------------------------
    pw_h = ug.canvas(N_COL, 2)[1]
    fig, rows = _new_strip(n_row=2)
    draw_reference_row(rows[0], drawn, xlim, ylim, marker_clip)
    ug.group_bands(fig, rows[0], REFERENCE_BANDS, FIG_CFG, colour=_OI["black"])
    draw_score_row(rows[1], drawn, xlim, ylim, arm_panels, a_clips)
    ug.group_bands(fig, rows[1], ARM_BANDS, FIG_CFG, colour=_OI["black"])
    ug.dress(fig, "The Treg-only map: the reference layout, and the signatures projected onto it",
             f"All {len(drawn):,} sorted Treg cells, one frame shared by all twelve panels. Top "
             "row tissue, Treg identity and donor, the four genes pooled onto one expression "
             "bar; bottom row per-cell AUCell rescaled per panel onto one 0-to-1 bar.",
             [COORD_LINE, COUNTERPART_LINE, TIER_LINE], FIG_CFG, colour=_OI["black"])
    save_overview(
        fig, STAGE, "umap_treg_patchwork",
        table=pd.concat([ref_tbl, score_table(summary, ARM_STRIP_SETS)], ignore_index=True),
        finding=("The Treg-only reference layout and its signature colouring on one canvas, so "
                 "the tissue separation, the identity genes and the mouse-derived arms are read "
                 "against each other in one view: the arms brighten the synovial-fluid territory "
                 "that the tissue panel marks, while the identity genes stay uniform across it."),
        script=SCRIPT, fn="draw_reference_row",
        config_kv=(f"{cfg_common}, rows = 2 x {N_COL}, canvas = {fw:.1f} x {pw_h:.1f} in, "
                   "top row = " + ", ".join(["tissue"] + MARKER_CHANNELS + ["donor"])
                   + "; bottom row = " + ", ".join(c for c, _ in arm_panels)
                   + ", rescaled to panel clip onto [0, 1] on one bar; AUCell limits rescaled "
                     "over: " + clip_note(arm_panels, a_clips)),
        input=f"{SUBSTRATE_REL}, {FULL_SUBSTRATE_REL}, {SUMMARY_REL}, {MANIFEST_REL}, "
              f"{MARKER_SUBSTRATE_REL}",
        how_to_read=(
            f"The two strips this stage ships separately, `umap_treg_reembedding` above "
            f"`umap_treg_arms`, on one canvas at identical panel size so a column reads top "
            f"to bottom. Nothing new is drawn, and both rows hold the identical frame of "
            f"cells at identical coordinates. The rows share cells and coordinates, and "
            f"their units differ: the top row is categorical annotation plus log-normalised "
            f"expression, the bottom row per-cell AUCell of a gene set rescaled per panel "
            f"onto a single 0-to-1 bar. Each row's own "
            f"caption carries its full reading, the bottom row's rules carry each "
            f"signature's provenance, and the source table stacks both rows' summaries with "
            f"a `metric` column separating them. {pair_line_short}"),
        width=fw, height=pw_h, config=FIG_CFG)
    plt.close(fig)

    # --- 6. captions for the same-stem source tables -------------------------
    write_caption(
        STAGE, "tables/_overview/umap_treg_reembedding.csv",
        finding=(f"Per (gene x tissue) summaries of the {len(MARKER_CHANNELS)} Treg identity "
                 "genes drawn on the Treg-only map, giving the numbers behind the marker "
                 "panels: all four sit close between the two tissues, so the gate holds across "
                 "the map."),
        script=SCRIPT, fn="channel_table",
        config_kv=(f"rows = {len(MARKER_CHANNELS)} genes x Treg x 2 tissues, "
                   "metric = log_normalised_expression"),
        input=f"{SUBSTRATE_REL}, {MARKER_SUBSTRATE_REL}",
        how_to_read=("One row per (`set_name` x `tissue`) over the drawn Treg cells, with the "
                     "mean, median and standard deviation of that gene's log-normalised "
                     "expression and the cell and donor counts behind it. `metric` reads "
                     "`log_normalised_expression`, so these values sit on a different scale "
                     "from the AUCell tables in this directory. Most cells carry zero for each "
                     "of these genes, which is ordinary for single-cell counts, so a median of "
                     "zero alongside a positive mean is the expected shape and the mean is the "
                     "column that ranks strata. Cells are pooled across donors, so the unit of "
                     "replication is the cell and every tissue difference here is "
                     "pseudoreplicated. Annotation tier, no test and no effect size."),
        config=FIG_CFG)

    for stem, sets in (("umap_treg_arms", ARM_STRIP_SETS),
                       ("umap_treg_programs", PROGRAM_SETS)):
        write_caption(
            STAGE, f"tables/_overview/{stem}.csv",
            finding=(f"Per-cell AUCell summaries of the {len(sets)} sets drawn in "
                     f"`figures/_overview/{stem}.png` — {', '.join(sets)} — restricted to the "
                     "Treg gate, one row per tissue, so the colouring reads as numbers."),
            script=SCRIPT, fn="score_table",
            config_kv=f"rows = {len(sets)} sets x Treg x 2 tissues, metric = AUCell",
            input=SUMMARY_REL,
            how_to_read=("A restriction of the narrative scoring summary table to the Treg "
                         "gate and the sets this figure draws. One row per (`set_name` x "
                         "`tissue`) with the mean, median and standard deviation of the "
                         "per-cell AUCell score and the cell and donor counts behind it. "
                         "These are the values the full-object figure's Treg rows carry too, "
                         "because the scores are joined on barcode and not recomputed, so a "
                         "difference between the paired figures is a difference of layout. "
                         "AUCell is bounded in [0, 1] and its scale depends on set size, so "
                         "values compare across tissue within a `set_name`. Cells are pooled "
                         "across donors, so the unit of replication is the cell and every "
                         "tissue difference is pseudoreplicated. `evidence_tier` reads "
                         "`secondary_percell`."),
            config=FIG_CFG)

    write_caption(
        STAGE, "tables/_overview/umap_treg_signatures.csv",
        finding=("Per (channel x tissue) summaries of the three candidate harvest signatures "
                 "within the Treg gate, giving the numbers behind the counterpart figure's "
                 "colouring: all three channels sit higher in synovial fluid than in paired "
                 "blood."),
        script=SCRIPT, fn="channel_table",
        config_kv=(f"rows = {len(SIGNATURE_CHANNELS)} channels x Treg x 2 tissues, "
                   "metric = scanpy score_genes module score"),
        input=MARKER_SUBSTRATE_REL,
        how_to_read=(f"One row per (`set_name` x `tissue`) over the {len(drawn):,} Treg "
                     "cells, with the mean, median and standard deviation of the module score "
                     "and the cell and donor counts behind it. `metric` reads "
                     "`scanpy_score_genes_module_score`. These values are mean-centred against "
                     "a sampled background and signed, so zero means at background and a "
                     "negative mean records a position on that scale; they compare across "
                     "tissue within a `set_name`. `WT_heat_up` is the mouse anchor arm carried "
                     "as annotation only and never as a selection predicate. Cells are pooled "
                     "across donors, so the unit of replication is the cell and every tissue "
                     "difference is pseudoreplicated. Annotation tier; no test, no effect "
                     "size."),
        config=FIG_CFG)

    write_caption(
        STAGE, "tables/_overview/umap_treg_patchwork.csv",
        finding=("Both rows of the stacked Treg-only layout in one table: the four identity "
                 "genes in log-normalised expression and the six gene sets in per-cell AUCell, "
                 "one row per channel and tissue, with `metric` keeping the two scales apart."),
        script=SCRIPT, fn="channel_table",
        config_kv=(f"rows = {len(MARKER_CHANNELS)} genes + {len(ARM_STRIP_SETS)} sets, "
                   "x Treg x 2 tissues; metrics = log_normalised_expression and AUCell"),
        input=f"{SUBSTRATE_REL}, {MARKER_SUBSTRATE_REL}, {SUMMARY_REL}",
        how_to_read=("The marker table and the arm score table stacked, one row per "
                     "(`set_name` x `tissue`). Read `metric` first: "
                     "`log_normalised_expression` rows are the top row's marker panels and "
                     "`AUCell` rows the bottom row's score panels, and the two scales are "
                     "unrelated, so a comparison stays inside one metric. Cells are pooled "
                     "across donors throughout, so the unit of replication is the cell and "
                     "every tissue difference is pseudoreplicated. Annotation tier, no test "
                     "and no effect size."),
        config=FIG_CFG)

    print("[17_treg_reembedding_viz] wrote 5 overviews (umap_treg_reembedding, umap_treg_arms, "
          "umap_treg_programs, umap_treg_signatures, umap_treg_patchwork) + 5 table captions")


if __name__ == "__main__":
    main()
