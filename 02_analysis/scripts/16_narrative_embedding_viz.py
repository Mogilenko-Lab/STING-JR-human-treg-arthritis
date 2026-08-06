#!/usr/bin/env python
"""
16_narrative_embedding_viz.py: VIZ ONLY (no statistics).
=============================================================================
Six-panel UMAP strips of the frozen 99,915-cell sorted JIA T-cell annotation (GSE160097),
drawn from `03_results/interactive/16_narrative_embedding.parquet`. One sampled frame per
strip, so a cell holds its place across every panel:

  umap_full_reference   tissue, the frozen FACS sort gate, three Treg identity genes, donor
  umap_full_arms        the three mouse-derived up arms, ruled off from three lenses
  umap_full_programs    six curated lenses: interferon/STING beside inflammation/activation
  umap_full_patchwork   the reference strip over the arm strip, for layout

GEOMETRY comes from `02_analysis/helpers/umap_grid.py`, shared with
`17_treg_reembedding_viz.py`, so strips from either stage stack at identical panel size.

THIS MAP CARRIES NO BATCH CORRECTION. The coordinates are `01_qc_filter.py`'s recipe carried
verbatim through `08_harvest_readout.parquet`. Every footer states the donor structure that
leaves, read from `17_treg_reembedding`'s mixing table.

TIER. An embedding places cells against their labels. Claims rest on donor-level pseudobulk
differential expression within the frozen cell states.

COLUMNS. Score panels colour on `*_AUCell`, sized from the manifest's
`n_genes_found_in_object`. Marker panels are log-normalised expression joined on barcode from
the `07_embedding` substrate. Two scored sets stay undrawn for want of a slot:
`Interaction_fdrOnly_up`, `HALLMARK_UNFOLDED_PROTEIN_RESPONSE`.

ONE BAR PER SCORE STRIP. Each panel keeps its own clip and is rescaled across it onto [0, 1],
so the strip carries one gradient and an assembled figure needs one key. The clip is unchanged,
so no panel's picture changes; what the bar stops carrying is LEVEL, which moves to the
caption's config line and the source table. The AUCell ranges run 0.002-0.032 to 0.073-0.225, so a bar in real
units would flatten half the row.

ONE BAR FOR THE MARKERS TOO, but in REAL units: FOXP3, IL2RA and CTLA4 are one quantity with
comparable ranges, so pooling their clip keeps level readable between the genes. They stay off
the score bar — expression and AUCell do not share a key.

Run in-container from the compartment root, after 16_narrative_scoring.py and
17_treg_reembedding.py:
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

COMPARTMENT_ROOT = Path(__file__).resolve().parents[2]

sys.path.insert(0, str(COMPARTMENT_ROOT))
sys.path.insert(0, str(COMPARTMENT_ROOT / "02_analysis"))
os.chdir(COMPARTMENT_ROOT)

from config import (  # noqa: E402
    ARM_STRIP_SETS,
    PATHS,
    POPULATION_COLORS,
    PROGRAM_STRIP_SETS,
    TISSUE_COLORS,
)
from helpers.figure_style import (  # noqa: E402
    FIG_CFG, save_overview, set_paper_style, write_caption,
)
from helpers import umap_grid as ug  # noqa: E402

STAGE = "16_narrative_scoring"
SCRIPT = "02_analysis/scripts/16_narrative_embedding_viz.py"
SUBSTRATE = "16_narrative_embedding.parquet"
SUBSTRATE_REL = f"03_results/interactive/{SUBSTRATE}"
SUMMARY_REL = "03_results/16_narrative_scoring/tables/narrative_score_summary.csv"
MANIFEST_REL = "03_results/16_narrative_scoring/tables/narrative_scoring_manifest.csv"
MARKER_SUBSTRATE_REL = "03_results/07_embedding/tables/hook_factor_substrate.parquet"
MIXING_REL = "03_results/17_treg_reembedding/tables/treg_reembedding_mixing.csv"

N_COL = 6
SAMPLE_N = 60_000       # cells drawn per strip; sampled ONCE and reused by every panel
SAMPLE_SEED = 0

_OI = (FIG_CFG.get("colors", {}) or {}).get("okabe_ito", {}) or {}
TISSUE_COL = TISSUE_COLORS
STATE_COL = POPULATION_COLORS
# Seven donors, seven Okabe-Ito hues. Yellow stays and orange sits out: orange and
# vermillion are the one confusable pair at scatter-point size.
DONOR_HUES = ["vermillion", "sky_blue", "bluish_green", "yellow",
              "blue", "reddish_purple", "black"]

TISSUE_LABEL = {"synovial_fluid": "synovial fluid", "peripheral_blood": "paired blood"}

TIER_LINE = ("Annotation tier. Claims in this compartment rest on donor-level pseudobulk "
             "differential expression within these frozen cell states.")

# --- the reference strip ------------------------------------------------------
# Three canonical Treg identity genes beside the sort gate they corroborate, so a reader
# checks the gate against expression.
MARKER_PANELS = [
    ("FOXP3", "FOXP3\nTreg lineage TF"),
    ("IL2RA", "IL2RA\nCD25"),
    ("CTLA4", "CTLA4\nsuppressive effector"),
]
MARKER_CHANNELS = [c for c, _ in MARKER_PANELS]

REFERENCE_BANDS = [
    ("frozen annotation", 0, 1, True),
    ("core Treg identity genes", 2, 4, True),
    ("donor of origin", 5, 5, False),
]

# --- the arm strip: anchor-dependent columns ruled off from anchor-independent ones ---
# Panel list and order from `analysis_config.yaml::percell_map_panels`, which the Treg-only
# counterpart, the violin figure and the sweep-coverage audit read too.
# WT_heat_up and KO_heat_up share 182 genes and go on ONE colour scale. Every other panel
# spans its own range and keeps its own bar.
ARM_SETS = ARM_STRIP_SETS[:3]
ARM_LENS_SETS = ARM_STRIP_SETS[3:]
ARM_SHARED_SCALE = ["WT_heat_up_AUCell", "KO_heat_up_AUCell"]
# The third band names the accession, not "project-derived": a reader can check GSE161426.
ARM_BANDS = [
    ("mouse 39 °C-derived up arms — anchor-dependent", 0, 2, True),
    ("curated, anchor-independent", 3, 4, True),
    ("derived here from GSE161426", 5, 5, False),
]

# --- the curated-lens strip, from the same one declaration --------------------
PROGRAM_SETS = list(PROGRAM_STRIP_SETS)
PROGRAM_BANDS = [
    ("cGAS-STING and type-I interferon", 0, 2, True),
    ("inflammation and activation", 3, 5, False),
]


# =============================================================================
# Inputs
# =============================================================================
def load_marker_columns(frame: pd.DataFrame) -> pd.DataFrame:
    """Add the 07_embedding marker columns to a frame, matched on barcode.

    The values come from that file so these panels and its own marker figure colour by the
    same numbers. An incomplete join halts: a reindexed column would put one cell's
    expression at another cell's coordinates.
    """
    src = pd.read_parquet(COMPARTMENT_ROOT / MARKER_SUBSTRATE_REL)
    cols = src[MARKER_CHANNELS].copy()
    cols.index = src.index.astype(str)
    out = frame.join(cols, on="barcode")
    missing = int(out[MARKER_CHANNELS].isna().any(axis=1).sum())
    if missing:
        raise AssertionError(
            f"{missing} of {len(out)} cells carry no {MARKER_SUBSTRATE_REL} expression; the "
            "two substrates disagree on barcodes, so the marker panels cannot be drawn.")
    return out


def set_sizes(manifest: pd.DataFrame) -> dict:
    """{set_name: genes found in the object} — the size each score was really computed over."""
    return dict(zip(manifest["set_name"], manifest["n_genes_found_in_object"].astype(int)))


def score_panels(sets: list, sizes: dict) -> list:
    """(column, title) per scored set, titles of uniform height across the strip."""
    return list(zip((f"{s}_AUCell" for s in sets), ug.set_titles(sets, sizes)))


def full_object_mixing(mix: pd.DataFrame, key: str) -> tuple[float, float]:
    """Observed and chance same-`key` neighbour fraction on THIS map, from stage 17's table.

    Stage 17 is the one place donor structure on this map is measured, and a figure stating
    that the map is uncorrected carries the number with it.
    """
    hit = mix[(mix["embedding"] == "full_object_restricted") & (mix["grouping_key"] == key)
              & (mix["group"] == "_all_")]
    if hit.empty:
        raise ValueError(
            f"no full_object_restricted / {key} / _all_ row in {MIXING_REL}; run "
            "02_analysis/scripts/17_treg_reembedding.py first, which is where donor "
            "neighbourhood composition on this map is measured.")
    return float(hit["observed_same_frac"].iloc[0]), float(hit["expected_same_frac"].iloc[0])


def footer(mix: pd.DataFrame) -> list:
    """The two standing lines every strip in this stage carries."""
    obs, exp = full_object_mixing(mix, "donor")
    return [("Coordinates: UMAP of the frozen annotation, no batch correction — same-donor "
             f"neighbours {obs:.3f} at k = 30 against {exp:.3f} expected from the donor "
             "proportions."),
            TIER_LINE]


# =============================================================================
# Source tables. Aggregates only: a per-cell frame is not a readable source CSV.
# =============================================================================
def stratum_table(full: pd.DataFrame, drawn: pd.DataFrame, channels: list) -> pd.DataFrame:
    """Per (cell state x tissue x donor) counts, plus the mean of every channel drawn.

    One row per stratum, one column per panel channel, so every panel has its numbers beside
    the counts.
    """
    keys = ["coarse_label", "tissue", "donor"]
    tbl = full.groupby(keys, observed=True).size().rename("n_cells").reset_index()
    got = drawn.groupby(keys, observed=True).size().rename("n_cells_drawn").reset_index()
    tbl = tbl.merge(got, on=keys, how="left")
    tbl["n_cells_drawn"] = tbl["n_cells_drawn"].fillna(0).astype(int)
    tbl["frac_of_total"] = tbl["n_cells"] / len(full)
    means = (full.groupby(keys, observed=True)[channels].mean()
             .rename(columns={c: f"mean_{c}" for c in channels}).reset_index())
    tbl = tbl.merge(means, on=keys, how="left")
    return tbl.sort_values(keys).reset_index(drop=True)


def score_table(summary: pd.DataFrame, sets: list) -> pd.DataFrame:
    """The committed per (set x cell state x tissue) summary, restricted to the panels drawn."""
    tbl = summary[summary["set_name"].isin(sets)].copy()
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


def _gate_mean(full: pd.DataFrame, gene: str, label: str) -> float:
    return float(full.loc[full["coarse_label"] == label, gene].mean())


# =============================================================================
# Strips
# =============================================================================
def _new_strip(n_row: int = 1):
    fig = plt.figure(figsize=ug.canvas(N_COL, n_row))
    return fig, ug.grid(fig, N_COL, n_row)


def draw_reference_row(axes: list, drawn: pd.DataFrame, xlim, ylim, marker_clip: tuple):
    """Tissue, sort gate, three Treg identity genes, donor — over one shared frame.

    The markers are one quantity at comparable ranges (0-3.16, 0-2.75, 0-2.35), so they pool
    onto one clip and one bar in REAL units and brightness compares between the genes. The
    score rows are rescaled instead because their ranges differ tenfold.
    """
    donors = sorted(drawn["donor"].unique())
    donor_col = {d: _OI[h] for d, h in zip(donors, DONOR_HUES)}
    donor_label = {d: d.replace("JIA_patient_", "JIA ") for d in donors}

    ug.scatter_cat(axes[0], drawn, "tissue", TISSUE_COL, "Tissue of origin\nsampling site",
                   xlim, ylim, FIG_CFG, order=["synovial_fluid", "peripheral_blood"],
                   labels=TISSUE_LABEL, ncol=2)
    ug.scatter_cat(axes[1], drawn, "coarse_label", STATE_COL,
                   "Frozen cell state\nFACS sort gate", xlim, ylim, FIG_CFG,
                   order=["Treg", "Tcon", "CD8"], ncol=3)
    sc = None
    for ax, (col, title) in zip(axes[2:5], MARKER_PANELS):
        sc = ug.scatter_cont(ax, drawn, col, title, xlim, ylim, FIG_CFG,
                             vlim=marker_clip, unit=ug.UNIT_EXPR, colorbar=False)
    ug.row_colorbar(axes[4].figure, axes[4], sc, ug.UNIT_EXPR, FIG_CFG)
    ug.scatter_cat(axes[5], drawn, "donor", donor_col, "Donor\n7 JIA patients", xlim, ylim,
                   FIG_CFG, order=donors, labels=donor_label, ncol=4)


def draw_score_row(axes: list, drawn: pd.DataFrame, xlim, ylim, panels: list, clips: dict):
    """One row of continuous AUCell panels on ONE shared bar.

    Each panel keeps the clip it was handed and is rescaled across it onto [0, 1]. The ranges
    span 0.002-0.032 to 0.073-0.225, so a bar in real units would flatten the narrow panels;
    rescaling keeps every panel's structure and moves level to the caption and source table.
    """
    sc = None
    for ax, (col, title) in zip(axes, panels):
        sc = ug.scatter_cont(ax, drawn, col, title, xlim, ylim, FIG_CFG,
                             vlim=clips[col], unit=ug.UNIT_RESCALED,
                             rescale=True, colorbar=False)
    ug.row_colorbar(axes[-1].figure, axes[-1], sc, ug.UNIT_RESCALED, FIG_CFG)


def clip_note(panels: list, clips: dict) -> str:
    """The limits each panel was rescaled over. The bar reads 0 to 1, so without these the
    figure states positions within unstated ranges."""
    return "; ".join(f"{col} {clips[col][0]:.4f}-{clips[col][1]:.4f}" for col, _ in panels)


# =============================================================================
def main() -> None:
    set_paper_style(config=FIG_CFG)

    full = pd.read_parquet(PATHS.interactive / SUBSTRATE)
    summary = pd.read_csv(PATHS.tables(STAGE) / "narrative_score_summary.csv")
    manifest = pd.read_csv(PATHS.tables(STAGE) / "narrative_scoring_manifest.csv")
    mix = pd.read_csv(COMPARTMENT_ROOT / MIXING_REL)
    sizes = set_sizes(manifest)
    full = load_marker_columns(full)
    print(f"[16_narrative_embedding_viz] substrate {full.shape[0]:,} cells x {full.shape[1]} cols "
          f"(incl. {len(MARKER_CHANNELS)} joined marker columns), "
          f"{full['donor'].nunique()} donors, {full['coarse_label'].nunique()} frozen cell states")

    drawn = ug.sample_frame(full, SAMPLE_N, SAMPLE_SEED)
    xlim, ylim = ug.data_box(drawn)
    print(f"[16_narrative_embedding_viz] one shared frame of {len(drawn):,} cells, "
          f"box x {xlim[0]:.2f}..{xlim[1]:.2f}")

    foot = footer(mix)
    donor_obs, donor_exp = full_object_mixing(mix, "donor")

    # Colour limits: the arm pair pooled onto one scale, every other panel on its own.
    arm_shared = ug.clip_of(drawn, ARM_SHARED_SCALE)
    arm_panels = score_panels(ARM_STRIP_SETS, sizes)
    arm_clips = {col: (arm_shared if col in ARM_SHARED_SCALE else ug.clip_of(drawn, [col]))
                 for col, _ in arm_panels}
    program_panels = score_panels(PROGRAM_SETS, sizes)
    program_clips = {col: ug.clip_of(drawn, [col]) for col, _ in program_panels}
    marker_clip = ug.clip_of(drawn, MARKER_CHANNELS)   # one clip over all three, one bar
    print(f"[16_narrative_embedding_viz] WT_heat_up + KO_heat_up shared colour scale "
          f"{arm_shared[0]:.4f} to {arm_shared[1]:.4f}")

    fw, fh = ug.canvas(N_COL)
    cfg_common = (f"figures.dpi = 300, figures.rasterized_dpi = 600, sample_n = {SAMPLE_N}, "
                  f"sample_seed = {SAMPLE_SEED}, panel = {ug.PANEL_W_IN:.2f} x "
                  f"{ug.PANEL_H_IN:.2f} in x {N_COL} "
                  f"columns, point_size = {ug.point_size(len(drawn))}")

    # --- 1. reference strip --------------------------------------------------
    fig, rows = _new_strip()
    draw_reference_row(rows[0], drawn, xlim, ylim, marker_clip)
    ug.group_bands(fig, rows[0], REFERENCE_BANDS, FIG_CFG, colour=_OI["black"])
    ug.dress(fig, "Sorted JIA T cells: tissue, sort gate, Treg identity and donor",
             f"UMAP of the frozen annotation. {len(drawn):,} of {len(full):,} cells drawn, one "
             "sampled frame shared by all six panels.", foot, FIG_CFG, colour=_OI["black"])
    save_overview(
        fig, STAGE, "umap_full_reference",
        table=stratum_table(full, drawn, MARKER_CHANNELS),
        finding=("One sampled frame of the frozen 99,915-cell sorted JIA T-cell map coloured six "
                 "ways, establishing the layout every score colouring of this substrate is read "
                 "against: all seven JIA donors contribute cells to both the synovial-fluid and "
                 "the paired-blood side, the sort gates occupy largely distinct territory within "
                 "each tissue, and the three Treg identity genes corroborate the Treg gate from "
                 f"expression — FOXP3 mean {_gate_mean(full, 'FOXP3', 'Treg'):.3f} in the Treg "
                 f"gate against {_gate_mean(full, 'FOXP3', 'Tcon'):.3f} in Tcon and "
                 f"{_gate_mean(full, 'FOXP3', 'CD8'):.3f} in CD8. The map carries no batch "
                 f"correction, which leaves {donor_obs:.3f} same-donor neighbours at k = 30 "
                 f"against {donor_exp:.3f} expected."),
        script=SCRIPT, fn="draw_reference_row",
        config_kv=(f"{cfg_common}, colours = colors.okabe_ito (tissue vermillion/blue, state "
                   "green/orange/pink, donor 7 hues), markers = "
                   f"{', '.join(MARKER_CHANNELS)} on {ug.CMAP}, "
                   f"clip_percentiles = {list(ug.CLIP)}, marker_scale = one pooled clip at "
                   f"{marker_clip[0]:.4f}-{marker_clip[1]:.4f} over all "
                   f"{len(MARKER_CHANNELS)} genes, one bar"),
        input=f"{SUBSTRATE_REL}, {MARKER_SUBSTRATE_REL}, {MIXING_REL}",
        how_to_read=(
            f"Six panels over one sampled frame at identical coordinates, so a cell sits in "
            f"the same place in all six. Panels 1 and 2 are the annotation this compartment "
            f"is built on: tissue of origin, synovial fluid in vermillion and paired blood "
            f"in blue, then the frozen FACS sort gate. Panels 3 to 5 are log-normalised "
            f"expression of FOXP3, the lineage transcription factor, IL2RA, the CD25 chain, "
            f"and CTLA4, the suppressive effector, joined on barcode from the 07_embedding "
            f"substrate, so the gate is checked against expression; most cells read zero, "
            f"which is ordinary for single-cell counts, so the three pool onto ONE clip at "
            f"the 2nd and 98th percentile, {marker_clip[0]:.2f} to {marker_clip[1]:.2f}, "
            f"with the highest-expressing cells drawn last. One bar serves all three and it "
            f"is in real units, so brightness compares between the genes as well as within "
            f"one. Panel 6 is "
            f"donor. THESE COORDINATES CARRY NO BATCH CORRECTION, so the donor panel shows "
            f"real donor structure: {donor_obs:.3f} same-donor neighbours at k = 30 against "
            f"{donor_exp:.3f} expected. Donor is crossed with tissue by design, every "
            f"patient contributing to both arms, so that structure sits inside each tissue. "
            f"The source table gives the per-stratum counts and marker means. Annotation "
            f"tier: claims rest on donor-level pseudobulk differential expression within "
            f"these frozen cell states."),
        width=fw, height=fh, config=FIG_CFG)
    plt.close(fig)

    # --- 2. mouse-derived up arms beside the lenses they are read against ----
    fig, rows = _new_strip()
    draw_score_row(rows[0], drawn, xlim, ylim, arm_panels, arm_clips)
    ug.group_bands(fig, rows[0], ARM_BANDS, FIG_CFG, colour=_OI["black"])
    ug.dress(fig, "Mouse 39 °C-derived up arms, and the lenses they are read against",
             "Colour is per-cell AUCell rescaled to each panel's own 2nd-98th percentile, so "
             "all six read on one 0-to-1 bar. WT_heat_up and KO_heat_up are rescaled together "
             f"({arm_shared[0]:.3f} to {arm_shared[1]:.3f}); AUCell values are in the source "
             "table.",
             foot, FIG_CFG, colour=_OI["black"])
    save_overview(
        fig, STAGE, "umap_full_arms",
        table=score_table(summary, ARM_STRIP_SETS),
        finding=("All three mouse 39 °C-derived up arms colour the synovial-fluid side of the map "
                 "brighter than the paired-blood side in every frozen sort gate, and so do all "
                 "three anchor-independent lenses ruled off beside them, which is the point of "
                 "putting them in one row: the per-cell AUCell means put the WT_heat_up shift at "
                 f"{_shift(summary, 'WT_heat_up', 'Treg')} in Treg, "
                 f"{_shift(summary, 'WT_heat_up', 'Tcon')} in Tcon and "
                 f"{_shift(summary, 'WT_heat_up', 'CD8')} in CD8, while the curated hypoxia lens "
                 f"runs {_shift(summary, 'HALLMARK_HYPOXIA', 'Treg')} and the curated "
                 f"proteostasis core {_shift(summary, 'HSR_core', 'Treg')} over the same Treg "
                 "cells, so the anchor arm's tissue colouring is not distinctive of the anchor."),
        script=SCRIPT, fn="draw_score_row",
        config_kv=(f"{cfg_common}, cmap = {ug.CMAP}, clip_percentiles = {list(ug.CLIP)}, "
                   "columns = " + ", ".join(c for c, _ in arm_panels)
                   + f", shared_scale = WT_heat_up + KO_heat_up at "
                     f"{arm_shared[0]:.4f}-{arm_shared[1]:.4f}"
                   + ", colour = rescaled to panel clip onto [0, 1], one bar for the row; "
                     "AUCell limits rescaled over: " + clip_note(arm_panels, arm_clips)),
        input=f"{SUBSTRATE_REL}, {SUMMARY_REL}, {MANIFEST_REL}, {MIXING_REL}",
        how_to_read=(
            f"Six panels over the reference strip's frame and bounding box. THE VERTICAL "
            f"RULES CARRY THE PROVENANCE. The left three are anchor-dependent: WT_heat_up is "
            f"the mouse WT iTreg 39-versus-37 °C up arm in human projection, KO_heat_up the "
            f"same contrast in cGAS-knockout iTregs, Interaction_up the genotype-by- "
            f"temperature arm at 7 genes, small enough that one gene moves the score, so "
            f"read it for location and treat its spread as noise. The middle two are curated "
            f"and versioned: HALLMARK_HYPOXIA and the activation-free HSR_core proteostasis "
            f"lens. The right one is eTreg_up, this compartment's own GSE161426 effector- "
            f"Treg contrast over 4 synovial against 14 blood donors, derived for exploration "
            f"and ruled off for that reason. ONE BAR SERVES THE ROW, and it reads 0 to 1, "
            f"not AUCell: each panel is clipped to its own 2nd and 98th percentile as "
            f"before and then rescaled across that clip, so the picture is unchanged and "
            f"the panels merely stop needing six separate keys. Brightness therefore "
            f"compares tissue WITHIN a panel and says nothing about level BETWEEN panels — "
            f"the AUCell limits each panel was rescaled over are in the config line below "
            f"and the values themselves in the source table. WT_heat_up and KO_heat_up "
            f"share 182 genes and one clip, {arm_shared[0]:.4f} to {arm_shared[1]:.4f}, so "
            f"those two do compare pixel for pixel. Titles give the symbols scored. "
            f"Cells pool across donors, so a tissue difference here is pseudoreplicated. "
            f"Temperature and hypoxia are both imposed by the inflamed joint and stay "
            f"entangled in cross-sectional human data, so the hypoxia panel carries no HIF "
            f"claim, and similarity across a rule is a reason to test. Claims rest on donor- "
            f"level pseudobulk differential expression within the frozen cell states."),
        width=fw, height=fh, config=FIG_CFG)
    plt.close(fig)

    # --- 3. curated lenses: interferon/STING against inflammation/activation --
    fig, rows = _new_strip()
    draw_score_row(rows[0], drawn, xlim, ylim, program_panels, program_clips)
    ug.group_bands(fig, rows[0], PROGRAM_BANDS, FIG_CFG, colour=_OI["black"])
    ug.dress(fig, "Curated program lenses on the JIA T-cell map",
             "Colour is per-cell AUCell rescaled to each panel's own 2nd-98th percentile, so all "
             "six read on one 0-to-1 bar; AUCell values are in the source table. Every set here "
             "is curated, versioned and independent of the mouse anchor.",
             foot, FIG_CFG, colour=_OI["black"])
    save_overview(
        fig, STAGE, "umap_full_programs",
        table=score_table(summary, PROGRAM_SETS),
        finding=("Six curated lenses on one map separate two readings that a single interferon "
                 "panel would merge: the 21 published IFN-independent STING-activation genes sit "
                 "far lower in the Treg gate than in Tcon or CD8 in both tissues (per-cell AUCell "
                 f"mean {_mean_of(summary, 'sting_specific_published', 'Treg', 'peripheral_blood'):.4f} "
                 "in Treg blood against "
                 f"{_mean_of(summary, 'sting_specific_published', 'Tcon', 'peripheral_blood'):.4f} "
                 "in Tcon and "
                 f"{_mean_of(summary, 'sting_specific_published', 'CD8', 'peripheral_blood'):.4f} "
                 "in CD8), so that panel reports a sort-gate difference alongside a tissue one, "
                 "while the generic type-I interferon axis and the three inflammation and "
                 "activation programs brighten the synovial-fluid side across all three gates "
                 f"({_shift(summary, 'ifn_generic_axis', 'Treg')} and "
                 f"{_shift(summary, 'HALLMARK_TNFA_SIGNALING_VIA_NFKB', 'Treg')} over Treg "
                 "cells)."),
        script=SCRIPT, fn="draw_score_row",
        config_kv=(f"{cfg_common}, cmap = {ug.CMAP}, clip_percentiles = {list(ug.CLIP)}, "
                   "columns = " + ", ".join(c for c, _ in program_panels)
                   + "; per-panel limits, no pooling across sets"
                   + ", colour = rescaled to panel clip onto [0, 1], one bar for the row; "
                     "AUCell limits rescaled over: " + clip_note(program_panels, program_clips)),
        input=f"{SUBSTRATE_REL}, {SUMMARY_REL}, {MANIFEST_REL}, {MIXING_REL}",
        how_to_read=(
            f"Six panels on the arm strip's frame and colormap; titles give the symbols "
            f"scored. The rule splits two families. Left of it: sting_specific_published, "
            f"the 21 published IFN-independent STING-activation genes; ifn_generic_axis, a "
            f"200-gene generic type-I interferon axis carrying the thinnest intersection in "
            f"the strip; and HALLMARK_INTERFERON_ALPHA_RESPONSE. Right of it: "
            f"HALLMARK_TNFA_SIGNALING_VIA_NFKB, HALLMARK_INFLAMMATORY_RESPONSE and "
            f"HALLMARK_IL2_STAT5_SIGNALING, the programs the first family has to be "
            f"distinguished from. A synovial-high colouring shared by both families is "
            f"generic inflammation; only a pattern the left family carries and the right one "
            f"lacks would be specific to STING or interferon. ONE BAR SERVES THE ROW, and it "
            f"reads 0 to 1, not AUCell: each panel is clipped per set to the 2nd and 98th "
            f"percentile with the highest drawn last, exactly as before, and then rescaled "
            f"across that clip. Brightness therefore compares tissue WITHIN a panel and says "
            f"nothing about level BETWEEN panels; the limits each panel was rescaled over are "
            f"in the config line below and the values in the source table. The six sets span "
            f"0.001-0.068 to 0.073-0.225 in AUCell, which is why a bar in real units would "
            f"flatten half of them. The published STING set is "
            f"21 genes and its IFN-β validation in the positive-control compartment is "
            f"underpowered at three donors, so a bright or dim panel there is consistent "
            f"with STING pathway activity and never proof of it. Cells pool across donors, "
            f"so tissue differences are pseudoreplicated. Claims rest on donor-level "
            f"pseudobulk differential expression within the frozen cell states."),
        width=fw, height=fh, config=FIG_CFG)
    plt.close(fig)

    # --- 4. the two strips on one canvas, for layout -------------------------
    pw_h = ug.canvas(N_COL, 2)[1]
    fig, rows = _new_strip(n_row=2)
    draw_reference_row(rows[0], drawn, xlim, ylim, marker_clip)
    ug.group_bands(fig, rows[0], REFERENCE_BANDS, FIG_CFG, colour=_OI["black"])
    draw_score_row(rows[1], drawn, xlim, ylim, arm_panels, arm_clips)
    ug.group_bands(fig, rows[1], ARM_BANDS, FIG_CFG, colour=_OI["black"])
    ug.dress(fig, "Sorted JIA T cells: the reference layout, and the signatures projected onto it",
             f"One sampled frame of {len(drawn):,} of {len(full):,} cells, shared by all twelve "
             "panels. Top row annotation and Treg identity, the three genes pooled onto one "
             "expression bar; bottom row per-cell AUCell rescaled per panel onto one 0-to-1 bar.",
             foot, FIG_CFG, colour=_OI["black"])
    save_overview(
        fig, STAGE, "umap_full_patchwork",
        table=stratum_table(full, drawn, MARKER_CHANNELS + [f"{s}_AUCell" for s in ARM_STRIP_SETS]),
        finding=("The reference layout and the signature colouring of the same 60,000-cell frame "
                 "on one canvas, so the sort gate, the Treg identity genes and the mouse-derived "
                 "arms are read against each other without turning a page: the arms brighten the "
                 "synovial-fluid side of every gate, and the gate that FOXP3 marks is not the "
                 "gate where they brighten most."),
        script=SCRIPT, fn="draw_reference_row",
        config_kv=(f"{cfg_common}, rows = 2 x {N_COL}, canvas = {fw:.1f} x {pw_h:.1f} in, "
                   "top row = " + ", ".join(["tissue", "coarse_label"] + MARKER_CHANNELS
                                            + ["donor"])
                   + "; bottom row = " + ", ".join(c for c, _ in arm_panels)
                   + ", rescaled to panel clip onto [0, 1] on one bar; AUCell limits rescaled "
                     "over: " + clip_note(arm_panels, arm_clips)),
        input=f"{SUBSTRATE_REL}, {SUMMARY_REL}, {MANIFEST_REL}, {MARKER_SUBSTRATE_REL}, "
              f"{MIXING_REL}",
        how_to_read=(
            f"The two strips this stage ships separately, `umap_full_reference` above "
            f"`umap_full_arms`, on one canvas at identical panel size so a column reads top "
            f"to bottom. Nothing new is drawn, and both rows hold the identical frame of "
            f"cells at identical coordinates. The rows share cells and coordinates, and "
            f"their units differ: the top row is categorical annotation plus log-normalised "
            f"expression, the bottom row per-cell AUCell of a gene set rescaled per panel "
            f"onto a single 0-to-1 bar. Each row's own "
            f"caption carries its full reading, the bottom row's rules carry each "
            f"signature's provenance, and the source table gives every panel's channel as a "
            f"column over the same strata. The coordinates carry no batch correction, "
            f"leaving {donor_obs:.3f} same-donor neighbours at k = 30 against "
            f"{donor_exp:.3f} expected. Annotation tier throughout; claims rest on donor- "
            f"level pseudobulk differential expression within the frozen cell states."),
        width=fw, height=pw_h, config=FIG_CFG)
    plt.close(fig)

    # --- 5. captions for the same-stem source tables -------------------------
    write_caption(
        STAGE, "tables/_overview/umap_full_reference.csv",
        finding=("The populated (cell state x tissue x donor) strata behind the reference map, "
                 "with each stratum's mean FOXP3, IL2RA and CTLA4: all seven donors appear in "
                 "both tissues, and the three Treg identity genes are highest in the Treg strata, "
                 "so the sort gate and the expression agree."),
        script=SCRIPT, fn="stratum_table",
        config_kv=f"sample_n = {SAMPLE_N}, sample_seed = {SAMPLE_SEED}",
        input=f"{SUBSTRATE_REL}, {MARKER_SUBSTRATE_REL}",
        how_to_read=("One row per (`coarse_label` x `tissue` x `donor`) stratum present in the "
                     "frozen annotation. `n_cells` is the full-object count and `n_cells_drawn` "
                     "the count in the sampled frame the figure draws, so the two say how "
                     "faithful the drawn frame is to the object. `frac_of_total` is `n_cells` "
                     f"over {len(full):,}. Each `mean_<gene>` column is that stratum's mean "
                     "log-normalised expression, the numbers behind the three marker panels; "
                     "they are means over zero-inflated per-cell values, so they rank strata and "
                     "do not estimate a per-cell level. Absent combinations carry no row. "
                     "Annotation tier, no test and no effect size."),
        config=FIG_CFG)

    for stem, sets in (("umap_full_arms", ARM_STRIP_SETS),
                       ("umap_full_programs", PROGRAM_SETS)):
        write_caption(
            STAGE, f"tables/_overview/{stem}.csv",
            finding=(f"Per-cell AUCell summaries of the {len(sets)} sets drawn in "
                     f"`figures/_overview/{stem}.png` — {', '.join(sets)} — one row per cell "
                     "state and tissue, so the colouring reads as numbers."),
            script=SCRIPT, fn="score_table",
            config_kv=f"rows = {len(sets)} sets x 3 frozen cell states x 2 tissues, metric = AUCell",
            input=SUMMARY_REL,
            how_to_read=("A restriction of the stage summary table to the sets this figure draws. "
                         "One row per (`set_name` x `coarse_label` x `tissue`) with the mean, "
                         "median and standard deviation of the per-cell AUCell score and the cell "
                         "and donor counts behind it. AUCell is bounded in [0, 1] and its scale "
                         "depends on set size, so values compare across tissue within a "
                         "`set_name`; a cross-set comparison takes the gene lists themselves. "
                         "Cells are pooled across donors, so the unit of replication is the cell "
                         "and every tissue difference is pseudoreplicated. `evidence_tier` reads "
                         "`secondary_percell` throughout."),
            config=FIG_CFG)

    write_caption(
        STAGE, "tables/_overview/umap_full_patchwork.csv",
        finding=("Every channel the stacked layout draws, as one row per (cell state x tissue x "
                 "donor) stratum: the cell counts, the three Treg identity gene means and the six "
                 "per-cell AUCell means, so all twelve panels are readable as numbers from one "
                 "table."),
        script=SCRIPT, fn="stratum_table",
        config_kv=(f"sample_n = {SAMPLE_N}, sample_seed = {SAMPLE_SEED}, channels = "
                   f"{len(MARKER_CHANNELS)} marker genes + {len(ARM_STRIP_SETS)} AUCell sets"),
        input=f"{SUBSTRATE_REL}, {MARKER_SUBSTRATE_REL}",
        how_to_read=("One row per (`coarse_label` x `tissue` x `donor`) stratum, with the "
                     "full-object and drawn cell counts and then one `mean_<channel>` column per "
                     "panel of the stacked figure. The `mean_<gene>` columns are log-normalised "
                     "expression and the `mean_<set>_AUCell` columns are rank-based scores in "
                     "[0, 1] whose scale depends on set size, so a value compares across strata "
                     "within its own column and never across columns. Means over cells within a "
                     "stratum: the unit is the cell, so nothing here is a donor-level effect. "
                     "Annotation tier, no test and no effect size."),
        config=FIG_CFG)

    print("[16_narrative_embedding_viz] wrote 4 overviews (umap_full_reference, umap_full_arms, "
          "umap_full_programs, umap_full_patchwork) + 4 table captions")


if __name__ == "__main__":
    main()
