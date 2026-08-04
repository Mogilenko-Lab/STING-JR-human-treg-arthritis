#!/usr/bin/env python
"""
13_arm_decomposition_viz.py: VIZ ONLY. Which curated programs contain each mouse-derived
UP arm, and how much of each arm is left over?
=========================================================================================
One deliberately plain panel: a horizontal stacked bar per mouse-derived up arm, band per
curated lens, plus the remainder no lens claims. No clustering, no ordination, no novel
method. The only arithmetic on this face is a sum of the fractional weights already
committed by 13_arm_decomposition.py.

Band widths use the FRACTIONAL accounting, which is the only accounting under which
"fraction of arm" is a true statement. A gene claimed by k lenses gives 1/k to each, so an
arm's bands total exactly 1.0 and the remainder is a share a reader can trust. The number
printed inside a band is the DUPLICATED count, how many of the arm's genes that lens
contains, and those counts sum to more than the arm. Both readings are on the face because
either one alone invites the wrong conclusion: the widths would hide how many genes a lens
actually contains, and the counts would imply a partition that does not exist.

The canvas carries the geometry and one line of key. The three readings a reader needs
before interpreting the geometry, that this is containment and carries no enrichment
statistic, that the lenses overlap so the counts over-count, and that the four arms share
structure by construction, live in `how_to_read=` and therefore in the stage README, with
every number kept. Long prose on the canvas was what made this supplement unreadable when
shrunk into a column.

Reads only committed 13_arm_decomposition tables. Run from the compartment root, AFTER
13_arm_decomposition.py:
  python 02_analysis/scripts/13_arm_decomposition_viz.py
"""
from __future__ import annotations

import os
import sys
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

COMPARTMENT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(COMPARTMENT_ROOT))
sys.path.insert(0, str(COMPARTMENT_ROOT / "02_analysis"))
os.chdir(COMPARTMENT_ROOT)

from config import PATHS  # noqa: E402
from helpers.figure_style import (  # noqa: E402
    FIG_CFG,
    purge_figures,
    round_numeric_cols,
    save_overview,
    set_paper_style,
)

STAGE = "13_arm_decomposition"
SCRIPT = "02_analysis/scripts/13_arm_decomposition_viz.py"
TIER = "secondary_annotation"
RESIDUAL = "unassigned"

ARM_ORDER = ["WT_heat_up", "KO_heat_up", "Interaction_up", "Interaction_up_fdrOnly"]

# --- declared palette --------------------------------------------------------------------
# Every hue is read from the project config (`colors.okabe_ito` for the eight categorical
# slots, `colors.diverging.up` for the ninth) so this panel shares the compartment's
# colourblind-safe palette instead of inventing one. The remainder gets a neutral grey on
# purpose: it names no program, so it must not look like one.
_OKABE = (FIG_CFG.get("colors", {}) or {}).get("okabe_ito", {}) or {}
_DIVERGING = (FIG_CFG.get("colors", {}) or {}).get("diverging", {}) or {}
PROGRAM_COL = {
    "nfkb_tnfa": _OKABE["vermillion"],
    "inflammatory": _OKABE["orange"],
    "hypoxia": _OKABE["blue"],
    "t_activation": _OKABE["bluish_green"],
    "ifn_type_i": _OKABE["sky_blue"],
    "ifn_generic_axis": _OKABE["reddish_purple"],
    "sting_specific_published": _OKABE["black"],
    "upr_er": _OKABE["yellow"],
    "hsr_curated": _DIVERGING["up"],
    RESIDUAL: "#BFBFBF",
}

# Human-readable band labels. Each names the lens by what it IS, so no band name can be
# read as a mechanism the membership count does not license.
PROGRAM_LABEL = {
    "nfkb_tnfa": "TNFA / NF-kB signalling",
    "inflammatory": "inflammatory response",
    "hypoxia": "hypoxia",
    "t_activation": "IL2-STAT5 activation",
    "ifn_type_i": "type-I interferon (Hallmark)",
    "ifn_generic_axis": "generic type-I interferon axis",
    "sting_specific_published": "published IFN-independent STING genes",
    "upr_er": "unfolded-protein response",
    "hsr_curated": "curated HSR core (Reactome/GO)",
    RESIDUAL: "no curated lens contains it",
}

# Short forms used ONLY in the right-hand hairline-lens notes, where the full label would
# wrap past the row. Each resolves unambiguously against the legend, which always carries
# the full name and the lens size.
PROGRAM_SHORT = {
    "ifn_generic_axis": "generic IFN axis",
    "ifn_type_i": "Hallmark IFN-alpha",
    "hsr_curated": "curated HSR core",
    "sting_specific_published": "published STING genes",
    "upr_er": "unfolded-protein response",
    "nfkb_tnfa": "TNFA / NF-kB",
    "inflammatory": "inflammatory response",
    "hypoxia": "hypoxia",
    "t_activation": "IL2-STAT5",
}

_F = FIG_CFG.get("figures", {}) or {}
ANNOT_SIZE = float(_F["axis_text_size"])
LEGEND_SIZE = float(_F["legend_text_size"])
# A band narrower than this carries no in-band number: the digits would collide with the
# band edges. Nothing is dropped: every count is in the sibling table and the caption names
# the small ones explicitly.
LABEL_FLOOR = 0.045
# The canvas this panel is exported on. Read from the config's two-column preset because that
# is what `save_overview(..., wide=True)` will set AFTER the figure is built. At draw time
# `fig.get_figwidth()` still returns the rcParams default, so deriving the note's wrap column from it
# silently under-wraps to the floor. Wrapping to the wrong width is not cosmetic: `save_figure`
# exports with `bbox_inches="tight"`, so a line past the right edge WIDENS the emitted canvas
# and shrinks every font relative to the page.
CANVAS_W_IN = float(_F["width_wide"])
AX_LEFT, AX_RIGHT = 0.155, 0.985


# ===========================================================================
# 1. The plotted table: a reshape of the committed membership tables
# ===========================================================================
def composition_table() -> pd.DataFrame:
    """One row per plotted band: the duplicated count and the fractional band width.

    Both columns come from tables 13_arm_decomposition.py already wrote. `n_intersect` is
    copied from `arm_program_summary.csv`, and the band width is the sum of the
    `weight_fractional` column of `arm_program_gene.csv`. No statistic is computed here.
    """
    tdir = PATHS.tables(STAGE)
    summary = pd.read_csv(tdir / "arm_program_summary.csv")
    gene = pd.read_csv(tdir / "arm_program_gene.csv")

    frac = (gene.groupby(["arm", "program"])["weight_fractional"].sum()
            .rename("weight_fractional_sum").reset_index())
    mult = pd.read_csv(tdir / "arm_program_multiplicity.csv")

    df = summary.merge(frac, on=["arm", "program"], how="left")
    df["weight_fractional_sum"] = df["weight_fractional_sum"].fillna(0.0)
    df["frac_of_arm_fractional"] = df["weight_fractional_sum"] / df["n_arm"]

    # Per-arm not-a-partition counts, repeated on every row so a caption or an in-figure
    # note can quote them without leaving this file.
    for arm in ARM_ORDER:
        m = mult[mult["arm"] == arm]["n_programs"].astype(int)
        sel = df["arm"] == arm
        df.loc[sel, "arm_n_claimed"] = int((m > 0).sum())
        df.loc[sel, "arm_n_unassigned"] = int((m == 0).sum())
        df.loc[sel, "arm_n_claims_total"] = int(m.sum())
        df.loc[sel, "arm_n_excess_claims"] = int(m.sum() - (m > 0).sum())
        df.loc[sel, "arm_max_lenses_per_gene"] = int(m.max())

    # Band order: named lenses by total fractional share across all arms (widest first),
    # remainder always last so it reads as the tail rather than as another program.
    share = (df[df["program"] != RESIDUAL].groupby("program")["weight_fractional_sum"]
             .sum().sort_values(ascending=False))
    order = list(share.index) + [RESIDUAL]
    df["band_order"] = df["program"].map(order.index)
    df["_arm"] = df["arm"].map(ARM_ORDER.index)
    df["is_partition"] = False
    df["evidence_tier"] = TIER
    df["measurement"] = "membership_not_enrichment"

    cols = ["arm", "gate", "n_arm", "program", "curated_set", "n_curated_set",
            "n_intersect", "frac_of_arm", "weight_fractional_sum",
            "frac_of_arm_fractional", "band_order", "arm_n_claimed",
            "arm_n_unassigned", "arm_n_claims_total", "arm_n_excess_claims",
            "arm_max_lenses_per_gene", "is_partition", "measurement", "evidence_tier"]
    return (df.sort_values(["_arm", "band_order"])[cols].reset_index(drop=True))


def band_order(df: pd.DataFrame) -> list[str]:
    return list(df.sort_values("band_order")["program"].drop_duplicates())


# ===========================================================================
# 2. The panel
# ===========================================================================
def sub_floor_note(sub: pd.DataFrame, bands: list[str]) -> str:
    """The lenses whose band is too narrow to carry a printed number, with their counts.

    A band under `LABEL_FLOOR` is a hairline, so its number cannot go inside it. Dropping the
    number instead would make the panel read as if those lenses contained nothing, the exact
    silent truncation the remainder is reported to avoid. So they go to the right of the bar,
    outside the 0-1 share range, where they cost no width and distort no band.
    """
    small = [(p, int(sub.loc[p, "n_intersect"])) for p in bands
             if p != RESIDUAL and p in sub.index
             and 0 < float(sub.loc[p, "frac_of_arm_fractional"]) < LABEL_FLOOR]
    if not small:
        return ""
    body = " · ".join(f"{PROGRAM_SHORT[p]} {n}" for p, n in small)
    return "\n".join(textwrap.wrap("also contains: " + body, width=46))


def plot_composition(df: pd.DataFrame):
    fig, ax = plt.subplots()
    fig.subplots_adjust(left=AX_LEFT, right=AX_RIGHT, top=0.94, bottom=0.34)
    bands = band_order(df)
    n = len(ARM_ORDER)

    for i, arm in enumerate(ARM_ORDER):
        y = n - 1 - i
        sub = df[df["arm"] == arm].set_index("program")
        left = 0.0
        for program in bands:
            if program not in sub.index:
                continue
            w = float(sub.loc[program, "frac_of_arm_fractional"])
            if w <= 0:
                continue
            ax.barh(y, w, left=left, height=0.62, color=PROGRAM_COL[program],
                    edgecolor="white", linewidth=0.8, zorder=2)
            if w >= LABEL_FLOOR:
                # The DUPLICATED gene count, not the band width. The two differ wherever a
                # lens shares genes with another, and the reader needs the count.
                txt = "white" if program in ("sting_specific_published", "hypoxia") else "black"
                ax.text(left + w / 2, y, f"{int(sub.loc[program, 'n_intersect'])}",
                        va="center", ha="center", fontsize=ANNOT_SIZE, color=txt, zorder=3)
            left += w
        note = sub_floor_note(sub, bands)
        if note:
            ax.text(1.03, y, note, va="center", ha="left", fontsize=ANNOT_SIZE)

    labels = []
    for arm in ARM_ORDER:
        r = df[df["arm"] == arm].iloc[0]
        labels.append(f"{arm}\n{int(r['n_arm'])} genes · gate {r['gate']}")
    ax.set_yticks(range(n))
    ax.set_yticklabels(list(reversed(labels)))
    ax.set_ylim(-0.55, n - 0.45)
    # Bands occupy 0-1; the space beyond it is the right-hand home for the hairline lenses,
    # so no band is widened to make room for its own label.
    ax.set_xlim(0, 1.75)
    ax.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xlabel("Share of the arm's genes (a gene in k lenses gives 1/k to each)")
    ax.set_title("Curated-lens membership of the mouse 39 °C-derived up arms")
    ax.spines[["top", "right"]].set_visible(False)

    handles = [Patch(facecolor=PROGRAM_COL[p], edgecolor="white",
                     label=f"{PROGRAM_LABEL[p]}"
                           + ("" if p == RESIDUAL else
                              f"  ({int(df[df['program'] == p]['n_curated_set'].iloc[0])} genes)"))
               for p in bands]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(AX_LEFT, 0.28),
               ncol=2, frameon=False, fontsize=LEGEND_SIZE)

    # ONE line of key, and only what a reader needs to interpret the geometry. The readings
    # that used to sit here as three paragraphs of capitals now live in `how_to_read=`, and so
    # in the stage README, with every number kept.
    #
    # Wrapped to the canvas, not to the sentence. `save_figure` exports with
    # `bbox_inches="tight"`, so a line running past the right edge silently WIDENS the emitted
    # canvas, the figure stops being the declared geometry, and every font shrinks relative to
    # the page. The wrap column is derived from the axes' own width so it tracks the config
    # geometry rather than a number guessed once.
    wrap_col = max(60, int(CANVAS_W_IN * (AX_RIGHT - AX_LEFT) / (ANNOT_SIZE * 0.0088)))
    key = ("Band widths use 1/k weights and total 1.0; a band's number is how many of the arm's "
           "genes that lens contains.")
    fig.text(AX_LEFT, 0.055, "\n".join(textwrap.wrap(key, width=wrap_col)),
             ha="left", va="top", fontsize=ANNOT_SIZE)
    return fig


# ===========================================================================
def main() -> None:
    set_paper_style(config=FIG_CFG)
    purge_figures(STAGE, "arm_program", overview=True, config=FIG_CFG)

    df = composition_table()
    fig = plot_composition(df)

    wt = df[df["arm"] == "WT_heat_up"].set_index("program")
    ko = df[df["arm"] == "KO_heat_up"].set_index("program")
    inter = df[df["arm"] == "Interaction_up"].set_index("program")
    fdr_only = df[df["arm"] == "Interaction_up_fdrOnly"].set_index("program")

    save_overview(
        fig, STAGE, "arm_program_composition",
        table=round_numeric_cols(df),
        finding=(
            "Nine curated anchor-independent lenses contain "
            f"{int(wt.loc[RESIDUAL, 'arm_n_claimed'])} of the 199 WT_heat_up genes and "
            f"{int(ko.loc[RESIDUAL, 'arm_n_claimed'])} of the 218 KO_heat_up genes, leaving "
            f"remainders of {int(wt.loc[RESIDUAL, 'n_intersect'])} and "
            f"{int(ko.loc[RESIDUAL, 'n_intersect'])} genes as the largest single part of each "
            "large arm, while what the lenses do claim is dominated by inflammatory gene "
            f"content ({int(wt.loc['nfkb_tnfa', 'n_intersect'])} TNFA/NF-kB and "
            f"{int(wt.loc['inflammatory', 'n_intersect'])} inflammatory-response genes in "
            f"WT_heat_up) against {int(wt.loc['hsr_curated', 'n_intersect'])} in the curated HSR "
            f"core and {int(wt.loc['sting_specific_published', 'n_intersect'])} of the 21 "
            "published IFN-independent STING genes. The thin Interaction arms invert that shape: "
            f"Hallmark type-I interferon contains {int(inter.loc['ifn_type_i', 'n_intersect'])} "
            f"of the 7 genes at the fdr_logfc gate and "
            f"{int(fdr_only.loc['ifn_type_i', 'n_intersect'])} of the 18 at the relaxed fdr_only "
            "gate."),
        script=SCRIPT, fn="plot_composition",
        config_kv=(f"colors.okabe_ito + colors.diverging.up; label_floor={LABEL_FLOOR}; "
                   f"accounting=fractional (1/n_programs_for_gene); evidence_tier={TIER}"),
        input=("03_results/13_arm_decomposition/tables/arm_program_summary.csv, "
               "03_results/13_arm_decomposition/tables/arm_program_gene.csv, "
               "03_results/13_arm_decomposition/tables/arm_program_multiplicity.csv"),
        how_to_read=(
            "Each band counts how many of one arm's genes a frozen curated lens contains. That "
            "is set arithmetic over committed gene lists, so no NES, FDR, direction or effect "
            "size appears here or anywhere in these tables. One row per arm, named by how it was "
            "derived and labelled with its gene count and the mouse anchor's gate (fdr_logfc, "
            "with fdr_only the relaxed Interaction variant frozen as Interaction_fdrOnly_up.txt). "
            "Band width is the fractional share: a gene in k lenses gives 1/k to each, so widths "
            "total 1.0. The number in a band is the duplicated count of that arm's genes the lens "
            "contains, so numbers and widths disagree. A band too narrow for a digit carries its "
            "count to the right of the bar, so every lens keeps its number. Grey is the "
            "remainder, the genes no lens contains, left unnamed on purpose. "
            "The lenses overlap, so the bands are a containment tally: "
            f"{int(wt.loc[RESIDUAL, 'arm_n_claimed'])} of the "
            f"{int(wt.loc[RESIDUAL, 'n_arm'])} WT_heat_up genes are contained by at least one "
            f"lens, and those {int(wt.loc[RESIDUAL, 'arm_n_claimed'])} carry "
            f"{int(wt.loc[RESIDUAL, 'arm_n_claims_total'])} memberships, with up to "
            f"{int(wt.loc[RESIDUAL, 'arm_max_lenses_per_gene'])} lenses on one gene, so the "
            "printed counts exceed the claimed genes by "
            f"{int(wt.loc[RESIDUAL, 'arm_n_excess_claims'])}. Per-gene multiplicity is in "
            "arm_program_multiplicity.csv. The four rows share structure by construction: the "
            "mouse contrasts are linearly dependent as model coefficients (WT_heat = KO_heat + "
            "Interaction), and the two Interaction rows are one contrast at two gates, so "
            "agreement between rows is expected. That algebra holds for the coefficients and "
            "stops at the thresholded lists: WT_heat_up and KO_heat_up share 182 genes, "
            "Interaction_up shares none with either, and Interaction_up_fdrOnly holds all 7 "
            "Interaction_up genes among its 18. Annotation tier, firewalled from the "
            "donor-pseudobulk claim spine; no effect-size row."),
        config=FIG_CFG, wide=True, height=7.4,
    )
    plt.close(fig)

    print(f"[{STAGE}_viz] wrote 1 overview (arm_program_composition) over "
          f"{len(band_order(df))} bands x {len(ARM_ORDER)} arms")
    print(df[["arm", "program", "n_intersect", "frac_of_arm",
              "frac_of_arm_fractional"]].to_string(index=False))


if __name__ == "__main__":
    main()
