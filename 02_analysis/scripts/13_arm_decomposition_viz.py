#!/usr/bin/env python
"""
13_arm_decomposition_viz.py: VIZ ONLY. Which curated programs contain each mouse-derived
UP arm, and how much of each arm is left over?
=========================================================================================
Two panels of set arithmetic over committed gene lists.

`arm_program_composition` — a horizontal stacked bar per mouse-derived up arm, band per
curated lens, plus the remainder no lens claims. No clustering, no ordination, no novel
method. The only arithmetic on this face is a sum of the fractional weights already
committed by 13_arm_decomposition.py.

`arm_hypoxia_euler` — the WT_heat_up-versus-HALLMARK_HYPOXIA overlap drawn instead of
quoted, area-proportional, in TWO vocabularies side by side. The composition bar reports
that overlap as a single number against the frozen 200-gene lens; that number is not the
one a GSEA of hypoxia in the Treg contrast is computed on, because a ranked list carries
only the genes that survived detection and `filterByExpr`. Drawing one vocabulary alone
would contradict the other and read as an error, so both are drawn on one shared
area-per-gene scale and the loss between them becomes the visible quantity.

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

import math
import os
import sys
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import Circle, Patch  # noqa: E402
from scipy.optimize import brentq  # noqa: E402

COMPARTMENT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(COMPARTMENT_ROOT))
sys.path.insert(0, str(COMPARTMENT_ROOT / "02_analysis"))
os.chdir(COMPARTMENT_ROOT)

from config import CONFIG, PATHS  # noqa: E402
from helpers.figure_style import (  # noqa: E402
    FIG_CFG,
    purge_figures,
    round_numeric_cols,
    save_overview,
    set_paper_style,
    write_caption,
)
from helpers.geneset_utils import load_alias_map, resolve_symbols  # noqa: E402

STAGE = "13_arm_decomposition"
SCRIPT = "02_analysis/scripts/13_arm_decomposition_viz.py"
TIER = "secondary_annotation"
RESIDUAL = "unassigned"

ARM_ORDER = ["WT_heat_up", "KO_heat_up", "Interaction_up", "Interaction_up_fdrOnly"]

# --- the one pair the Euler draws, and the ranking it is restricted to -------------------
# The Treg gate is the compartment's most directly Treg-relevant ranking and the one whose
# hypoxia enrichment is quoted downstream, so it is the vocabulary the second panel uses.
EULER_ARM = "WT_heat_up"
EULER_LENS = "hypoxia"
EULER_LENS_SET = "HALLMARK_HYPOXIA"
EULER_POPULATION = "treg"

# --- declared palette --------------------------------------------------------------------
# Every hue is read from the project config (`colors.okabe_ito` for the eight categorical
# slots, `colors.diverging.up` for the ninth) so this panel shares the compartment's
# colourblind-safe palette instead of inventing one. The remainder gets a neutral grey on
# purpose: it names no program, so it must not look like one.
_OKABE = (FIG_CFG.get("colors", {}) or {}).get("okabe_ito", {}) or {}
_DIVERGING = (FIG_CFG.get("colors", {}) or {}).get("diverging", {}) or {}
# The neutral grey config reserves for a reference annotation drawn over data, so a leader
# line can never be mistaken for a data series.
_REFERENCE_LINE = (FIG_CFG.get("colors", {}) or {})["reference_line"]
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


def arm_gene_sets() -> dict:
    """The gene set of each arm, read from `arm_program_gene.csv`.

    The overlaps between arms are stated in the caption; reading them here keeps a printed
    overlap from outliving the frozen lists the arms were thresholded from.
    """
    gene = pd.read_csv(PATHS.tables(STAGE) / "arm_program_gene.csv")
    return {str(a): set(d["gene"]) for a, d in gene.groupby("arm")}


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
# 3. The Euler: one arm against one lens, in two vocabularies
# ===========================================================================
# Set membership is NOT recomputed here. The arm's genes and its intersection with the
# hypoxia lens are read from the tables 13_arm_decomposition.py committed, the frozen lens
# file is hash-verified against that stage's own pin before it is read, and the read is
# asserted to reproduce the committed intersection gene for gene. What this section adds is
# a RESTRICTION of that membership to a declared vocabulary, which is a different question
# from who is a member.
EULER_COL_ARM = _OKABE["orange"]
EULER_COL_LENS = PROGRAM_COL["hypoxia"]
EULER_FILL_ALPHA = 0.45
# Room above the taller circle for the shared-count leader, in gene-area units.
EULER_HEAD_ROOM = 2.6
EULER_MARGIN = 0.7
# The band below the panels, wide enough that the legend and the key never overlap. Both
# are placed by hand for that reason: a collision there is invisible in a text report.
EULER_BOTTOM = 0.29
EULER_LEGEND_Y = 0.255
EULER_KEY_Y = 0.145


def frozen_lens_path(set_id: str) -> Path:
    """The frozen lens file this compartment declares for `set_id`, read from config.

    Selected out of the declared `project_frozen` file list rather than spelled out, so a
    config that stops declaring the set stops this figure instead of silently reading a
    path that happens to still exist on disk.
    """
    files = (((CONFIG.get("unbiased_enrichment", {}) or {}).get("project_frozen", {}) or {})
             .get("files", []) or [])
    hits = [f for f in files if Path(f).stem == set_id]
    if len(hits) != 1:
        raise ValueError(
            f"analysis_config.yaml::unbiased_enrichment.project_frozen.files declares "
            f"{len(hits)} file(s) with stem {set_id}; expected exactly one.")
    return COMPARTMENT_ROOT / hits[0]


def read_gene_list(path: Path) -> list[str]:
    """Newline-delimited HGNC symbols, order preserved, blank lines dropped."""
    return [ln.strip() for ln in Path(path).read_text().splitlines() if ln.strip()]


def _symbol_column(rel_path: str) -> list[str]:
    """The `gene_symbol` column of one committed vocabulary table."""
    return pd.read_csv(COMPARTMENT_ROOT / rel_path)["gene_symbol"].astype(str).tolist()


def vocabulary_layers() -> dict:
    """The three nested vocabulary layers plus the alias map, all from config.

    `ranked` is the layer the Euler's second panel restricts to; `matrix` and `reference`
    exist so a symbol absent from `ranked` can be attributed to the expression filter, to
    never having been detected, or to the reference build, rather than all three at once.
    """
    sa = CONFIG.get("symbol_alias") or {}
    for key in ("map_path", "matrix_vocabulary", "reference_feature_symbols", "ranked_list"):
        if not sa.get(key):
            raise ValueError(f"analysis_config.yaml::symbol_alias has no `{key}`")
    ranked_rel = str(sa["ranked_list"]).replace("{population}", EULER_POPULATION)
    ranked = pd.read_csv(COMPARTMENT_ROOT / ranked_rel, sep="\t", header=None,
                         names=["symbol", "stat"])
    return {
        "alias_map": load_alias_map(COMPARTMENT_ROOT / sa["map_path"]),
        "matrix": set(_symbol_column(sa["matrix_vocabulary"])),
        "reference": set(_symbol_column(sa["reference_feature_symbols"])),
        "ranked": set(ranked["symbol"].dropna().astype(str)),
        "ranked_rel": ranked_rel,
        "map_rel": str(sa["map_path"]),
    }


def restrict_to_ranked(genes: list[str], voc: dict) -> dict:
    """Resolve `genes` into the matrix vintage, then keep what the ranked list carries.

    Returns the surviving symbols plus the three-way ledger the compartment requires:
    matched, matched-via-alias, and genuinely absent split by WHY. Alias resolution only
    ever adds, so `n_exact` cannot move and the recovery is reported rather than absorbed.
    """
    nominal = list(dict.fromkeys(genes))
    resolved, applied = resolve_symbols(nominal, voc["alias_map"], voc["matrix"])
    kept = [g for g in resolved if g in voc["ranked"]]
    exact = [g for g in nominal if g in voc["ranked"]]
    via_alias = [(ref, tgt) for ref, tgt in applied
                 if tgt in voc["ranked"] and ref not in voc["ranked"]]
    lost = [g for g in nominal if g not in voc["ranked"]
            and g not in {ref for ref, _ in via_alias}]
    return {
        "kept": set(kept),
        "n_nominal": len(nominal),
        "n_exact": len(set(exact)),
        "n_via_alias": len(via_alias),
        "alias_pairs": [f"{ref}->{tgt}" for ref, tgt in via_alias],
        # A symbol the matrix carries but the ranking does not was dropped by filterByExpr;
        # one the reference build carries but the matrix does not was never detected in
        # sorted T cells; one the reference build lacks is a vocabulary miss outright.
        "n_expression_filtered": len([g for g in lost if g in voc["matrix"]]),
        "n_undetected": len([g for g in lost
                             if g not in voc["matrix"] and g in voc["reference"]]),
        "n_absent_from_reference": len([g for g in lost
                                        if g not in voc["matrix"]
                                        and g not in voc["reference"]]),
    }


def lens_area(r_a: float, r_b: float, d: float) -> float:
    """Area shared by two circles of radii `r_a`, `r_b` whose centres are `d` apart."""
    if d >= r_a + r_b:
        return 0.0
    if d <= abs(r_a - r_b):
        return math.pi * min(r_a, r_b) ** 2
    seg_a = r_a ** 2 * math.acos((d * d + r_a * r_a - r_b * r_b) / (2 * d * r_a))
    seg_b = r_b ** 2 * math.acos((d * d + r_b * r_b - r_a * r_a) / (2 * d * r_b))
    kite = 0.5 * math.sqrt((r_a + r_b - d) * (d + r_a - r_b)
                           * (d - r_a + r_b) * (d + r_a + r_b))
    return seg_a + seg_b - kite


def solve_euler(n_a: int, n_b: int, n_both: int) -> dict:
    """Circle geometry whose three areas EQUAL the three region counts, in gene units.

    A two-set Euler is exactly solvable for every valid configuration, which is why this
    figure can be area-proportional rather than merely suggestive: the shared area falls
    continuously and monotonically from min(n_a, n_b) when one circle sits inside the other
    to zero when they part, so one distance hits any admissible overlap exactly. The
    configurations with no exact solution begin at three sets. Areas are in gene units, so
    the same scale serves every panel and the residual is a gene count.
    """
    if n_both > min(n_a, n_b):
        raise ValueError(f"overlap {n_both} exceeds the smaller set ({min(n_a, n_b)})")
    r_a, r_b = math.sqrt(n_a / math.pi), math.sqrt(n_b / math.pi)
    lo, hi = abs(r_a - r_b), r_a + r_b
    if n_both == min(n_a, n_b):
        d = lo
    elif n_both == 0:
        d = hi
    else:
        d = brentq(lambda x: lens_area(r_a, r_b, x) - n_both, lo, hi,
                   xtol=1e-12, rtol=8.9e-16)
    solved = lens_area(r_a, r_b, d)
    y_cross = 0.0 if d >= r_a + r_b else math.sqrt(
        max(0.0, r_a ** 2 - ((d * d + r_a * r_a - r_b * r_b) / (2 * d)) ** 2))
    return {"r_a": r_a, "r_b": r_b, "d": d, "solved_area": solved,
            "residual": solved - n_both, "y_cross": y_cross,
            "cx_a": -d / 2.0, "cx_b": d / 2.0,
            "x_lens": (r_a - r_b) / 2.0,
            "x_a_only": -(r_a + r_b) / 2.0, "x_b_only": (r_a + r_b) / 2.0}


def euler_universes() -> list[dict]:
    """The two vocabularies the Euler draws, each with its regions solved and its ledger.

    Panel one is the frozen sets exactly as curated, which is the universe the composition
    bar and the committed membership tables report. Panel two restricts both sets to the
    Treg donor-pseudobulk ranked list, which is the universe an enrichment statistic on
    that ranking is actually computed over.
    """
    tdir = PATHS.tables(STAGE)
    summary = pd.read_csv(tdir / "arm_program_summary.csv")
    gene = pd.read_csv(tdir / "arm_program_gene.csv")

    arm = sorted(set(gene[gene["arm"] == EULER_ARM]["gene"].astype(str)))
    row = summary[(summary["arm"] == EULER_ARM) & (summary["program"] == EULER_LENS)]
    if len(row) != 1:
        raise ValueError(f"arm_program_summary.csv has no single {EULER_ARM}/{EULER_LENS} row")
    row = row.iloc[0]
    committed_both = sorted(str(row["genes"]).split(";")) if str(row["genes"]) else []

    lens_path = frozen_lens_path(EULER_LENS_SET)
    lens = read_gene_list(lens_path)
    if len(lens) != int(row["n_curated_set"]):
        raise AssertionError(
            f"{EULER_LENS_SET} holds {len(lens)} genes, the committed membership table "
            f"declares {int(row['n_curated_set'])}")
    if len(arm) != int(row["n_arm"]):
        raise AssertionError(
            f"{EULER_ARM} holds {len(arm)} genes in arm_program_gene.csv, "
            f"arm_program_summary.csv declares {int(row['n_arm'])}")
    both = sorted(set(arm) & set(lens))
    if both != committed_both:
        raise AssertionError(
            "the frozen lens read here does not reproduce the committed intersection:\n"
            f"  this read: {both}\n  arm_program_summary.csv: {committed_both}")

    voc = vocabulary_layers()
    arm_r = restrict_to_ranked(arm, voc)
    lens_r = restrict_to_ranked(lens, voc)
    print(f"[{STAGE}_viz] {EULER_LENS_SET} in the {EULER_POPULATION} ranking: "
          f"{lens_r['n_exact']} exact + {lens_r['n_via_alias']} via alias "
          f"= {len(lens_r['kept'])} of {lens_r['n_nominal']} "
          f"({', '.join(lens_r['alias_pairs']) or 'no alias recovery'})")
    print(f"[{STAGE}_viz] {EULER_ARM} in the {EULER_POPULATION} ranking: "
          f"{arm_r['n_exact']} exact + {arm_r['n_via_alias']} via alias "
          f"= {len(arm_r['kept'])} of {arm_r['n_nominal']}")

    return [
        {"universe": "frozen_sets",
         "universe_label": "Frozen sets as curated",
         "universe_note": "every gene of both committed lists",
         "vocabulary": "(none — the curated lists in full)",
         "arm": set(arm), "lens": set(lens), "ledger_arm": None, "ledger_lens": None},
        {"universe": "treg_ranked_list",
         "universe_label": f"Restricted to the {EULER_POPULATION.capitalize()} ranked list",
         "universe_note": f"only genes {Path(voc['ranked_rel']).name} carries "
                          f"({len(voc['ranked']):,} symbols)",
         "vocabulary": voc["ranked_rel"],
         "arm": arm_r["kept"], "lens": lens_r["kept"],
         "ledger_arm": arm_r, "ledger_lens": lens_r},
    ]


def euler_table(universes: list[dict]) -> pd.DataFrame:
    """One row per (universe, region): the count drawn, and the geometry that drew it."""
    rows = []
    for u in universes:
        a, b = u["arm"], u["lens"]
        geo = solve_euler(len(a), len(b), len(a & b))
        u["geometry"] = geo
        regions = [
            ("arm_only", f"in {EULER_ARM} only", sorted(a - b)),
            ("shared", f"in both {EULER_ARM} and {EULER_LENS_SET}", sorted(a & b)),
            ("lens_only", f"in {EULER_LENS_SET} only", sorted(b - a)),
        ]
        led_a, led_l = u["ledger_arm"], u["ledger_lens"]
        for region, label, genes in regions:
            rows.append({
                "universe": u["universe"],
                "universe_label": u["universe_label"],
                "vocabulary": u["vocabulary"],
                "region": region,
                "region_label": label,
                "n_genes": len(genes),
                "n_arm": len(a),
                "n_lens": len(b),
                "frac_of_arm": len(genes) / len(a),
                "arm_n_nominal": led_a["n_nominal"] if led_a else len(a),
                "arm_n_exact_match": led_a["n_exact"] if led_a else len(a),
                "arm_n_via_alias": led_a["n_via_alias"] if led_a else 0,
                "arm_n_expression_filtered": led_a["n_expression_filtered"] if led_a else 0,
                "arm_n_undetected": led_a["n_undetected"] if led_a else 0,
                "arm_n_absent_from_reference": (led_a["n_absent_from_reference"]
                                                if led_a else 0),
                "lens_n_nominal": led_l["n_nominal"] if led_l else len(b),
                "lens_n_exact_match": led_l["n_exact"] if led_l else len(b),
                "lens_n_via_alias": led_l["n_via_alias"] if led_l else 0,
                "lens_n_expression_filtered": (led_l["n_expression_filtered"]
                                               if led_l else 0),
                "lens_n_undetected": led_l["n_undetected"] if led_l else 0,
                "lens_n_absent_from_reference": (led_l["n_absent_from_reference"]
                                                 if led_l else 0),
                "lens_alias_pairs_applied": ";".join(led_l["alias_pairs"]) if led_l else "",
                "arm_alias_pairs_applied": ";".join(led_a["alias_pairs"]) if led_a else "",
                "circle_radius_arm": geo["r_a"],
                "circle_radius_lens": geo["r_b"],
                "centre_distance": geo["d"],
                "shared_area_solved": geo["solved_area"],
                "shared_area_residual_genes": geo["residual"],
                "is_area_proportional": True,
                "measurement": "membership_not_enrichment",
                "evidence_tier": TIER,
                "genes": ";".join(genes),
            })
    return pd.DataFrame(rows)


def plot_arm_hypoxia_euler(universes: list[dict]):
    """Two exactly area-proportional two-circle Eulers on ONE shared area-per-gene scale."""
    geos = [u["geometry"] for u in universes]
    half_x = max(g["d"] / 2 + max(g["r_a"], g["r_b"]) for g in geos) + EULER_MARGIN
    r_max = max(max(g["r_a"], g["r_b"]) for g in geos)
    ylim = (-(r_max + EULER_MARGIN), r_max + EULER_HEAD_ROOM)

    fig, axes = plt.subplots(1, len(universes))
    fig.subplots_adjust(left=0.015, right=0.985, top=0.82, bottom=EULER_BOTTOM, wspace=0.04)

    for ax, u in zip(axes, universes):
        g = u["geometry"]
        a, b = u["arm"], u["lens"]
        for cx, r, colour in ((g["cx_a"], g["r_a"], EULER_COL_ARM),
                              (g["cx_b"], g["r_b"], EULER_COL_LENS)):
            ax.add_patch(Circle((cx, 0), r, facecolor=colour, alpha=EULER_FILL_ALPHA,
                                edgecolor=colour, linewidth=1.6, zorder=2))
        # Counts sit in the region they measure. The shared count goes inside the lens and
        # its name above the diagram, because the lens is too narrow to hold both.
        ax.text(g["x_a_only"], 0.9, f"{len(a - b)}", ha="center", va="center",
                fontsize=ANNOT_SIZE + 3, zorder=4)
        ax.text(g["x_a_only"], -0.9, f"{EULER_ARM} only", ha="center", va="center",
                fontsize=ANNOT_SIZE, zorder=4)
        ax.text(g["x_b_only"], 0.9, f"{len(b - a)}", ha="center", va="center",
                fontsize=ANNOT_SIZE + 3, zorder=4)
        ax.text(g["x_b_only"], -0.9, f"{EULER_LENS_SET}\nonly", ha="center", va="center",
                fontsize=ANNOT_SIZE, zorder=4)
        ax.text(g["x_lens"], 0.0, f"{len(a & b)}", ha="center", va="center",
                fontsize=ANNOT_SIZE + 3, zorder=4)
        head = r_max + EULER_HEAD_ROOM - 1.1
        ax.plot([g["x_lens"], g["x_lens"]], [g["y_cross"] + 0.35, head - 0.5],
                color=_REFERENCE_LINE, linewidth=0.9, zorder=1)
        ax.text(g["x_lens"], head, "in both", ha="center", va="bottom",
                fontsize=ANNOT_SIZE, color=_REFERENCE_LINE, zorder=4)

        # Two lines only. The vocabulary's fuller description would run into the
        # neighbouring panel's title at this canvas width, so it lives in the caption.
        ax.set_title(f"{u['universe_label']}\n{EULER_ARM} {len(a)} · "
                     f"{EULER_LENS_SET} {len(b)}", fontsize=ANNOT_SIZE + 1)
        ax.set_xlim(-half_x, half_x)
        ax.set_ylim(*ylim)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)

    handles = [Patch(facecolor=EULER_COL_ARM, alpha=EULER_FILL_ALPHA,
                     edgecolor=EULER_COL_ARM,
                     label=f"{EULER_ARM} — the mouse WT iTreg 39-versus-37 °C up arm"),
               Patch(facecolor=EULER_COL_LENS, alpha=EULER_FILL_ALPHA,
                     edgecolor=EULER_COL_LENS,
                     label=f"{EULER_LENS_SET} — frozen MSigDB Hallmark, "
                           "anchor-independent")]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, EULER_LEGEND_Y),
               ncol=2, frameon=False, fontsize=LEGEND_SIZE, columnspacing=2.4)
    fig.suptitle("The mouse 39 °C-derived up arm against curated hypoxia, "
                 "drawn to scale in two vocabularies", fontsize=float(_F["title_size"]))

    led = universes[-1]["ledger_lens"]
    led_a = universes[-1]["ledger_arm"]
    # Both recovery counts are read from the ledgers rather than asserted, so a change in
    # either frozen list cannot leave this sentence stating a number the figure disproves.
    arm_recovery = (
        "and the arm recovers none" if not led_a["n_via_alias"] else
        f"and the arm recovers {led_a['n_via_alias']} of its own "
        f"({', '.join(p.replace('->', ' as ') for p in led_a['alias_pairs'])})"
    )
    key = ("Every area is exactly proportional to its gene count and both panels share one "
           "area-per-gene scale, so the right panel is smaller because fewer genes are "
           "testable in this contrast. Right is the "
           f"vocabulary an enrichment statistic on the {EULER_POPULATION.capitalize()} "
           f"ranking is computed over; {led['n_via_alias']} of its "
           f"{len(universes[-1]['lens'])} hypoxia genes match only after alias resolution "
           f"({', '.join(p.replace('->', ' as ') for p in led['alias_pairs'])}), "
           f"{arm_recovery}. This is membership: no NES, FDR or effect size here.")
    wrap_col = max(60, int(CANVAS_W_IN * 0.97 / (ANNOT_SIZE * 0.0088)))
    fig.text(0.015, EULER_KEY_Y, "\n".join(textwrap.wrap(key, width=wrap_col)),
             ha="left", va="top", fontsize=ANNOT_SIZE)
    return fig


# ===========================================================================
def main() -> None:
    set_paper_style(config=FIG_CFG)
    purge_figures(STAGE, "arm_program", overview=True, config=FIG_CFG)
    purge_figures(STAGE, "arm_hypoxia", overview=True, config=FIG_CFG)

    df = composition_table()
    gsets = arm_gene_sets()
    fig = plot_composition(df)

    wt = df[df["arm"] == "WT_heat_up"].set_index("program")
    ko = df[df["arm"] == "KO_heat_up"].set_index("program")
    inter = df[df["arm"] == "Interaction_up"].set_index("program")
    fdr_only = df[df["arm"] == "Interaction_up_fdrOnly"].set_index("program")

    save_overview(
        fig, STAGE, "arm_program_composition",
        table=round_numeric_cols(df),
        finding=(
            f"{len(band_order(df)) - 1} curated anchor-independent lenses contain "
            f"{int(wt.loc[RESIDUAL, 'arm_n_claimed'])} of the "
            f"{int(wt.loc[RESIDUAL, 'n_arm'])} WT_heat_up genes and "
            f"{int(ko.loc[RESIDUAL, 'arm_n_claimed'])} of the "
            f"{int(ko.loc[RESIDUAL, 'n_arm'])} KO_heat_up genes, leaving remainders of "
            f"{int(wt.loc[RESIDUAL, 'n_intersect'])} and "
            f"{int(ko.loc[RESIDUAL, 'n_intersect'])} genes as the largest single part of each "
            "large arm. Inflammatory gene content dominates what the lenses do claim: "
            f"{int(wt.loc['nfkb_tnfa', 'n_intersect'])} TNFA/NF-kB and "
            f"{int(wt.loc['inflammatory', 'n_intersect'])} inflammatory-response genes in "
            f"WT_heat_up, against {int(wt.loc['hsr_curated', 'n_intersect'])} in the curated HSR "
            f"core and {int(wt.loc['sting_specific_published', 'n_intersect'])} of the "
            f"{int(wt.loc['sting_specific_published', 'n_curated_set'])} published "
            "IFN-independent STING genes. The thin Interaction arms invert that shape: "
            f"Hallmark type-I interferon contains {int(inter.loc['ifn_type_i', 'n_intersect'])} "
            f"of the {int(inter.loc[RESIDUAL, 'n_arm'])} genes at the fdr_logfc gate and "
            f"{int(fdr_only.loc['ifn_type_i', 'n_intersect'])} of the "
            f"{int(fdr_only.loc[RESIDUAL, 'n_arm'])} at the relaxed fdr_only gate."),
        script=SCRIPT, fn="plot_composition",
        config_kv=(f"colors.okabe_ito + colors.diverging.up; label_floor={LABEL_FLOOR}; "
                   f"accounting=fractional (1/n_programs_for_gene); evidence_tier={TIER}"),
        input=("03_results/13_arm_decomposition/tables/arm_program_summary.csv, "
               "03_results/13_arm_decomposition/tables/arm_program_gene.csv, "
               "03_results/13_arm_decomposition/tables/arm_program_multiplicity.csv"),
        how_to_read=(
            "Each band counts how many of one arm's genes a frozen curated lens contains. That "
            "is set arithmetic over committed gene lists, so no NES, FDR, direction or effect "
            "size appears here or in these tables. One row per arm, named by how it was "
            "derived and labelled with its gene count and the mouse anchor's gate (fdr_logfc, "
            "with fdr_only the relaxed Interaction variant frozen as Interaction_fdrOnly_up.txt). "
            "Band width is the fractional share: a gene in k lenses gives 1/k to each, so widths "
            "total 1.0. The number in a band is the duplicated count of that arm's genes the lens "
            "contains, so numbers and widths measure different things. A band too narrow for a "
            "digit carries its count to the right of the bar, so every lens keeps its number. "
            "Grey is the remainder, the genes no lens contains, left unnamed on purpose. "
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
            f"stops at the thresholded lists: WT_heat_up and KO_heat_up share "
            f"{len(gsets['WT_heat_up'] & gsets['KO_heat_up'])} genes, Interaction_up shares "
            f"{len(gsets['Interaction_up'] & (gsets['WT_heat_up'] | gsets['KO_heat_up']))} with "
            f"either, and Interaction_up_fdrOnly holds all "
            f"{len(gsets['Interaction_up'] & gsets['Interaction_up_fdrOnly'])} Interaction_up "
            f"genes among its {len(gsets['Interaction_up_fdrOnly'])}. Annotation tier, "
            f"firewalled from the donor-pseudobulk claim spine; no effect-size row."),
        config=FIG_CFG, wide=True, height=7.4,
    )
    plt.close(fig)

    print(f"[{STAGE}_viz] wrote 1 overview (arm_program_composition) over "
          f"{len(band_order(df))} bands x {len(ARM_ORDER)} arms")
    print(df[["arm", "program", "n_intersect", "frac_of_arm",
              "frac_of_arm_fractional"]].to_string(index=False))

    # --- the arm-versus-hypoxia Euler, in both vocabularies ------------------
    universes = euler_universes()
    euler = euler_table(universes)
    fig = plot_arm_hypoxia_euler(universes)

    def _n(universe: str, region: str) -> int:
        return int(euler[(euler["universe"] == universe)
                         & (euler["region"] == region)]["n_genes"].iloc[0])

    nominal, ranked = "frozen_sets", "treg_ranked_list"
    n_arm_nom = _n(nominal, "arm_only") + _n(nominal, "shared")
    n_arm_rank = _n(ranked, "arm_only") + _n(ranked, "shared")
    n_lens_rank = _n(ranked, "lens_only") + _n(ranked, "shared")
    led = universes[-1]["ledger_lens"]
    led_a = universes[-1]["ledger_arm"]
    arm_recovery_long = (
        "The arm recovers none of them." if not led_a["n_via_alias"] else
        f"The arm recovers {led_a['n_via_alias']} of its own the same way "
        f"({', '.join(led_a['alias_pairs'])}), which is what lifts it from "
        f"{led_a['n_exact']} to {n_arm_rank}."
    )
    lost_shared = sorted(set(str(euler[(euler["universe"] == nominal)
                                       & (euler["region"] == "shared")]["genes"].iloc[0]).split(";"))
                         - set(str(euler[(euler["universe"] == ranked)
                                         & (euler["region"] == "shared")]["genes"].iloc[0]).split(";")))
    max_resid = float(euler["shared_area_residual_genes"].abs().max())

    save_overview(
        fig, STAGE, "arm_hypoxia_euler",
        table=round_numeric_cols(euler),
        finding=(
            f"Curated hypoxia accounts for a small minority of the mouse 39 °C-derived up arm "
            f"in either vocabulary, and the vocabulary decides how small: the frozen lists "
            f"share {_n(nominal, 'shared')} genes of the arm's {n_arm_nom} against "
            f"{_n(nominal, 'lens_only')} hypoxia genes the arm does not carry, while "
            f"restricting both to the {EULER_POPULATION.capitalize()} donor-pseudobulk ranked "
            f"list leaves only {n_arm_rank} of the arm's {n_arm_nom} genes testable and drops "
            f"the shared count to {_n(ranked, 'shared')} — so the arm's hypoxia content is "
            f"{_n(nominal, 'shared') / n_arm_nom:.0%} of the curated arm but "
            f"{_n(ranked, 'shared') / n_arm_rank:.0%} of the part of it this contrast can "
            f"actually test, and {led['n_via_alias']} of the {n_lens_rank} testable hypoxia "
            "genes are visible only because alias resolution recovered them."),
        script=SCRIPT, fn="plot_arm_hypoxia_euler",
        config_kv=(f"arm = {EULER_ARM}; lens = {EULER_LENS_SET} "
                   f"(gene_sets.project_frozen); vocabulary = symbol_alias.ranked_list at "
                   f"population={EULER_POPULATION}; alias pairs from symbol_alias.map_path, "
                   f"accepted only; colours = colors.okabe_ito.orange + the hypoxia band hue "
                   f"of this stage's program palette; fill_alpha={EULER_FILL_ALPHA}; "
                   f"evidence_tier={TIER}"),
        input=("03_results/13_arm_decomposition/tables/arm_program_summary.csv, "
               "03_results/13_arm_decomposition/tables/arm_program_gene.csv, "
               "00_data/references/msigdb_hallmark/HALLMARK_HYPOXIA.txt, "
               "00_data/references/symbol_alias/symbol_alias_map.csv, "
               "03_results/03_pseudobulk/tables/ranked_treg.tsv, "
               "03_results/03_pseudobulk/tables/gene_symbols.csv, "
               "03_results/00_build/tables/reference_feature_symbols.csv"),
        how_to_read=(
            "The same two gene lists, read in two vocabularies. Every area equals its gene "
            "count, solved numerically; the largest residual across "
            f"both panels is {max_resid:.1e} genes. A two-set Euler is exactly solvable for "
            "every valid configuration, because the shared area falls continuously from the "
            "smaller set's size to zero as the circles part, so every area here is exact; the "
            "configurations with no exact solution begin at three sets. "
            "Both panels share one area-per-gene scale and one bounding box, so the right "
            "panel is smaller because it holds fewer genes. Orange is the mouse WT iTreg "
            "39-versus-37 °C up arm in human projection; blue is frozen MSigDB Hallmark "
            "hypoxia, curated independently of the anchor, so the overlap is an independent "
            "measurement. Each region carries its count, the "
            "shared one inside the lens with its name above on a grey leader because the lens "
            f"is too narrow for both. Left is the frozen lists in full, {n_arm_nom} and "
            f"{_n(nominal, 'lens_only') + _n(nominal, 'shared')} genes, the universe the "
            "composition bar and the committed membership tables report. Right keeps only what "
            f"the {EULER_POPULATION.capitalize()} donor-pseudobulk ranked list carries, the "
            "universe an enrichment statistic on that ranking is computed over: the arm falls "
            f"to {n_arm_rank} and hypoxia to {n_lens_rank}. So the panels report "
            f"{_n(nominal, 'shared')} and {_n(ranked, 'shared')} shared genes and both are "
            "correct, each within its own vocabulary; quoting either without its vocabulary is "
            f"the misreading this figure exists to prevent. The {len(lost_shared)} shared genes "
            f"the left panel carries and the right panel drops are {', '.join(lost_shared)}. "
            "Absence has three causes and the source table splits them: a symbol in the count "
            "matrix and outside the ranking was dropped by filterByExpr, a symbol in the "
            "CellRanger reference and outside the matrix was never detected in sorted T cells, "
            "and a symbol outside the reference is a vocabulary miss. "
            f"Alias resolution runs first and only ever adds, so {led['n_exact']} hypoxia "
            f"genes match exactly and {led['n_via_alias']} more only once their current "
            "symbols resolve into the hg19-vintage vocabulary this matrix carries "
            f"({', '.join(led['alias_pairs'])}) — which is what lifts the testable size from "
            f"{led['n_exact']} to {n_lens_rank}. {arm_recovery_long} This is "
            "membership: the figure and its table carry no NES, FDR, direction or effect size, "
            "and no row reaches effect_sizes_treg_arthritis.csv or any 03_results/master/ "
            "accumulator. A small overlap bounds how much of the arm is hypoxia gene content. "
            "Whether temperature and hypoxia are separable in this niche is undecidable from "
            "cross-sectional human data. Annotation tier."),
        config=FIG_CFG, wide=True, height=6.6,
    )
    plt.close(fig)

    write_caption(
        STAGE, "tables/_overview/arm_hypoxia_euler.csv",
        finding=(f"The plotted regions with their gene names, and the ledger behind the two "
                 f"vocabularies: {led['n_exact']} of the {int(led['n_nominal'])} frozen "
                 f"hypoxia genes match the "
                 f"{EULER_POPULATION.capitalize()} ranked list by exact symbol and "
                 f"{led['n_via_alias']} more only after alias resolution, so what looks like a "
                 f"{led['n_exact']}-gene set is a {n_lens_rank}-gene one."),
        script=SCRIPT, fn="euler_table",
        config_kv=(f"rows = {len(universes)} vocabularies x 3 regions; "
                   f"arm = {EULER_ARM}; lens = {EULER_LENS_SET}; "
                   f"vocabulary = symbol_alias.ranked_list at population={EULER_POPULATION}"),
        input=("03_results/13_arm_decomposition/tables/arm_program_summary.csv, "
               "03_results/13_arm_decomposition/tables/arm_program_gene.csv, "
               "00_data/references/msigdb_hallmark/HALLMARK_HYPOXIA.txt, "
               "00_data/references/symbol_alias/symbol_alias_map.csv, "
               "03_results/03_pseudobulk/tables/ranked_treg.tsv, "
               "03_results/03_pseudobulk/tables/gene_symbols.csv, "
               "03_results/00_build/tables/reference_feature_symbols.csv"),
        how_to_read=(
            "One row per (`universe` x `region`), six rows. `region` is `arm_only`, `shared` or "
            "`lens_only` and `n_genes` is the count the figure draws that area to; `genes` "
            "names them, semicolon-delimited and sorted, so any region can be checked gene by "
            "gene. `universe` is `frozen_sets` for the lists as curated and "
            "`treg_ranked_list` for the restriction, and `vocabulary` records the file the "
            "restriction was made against. The `arm_*` and `lens_*` columns repeat that "
            "universe's ledger on every row: `n_nominal` is the curated size, `n_exact_match` "
            "how many symbols match the ranked list verbatim, `n_via_alias` how many more are "
            "recovered by resolving a current symbol into this matrix's hg19 vintage (named in "
            "`*_alias_pairs_applied`), and the three `n_absent_*` columns split what is left: "
            "`expression_filtered` for symbols the count matrix carries but filterByExpr "
            "dropped, `undetected` for symbols the CellRanger reference carries but the matrix "
            "does not, `absent_from_reference` for a vocabulary miss outright. Those three are "
            "reported separately because collapsing them reads a power fact and a nomenclature "
            "fact as the same biological absence. On `frozen_sets` rows the ledger columns are "
            "trivial by construction, since that universe applies no restriction. "
            "`circle_radius_*`, `centre_distance` and `shared_area_solved` are the geometry the "
            "figure drew, and `shared_area_residual_genes` is how far the drawn shared area "
            "misses the count in gene units — the proof behind `is_area_proportional`. This is "
            "membership: the file carries no NES, p-value or effect size. Annotation tier."),
        config=FIG_CFG)

    print(f"[{STAGE}_viz] wrote 1 overview (arm_hypoxia_euler) over "
          f"{len(universes)} vocabularies; max area residual {max_resid:.2e} genes")
    print(euler[["universe", "region", "n_genes", "n_arm", "n_lens",
                 "shared_area_residual_genes"]].to_string(index=False))


if __name__ == "__main__":
    main()
