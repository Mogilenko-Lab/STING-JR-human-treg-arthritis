#!/usr/bin/env python
"""
14_unbiased_enrichment_viz.py: VIZ ONLY (no statistics).
=============================================================================
Draws one panel: the pre-ranked enrichment score of each mouse-derived up arm on the sorted
JIA synovial-fluid-versus-paired-blood ranked list, in each of the three frozen sort labels.
Every number on the panel is read from a committed table written by
02_analysis/scripts/14_unbiased_enrichment.R. The work here is formatting plus a line count
of the frozen arm files, which supplies the denominator each effective set size is read
against.

The panel sits on the confirmatory tier: donor-level pseudobulk within frozen sort labels,
limma-voom moderated t, then pre-ranked fgsea. That is the only tier in this compartment that
may support a claim, so the tier is stated on the figure face.

Two design constraints are load-bearing.

NO ERROR BARS. A pre-ranked enrichment score carries no standard error and no interval, so
whiskers here would be drawn from nothing. The geometry is an ordered dot plot with an
aligned annotation column, and that column carries the two numbers a reader needs to size
each point: how many of the arm's genes reached that population's ranked list, and the
adjusted p.

ONE MULTIPLE-TESTING FAMILY, NAMED. `gsea_all.csv` carries both a per-collection `padj`
(three tests, the three mouse-derived arms) and a `padj_pooled` corrected across every set
that population's sweep asked about (11,236 in Treg, 11,459 in Tcon, 11,242 in CD8). This
panel uses `padj_pooled` throughout, the more conservative of the two, and names it in the
column header. Both columns travel in the same-stem source table, keeping the other family
checkable.

Input  (03_results/14_unbiased_enrichment/tables/):
  gsea_all.csv                                     NES, both FDR families, set sizes
  ../mouse_anchor/03_results/human_projection/signatures/<arm>/<arm>_up.txt   nominal sizes

Output (03_results/14_unbiased_enrichment/):
  figures/_overview/arm_nes_by_cell_state.{pdf,png}
  tables/_overview/arm_nes_by_cell_state.csv       the nine plotted rows
  README.md                                        caption (via save_overview)

Run in-container from the compartment root, AFTER 14_unbiased_enrichment.R:
  python 02_analysis/scripts/14_unbiased_enrichment_viz.py
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
from matplotlib.lines import Line2D  # noqa: E402

COMPARTMENT_ROOT = Path(__file__).resolve().parents[2]

sys.path.insert(0, str(COMPARTMENT_ROOT))
sys.path.insert(0, str(COMPARTMENT_ROOT / "02_analysis"))
os.chdir(COMPARTMENT_ROOT)

from config import PARAMS, PATHS, POPULATION_COLORS  # noqa: E402
from helpers.figure_style import (  # noqa: E402
    FIG_CFG,
    purge_figures,
    save_overview,
    set_paper_style,
)

STAGE = "14_unbiased_enrichment"
SCRIPT = "02_analysis/scripts/14_unbiased_enrichment_viz.py"
STEM = "arm_nes_by_cell_state"
GATE_DB = "mouse_projection"

# Rows run top to bottom in the order given; the arm subdirectory of the frozen
# mouse-to-human projection contract supplies the nominal size of each arm.
ARMS = [
    ("WT_heat_up", "WT_heat"),
    ("KO_heat_up", "KO_heat"),
    ("Interaction_up", "Interaction"),
]
# Vertical offset of each cell state inside its arm row, so three points on one
# row stay separately readable and each keeps its own annotation line.
CELL_STATES = [("Treg", 0.26), ("Tcon", 0.0), ("CD8", -0.26)]

_F = FIG_CFG["figures"]
_OI = FIG_CFG["colors"]["okabe_ito"]

SZ_TITLE = float(_F["title_size"])
SZ_AXIS_TITLE = float(_F["axis_title_size"])
SZ_AXIS_TEXT = float(_F["axis_text_size"])
SZ_LEGEND = float(_F["legend_text_size"])
SZ_CAPTION = float(_F["caption_size"])
# scatter takes an AREA; the config carries a marker diameter in points.
MARKER_AREA = (float(_F["point_size"]) * 5.0) ** 2
LINE_W = float(_F["line_width"])
# Footnote wrap width in characters, so a footnote never runs past the canvas edge.
FOOT_WRAP = 138

# The one population palette, read from `colors.populations`.
STATE_COLOR = POPULATION_COLORS
GUIDE = _OI["black"]


def nominal_arm_sizes() -> dict:
    """Line count of each frozen mouse-derived up arm, the denominator of a set size.

    An effective set size only means something against the nominal one, so both
    travel with the panel. Read from the same frozen projection contract the sweep
    scored, so numerator and denominator cannot drift apart.
    """
    out = {}
    for arm, subdir in ARMS:
        path = PATHS.signature_contract / "signatures" / subdir / f"{arm}.txt"
        if not path.exists():
            raise FileNotFoundError(f"[14_viz] frozen arm file not found: {path}")
        out[arm] = len({ln.strip() for ln in path.read_text().splitlines() if ln.strip()})
    return out


def fmt_p(p: float) -> str:
    """An adjusted p for an on-face label: fixed above 0.001, else one-digit scientific."""
    if pd.isna(p):
        return "n/a"
    return f"{p:.3f}" if p >= 0.001 else f"{p:.0e}"


def plotted_rows() -> pd.DataFrame:
    """The nine rows behind the panel, read from the committed sweep table."""
    path = PATHS.tables(STAGE) / "gsea_all.csv"
    if not path.exists():
        raise FileNotFoundError(
            f"[14_viz] {path} not found. Run 02_analysis/scripts/14_unbiased_enrichment.R first.")
    sweep = pd.read_csv(path)
    arms = [a for a, _ in ARMS]
    states = [s for s, _ in CELL_STATES]
    sub = sweep[(sweep["database"] == GATE_DB) & sweep["pathway_id"].isin(arms)].copy()

    nominal = nominal_arm_sizes()
    rows = []
    for arm, _ in ARMS:
        for state, _off in CELL_STATES:
            hit = sub[(sub["pathway_id"] == arm) & (sub["population"] == state)]
            if len(hit) != 1:
                raise ValueError(
                    f"[14_viz] expected exactly one {GATE_DB} row for {arm} in {state}, "
                    f"found {len(hit)}")
            r = hit.iloc[0]
            rows.append({
                "arm": arm,
                "cell_state": state,
                "contrast": r["contrast"],
                "nes": float(r["nes"]),
                "pvalue": float(r["pvalue"]),
                "padj_per_collection": float(r["padj"]),
                "padj_pooled": float(r["padj_pooled"]),
                "genes_in_ranked_list": int(r["set_size"]),
                "genes_in_arm": int(nominal[arm]),
                "leading_edge_size": int(r["leading_edge_size"]),
                "n_sets_in_pooled_family": int(r["n_tests_pooled"]),
            })
    out = pd.DataFrame(rows)
    out["arm"] = pd.Categorical(out["arm"], categories=arms, ordered=True)
    out["cell_state"] = pd.Categorical(out["cell_state"], categories=states, ordered=True)
    return out


def build_figure(rows: pd.DataFrame, fdr: float, width: float, height: float):
    """Ordered NES dot plot with an aligned annotation column. No intervals are drawn."""
    fig = plt.figure(figsize=(width, height))
    # Explicit rectangles, in place of a layout engine: the annotation column has to
    # keep its rows aligned with the dots after the exporter fixes the canvas size.
    ax = fig.add_axes((0.075, 0.275, 0.415, 0.585))
    axt = fig.add_axes((0.535, 0.275, 0.455, 0.585), sharey=ax)

    n_arms = len(ARMS)
    y_of = {}
    for i, (arm, _) in enumerate(ARMS):
        base = n_arms - 1 - i                      # first arm on the top row
        for state, off in CELL_STATES:
            y_of[(arm, state)] = base + off

    nes_max = float(rows["nes"].max())
    x_hi = max(3.0, nes_max + 0.35)

    for _, r in rows.iterrows():
        y = y_of[(r["arm"], r["cell_state"])]
        col = STATE_COLOR[r["cell_state"]]
        sig = bool(r["padj_pooled"] < fdr)
        # A light guide runs the width of the panel so the eye carries each dot
        # across to its own annotation line.
        ax.plot([0, x_hi], [y, y], lw=LINE_W * 0.4, color=GUIDE, alpha=0.14, zorder=1)
        ax.scatter(r["nes"], y, s=MARKER_AREA,
                   facecolor=col if sig else "white",
                   edgecolor=col,
                   linewidths=LINE_W * (1.4 if sig else 2.2), zorder=3)

    ax.set_yticks(range(n_arms))
    ax.set_yticklabels([arm for arm, _ in ARMS][::-1], fontsize=SZ_AXIS_TEXT)
    ax.set_ylim(-0.62, n_arms - 1 + 0.80)
    ax.set_xlim(0, x_hi)
    ax.set_xlabel("NES, synovial fluid over paired blood\n(positive = higher in synovial fluid)",
                  fontsize=SZ_AXIS_TITLE)
    ax.tick_params(axis="x", labelsize=SZ_AXIS_TEXT)
    ax.spines["left"].set_visible(True)

    handles = [
        Line2D([0], [0], marker="o", linestyle="", markerfacecolor=STATE_COLOR[s],
               markeredgecolor=STATE_COLOR[s], markersize=SZ_LEGEND * 0.8, label=s)
        for s, _ in CELL_STATES
    ]
    # Upper left: every plotted score sits at NES 1.4 or above, so the low-NES side of
    # the top row is the one region where a box crosses no dot and no guide segment
    # that carries a dot across to its annotation line. The filled-versus-open
    # convention is stated in the footnote instead, which keeps this box to three keys.
    ax.legend(handles=handles, loc="upper left", frameon=True, framealpha=0.92,
              fontsize=SZ_LEGEND, title="cell state", title_fontsize=SZ_LEGEND,
              borderaxespad=0.4, handletextpad=0.4)

    # ---- annotation column, one line per plotted dot ----
    axt.set_axis_off()
    axt.set_xlim(0, 1)
    col_x = {"state": 0.02, "genes": 0.34, "fdr": 0.70}
    head_y = n_arms - 1 + 0.58
    for key, label in (("state", "cell state"),
                       ("genes", "arm genes in list"),
                       ("fdr", "FDR, pooled")):
        axt.text(col_x[key], head_y, label, ha="left", va="center",
                 fontsize=SZ_AXIS_TEXT, fontweight="bold")
    for _, r in rows.iterrows():
        y = y_of[(r["arm"], r["cell_state"])]
        col = STATE_COLOR[r["cell_state"]]
        # The guide continues across the gap so a dot and its annotation line stay tied.
        axt.plot([0, 1], [y, y], lw=LINE_W * 0.4, color=GUIDE, alpha=0.14, zorder=1)
        axt.text(col_x["state"], y, str(r["cell_state"]), ha="left", va="center",
                 fontsize=SZ_AXIS_TEXT, color=col)
        axt.text(col_x["genes"], y,
                 f"{int(r['genes_in_ranked_list'])} of {int(r['genes_in_arm'])}",
                 ha="left", va="center", fontsize=SZ_AXIS_TEXT)
        axt.text(col_x["fdr"], y, fmt_p(float(r["padj_pooled"])), ha="left", va="center",
                 fontsize=SZ_AXIS_TEXT)

    fig.text(0.5, 0.945,
             "Donor pseudobulk NES by cell state, mouse-derived up arms",
             ha="center", va="center", fontsize=SZ_TITLE, fontweight="bold")
    # Hard-wrapped: an unwrapped footnote wider than the canvas makes the tight
    # bounding box grow at export and the panel shrinks inside it.
    foot = "\n".join([
        textwrap.fill("Confirmatory tier: donor-level pseudobulk within frozen sort labels, "
                      "limma-voom moderated t, then pre-ranked fgsea.", FOOT_WRAP),
        textwrap.fill(f"A filled dot clears FDR {fdr:g} and an open dot sits above it. FDR is "
                      "Benjamini-Hochberg pooled across every set the population's sweep tested "
                      f"({rows['n_sets_in_pooled_family'].min():,} to "
                      f"{rows['n_sets_in_pooled_family'].max():,}). A score of this kind carries "
                      "no interval, so none is drawn.", FOOT_WRAP),
    ])
    fig.text(0.075, 0.135, foot, ha="left", va="top", fontsize=SZ_CAPTION, linespacing=1.6)
    return fig


def main() -> None:
    set_paper_style(config=FIG_CFG)
    purge_figures(STAGE, STEM, overview=True, config=FIG_CFG)

    fdr = float(PARAMS.gsea_fdr)
    rows = plotted_rows()
    width, height = 11.0, 6.8
    fig = build_figure(rows, fdr, width, height)

    wt = rows[rows["arm"].eq("WT_heat_up")].set_index("cell_state")
    inter = rows[rows["arm"].eq("Interaction_up")]
    ko = rows[rows["arm"].eq("KO_heat_up")].set_index("cell_state")
    n_inter_sig = int((inter["padj_pooled"] < fdr).sum())

    save_overview(
        fig, STAGE, STEM,
        table=rows,
        finding=(
            "The mouse 39 °C-derived up arm rises on the synovial-fluid side of the paired "
            f"contrast in all three sorted cell states, at NES {wt.loc['Treg', 'nes']:.4f} in "
            f"Treg ({int(wt.loc['Treg', 'genes_in_ranked_list'])} of "
            f"{int(wt.loc['Treg', 'genes_in_arm'])} arm genes reaching the ranked list), "
            f"{wt.loc['Tcon', 'nes']:.4f} in Tcon "
            f"({int(wt.loc['Tcon', 'genes_in_ranked_list'])}) and "
            f"{wt.loc['CD8', 'nes']:.4f} in CD8 "
            f"({int(wt.loc['CD8', 'genes_in_ranked_list'])}), every one below pooled FDR 1e-4. "
            "The Treg score sits between the Tcon and CD8 scores, so the separation reads as "
            f"pan-T. KO_heat_up tracks it row for row (NES {ko.loc['Treg', 'nes']:.4f} in Treg "
            f"on {int(ko.loc['Treg', 'genes_in_ranked_list'])} of "
            f"{int(ko.loc['Treg', 'genes_in_arm'])} genes). The 7-gene interaction arm reaches "
            f"NES {float(inter['nes'].min()):.4f} to {float(inter['nes'].max()):.4f} on "
            f"{int(inter['genes_in_ranked_list'].max())} testable genes and clears pooled FDR "
            f"{fdr:g} in {'none' if n_inter_sig == 0 else n_inter_sig} of the three cell states, "
            f"so at that size it carries no direction here."
        ),
        script=SCRIPT, fn="build_figure",
        config_kv=(f"thresholds.gsea_fdr = {fdr}; gsea_min_size = {PARAMS.gsea_min_size}; "
                   f"gsea_max_size = {PARAMS.gsea_max_size}"),
        input="03_results/14_unbiased_enrichment/tables/gsea_all.csv",
        how_to_read=(
            "One dot per mouse-derived up arm and cell state, at the confirmatory tier: "
            "donor-level pseudobulk within frozen sort labels, limma-voom moderated t, then "
            "pre-ranked fgsea. Rows are the three arms; inside a row the three cell states are "
            "offset vertically and coloured, each with its own annotation line. The x position is "
            "the normalised enrichment score for synovial fluid over paired blood. A filled dot "
            f"clears the config FDR threshold of {fdr:g} and an open dot sits above it. The "
            "annotation column gives how many of the arm's genes reached that population's "
            "ranked list, against how many the frozen arm holds, then the adjusted p. Read that "
            f"count with the score: resolution scales with it, and the arms here span "
            f"{int(rows['genes_in_ranked_list'].min())} to "
            f"{int(rows['genes_in_ranked_list'].max())} testable genes. The adjusted p is "
            "Benjamini-Hochberg pooled across every set that "
            "population's sweep tested; the same-stem source table also carries the "
            "per-collection value over the three arms alone; the two agree on every row but "
            "Interaction_up in CD8, where they read 0.035 per-collection and 0.172 pooled. A "
            "score of this "
            "kind has no interval, so none is drawn. An arm rising here means its gene content "
            "moves with the synovial-fluid side of this ranking; naming follows how the arm was "
            "derived, from mouse iTreg 37 versus 39 °C contrasts, and the reading stays "
            "correlative."
        ),
        config=FIG_CFG, width=width, height=height)

    print(f"[14_unbiased_enrichment_viz] wrote {STEM} from {len(rows)} plotted rows")


if __name__ == "__main__":
    main()
