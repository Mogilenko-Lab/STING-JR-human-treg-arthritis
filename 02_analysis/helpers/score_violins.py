"""
02_analysis/helpers/score_violins.py — the per-cell score violin panel, shared.
=============================================================================
One layout for both violin figures sitting under the per-cell score maps:
`arm_score_violins` and `program_score_violins`. They differ by gene set and nothing
else, so geometry, glyph key and summary arithmetic live here once.

A panel is one gene set: six violins (three frozen sort labels x two tissues), each
panel on its own y axis because AUCell depends on set size, so a level compares inside
a panel only.

Under each panel sits a companion row giving the tissue separation inside each sort
label as Cliff's delta — bounded, unit-free, one shared axis across every panel. It is
there because a label-selective tissue effect is a difference of differences, which
violins alone read worst. Descriptive only: cells pool across donors of unequal yield, so
each pooled delta is drawn over the per-donor deltas behind it and carries no interval.
Testing belongs to the donor-level panels under
03_results/14_unbiased_enrichment/figures/_overview/.

Consumed by:
  02_analysis/scripts/16_narrative_scoring_arms_viz.py
  02_analysis/scripts/16_narrative_scoring_programs_viz.py
"""
from __future__ import annotations

import textwrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.stats import rankdata

CELL_STATES = ["Treg", "Tcon", "CD8"]

DODGE = 0.21
VIOLIN_W = 0.36
# Every set reaches 0.0, and a limit pinned at 0 puts the violin's closing edge on the spine.
Y_PAD_FRAC = 0.04
DONOR_SPREAD = 0.085
# Pinned to the statistic's own bound, so a marker in one figure compares with a marker in
# the other.
DELTA_LIMIT = 1.0
TITLE_CHARS = 27


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    """Cliff's delta of `a` against `b`: 2 * P(a > b) + P(a == b) - 1, in [-1, 1].

    Mid-ranks give the tie-corrected value exactly and stay linear in n, where the n*m
    pairwise form is 200 million comparisons at 15,000 cells a side. Ties are common here:
    a thin set leaves many cells at zero.
    """
    n_a, n_b = len(a), len(b)
    if n_a == 0 or n_b == 0:
        return float("nan")
    r = rankdata(np.concatenate([a, b]))
    u_a = float(r[:n_a].sum()) - n_a * (n_a + 1) / 2.0
    return 2.0 * (u_a / (n_a * n_b)) - 1.0


def summary_table(df: pd.DataFrame, panels: list, tissues: list,
                  donor_key: str = "donor", label_key: str = "coarse_label",
                  tissue_key: str = "tissue") -> pd.DataFrame:
    """The rows behind one figure: one per gene set, sort label and tissue.

    Five-number summary per violin, which is what a reader checking one wants. The
    `delta_*` and `median_shift_*` columns describe the PAIR, so they repeat on both rows
    of a pair.
    """
    num, den = tissues[0][0], tissues[1][0]
    rows = []
    for col, set_name, n_genes in panels:
        for state in CELL_STATES:
            in_state = df[label_key].eq(state)
            sf = df.loc[in_state & df[tissue_key].eq(num), col].to_numpy()
            pb = df.loc[in_state & df[tissue_key].eq(den), col].to_numpy()
            delta = cliffs_delta(sf, pb)
            per_donor = donor_deltas(df, col, state, tissues, donor_key,
                                     label_key, tissue_key)
            shift = float(np.median(sf) - np.median(pb))
            for tissue, tissue_label in tissues:
                v = df.loc[in_state & df[tissue_key].eq(tissue), col].to_numpy()
                rows.append({
                    "gene_set": set_name,
                    "genes_scored": n_genes,
                    "cell_state": state,
                    "tissue": tissue_label,
                    "n_cells": int(v.size),
                    "n_donors": int(df.loc[in_state & df[tissue_key].eq(tissue),
                                           donor_key].nunique()),
                    "mean": float(np.mean(v)),
                    "median": float(np.median(v)),
                    "q25": float(np.percentile(v, 25)),
                    "q75": float(np.percentile(v, 75)),
                    "min": float(np.min(v)),
                    "max": float(np.max(v)),
                    "frac_at_zero": float(np.mean(v == 0.0)),
                    "median_shift_sf_minus_pb": shift,
                    "delta_cliffs_sf_vs_pb": delta,
                    "delta_n_donors_paired": len(per_donor),
                    "delta_donor_min": float(min(per_donor)) if per_donor else np.nan,
                    "delta_donor_max": float(max(per_donor)) if per_donor else np.nan,
                    "delta_donors_agreeing_in_sign":
                        int(sum(np.sign(d) == np.sign(delta) for d in per_donor)),
                })
    return pd.DataFrame(rows)


def donor_deltas(df: pd.DataFrame, col: str, state: str, tissues: list,
                 donor_key: str = "donor", label_key: str = "coarse_label",
                 tissue_key: str = "tissue") -> list:
    """Cliff's delta inside each donor carrying BOTH tissues in this sort label.

    The contrast is paired, so a donor carrying one arm alone sits the comparison out and
    the count per sort label is reported.
    """
    num, den = tissues[0][0], tissues[1][0]
    sub = df[df[label_key].eq(state)]
    out = []
    for _donor, g in sub.groupby(donor_key, observed=True):
        sf = g.loc[g[tissue_key].eq(num), col].to_numpy()
        pb = g.loc[g[tissue_key].eq(den), col].to_numpy()
        if len(sf) and len(pb):
            out.append(cliffs_delta(sf, pb))
    return out


def panel_titles(panels: list) -> list:
    """Titles wrapped on underscores, padded to uniform depth. The identifier stays whole,
    so a panel is findable in the scoring manifest and in the sweep."""
    wrapped = []
    for _col, name, _n in panels:
        lines, cur = [], ""
        for part in name.split("_"):
            candidate = part if not cur else f"{cur}_{part}"
            if len(candidate) > TITLE_CHARS and cur:
                lines.append(cur)
                cur = f"_{part}"
            else:
                cur = candidate
        lines.append(cur)
        wrapped.append(lines)
    depth = max(len(w) for w in wrapped)
    return ["\n" * (depth - len(w)) + "\n".join(w) + f"\n{n} genes scored"
            for w, (_c, _s, n) in zip(wrapped, panels)]


def build_figure(df: pd.DataFrame, panels: list, tissues: list, tissue_colors: dict,
                 state_colors: dict, summary: pd.DataFrame, cfg: dict,
                 title: str, subtitle: str, width: float, height: float, n_col: int = 3,
                 donor_key: str = "donor", label_key: str = "coarse_label",
                 tissue_key: str = "tissue"):
    """Draw the whole figure: `n_col` columns of violin panel over delta row.

    `panels` is [(column, set_name, genes_scored)] in draw order, matching the panel order
    of the map this figure sits under.
    """
    f = cfg["figures"]
    ink = cfg["colors"]["okabe_ito"]["black"]
    sz_title = float(f["title_size"])
    sz_sub = float(f["subtitle_size"])
    sz_axis_title = float(f["axis_title_size"])
    sz_axis_text = float(f["axis_text_size"])
    sz_legend = float(f["legend_text_size"])
    sz_strip = float(f["strip_size"])
    sz_caption = float(f["caption_size"])
    lw = float(f["line_width"])

    n_row = int(np.ceil(len(panels) / n_col))
    fig = plt.figure(figsize=(width, height))

    # Explicit rectangles: the delta row must stay locked under its own violins after the
    # exporter fixes the canvas size.
    left, gap_x = 0.062, 0.052
    span = (0.975 - left - (n_col - 1) * gap_x) / n_col
    row_base = [0.535, 0.140][:n_row] if n_row == 2 else [0.140]
    h_delta, h_gap, h_violin = 0.074, 0.017, 0.216

    donor_values = {}
    for col, set_name, _n in panels:
        for state in CELL_STATES:
            donor_values[(set_name, state)] = donor_deltas(
                df, col, state, tissues, donor_key, label_key, tissue_key)
    dlim = DELTA_LIMIT
    titles = panel_titles(panels)

    x_pos = list(range(len(CELL_STATES)))
    for idx, (col, set_name, n_genes) in enumerate(panels):
        r, c = divmod(idx, n_col)
        x0 = left + c * (span + gap_x)
        base = row_base[r]
        axd = fig.add_axes((x0, base, span, h_delta))
        axv = fig.add_axes((x0, base + h_delta + h_gap, span, h_violin))

        drawn_max = 0.0
        for xi, state in enumerate(CELL_STATES):
            in_state = df[label_key].eq(state)
            for sgn, (tissue, _label) in zip((-1, 1), tissues):
                v = df.loc[in_state & df[tissue_key].eq(tissue), col].to_numpy()
                drawn_max = max(drawn_max, float(np.max(v)))
                parts = axv.violinplot([v], positions=[xi + sgn * DODGE], widths=VIOLIN_W,
                                       showextrema=False, showmedians=True)
                for body in parts["bodies"]:
                    body.set_facecolor(tissue_colors[tissue])
                    body.set_edgecolor(tissue_colors[tissue])
                    body.set_alpha(0.72)
                    body.set_linewidth(lw * 0.8)
                med = parts["cmedians"]
                med.set_color(ink)
                med.set_linewidth(lw * 1.3)

        axv.set_xticks(x_pos)
        axv.set_xticklabels([])
        axv.set_xlim(-0.55, len(CELL_STATES) - 0.45)
        axv.tick_params(axis="y", labelsize=sz_axis_text)
        axv.set_title(titles[idx], fontsize=sz_strip, fontweight="bold", linespacing=1.15)
        # Headroom both ends; the tick the lower pad opens below zero is dropped, since
        # AUCell is bounded in [0, 1].
        pad = Y_PAD_FRAC * drawn_max
        axv.set_ylim(-pad, drawn_max + pad)
        axv.set_yticks([t for t in axv.get_yticks() if 0.0 <= t <= drawn_max + pad])

        # ---- the delta row ----
        axd.axhline(0.0, color=ink, lw=lw * 0.8, alpha=0.5, zorder=1)
        for xi, state in enumerate(CELL_STATES):
            colour = state_colors[state]
            hit = summary[summary["gene_set"].eq(set_name)
                          & summary["cell_state"].eq(state)].iloc[0]
            d = float(hit["delta_cliffs_sf_vs_pb"])
            per_donor = donor_values[(set_name, state)]
            # Fixed lattice, so the figure is identical on every run.
            if per_donor:
                offs = np.linspace(-DONOR_SPREAD, DONOR_SPREAD, len(per_donor))
                axd.scatter(xi + offs, per_donor, s=(sz_legend * 0.42) ** 2,
                            facecolor="white", edgecolor=colour,
                            linewidths=lw * 0.8, zorder=2)
            axd.plot([xi, xi], [0.0, d], color=colour, lw=lw * 1.6, zorder=3,
                     solid_capstyle="butt")
            axd.scatter([xi], [d], s=(sz_legend * 0.78) ** 2, marker="D",
                        facecolor=colour, edgecolor=colour, zorder=4)

        axd.set_xticks(x_pos)
        axd.set_xticklabels(CELL_STATES, fontsize=sz_axis_text)
        axd.set_xlim(-0.55, len(CELL_STATES) - 0.45)
        axd.set_ylim(-dlim, dlim)
        # Half-scale ticks: most sets sit well inside the axis, and comparing two markers a
        # panel apart needs an interior reference.
        axd.set_yticks([-dlim, -dlim / 2, 0.0, dlim / 2, dlim])
        axd.set_yticklabels([f"{-dlim:g}", "", "0", "", f"{dlim:g}"])
        axd.axhline(dlim / 2, color=ink, lw=lw * 0.4, alpha=0.18, zorder=1)
        axd.axhline(-dlim / 2, color=ink, lw=lw * 0.4, alpha=0.18, zorder=1)
        axd.tick_params(axis="y", labelsize=sz_axis_text * 0.86, length=2.5)
        axd.set_facecolor("#F7F7F7")

        if c == 0:
            axv.set_ylabel("AUCell, per cell", fontsize=sz_axis_title)
            axd.set_ylabel("Cliff's δ", fontsize=sz_axis_title * 0.86)

    handles = [Patch(facecolor=tissue_colors[t], edgecolor=tissue_colors[t], alpha=0.72,
                     label=lab) for t, lab in tissues]
    handles += [
        Line2D([0], [0], marker="D", linestyle="", markerfacecolor=ink,
               markeredgecolor=ink, markersize=sz_legend * 0.62,
               label="Cliff's δ, all cells of that sort label"),
        Line2D([0], [0], marker="o", linestyle="", markerfacecolor="white",
               markeredgecolor=ink, markersize=sz_legend * 0.52,
               label="Cliff's δ within one donor"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.062),
               ncol=4, frameon=False, fontsize=sz_legend)

    fig.text(0.5, 0.962, title, ha="center", va="center",
             fontsize=sz_title, fontweight="bold")
    fig.text(0.5, 0.929, subtitle, ha="center", va="center", fontsize=sz_sub)
    # Glyph key only; how to read it lives in `how_to_read` and lands in the README.
    fig.text(0.020, 0.036, face_note(df, summary, dlim, donor_key),
             ha="left", va="top", fontsize=sz_caption, linespacing=1.5)
    return fig, dlim


def face_note(df: pd.DataFrame, summary: pd.DataFrame, dlim: float,
              donor_key: str = "donor") -> str:
    """The glyph key drawn on the figure face."""
    n_cells = int(len(df))
    n_donors = int(df[donor_key].nunique())
    paired = sorted(set(int(v) for v in summary["delta_n_donors_paired"]))
    paired_txt = (f"{paired[0]}" if len(paired) == 1
                  else f"{paired[0]} to {paired[-1]}")
    return textwrap.fill(
        f"Black line marks the median of a violin. The grey row under each panel gives "
        f"Cliff's δ of synovial fluid against paired blood inside one sort label — the "
        f"chance a synovial cell outscores a blood cell of that label, ties counted as "
        f"half, on one shared ±{dlim:g} axis across all panels. Open circles behind each "
        f"δ are the same quantity computed inside one donor ({paired_txt} donors carry "
        f"both tissues in a sort label). Descriptive: {n_cells:,} sorted cells from "
        f"{n_donors} donors, one vote per cell, no interval and no test.", 168)


def how_to_read(df: pd.DataFrame, summary: pd.DataFrame, panels: list, dlim: float,
                map_figure: str, confirmatory_figure: str,
                donor_key: str = "donor") -> str:
    """The README `how to read`, shared by both figures bar their two names."""
    n_cells = int(len(df))
    n_donors = int(df[donor_key].nunique())
    per_donor = df[donor_key].value_counts()
    lo, hi = int(per_donor.min()), int(per_donor.max())
    return (
        f"Annotation tier. One panel per gene set, in the panel order of {map_figure}, so "
        "the two figures lay side by side and a column of that map has a distribution "
        "here. Inside a panel the x axis is the three frozen sort labels and the two "
        "violins of a label are the two tissues: warm is synovial fluid, cool is paired "
        "peripheral blood, black line at the median. AUCell is a rank-based score in 0 to "
        "1, the area under a cell's gene-recovery curve for that set, so it is robust to "
        "library size and composition. A panel title counts the genes the score was really "
        "computed over — the set's genes present in this object after symbol resolution, "
        "which is the same count the map's own panel title carries and is smaller than the "
        "set's nominal size. Each panel keeps its own y axis because the sets "
        f"range from {min(n for _c, _s, n in panels)} to "
        f"{max(n for _c, _s, n in panels)} genes scored and AUCell is computed against "
        "each cell's own ranking, so a level compares within a panel and a shape compares "
        "anywhere. Both ends of every y axis carry headroom, and the score is bounded in "
        "[0, 1]. "
        "The grey row under each panel answers the question the violins are worst at, "
        "which is whether a tissue difference is bigger in one sort label than another. It "
        "gives Cliff's δ of synovial fluid against paired blood inside one sort label: the "
        "probability that a randomly drawn synovial cell outscores a randomly drawn blood "
        "cell of the same label, ties counted as half, rescaled onto -1 to +1. It is "
        f"unit-free and bounded, so ONE ±{dlim:g} axis serves every panel and a difference "
        "between sort labels in one panel is legible against the same difference in any "
        "other. A label-selective tissue effect is a δ that stands away from the other two "
        "in its own panel and does not do so in the panels beside it. "
        "Behind each pooled δ sit the per-donor δ values, each computed inside one donor's "
        "own paired cells; a donor missing one tissue in a sort label contributes none, "
        "which is why the count differs between labels. Read the donor cloud first: a "
        "pooled marker whose donors straddle zero is one donor's result. "
        f"Every one of the {n_cells:,} cells casts one vote and the {n_donors} donors "
        f"contributed {lo:,} to {hi:,} cells each, so both rows are pseudoreplicated and "
        "neither carries a p-value or an interval. Ranking the sort labels, and testing "
        f"any of this, is the job of the donor-level panel {confirmatory_figure}, where "
        "each donor carries one vote inside a frozen label. The same-stem source table "
        f"gives the cell count, mean, median, quartiles, range and zero fraction of all "
        f"{len(summary)} violins together with each pair's δ, its median shift and its "
        "donor range. Naming follows how each set was derived and the reading stays "
        "correlative."
    )
