#!/usr/bin/env python
"""
14_program_nes_by_cell_state_viz.py: VIZ ONLY (no statistics).
=============================================================================
The confirmatory companion to `arm_nes_by_cell_state`. Same geometry, same cell-state colours,
same filled-versus-open significance convention, same annotation column, same tier line, so
the two panels lie side by side and read identically. Where that panel carries the three
mouse-derived up arms, this one carries the two program families a reader will ask about next,
in three labelled groups (see GROUPS): the four oxygen-named sets, the six cGAS-STING and
interferon sets, and the five remaining panels of the per-cell score maps.

EVERY SET THE MAPS COLOUR BY HAS A ROW ON THIS GEOMETRY. `assert_map_panels_covered` checks
this panel plus its sibling against `analysis_config.yaml::percell_map_panels`, the
declaration those maps read, so a lens met on a map always has a donor-level row to turn to.

WHY THE FOUR OXYGEN SETS TRAVEL TOGETHER. They are four differently-built sets for one named
biology, and they disagree with each other here and in the mouse anchor. Drawing
HALLMARK_HYPOXIA on its own invites a reader to take one set's score as the family's, so the
family is drawn.

ONE DERIVATION PER NUMBER. A row the committed `tables/_overview/named_sets_in_sweep.csv`
already carries is read from there, with its verified `nes`, `padj_pooled`, `set_size` and
rank columns; the rest come from `tables/gsea_all.csv`, the sweep output that table was
built from. Each output row records which table it came from, which keeps the two committed
tables in agreement.

SET SIZE IS PART OF THE COMPARISON. A large set clears pooled significance far more often
than a small one in this sweep, so the gene count reaching each population's ranked list is
printed beside every dot and the sweep's own per-band baseline rates are quoted in the
caption.

ABSENCE IS DRAWN. A set that fell below `gsea_min_size` in a population has no row in the
sweep. Those cells are labelled `not tested in <population>`, which keeps an untested cell
distinct from a null result.

Input  (03_results/14_unbiased_enrichment/tables/):
  _overview/named_sets_in_sweep.csv    verified rows for the sets that table carries
  gsea_all.csv                         NES, both FDR families, set sizes
  sweep_setsize_baseline.csv           per-band pooled-significance base rates

Output (03_results/14_unbiased_enrichment/):
  figures/_overview/program_nes_by_cell_state.{pdf,png}
  tables/_overview/program_nes_by_cell_state.csv    one row per plotted cell
  README.md                                          caption (via save_overview)

Run in-container from the compartment root, AFTER 14_unbiased_enrichment.R and
14_sweep_named_sets.R:
  python 02_analysis/scripts/14_program_nes_by_cell_state_viz.py
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

from config import MAP_PANEL_SWEEP_ID, PARAMS, PATHS, POPULATION_COLORS  # noqa: E402
from helpers.figure_style import (  # noqa: E402
    FIG_CFG,
    purge_figures,
    save_overview,
    set_paper_style,
)

STAGE = "14_unbiased_enrichment"
SCRIPT = "02_analysis/scripts/14_program_nes_by_cell_state_viz.py"
STEM = "program_nes_by_cell_state"

NAMED_SETS_CSV = "tables/_overview/named_sets_in_sweep.csv"
GSEA_ALL_CSV = "tables/gsea_all.csv"
BASELINE_CSV = "tables/sweep_setsize_baseline.csv"
ARM_NES_CSV = "tables/_overview/arm_nes_by_cell_state.csv"   # the sibling panel

# The three groups, in draw order. Rows inside a group are ordered by descending Treg
# NES, so the ordering belongs to one population and is stated in the caption.
GROUPS = [
    ("Oxygen and HIF response", [
        "HALLMARK_HYPOXIA",
        "GOBP_RESPONSE_TO_OXYGEN_LEVELS",
        "GOBP_CELLULAR_RESPONSE_TO_OXYGEN_LEVELS",
        "REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA",
    ]),
    ("cGAS-STING and interferon", [
        "sting_specific_up",
        "GOBP_CGAS_STING_SIGNALING_PATHWAY",
        "REACTOME_STING_MEDIATED_INDUCTION_OF_HOST_IMMUNE_RESPONSES",
        "WP_STING_PATHWAY_IN_KAWASAKILIKE_DISEASE_AND_COVID19",
        "ifn_only_up",
        "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    ]),
    # The remaining per-cell map panels, so every lens a reader meets on a map has a
    # donor-level row somewhere on this geometry.
    ("Proteostasis, inflammation and the effector-Treg reference", [
        "HSR_core",
        "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
        "HALLMARK_INFLAMMATORY_RESPONSE",
        "HALLMARK_IL2_STAT5_SIGNALING",
        "eTreg_up",
    ]),
]
ORDER_BY_POPULATION = "Treg"

CELL_STATES = [("Treg", 0.26), ("Tcon", 0.0), ("CD8", -0.26)]
GROUP_GAP = 1.05          # blank data units inserted between groups
ROW_LABEL_WRAP = 30       # characters per line of a wrapped set identifier

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
FOOT_WRAP = 150

# The one population palette from `colors.populations`, the same source
# arm_nes_by_cell_state reads, so the two panels read as one pair.
STATE_COLOR = POPULATION_COLORS
GUIDE = _OI["black"]


def fmt_p(p: float) -> str:
    """An adjusted p for an on-face label: fixed above 0.001, else one-digit scientific."""
    if pd.isna(p):
        return "n/a"
    return f"{p:.3f}" if p >= 0.001 else f"{p:.0e}"


def wrap_set_id(set_id: str) -> str:
    """Wrap a long set identifier on its underscores so the axis label stays intact.

    The identifier is never abbreviated: a reader has to be able to find the row in
    `gsea_all.csv` from what the axis says.
    """
    parts = set_id.split("_")
    lines, cur = [], ""
    for part in parts:
        candidate = part if not cur else f"{cur}_{part}"
        if len(candidate) > ROW_LABEL_WRAP and cur:
            # The separator opens the continuation line, leaving the one
            # above, so no line ends on a dangling underscore.
            lines.append(cur)
            cur = f"_{part}"
        else:
            cur = candidate
    if cur:
        lines.append(cur)
    return "\n".join(lines)


def assert_map_panels_covered() -> int:
    """Halt unless every set the per-cell score maps colour by has a donor-level row.

    A reader arrives here from a map, so a lens drawn there with no row on this geometry
    leaves that colouring untested on the page.
    This panel and `arm_nes_by_cell_state` are one surface drawn twice, so the three mouse
    arms count as covered by the sibling, read from its committed table.
    Returns the number of panels checked.
    """
    sibling = PATHS.results / STAGE / ARM_NES_CSV
    if not sibling.exists():
        raise FileNotFoundError(
            f"[14_program_viz] {sibling} not found; the coverage audit reads the sibling "
            "panel's own table. Run 14_unbiased_enrichment_viz.py first.")
    drawn = {s for _g, ids in GROUPS for s in ids}
    drawn |= set(pd.read_csv(sibling)["arm"].astype(str))
    missing = sorted(v for v in MAP_PANEL_SWEEP_ID.values() if v not in drawn)
    if missing:
        raise ValueError(
            f"[14_program_viz] {len(missing)} set(s) drawn on a per-cell score map have no "
            f"row on this geometry: {missing}. Add them to GROUPS, or drop them from "
            "analysis_config.yaml::percell_map_panels.")
    return len(MAP_PANEL_SWEEP_ID)


def plotted_rows() -> pd.DataFrame:
    """The cells behind the panel, each read from a committed sweep table.

    A set that `named_sets_in_sweep.csv` already verified is taken from there, and the
    rest from `gsea_all.csv`, the sweep output that table was itself built from. A
    (set, population) cell with no row in either is carried with `tested = False` so the
    absence can be drawn.
    """
    named_path = PATHS.results / STAGE / NAMED_SETS_CSV
    gsea_path = PATHS.results / STAGE / GSEA_ALL_CSV
    for path, producer in ((named_path, "14_sweep_named_sets.R"),
                           (gsea_path, "14_unbiased_enrichment.R")):
        if not path.exists():
            raise FileNotFoundError(f"[14_program_viz] {path} not found. Run {producer} first.")
    named = pd.read_csv(named_path)
    sweep = pd.read_csv(gsea_path)

    states = [s for s, _ in CELL_STATES]
    rows = []
    for group, set_ids in GROUPS:
        for set_id in set_ids:
            in_named = set_id in set(named["pathway_id"])
            for state in states:
                base = {"group": group, "pathway_id": set_id, "cell_state": state}
                if in_named:
                    hit = named[(named["pathway_id"] == set_id)
                                & (named["population"] == state)]
                    src = NAMED_SETS_CSV
                else:
                    hit = sweep[(sweep["pathway_id"] == set_id)
                                & (sweep["population"] == state)]
                    src = GSEA_ALL_CSV
                if len(hit) > 1:
                    raise ValueError(
                        f"[14_program_viz] expected at most one row for {set_id} in "
                        f"{state}, found {len(hit)} in {src}")
                if len(hit) == 0 or (in_named and not bool(hit.iloc[0]["tested"])):
                    # Below gsea_min_size in this population, so the sweep never
                    # tested it. Carried as an explicit absence.
                    rows.append({**base, "source_table": src, "tested": False,
                                 "contrast": pd.NA, "nes": pd.NA, "pvalue": pd.NA,
                                 "padj_per_collection": pd.NA, "padj_pooled": pd.NA,
                                 "significant_pooled": pd.NA,
                                 "genes_in_ranked_list": pd.NA,
                                 "leading_edge_size": pd.NA,
                                 "rank_padj_pooled": pd.NA, "rank_nes_signed": pd.NA,
                                 "n_sets_in_pooled_family": pd.NA})
                    continue
                r = hit.iloc[0]
                rows.append({
                    **base,
                    "source_table": src,
                    "tested": True,
                    # named_sets_in_sweep keys on `population` alone, gsea_all names the
                    # contrast per population. One spelling for both, so the column is uniform.
                    "contrast": (r["contrast"] if "contrast" in r.index
                                 else f"SF_vs_PB_{state}"),
                    "nes": float(r["nes"]),
                    "pvalue": float(r["pvalue"]),
                    "padj_per_collection": float(r["padj"]),
                    "padj_pooled": float(r["padj_pooled"]),
                    "significant_pooled": bool(float(r["padj_pooled"]) < float(PARAMS.gsea_fdr)),
                    "genes_in_ranked_list": int(r["set_size"]),
                    "leading_edge_size": int(r["leading_edge_size"]),
                    # Two different ranks travel with a named set: its place by pooled FDR
                    # and its place by signed NES. They are not the same number, so both
                    # are carried and neither is quoted without saying which.
                    "rank_padj_pooled": (int(r["rank_padj_pooled"])
                                         if "rank_padj_pooled" in r.index
                                         and pd.notna(r["rank_padj_pooled"]) else pd.NA),
                    "rank_nes_signed": (int(r["rank_nes_signed"])
                                        if "rank_nes_signed" in r.index
                                        and pd.notna(r["rank_nes_signed"]) else pd.NA),
                    "n_sets_in_pooled_family": int(r["n_tests_pooled"]),
                })
    out = pd.DataFrame(rows)
    out["cell_state"] = pd.Categorical(out["cell_state"], categories=states, ordered=True)
    return out


def assert_one_fdr_family(rows: pd.DataFrame) -> int:
    """Every plotted FDR is `padj_pooled`. Verify it here on every run.

    The rows come from two committed tables, and both carry a per-collection `padj`
    alongside the pooled one,
    so a column that silently mixed the two families would look plausible and be wrong.
    For every row taken from the named-sets table this re-reads the same cell out of
    `gsea_all.csv` and requires the two pooled values to agree exactly. Returns the
    number of cells cross-checked.
    """
    sweep = pd.read_csv(PATHS.results / STAGE / GSEA_ALL_CSV)
    checked = 0
    for _, r in rows[rows["tested"] & rows["source_table"].eq(NAMED_SETS_CSV)].iterrows():
        hit = sweep[(sweep["pathway_id"] == r["pathway_id"])
                    & (sweep["population"] == r["cell_state"])]
        if len(hit) != 1:
            raise ValueError(f"[14_program_viz] {r['pathway_id']} in {r['cell_state']} has "
                             f"{len(hit)} rows in {GSEA_ALL_CSV}, expected exactly one")
        ref = float(hit.iloc[0]["padj_pooled"])
        if abs(float(r["padj_pooled"]) - ref) > 1e-12 * max(1.0, abs(ref)):
            raise ValueError(
                f"[14_program_viz] pooled FDR disagrees between tables for "
                f"{r['pathway_id']} in {r['cell_state']}: {float(r['padj_pooled']):.6e} from "
                f"{NAMED_SETS_CSV} against {ref:.6e} from {GSEA_ALL_CSV}. The plotted column "
                "would be mixing two multiple-testing families.")
        checked += 1
    n_family = rows.loc[rows["tested"], "n_sets_in_pooled_family"].nunique()
    if n_family != len(CELL_STATES):
        raise ValueError(f"[14_program_viz] expected one pooled-family size per cell state, "
                         f"found {n_family} distinct values")
    return checked


def order_rows(rows: pd.DataFrame) -> list:
    """(group, set_id) in draw order: groups as declared, sets by descending Treg NES."""
    ordered = []
    for group, set_ids in GROUPS:
        key = rows[(rows["group"] == group) & (rows["cell_state"] == ORDER_BY_POPULATION)]
        key = key.set_index("pathway_id")["nes"]
        ranked = sorted(set_ids, key=lambda s: (-float(key[s]) if pd.notna(key.get(s))
                                                else float("inf")))
        ordered.extend((group, s) for s in ranked)
    return ordered


def baseline_band(baseline: pd.DataFrame, population: str, band: str) -> pd.Series:
    hit = baseline[(baseline["population"] == population) & (baseline["band_label"] == band)]
    if hit.empty:
        raise ValueError(f"[14_program_viz] no baseline row for {population} / {band}")
    return hit.iloc[0]


def build_figure(rows: pd.DataFrame, order: list, fdr: float, width: float, height: float):
    """Grouped NES dot plot with an aligned annotation column. No intervals are drawn."""
    fig = plt.figure(figsize=(width, height))
    # Explicit rectangles, in place of a layout engine: the annotation column has to
    # keep its rows aligned with the dots after the exporter fixes the canvas size.
    ax = fig.add_axes((0.235, 0.185, 0.375, 0.720))
    axt = fig.add_axes((0.632, 0.185, 0.330, 0.720), sharey=ax)

    # y positions, top row first, with a blank band between groups.
    y_of_row, group_top, cursor = {}, {}, 0.0
    prev_group = None
    for group, set_id in order:
        if prev_group is not None and group != prev_group:
            cursor += GROUP_GAP
        y_of_row[(group, set_id)] = cursor
        group_top.setdefault(group, cursor)
        cursor += 1.0
        prev_group = group
    span = cursor - 1.0
    # Flip so the first declared row sits at the top of the panel.
    y_of_row = {k: span - v for k, v in y_of_row.items()}
    group_top = {g: span - v for g, v in group_top.items()}

    def dot_y(group, set_id, state_off):
        return y_of_row[(group, set_id)] + state_off

    tested = rows[rows["tested"]]
    nes_lo = float(tested["nes"].min())
    nes_hi = float(tested["nes"].max())
    x_lo = min(-0.25, nes_lo - 0.35)
    x_hi = max(3.0, nes_hi + 0.35)

    for (group, set_id) in order:
        for state, off in CELL_STATES:
            y = dot_y(group, set_id, off)
            col = STATE_COLOR[state]
            # A light guide runs the width of the panel so the eye carries each dot
            # across to its own annotation line.
            ax.plot([x_lo, x_hi], [y, y], lw=LINE_W * 0.4, color=GUIDE, alpha=0.14, zorder=1)
            hit = rows[(rows["group"] == group) & (rows["pathway_id"] == set_id)
                       & (rows["cell_state"] == state)].iloc[0]
            if not bool(hit["tested"]):
                continue                      # absence is written in the annotation column
            sig = bool(hit["significant_pooled"])
            ax.scatter(float(hit["nes"]), y, s=MARKER_AREA,
                       facecolor=col if sig else "white",
                       edgecolor=col,
                       linewidths=LINE_W * (1.4 if sig else 2.2), zorder=3)

    # NES crosses zero in this panel, so the no-change position is marked.
    ax.axvline(0.0, color=GUIDE, lw=LINE_W * 0.9, alpha=0.55, zorder=2)

    ax.set_yticks([y_of_row[k] for k in order])
    ax.set_yticklabels([wrap_set_id(s) for _, s in order], fontsize=SZ_AXIS_TEXT)
    ax.set_ylim(-0.62, span + 0.80)
    ax.set_xlim(x_lo, x_hi)
    ax.set_xlabel("NES, synovial fluid over paired blood\n(positive = higher in synovial fluid)",
                  fontsize=SZ_AXIS_TITLE)
    ax.tick_params(axis="x", labelsize=SZ_AXIS_TEXT)
    ax.spines["left"].set_visible(True)

    handles = [
        Line2D([0], [0], marker="o", linestyle="", markerfacecolor=STATE_COLOR[s],
               markeredgecolor=STATE_COLOR[s], markersize=SZ_LEGEND * 0.8, label=s)
        for s, _ in CELL_STATES
    ]
    # OUTSIDE the panel, in the band between the title and the plot area, laid out in one
    # row. Ten rows of set identifiers and two group headers leave no interior region a
    # box can occupy without covering something, so the legend does not sit inside.
    # The filled-versus-open convention is stated in the footnote, which keeps this to
    # three keys.
    ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.0, 1.012), ncol=3,
              frameon=True, framealpha=0.92,
              fontsize=SZ_LEGEND, title="cell state", title_fontsize=SZ_LEGEND,
              borderaxespad=0.0, handletextpad=0.4, columnspacing=1.6)

    # ---- annotation column, one line per plotted dot ----
    axt.set_axis_off()
    axt.set_xlim(0, 1)
    col_x = {"state": 0.02, "genes": 0.37, "fdr": 0.75}
    head_y = span + 0.58
    for key, label in (("state", "cell state"),
                       ("genes", "genes in list"),
                       ("fdr", "FDR, pooled")):
        axt.text(col_x[key], head_y, label, ha="left", va="center",
                 fontsize=SZ_AXIS_TEXT, fontweight="bold")
    for (group, set_id) in order:
        for state, off in CELL_STATES:
            y = dot_y(group, set_id, off)
            col = STATE_COLOR[state]
            # The guide continues across the gap so a dot and its annotation line stay tied.
            axt.plot([0, 1], [y, y], lw=LINE_W * 0.4, color=GUIDE, alpha=0.14, zorder=1)
            axt.text(col_x["state"], y, state, ha="left", va="center",
                     fontsize=SZ_AXIS_TEXT, color=col)
            hit = rows[(rows["group"] == group) & (rows["pathway_id"] == set_id)
                       & (rows["cell_state"] == state)].iloc[0]
            if not bool(hit["tested"]):
                axt.text(col_x["genes"], y, f"not tested in {state}", ha="left", va="center",
                         fontsize=SZ_AXIS_TEXT, style="italic", color=GUIDE, alpha=0.75)
                continue
            axt.text(col_x["genes"], y, f"{int(hit['genes_in_ranked_list'])}",
                     ha="left", va="center", fontsize=SZ_AXIS_TEXT)
            axt.text(col_x["fdr"], y, fmt_p(float(hit["padj_pooled"])), ha="left", va="center",
                     fontsize=SZ_AXIS_TEXT)

    # ---- group labels and the separator between them ----
    for group, _ in GROUPS:
        ax.text(x_lo, group_top[group] + 0.62, group, ha="left", va="center",
                fontsize=SZ_AXIS_TITLE, fontweight="bold")
    seps = sorted(group_top.values())[:-1]
    for y in seps:
        for a in (ax, axt):
            a.axhline(y + 0.95, color=GUIDE, lw=LINE_W * 0.7, alpha=0.35,
                      linestyle=(0, (4, 3)), zorder=1)

    fig.text(0.5, 0.972,
             "Donor pseudobulk NES by cell state, every program the score maps colour by",
             ha="center", va="center", fontsize=SZ_TITLE, fontweight="bold")
    # Hard-wrapped: an unwrapped footnote wider than the canvas makes the tight
    # bounding box grow at export and the panel shrinks inside it.
    foot = "\n".join([
        textwrap.fill("Confirmatory tier: donor-level pseudobulk within frozen sort labels, "
                      "limma-voom moderated t, then pre-ranked fgsea.", FOOT_WRAP),
        textwrap.fill(f"A filled dot clears FDR {fdr:g} and an open dot sits above it. FDR is "
                      "Benjamini-Hochberg pooled across every set the population's sweep tested "
                      f"({int(tested['n_sets_in_pooled_family'].min()):,} to "
                      f"{int(tested['n_sets_in_pooled_family'].max()):,}). A score of this kind "
                      "carries no interval, so none is drawn.", FOOT_WRAP),
        textwrap.fill("Set size drives much of this comparison, so the number of the set's genes "
                      "reaching that population's ranked list is given beside every dot. Rows "
                      f"inside a group run by descending {ORDER_BY_POPULATION} NES.", FOOT_WRAP),
    ])
    fig.text(0.020, 0.108, foot, ha="left", va="top", fontsize=SZ_CAPTION, linespacing=1.6)
    return fig


def main() -> None:
    set_paper_style(config=FIG_CFG)
    purge_figures(STAGE, STEM, overview=True, config=FIG_CFG)

    fdr = float(PARAMS.gsea_fdr)
    rows = plotted_rows()
    order = order_rows(rows)
    baseline = pd.read_csv(PATHS.results / STAGE / BASELINE_CSV)

    n_panels = assert_map_panels_covered()
    print(f"[14_program_nes_by_cell_state_viz] map coverage: all {n_panels} per-cell map "
          "panels carry a row on this geometry")

    n_checked = assert_one_fdr_family(rows)
    print(f"[14_program_nes_by_cell_state_viz] pooled-FDR family gate: {n_checked} cells "
          f"cross-checked between {NAMED_SETS_CSV} and {GSEA_ALL_CSV}, all agree")

    n_absent = int((~rows["tested"]).sum())
    print(f"[14_program_nes_by_cell_state_viz] {len(rows)} cells, {n_absent} untested "
          f"({'none' if n_absent == 0 else 'drawn as an explicit absence'})")

    width, height = 13.0, 17.0
    fig = build_figure(rows, order, fdr, width, height)

    def cell(set_id: str, state: str) -> pd.Series:
        hit = rows[(rows["pathway_id"] == set_id) & (rows["cell_state"] == state)]
        if hit.empty:
            raise ValueError(f"[14_program_viz] no plotted cell for {set_id} / {state}")
        return hit.iloc[0]

    hyp_t = cell("HALLMARK_HYPOXIA", "Treg")
    cro_t = cell("GOBP_CELLULAR_RESPONSE_TO_OXYGEN_LEVELS", "Treg")
    sting_tcon = cell("sting_specific_up", "Tcon")
    sting_treg = cell("sting_specific_up", "Treg")
    sting_cd8 = cell("sting_specific_up", "CD8")
    ifn_cd8 = cell("ifn_only_up", "CD8")
    path_sting = rows[rows["pathway_id"].isin([
        "GOBP_CGAS_STING_SIGNALING_PATHWAY",
        "REACTOME_STING_MEDIATED_INDUCTION_OF_HOST_IMMUNE_RESPONSES",
        "WP_STING_PATHWAY_IN_KAWASAKILIKE_DISEASE_AND_COVID19"]) & rows["tested"]]

    etreg_t = cell("eTreg_up", "Treg")
    hsr_t = cell("HSR_core", "Treg")
    hsr_tcon = cell("HSR_core", "Tcon")
    hsr_cd8 = cell("HSR_core", "CD8")

    big = baseline_band(baseline, "Treg", "130 to 150 genes")
    small = baseline_band(baseline, "Treg", "10 to 22 genes")

    save_overview(
        fig, STAGE, STEM,
        table=rows,
        finding=(
            "All four oxygen-named sets rise on the synovial-fluid side of the paired contrast "
            f"in Treg, and they spread across the panel: HALLMARK_HYPOXIA reaches NES "
            f"{float(hyp_t['nes']):.4f} at pooled FDR {fmt_p(float(hyp_t['padj_pooled']))} on "
            f"{int(hyp_t['genes_in_ranked_list'])} genes while "
            f"GOBP_CELLULAR_RESPONSE_TO_OXYGEN_LEVELS reaches {float(cro_t['nes']):.4f} at "
            f"{fmt_p(float(cro_t['padj_pooled']))} on {int(cro_t['genes_in_ranked_list'])}, so a "
            "reading taken from one of the four holds for that set alone. In the second group "
            "the two interferon sets carry the strongest cGAS-STING-family rows, ifn_only_up "
            f"reaching NES {float(ifn_cd8['nes']):.4f} at "
            f"{fmt_p(float(ifn_cd8['padj_pooled']))} in CD8. sting_specific_up clears pooled FDR "
            f"{fdr:g} in Tcon at NES {float(sting_tcon['nes']):.4f} on "
            f"{int(sting_tcon['genes_in_ranked_list'])} genes and sits above the threshold in "
            f"Treg ({float(sting_treg['nes']):.4f}, {fmt_p(float(sting_treg['padj_pooled']))}) "
            f"and CD8 ({float(sting_cd8['nes']):.4f}, {fmt_p(float(sting_cd8['padj_pooled']))}), "
            f"and the three pathway-database STING terms stay between NES "
            f"{float(path_sting['nes'].min()):.2f} and {float(path_sting['nes'].max()):.2f} on "
            f"{int(path_sting['genes_in_ranked_list'].min())} to "
            f"{int(path_sting['genes_in_ranked_list'].max())} genes. "
            f"The third group holds the strongest and the weakest rows of the figure: "
            f"eTreg_up, this compartment's own synovial-versus-blood contrast on another "
            f"cohort, reaches NES {float(etreg_t['nes']):.4f} at "
            f"{fmt_p(float(etreg_t['padj_pooled']))} in Treg, which marks what a set built to "
            f"separate exactly these two tissues reaches here, while the curated proteostasis "
            f"core HSR_core changes sign between the sorted populations at NES "
            f"{float(hsr_t['nes']):.4f} in Treg against {float(hsr_tcon['nes']):.4f} in Tcon "
            f"and {float(hsr_cd8['nes']):.4f} in CD8, clearing pooled FDR {fdr:g} in none of "
            f"the three."
        ),
        script=SCRIPT, fn="build_figure",
        config_kv=(f"thresholds.gsea_fdr = {fdr}; gsea_min_size = {PARAMS.gsea_min_size}; "
                   f"gsea_max_size = {PARAMS.gsea_max_size}; row order = descending "
                   f"{ORDER_BY_POPULATION} NES within a group"),
        input=(f"03_results/{STAGE}/{NAMED_SETS_CSV}, 03_results/{STAGE}/{GSEA_ALL_CSV}, "
               f"03_results/{STAGE}/{BASELINE_CSV}"),
        how_to_read=(
            "The companion to arm_nes_by_cell_state under "
            "03_results/14_unbiased_enrichment/, drawn to the same geometry so the two can be "
            "laid side by side. One dot per gene set and cell state, at the confirmatory tier: "
            "donor-level pseudobulk within frozen sort labels, limma-voom moderated t, then "
            "pre-ranked fgsea. Rows are gene sets in three labelled groups, ordered inside a group "
            f"by descending {ORDER_BY_POPULATION} NES, and inside a row the three cell states are "
            "offset vertically and coloured, each with its own annotation line. The x position is "
            "the normalised enrichment score for synovial fluid over paired blood, with a "
            f"vertical rule at zero. A filled dot clears the config FDR threshold of {fdr:g} and "
            "an open dot sits above it. A cell reading not tested had fewer genes in that "
            "population's ranked list than gsea_min_size, so the sweep never scored it and the "
            "cell records an untested set. "
            "Read every score against the gene count beside it. In this sweep, size alone moves "
            f"the odds a long way: in {ORDER_BY_POPULATION} a set of "
            f"{big['band_label']} clears pooled FDR {fdr:g} in "
            f"{100 * float(big['frac_pooled_significant']):.1f}% of the "
            f"{int(big['n_sets_tested']):,} such sets tested, while a set of "
            f"{small['band_label']} clears it in "
            f"{100 * float(small['frac_pooled_significant']):.1f}% of "
            f"{int(small['n_sets_tested']):,}. The four oxygen-named sets are four differently "
            "built sets for one named biology, and their scores differ, so read the group as a "
            "family and each score as that set's own. The cGAS-STING group carries a positive "
            "result in Tcon, so read that group set by set too. "
            "The third group carries the remaining panels of the per-cell score maps, so every "
            "lens a reader meets on a map has a row on this geometry, the three mouse arms on "
            "the sibling panel arm_nes_by_cell_state and the rest here. "
            "The complete six-member cGAS-STING family of this sweep, including the two "
            "regulation-of terms whose signs disagree, is drawn in named_sets_in_sweep under the "
            "same stage. Temperature and hypoxia are both imposed by the inflamed joint and stay "
            "entangled in cross-sectional human data, so these rows give what the niche "
            "contrast contains and leave open which of the two drives it. A score of this kind "
            "has no interval, so none is drawn, and the reading stays correlative."
        ),
        config=FIG_CFG, width=width, height=height)

    print(f"[14_program_nes_by_cell_state_viz] wrote {STEM} from {len(rows)} plotted rows "
          f"({len(order)} sets x {len(CELL_STATES)} cell states)")


if __name__ == "__main__":
    main()
