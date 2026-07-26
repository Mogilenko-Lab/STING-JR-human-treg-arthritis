#!/usr/bin/env python
"""
09_heat_hypoxia_viz.py — VIZ ONLY. The heat-versus-hypoxia narrative in three beats.
===================================================================================
Synovial fluid is a hypoxic niche as well as an inflamed one, so the standing
objection to reading the mouse 39 °C enrichment as thermal is that it is really
hypoxia. Stage 09 answered by purging `HALLMARK_HYPOXIA` genes from the mouse
sets and re-running the same fgsea engine. These figures walk that answer:

  1. heat_purge_nes_paired        — how much of the enrichment survives the purge
  2. heat_hypoxia_colocalization — do heat-high and hypoxia-high mark the same cells
  3. heat_leadingedge_composition — what the enriching genes actually are

Reads only committed stage-09 tables and computes no statistic: every NES, FDR,
correlation and gene tally plotted here is read verbatim from a CSV written by
09_heat_hypoxia.py.

Run from the compartment root, AFTER 09_heat_hypoxia.py:
  python 02_analysis/scripts/09_heat_hypoxia_viz.py
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "02_analysis"))
os.chdir(ROOT)

from config import PARAMS, PATHS  # noqa: E402
from helpers.figure_style import (  # noqa: E402
    FIG_CFG,
    purge_figures,
    round_numeric_cols,
    save_overview,
    set_paper_style,
)

STAGE = "09_heat_hypoxia"
SCRIPT = "02_analysis/scripts/09_heat_hypoxia_viz.py"
POP_TAG = {"Treg": "treg", "Tcon": "tcon", "CD8": "cd8"}
ARM_ORDER = ["up", "down"]
ARM_SET = {"up": "WT_heat_up", "down": "WT_heat_down"}

# Declared palette constants -------------------------------------------------
# Mouse-arm diverging cue, IDENTICAL to 05_score_signatures_viz.R: heat-up =
# warm brown, heat-down = cool blue. Keyed by arm so the mapping can never come
# out of a positional vector in the wrong order.
ARM_COL = {"up": "#A6611A", "down": "#2166AC"}

# Leading-edge program colours come from the config Okabe-Ito palette, so the
# only literal here is the config key, not a hex string.
_OKABE = (FIG_CFG.get("colors", {}) or {}).get("okabe_ito", {}) or {}
LE_PROGRAMS = [
    ("heat_shock_proteostasis", "heat shock / proteostasis", _OKABE["vermillion"]),
    ("hypoxia_HIF", "hypoxia / HIF overlap", _OKABE["sky_blue"]),
    ("immediate_early_stress", "immediate-early stress", _OKABE["orange"]),
    ("effector_activation", "effector / activation", _OKABE["bluish_green"]),
    ("other", "other", _OKABE["reddish_purple"]),
]

_F = FIG_CFG.get("figures", {}) or {}
ANNOT_SIZE = float(_F["axis_text_size"])
LEGEND_SIZE = float(_F["legend_text_size"])
FDR = float(PARAMS.gsea_fdr)


def fmt_fdr(p: float) -> str:
    """Render an FDR for an in-figure label: fixed below 3 decimals, else scientific."""
    if pd.isna(p):
        return "FDR n/a"
    return f"FDR {p:.3f}" if p >= 0.001 else f"FDR {p:.0e}"


# ===========================================================================
# 1. Paired full-versus-purged NES — does the enrichment survive the purge?
# ===========================================================================
def purge_paired_table() -> pd.DataFrame:
    """Marshal the per-population full/purged fgsea rows into one row per arm.

    `gene_purge_nes_comparison.csv` carries only the up arm, so the down arm is read
    from the same `gsea_{full,purged}_*.csv` files that table was built from. No
    number is recomputed — NES, p, FDR and set_size are all read as written.
    """
    tdir = PATHS.tables(STAGE)
    rows = []
    for pop, tag in POP_TAG.items():
        full = pd.read_csv(tdir / f"gsea_full_{tag}.csv").set_index("pathway_id")
        purged = pd.read_csv(tdir / f"gsea_purged_{tag}.csv").set_index("pathway_id")
        for arm in ARM_ORDER:
            sid = ARM_SET[arm]
            if sid not in full.index or sid not in purged.index:
                continue
            rows.append({
                "population": pop,
                "arm": arm,
                "signature": sid,
                "contrast": "SF_vs_PB",
                "nes_full": float(full.loc[sid, "nes"]),
                "padj_full": float(full.loc[sid, "padj"]),
                "set_size_full": int(full.loc[sid, "set_size"]),
                "nes_purged": float(purged.loc[sid, "nes"]),
                "padj_purged": float(purged.loc[sid, "padj"]),
                "set_size_purged": int(purged.loc[sid, "set_size"]),
                "evidence_tier": "primary_pseudobulk",
            })
    return pd.DataFrame(rows)


def plot_purge_paired(df: pd.DataFrame):
    fig, ax = plt.subplots()
    n = len(df)
    ylabels = []
    for i, (_, r) in enumerate(df.iterrows()):
        y = n - 1 - i
        col = ARM_COL[r["arm"]]
        ax.plot([r["nes_full"], r["nes_purged"]], [y, y], color=col, lw=2.6, zorder=2,
                solid_capstyle="round")
        # Sizes differ so the two markers stay separately readable when the purge
        # changes nothing and they land on the same NES (the whole down arm).
        for nes, padj, marker, size in ((r["nes_full"], r["padj_full"], "D", 210),
                                        (r["nes_purged"], r["padj_purged"], "o", 95)):
            ax.scatter(nes, y, marker=marker, s=size, zorder=3, linewidths=1.8,
                       edgecolors=col, facecolors=col if padj < FDR else "white")
        ax.text(3.0, y, f"n {r['set_size_full']}→{r['set_size_purged']}  ·  "
                        f"{fmt_fdr(r['padj_purged'])}",
                va="center", ha="left", fontsize=ANNOT_SIZE)
        ylabels.append(f"{r['population']} · {r['arm']}")

    ax.axvline(0, color="black", lw=1)
    ax.set_yticks(range(n))
    ax.set_yticklabels(list(reversed(ylabels)))
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_xlim(0, 4.25)
    ax.set_xticks([0, 1, 2, 3])
    ax.set_xlabel("fgsea NES (synovial fluid vs paired blood)")
    ax.set_title("Purging hypoxia genes barely moves the mouse heat enrichment")
    handles = [
        Line2D([0], [0], color=ARM_COL["up"], lw=2.6, label="WT_heat up (warm)"),
        Line2D([0], [0], color=ARM_COL["down"], lw=2.6, label="WT_heat down (cool)"),
        Line2D([0], [0], marker="D", color="w", markerfacecolor="grey",
               markeredgecolor="grey", markersize=13, label="full mouse set"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="grey",
               markeredgecolor="grey", markersize=9, label="hypoxia-purged set"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="white",
               markeredgecolor="grey", markersize=9, markeredgewidth=1.8,
               label=f"open marker: FDR ≥ {FDR}"),
    ]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.14),
              ncol=3, frameon=False, fontsize=LEGEND_SIZE)
    fig.tight_layout()
    return fig


# ===========================================================================
# 2. Per-cell co-localization — are these the same cells?
# ===========================================================================
def colocalization_table() -> pd.DataFrame:
    df = pd.read_csv(PATHS.tables(STAGE) / "heat_hypoxia_colocalization.csv")
    cell = df[df["level"] == "cell"].copy()
    cell["population"] = pd.Categorical(cell["population"], list(POP_TAG), ordered=True)
    return cell.sort_values(["population", "method"]).reset_index(drop=True)


def plot_colocalization(df: pd.DataFrame):
    fig, ax = plt.subplots()
    methods = ["spearman", "pearson"]
    m_col = {"spearman": _OKABE["blue"], "pearson": _OKABE["sky_blue"]}
    width = 0.34
    pops = list(POP_TAG)
    for i, method in enumerate(methods):
        sub = df[df["method"] == method].set_index("population").reindex(pops)
        xpos = [x + (i - 0.5) * width for x in range(len(pops))]
        ax.bar(xpos, sub["r"].astype(float), width=width, color=m_col[method],
               label=method.capitalize())
        for xp, r in zip(xpos, sub["r"].astype(float)):
            ax.text(xp, r + 0.025, f"{r:.2f}", ha="center", va="bottom", fontsize=ANNOT_SIZE)

    ns = df[df["method"] == "spearman"].set_index("population").reindex(pops)["n"]
    ax.set_xticks(range(len(pops)))
    ax.set_xticklabels([f"{p}\nn = {int(v):,} cells" for p, v in zip(pops, ns)])
    ax.axhline(0, color="black", lw=1)
    ax.set_ylim(-0.05, 1.0)
    ax.set_ylabel("Correlation of heat and hypoxia score, within SF cells")
    ax.set_title("Heat-high and hypoxia-high are largely different cells")
    ax.legend(frameon=False, fontsize=LEGEND_SIZE, loc="upper right")
    ax.text(0.02, 0.94,
            "Full –0.05 to 1 scale, so the height of a bar is the whole story:\n"
            "co-localization is weak, and the niche's two stresses land on\n"
            "different cells. Donor-level SF means (n = 6–7) are unpowered and\n"
            "reported in the stage table only.",
            transform=ax.transAxes, ha="left", va="top", fontsize=ANNOT_SIZE)
    fig.tight_layout()
    return fig


# ===========================================================================
# 3. Leading-edge composition — what are the enriching genes?
# ===========================================================================
def leadingedge_table() -> pd.DataFrame:
    df = pd.read_csv(PATHS.tables(STAGE) / "leadingedge_composition.csv")
    keep = ["population", "signature", "n_leading_edge"]
    keep += [f"n_{k}" for k, _, _ in LE_PROGRAMS]
    keep += [f"frac_{k}" for k, _, _ in LE_PROGRAMS]
    keep += ["n_unclassified", "taxonomy_source", "evidence_tier"]
    df["population"] = pd.Categorical(df["population"], list(POP_TAG), ordered=True)
    return df.sort_values("population")[keep].reset_index(drop=True)


def plot_leadingedge(df: pd.DataFrame):
    fig, ax = plt.subplots()
    pops = list(POP_TAG)
    sub = df.set_index("population").reindex(pops)
    for i, pop in enumerate(pops):
        y = len(pops) - 1 - i
        left = 0.0
        for key, _, colour in LE_PROGRAMS:
            frac = float(sub.loc[pop, f"frac_{key}"])
            count = int(sub.loc[pop, f"n_{key}"])
            ax.barh(y, frac, left=left, height=0.62, color=colour, edgecolor="white",
                    linewidth=1.2)
            # Print the count inside the segment whenever it fits; the trace-level
            # heat-shock segment is the point of the figure, so its 2-3 genes must
            # stay legible rather than dropping out on a fixed width threshold.
            txt = str(count)
            if frac >= 0.025 * len(txt):
                ax.text(left + frac / 2, y, txt, ha="center", va="center",
                        fontsize=ANNOT_SIZE, color="black")
            left += frac
        ax.text(1.02, y, f"{int(sub.loc[pop, 'n_leading_edge'])} genes",
                va="center", ha="left", fontsize=ANNOT_SIZE)

    ax.set_yticks(range(len(pops)))
    ax.set_yticklabels(list(reversed(pops)))
    ax.set_xlim(0, 1.0)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xlabel("Fraction of the WT_heat up leading edge")
    ax.set_title("The surviving signal is carried mostly by activation genes")
    handles = [Patch(facecolor=c, label=lab) for _, lab, c in LE_PROGRAMS]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.2),
              ncol=3, frameon=False, fontsize=LEGEND_SIZE)
    fig.tight_layout()
    return fig


# ===========================================================================
def main() -> None:
    set_paper_style(config=FIG_CFG)
    purge_figures(STAGE, "heat_", overview=True, config=FIG_CFG)

    paired = purge_paired_table()
    fig = plot_purge_paired(paired)
    save_overview(
        fig, STAGE, "heat_purge_nes_paired",
        table=round_numeric_cols(paired),
        finding=("Removing every HALLMARK_HYPOXIA gene from the mouse 39 °C up-set costs "
                 "only 0.14 to 0.19 NES and leaves the synovial-fluid enrichment strong and "
                 "significant in all three sorted populations, so hypoxia does not explain it."),
        script=SCRIPT, fn="plot_purge_paired",
        config_kv=(f"thresholds.gsea_fdr={FDR}; gsea_min_size={PARAMS.gsea_min_size}; "
                   f"gsea_nperm={PARAMS.gsea_nperm}"),
        input="03_results/09_heat_hypoxia/tables/gsea_{full,purged}_{treg,tcon,cd8}.csv",
        how_to_read=(
            "One row per population and mouse arm; x is fgsea NES, positive = enriched toward "
            "the synovial-fluid end of the paired SF-vs-blood ranking. Each row pairs the full "
            "mouse set (diamond) with the hypoxia-purged set (circle); the connecting bar is "
            "what the purge cost. Warm brown = up arm, cool blue = down arm. A filled marker "
            f"means FDR < {FDR}, an open marker FDR at or above it — no other significance "
            "glyph is used. The right-hand text gives testable set size before and after the "
            "purge and the purged FDR; the down arm loses no genes because the hypoxia overlap "
            "sits entirely in the up arm. Primary donor-pseudobulk tier; correlative."),
        config=FIG_CFG, height=7.0,
    )
    plt.close(fig)

    coloc = colocalization_table()
    fig = plot_colocalization(coloc)
    save_overview(
        fig, STAGE, "heat_hypoxia_colocalization",
        table=round_numeric_cols(coloc),
        finding=("Within synovial-fluid cells the mouse heat score and the hypoxia score "
                 "correlate only weakly (Spearman 0.08 to 0.20), so the niche's thermal and "
                 "hypoxic readouts are carried by largely different cells rather than one "
                 "shared stress state."),
        script=SCRIPT, fn="plot_colocalization",
        config_kv="level=cell; tissue=synovial_fluid; evidence_tier=secondary_percell",
        input="03_results/09_heat_hypoxia/tables/heat_hypoxia_colocalization.csv",
        how_to_read=(
            "Bars are the within-SF, cell-level correlation between the per-cell WT_heat_up "
            "and HALLMARK_HYPOXIA AUCell scores, Spearman (dark) beside Pearson (light), with "
            "the cell count under each population. The y-axis deliberately runs the full "
            "-0.05 to 1 range: read the shortness of the bars, not their rank order. Positive "
            "r means a heat-high cell tends to be hypoxia-high. Donor-level SF means are "
            "unpowered at 6 to 7 donors and are left in the stage table rather than drawn. "
            "This is a secondary per-cell diagnostic of where the two scores sit, never "
            "pooled with the pseudobulk NES and never read as directional evidence."),
        config=FIG_CFG,
    )
    plt.close(fig)

    le = leadingedge_table()
    fig = plot_leadingedge(le)
    save_overview(
        fig, STAGE, "heat_leadingedge_composition",
        table=round_numeric_cols(le),
        finding=("Half the WT_heat up leading edge in synovial-fluid T cells is effector and "
                 "activation genes, with a hypoxia-overlap minority and only two or three "
                 "classic heat-shock genes, so surviving the hypoxia purge does not make the "
                 "enrichment thermally specific."),
        script=SCRIPT, fn="plot_leadingedge",
        config_kv="taxonomy=00_data/references/heat_leadingedge_taxonomy; evidence_tier=secondary_exploratory",
        input="03_results/09_heat_hypoxia/tables/leadingedge_composition.csv",
        how_to_read=(
            "One stacked bar per population, spanning that population's full WT_heat up "
            "leading edge; segment width is the fraction of leading-edge genes in each "
            "program and the number inside a segment is its gene count (printed where the "
            "segment is wide enough to hold it; every count is in the source table). Program "
            "assignment comes from a frozen external-model gene taxonomy, not from this run. "
            "The narrow heat-shock segment is the point, and it is why an activation-free "
            "proteostasis lens was built next. Exploratory secondary tier, never pooled with "
            "the pseudobulk NES."),
        config=FIG_CFG, height=7.0,
    )
    plt.close(fig)
    print("[09_heat_hypoxia_viz] wrote 3 overviews (purge pairing, co-localization, leading edge)")


if __name__ == "__main__":
    main()
