#!/usr/bin/env python
"""
09_heat_hypoxia_viz.py — VIZ ONLY. One bounded question, one panel.
===================================================================
The inflamed synovial niche is hypoxic as well as inflamed, so one bounded
question can be asked of the mouse 39 °C-derived enrichment: is it reducible to
the set's own HALLMARK_HYPOXIA-overlap gene content? That is a membership
question and it is answered by deleting those genes from the mouse sets and
re-running the same fgsea engine. It is not a question about temperature, and it
is not a question about whether hypoxia is a confound or a co-exposure — those
are not separable in cross-sectional human data, and nothing here licenses a
statement about either.

  heat_purge_nes_paired — how much of the enrichment survives the purge

That panel is the stage's published overview. The per-cell co-localization
correlations and the model-assigned leading-edge composition remain compute
resources in this stage's tables, read by the reactive review notebook and by
the cross-dataset layer; neither is drawn here. Whole-arm composition is carried
by 11_heat_decomposition/heatdecomp_arm_coverage.

Computes no statistic. Every NES, FDR and set size is read verbatim from a
committed CSV; the only derived quantities are plain counts over frozen gene
lists and differences of those counts.

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
from matplotlib.colors import to_rgba  # noqa: E402
from matplotlib.legend_handler import HandlerTuple  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

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

# The one contrasting edge in this module comes from the config Okabe-Ito
# palette, so the only literal here is the config key, not a hex string.
_OKABE = (FIG_CFG.get("colors", {}) or {}).get("okabe_ito", {}) or {}
EDGE_CONTRAST = _OKABE["black"]

_F = FIG_CFG.get("figures", {}) or {}
ANNOT_SIZE = float(_F["axis_text_size"])
LEGEND_SIZE = float(_F["legend_text_size"])
FDR = float(PARAMS.gsea_fdr)
# One translucency for every paired/overlapping marker in this module, so
# coincident markers read through each other instead of hiding one another.
MARKER_ALPHA = 0.65


def fmt_fdr(p: float) -> str:
    """Render an FDR for an in-figure label: fixed below 3 decimals, else scientific."""
    if pd.isna(p):
        return "FDR n/a"
    return f"FDR {p:.3f}" if p >= 0.001 else f"FDR {p:.0e}"


# ===========================================================================
# 1. Paired full-versus-purged NES — does the enrichment survive the purge?
# ===========================================================================
def _count_symbols(path: Path) -> int:
    """Plain line count of a frozen newline-delimited symbol file. No statistic."""
    if not path.exists():
        raise FileNotFoundError(f"[09_viz] frozen signature file missing: {path}")
    return len({ln.strip() for ln in path.read_text().splitlines() if ln.strip()})


def purge_paired_table() -> pd.DataFrame:
    """Marshal the per-population full/purged fgsea rows into one row per arm.

    `gene_purge_nes_comparison.csv` carries only the up arm, so the down arm is read
    from the same `gsea_{full,purged}_*.csv` files that table was built from. No
    number is recomputed — NES, p, FDR and set_size are all read as written.

    Two membership counts travel with every row so the caption can quote the
    nominal-versus-effective distinction from its own same-stem table rather than
    from a neighbour:

      `n_removed_nominal`  — genes the purge deletes from the frozen set FILE
                             (a line-count difference over two committed files);
      `n_removed_testable` — how many of those were in that population's ranked
                             list at all, i.e. the number the purge could actually
                             move (`set_size_full - set_size_purged`).

    They differ, and quoting the nominal count as "genes removed" overstates the
    purge. `delta_nes` is likewise carried so the caption never has to subtract
    two columns in prose.
    """
    tdir = PATHS.tables(STAGE)
    n_nominal = {arm: _count_symbols(tdir / "_signatures_full" / f"{ARM_SET[arm]}.txt")
                 for arm in ARM_ORDER}
    nominal = {
        arm: n_nominal[arm]
        - _count_symbols(tdir / "_signatures_purged" / f"{ARM_SET[arm]}.txt")
        for arm in ARM_ORDER
    }
    rows = []
    for pop, tag in POP_TAG.items():
        full = pd.read_csv(tdir / f"gsea_full_{tag}.csv").set_index("pathway_id")
        purged = pd.read_csv(tdir / f"gsea_purged_{tag}.csv").set_index("pathway_id")
        for arm in ARM_ORDER:
            sid = ARM_SET[arm]
            if sid not in full.index or sid not in purged.index:
                continue
            size_full = int(full.loc[sid, "set_size"])
            size_purged = int(purged.loc[sid, "set_size"])
            rows.append({
                "population": pop,
                "arm": arm,
                "signature": sid,
                "contrast": "SF_vs_PB",
                "nes_full": float(full.loc[sid, "nes"]),
                "padj_full": float(full.loc[sid, "padj"]),
                "set_size_full": size_full,
                "nes_purged": float(purged.loc[sid, "nes"]),
                "padj_purged": float(purged.loc[sid, "padj"]),
                "set_size_purged": size_purged,
                "delta_nes": float(purged.loc[sid, "nes"]) - float(full.loc[sid, "nes"]),
                "n_nominal": n_nominal[arm],
                "n_removed_nominal": nominal[arm],
                "n_removed_testable": size_full - size_purged,
                "evidence_tier": "primary_pseudobulk",
            })
    return pd.DataFrame(rows)


def _pair_handles(colour: str) -> tuple:
    """Legend handle for one arm: the full-set diamond beside the purged circle.

    Both are drawn at the plotted translucency and with no contrasting edge, so
    the key shows the default (FDR-not-passing) marker state and the outline key
    below it is the only thing that carries significance.
    """
    face = to_rgba(colour, MARKER_ALPHA)
    return (Line2D([0], [0], linestyle="none", marker="D", markerfacecolor=face,
                   markeredgecolor=face, markersize=15),
            Line2D([0], [0], linestyle="none", marker="o", markerfacecolor=face,
                   markeredgecolor=face, markersize=10))


def plot_purge_paired(df: pd.DataFrame):
    fig, ax = plt.subplots()
    n = len(df)
    ylabels = []
    for i, (_, r) in enumerate(df.iterrows()):
        y = n - 1 - i
        col = ARM_COL[r["arm"]]
        ax.plot([r["nes_full"], r["nes_purged"]], [y, y], color=col, lw=2.6, zorder=2,
                solid_capstyle="round")
        # Every marker is filled and translucent, larger full-set diamond first and
        # smaller purged circle on top: where the purge removes nothing and both land
        # on the same NES (the whole down arm) the pair reads as a darker circle
        # inside a lighter diamond rather than one marker hiding the other.
        for nes, padj, marker, size in ((r["nes_full"], r["padj_full"], "D", 300),
                                        (r["nes_purged"], r["padj_purged"], "o", 120)):
            passes = padj < FDR
            ax.scatter(nes, y, marker=marker, s=size, zorder=3,
                       facecolors=to_rgba(col, MARKER_ALPHA),
                       edgecolors=EDGE_CONTRAST if passes else to_rgba(col, MARKER_ALPHA),
                       linewidths=2.2 if passes else 0.0)
        # Effective set size against this arm's OWN nominal size, then what the
        # purge cost. Where the purge takes nothing there is no NES to quote, so
        # the row says that instead of printing a rounded zero.
        if int(r["n_removed_testable"]) == 0:
            cost = "purge removes no gene"
        else:
            cost = f"ΔNES {r['delta_nes']:+.3f}"
        ax.text(2.95, y,
                f"n {r['set_size_full']}→{r['set_size_purged']} of {r['n_nominal']}  ·  "
                f"{cost}  ·  {fmt_fdr(r['padj_purged'])}",
                va="center", ha="left", fontsize=ANNOT_SIZE)
        ylabels.append(f"{r['population']} · {r['arm']}")

    ax.axvline(0, color="black", lw=1)
    ax.set_yticks(range(n))
    ax.set_yticklabels(list(reversed(ylabels)))
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_xlim(0, 4.6)
    ax.set_xticks([0, 1, 2, 3])
    ax.set_xlabel("fgsea NES (synovial fluid vs paired blood)")
    # The title states the bounded answer in membership terms. It may not say
    # anything about temperature, and it may not call hypoxia a confound or a
    # co-exposure: the purge tests gene content and licenses nothing else.
    ax.set_title("The synovial-fluid enrichment is not reducible to the set's\n"
                 "HALLMARK_HYPOXIA-overlap gene content")
    # Q -> A marker, placed under the axis so it cannot collide with a marker:
    # this panel answers, and it says on its face what it answers, at what tier,
    # and where the answer stops. `bbox_inches="tight"` keeps below-axis text.
    up = df[df["arm"].eq("up")]
    ax.text(0.0, -0.135,
            "ANSWER — confirmatory tier (donor-level pseudobulk, limma-voom → fgsea).\n"
            f"Deleting the {int(up['n_removed_nominal'].iloc[0])} HALLMARK_HYPOXIA-overlap genes "
            f"takes {int(up['n_removed_testable'].min())}–{int(up['n_removed_testable'].max())} "
            "testable genes out of the up arm and costs "
            f"{abs(float(up['delta_nes'].max())):.3f}–{abs(float(up['delta_nes'].min())):.3f} NES;\n"
            "all three up arms stay significant, so the enrichment is not reducible to that gene "
            "content. That is the whole claim. It says nothing\n"
            "about temperature, and nothing about whether hypoxia is a confound or a co-exposure — "
            "those are not separable in cross-sectional human data.",
            transform=ax.transAxes, ha="left", va="top", fontsize=ANNOT_SIZE)
    # Three keys: each arm key carries its colour AND its two shapes, and the one
    # remaining key explains the outline. No open/filled convention anywhere.
    handles = [_pair_handles(ARM_COL["up"]), _pair_handles(ARM_COL["down"]),
               Line2D([0], [0], linestyle="none", marker="o",
                      markerfacecolor=to_rgba("grey", MARKER_ALPHA),
                      markeredgecolor=EDGE_CONTRAST, markeredgewidth=2.2, markersize=11)]
    # The down-arm key names where that arm is significant, read from the table
    # rather than asserted, so the panel cannot outlive "the up arm is the only
    # informative arm" by carrying it in a legend.
    down_sig = df[df["arm"].eq("down") & df["padj_full"].lt(FDR)]["population"].tolist()
    down_where = (", ".join(down_sig) + " only") if down_sig else "no population"
    labels = ["WT_heat up arm — diamond = full mouse set, circle = hypoxia-purged set",
              "WT_heat down arm — same pair, and the purge removes no gene, so they coincide; "
              f"this arm is not silent, reaching FDR below {FDR} in {down_where}",
              f"dark outline = FDR below {FDR} (every FDR is printed at right)"]
    ax.legend(handles, labels, handler_map={tuple: HandlerTuple(ndivide=None, pad=0.7)},
              loc="upper center", bbox_to_anchor=(0.5, -0.30), ncol=1, frameon=False,
              fontsize=LEGEND_SIZE, handlelength=3.2, handletextpad=1.0, labelspacing=0.7)
    fig.tight_layout()
    return fig


# ===========================================================================
def main() -> None:
    set_paper_style(config=FIG_CFG)
    purge_figures(STAGE, "heat_", overview=True, config=FIG_CFG)

    paired = purge_paired_table()
    fig = plot_purge_paired(paired)
    # Every NES, delta and gene count below is read from `paired` in this run. The mouse arms
    # are re-derived upstream in mouse_anchor, so a typed value here would go stale.
    up = paired[paired["arm"].eq("up")]
    up_moves = ", ".join(f"{float(r['nes_full']):.4f} to {float(r['nes_purged']):.4f} in "
                         f"{r['population']}" for _, r in up.iterrows())
    save_overview(
        fig, STAGE, "heat_purge_nes_paired",
        table=round_numeric_cols(paired),
        finding=(f"Deleting the {int(up['n_removed_nominal'].max())} HALLMARK_HYPOXIA-overlap "
                 f"genes from the mouse 39 °C-derived up-set takes "
                 f"{int(up['n_removed_testable'].min())} to "
                 f"{int(up['n_removed_testable'].max())} testable genes out of the arm and costs "
                 f"{up['delta_nes'].abs().min():.3f} to {up['delta_nes'].abs().max():.3f} NES — "
                 f"{up_moves} — leaving all three significant. The synovial-fluid enrichment "
                 f"therefore survives the removal of its HALLMARK_HYPOXIA-overlap gene content. "
                 f"This is a statement about gene content. Temperature is untested here, and "
                 f"cross-sectional human data leave hypoxia's status as confound or co-exposure "
                 f"undetermined."),
        script=SCRIPT, fn="plot_purge_paired",
        config_kv=(f"thresholds.gsea_fdr={FDR}; gsea_min_size={PARAMS.gsea_min_size}; "
                   f"gsea_nperm={PARAMS.gsea_nperm}"),
        input="03_results/09_heat_hypoxia/tables/gsea_{full,purged}_{treg,tcon,cd8}.csv",
        how_to_read=(
            "This is the confirmatory tier: donor-level pseudobulk within frozen sort labels, "
            "limma-voom then fgsea. Positive NES points toward synovial fluid. Each row pairs "
            "the full set (large diamond) with its purged form (small circle); the connecting "
            "bar is the NES cost. Warm brown is the up arm and cool blue the down arm. A dark "
            f"outline marks FDR below {FDR}. Right-hand text reports effective and nominal set "
            f"sizes, the NES cost, and purged FDR. Two gene counts differ and both are given: "
            f"{int(up['n_removed_nominal'].max())} genes come out of the frozen set, of which "
            f"{int(up['n_removed_testable'].min())} to {int(up['n_removed_testable'].max())} "
            f"were present in a ranked list. The Tcon down arm stays significant at the up "
            f"arm's sign. This licenses a membership statement. Correlative."),
        config=FIG_CFG, height=7.6,
    )
    plt.close(fig)
    print("[09_heat_hypoxia_viz] wrote 1 overview (purge pairing)")


if __name__ == "__main__":
    main()
