#!/usr/bin/env python
"""
10_hsr_lens_viz.py — VIZ ONLY. The curated HSR lens, drawn as a trend.
==========================================================================
Heat-shock readouts in inflamed tissue are confounded by activation, so a
curated, activation-free proteostasis lens was built to ask what is left. The
answer is a sign flip at trend level — HSR core points toward synovial fluid in
Treg and away from it in Tcon and CD8, while no population clears FDR 0.05.
Both halves of that sentence have to be visible in the figure, so nothing here
carries a star or any other glyph implying significance; the FDRs are printed.

  1. hsr_core_running_sum — the sign flip read off the ranked lists themselves

Reads committed stage-10 tables and computes no statistic.

Run from the compartment root, AFTER 10_hsr_lens.py:
  python 02_analysis/scripts/10_hsr_lens_viz.py
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "02_analysis"))
os.chdir(ROOT)

from config import PARAMS, PATHS, POPULATION_COLORS  # noqa: E402
from helpers.figure_style import (  # noqa: E402
    FIG_CFG,
    purge_figures,
    round_numeric_cols,
    save_overview,
    set_paper_style,
)

STAGE = "10_hsr_lens"
SCRIPT = "02_analysis/scripts/10_hsr_lens_viz.py"
POP_TAG = {"Treg": "treg", "Tcon": "tcon", "CD8": "cd8"}
POP_ORDER = list(POP_TAG)

# Declared palette constants, all resolved from the config so no hex literal appears in
# this script. The one population palette comes from `colors.populations`, so a population
# keeps one colour across the compartment and no figure declares its own assignment.
POP_COL = POPULATION_COLORS

_F = FIG_CFG.get("figures", {}) or {}
ANNOT_SIZE = float(_F["axis_text_size"])
LEGEND_SIZE = float(_F["legend_text_size"])
RS_HEIGHTS = [float(h) for h in _F["running_sum_heights"]]
RS_YLIM = [float(v) for v in _F["running_sum_ylim"]]
FDR = float(PARAMS.gsea_fdr)

# The canvas every running sum in this compartment is drawn on, so this panel and the
# runsum_<set> family are one geometry and one scale.
RS_W, RS_H = 11.0, 8.5


def fmt_fdr(p: float) -> str:
    if pd.isna(p):
        return "FDR n/a"
    return f"FDR {p:.3f}" if p >= 0.001 else f"FDR {p:.0e}"


# ===========================================================================
# HSR_core running sum — the sign flip, read off the ranked lists
# ===========================================================================
def _bool_col(s: pd.Series) -> pd.Series:
    """Robust to R writing a real bool or the strings TRUE/FALSE."""
    if s.dtype == bool:
        return s
    return s.astype(str).str.strip().str.upper().isin({"TRUE", "T", "1"})


def running_sum_traces() -> tuple[dict[str, pd.DataFrame], pd.DataFrame]:
    """Load the committed HSR_core running-sum traces plus their NES summary rows."""
    tdir = PATHS.tables(STAGE)
    nes = pd.read_csv(tdir / "hsr_lens_nes.csv")
    nes = nes[nes["signature"] == "HSR_core"].set_index("population")
    nominal_path = tdir / "_signatures_hsr" / "HSR_core.txt"
    if not nominal_path.exists():
        raise FileNotFoundError(f"[10_hsr_lens_viz] missing stage-local HSR core: {nominal_path}")
    n_nominal = len({g.strip() for g in nominal_path.read_text().splitlines() if g.strip()})
    traces, rows = {}, []
    for pop, tag in POP_TAG.items():
        path = tdir / f"runsum_interactive_hsr_gsea_{tag}_HSR_core.csv"
        if not path.exists():
            print(f"[10_hsr_lens_viz] {pop}: no HSR_core trace at {path} — skipping")
            continue
        tr = pd.read_csv(path, usecols=["rank", "stat", "running_es", "hit"])
        tr["hit"] = _bool_col(tr["hit"])
        # A fraction, so three rankings of different length share one axis.
        tr["rank_fraction"] = tr["rank"] / len(tr)
        traces[pop] = tr
        r = nes.loc[pop]
        n_effective = int(r["set_size"])
        testability = (
            "testable" if n_effective >= 15
            else "underpowered_reported" if n_effective >= 5
            else "untestable"
        )
        rows.append({
            "population": pop,
            "signature": "HSR_core",
            "contrast": "SF_vs_PB",
            "nes": float(r["nes"]),
            "pvalue": float(r["pvalue"]),
            "padj": float(r["padj"]),
            "set_size": n_effective,
            "n_nominal": n_nominal,
            "testability": testability,
            "n_ranked_genes": int(len(tr)),
            "evidence_tier": "secondary_annotation",
        })
    return traces, pd.DataFrame(rows)


def plot_running_sum(traces: dict[str, pd.DataFrame], summary: pd.DataFrame):
    # The three-panel running sum this compartment draws every other one in: the enrichment
    # traces, the gene-hit rug with one named row per population, and the three ranked
    # metrics, in the config's running_sum_heights proportions on one fractional-rank axis.
    fig, (ax, rug, met) = plt.subplots(
        3, 1, sharex=True, layout="constrained", height_ratios=RS_HEIGHTS)

    pops = [p for p in POP_ORDER if p in traces]
    for pop in pops:
        tr = traces[pop]
        row = summary[summary["population"] == pop].iloc[0]
        ax.plot(tr["rank_fraction"], tr["running_es"], color=POP_COL[pop], lw=2.0,
                label=(f"{pop}   NES {row['nes']:+.2f}, {fmt_fdr(row['padj'])}, "
                       f"{int(row['set_size'])} of {int(row['n_nominal'])} genes "
                       f"({row['testability']})"))
    ax.axhline(0, color="black", lw=1)
    # One enrichment-score range for every running sum in the project, mouse and human alike,
    # so the height of a curve carries the same meaning wherever it is drawn. A range fitted
    # to this panel's own data would make this modest excursion look like a strong one.
    ax.set_ylim(*RS_YLIM)
    ax.set_ylabel("Running enrichment score")
    ax.set_title(
        "Curated HSR-core sign flip is a trend, not a significant Treg effect\n"
        "SECONDARY ANNOTATION TIER — Treg positive; Tcon and CD8 negative"
    )
    ax.legend(frameon=False, fontsize=LEGEND_SIZE, loc="lower left")

    for i, pop in enumerate(pops):
        tr = traces[pop]
        hits = tr.loc[tr["hit"], "rank_fraction"].to_numpy()
        y = len(pops) - 1 - i
        rug.vlines(hits, y, y + 0.86, color=POP_COL[pop], lw=1.0)
    rug.set_ylim(-0.1, len(pops))
    rug.set_yticks([len(pops) - 1 - i + 0.43 for i in range(len(pops))])
    rug.set_yticklabels(pops, fontsize=ANNOT_SIZE)

    # The metric panel earns its space by showing the three rankings crossing zero at
    # comparable fractions, which is the assumption the shared x axis rests on.
    for pop in pops:
        tr = traces[pop]
        met.plot(tr["rank_fraction"], tr["stat"], color=POP_COL[pop], lw=1.0, alpha=0.75)
    met.axhline(0, color="black", lw=0.8)
    met.set_ylabel("Moderated t")
    met.set_xlabel("Position in the ranked list, as a fraction of its length "
                   "(0 = most synovial-fluid-up, 1 = most blood-up)")
    met.set_xlim(0, 1)
    return fig


# ===========================================================================
def main() -> None:
    set_paper_style(config=FIG_CFG)
    purge_figures(STAGE, "hsr_", overview=True, config=FIG_CFG)

    traces, summary = running_sum_traces()
    if not traces:
        print("[10_hsr_lens_viz] no HSR_core traces found — nothing drawn")
        return

    fig = plot_running_sum(traces, summary)
    save_overview(
        fig, STAGE, "hsr_core_running_sum",
        table=round_numeric_cols(summary),
        finding=(f"The curated HSR core changes sign at trend level: "
                 + ", ".join(f"{r['population']} NES {float(r['nes']):+.4f} at FDR "
                             f"{float(r['padj']):.4f}" for _, r in summary.iterrows())
                 + f", with {int(summary['set_size'].min())} of "
                 f"{int(summary['n_nominal'].iloc[0])} genes testable in every ranking. Every "
                 f"population sits above FDR {FDR}, so this secondary annotation supplies "
                 f"directional context. A Treg-selective effect is untested here."),
        script=SCRIPT, fn="plot_running_sum",
        config_kv=(f"figures.running_sum_heights={RS_HEIGHTS}; "
                   f"figures.running_sum_ylim={RS_YLIM}; running_sum_x=rank/n_ranked; "
                   f"thresholds.gsea_fdr={FDR}; evidence_tier=secondary_annotation"),
        input=("03_results/10_hsr_lens/tables/runsum_interactive_hsr_gsea_"
               "{treg,tcon,cd8}_HSR_core.csv, 03_results/10_hsr_lens/tables/hsr_lens_nes.csv"),
        how_to_read=(
            "Three stacked panels sharing one x axis, which is each gene's position in its own "
            "population's ranked list as a FRACTION of that list's length, most "
            "synovial-fluid-up at 0 and most blood-up at 1, because the three rankings differ "
            "in length. Top panel: the weighted running enrichment score as each ranked list is "
            "walked left to right; a positive, left-shifted excursion is synovial-fluid "
            "enrichment, a negative trace the opposite. Its y range is pinned to "
            f"[{RS_YLIM[0]}, {RS_YLIM[1]}], the one range every running sum in this project is "
            "drawn on, so the height of a curve means the same thing here as anywhere else. "
            "Middle panel: where each population's HSR core genes sit in its ranking, one "
            "labelled row per population in matching colour. Bottom panel: the ranked moderated "
            "t each curve above was computed on, which shows how much signal each rank carries "
            "and where the three rankings cross zero — the assumption the shared fractional axis "
            "rests on. Legend labels carry each NES and FDR; the Treg trace is a trend at FDR "
            f"{float(summary.loc[summary['population'].eq('Treg'), 'padj'].iloc[0]):.3f}. "
            "The legend also gives effective size against the "
            f"{int(summary['n_nominal'].iloc[0])}-gene nominal set and its testability band. "
            f"Secondary annotation tier; the confirmatory spine carries any Treg-selective "
            f"claim."),
        config=FIG_CFG, width=RS_W, height=RS_H,
    )
    plt.close(fig)
    print("[10_hsr_lens_viz] wrote 1 overview (HSR core running sum)")


if __name__ == "__main__":
    main()
