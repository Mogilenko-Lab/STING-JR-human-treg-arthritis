#!/usr/bin/env python
"""
10_hsr_lens_viz.py — VIZ ONLY. The activation-free HSR lens, drawn as a trend.
=============================================================================
Heat-shock readouts in inflamed tissue are confounded by activation, so a
curated, activation-free proteostasis lens was built to ask what is left. The
answer is selective in SIGN — HSR core enriches toward synovial fluid in Treg
and away from it in Tcon and CD8 — while sitting just outside significance in
Treg. Both halves of that sentence have to be visible in the figures, so no
figure here carries a star or any other glyph implying significance the numbers
do not have.

  1. hsr_nes_by_population  — the sign flip, with every FDR on the face of it
  2. hsr_core_running_sum   — the same flip read off the ranked lists themselves
  3. hsr_wtheatup_colocalization — do the empirical and curated lenses mark the same cells

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
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
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

STAGE = "10_hsr_lens"
SCRIPT = "02_analysis/scripts/10_hsr_lens_viz.py"
POP_TAG = {"Treg": "treg", "Tcon": "tcon", "CD8": "cd8"}
POP_ORDER = list(POP_TAG)
TERM_ORDER = ["HSR_core", "HSR_sensitivity"]

# Declared palette constants, all resolved from `colors.okabe_ito` in the config
# so no hex literal appears in this script.
_OKABE = (FIG_CFG.get("colors", {}) or {}).get("okabe_ito", {}) or {}
# HSR_core is the 56-gene subset of the 176-gene sensitivity lens, so the two get
# one hue in two depths rather than two unrelated colours.
TERM_COL = {"HSR_core": _OKABE["blue"], "HSR_sensitivity": _OKABE["sky_blue"]}
# Population colours match 05_score_signatures_viz.py, so a population keeps one
# colour across the compartment.
POP_COL = {"Treg": _OKABE["bluish_green"], "Tcon": _OKABE["orange"],
           "CD8": _OKABE["reddish_purple"]}
COLOC_COL = _OKABE["vermillion"]

_F = FIG_CFG.get("figures", {}) or {}
ANNOT_SIZE = float(_F["axis_text_size"])
LEGEND_SIZE = float(_F["legend_text_size"])
RS_HEIGHTS = [float(h) for h in _F["running_sum_heights"]]
FDR = float(PARAMS.gsea_fdr)


def fmt_fdr(p: float) -> str:
    if pd.isna(p):
        return "FDR n/a"
    return f"FDR {p:.3f}" if p >= 0.001 else f"FDR {p:.0e}"


def bar_style(padj: float, colour: str) -> dict:
    """Solid fill only when the bar clears FDR; otherwise open and hatched.

    The stage's whole point is a sign-selective result that is NOT significant, so
    significance is carried by the bar's own surface (legend-keyed) rather than by a
    star that would overstate it.
    """
    if pd.notna(padj) and padj < FDR:
        return {"color": colour, "edgecolor": colour, "linewidth": 1.4, "hatch": None}
    return {"color": "white", "edgecolor": colour, "linewidth": 2.0, "hatch": "//"}


# ===========================================================================
# 1. NES by population — the sign flip, with the FDRs on the face of it
# ===========================================================================
def hsr_nes_by_population() -> pd.DataFrame:
    df = pd.read_csv(PATHS.tables(STAGE) / "hsr_lens_nes.csv")
    plot_df = df[df["signature"].isin(TERM_ORDER)].copy()
    plot_df["population"] = pd.Categorical(plot_df["population"], POP_ORDER, ordered=True)
    plot_df["signature"] = pd.Categorical(plot_df["signature"], TERM_ORDER, ordered=True)
    return plot_df.sort_values(["population", "signature"]).reset_index(drop=True)


def plot_hsr_nes(plot_df: pd.DataFrame):
    x = np.arange(len(POP_ORDER))
    width = 0.34
    fig, ax = plt.subplots()
    for i, term in enumerate(TERM_ORDER):
        sub = plot_df[plot_df["signature"] == term].set_index("population").reindex(POP_ORDER)
        xpos = x + (i - 0.5) * width
        for xp, nes, padj in zip(xpos, sub["nes"].astype(float), sub["padj"].astype(float)):
            if not np.isfinite(nes):
                continue
            ax.bar(xp, nes, width=width, zorder=2, **bar_style(padj, TERM_COL[term]))
            va, off = ("bottom", 0.06) if nes >= 0 else ("top", -0.06)
            ax.text(xp, nes + off, f"{nes:+.2f}\n{fmt_fdr(padj)}", ha="center", va=va,
                    fontsize=ANNOT_SIZE)

    ax.axhline(0, color="black", lw=1, zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels(POP_ORDER)
    ax.set_ylim(-2.15, 2.35)
    ax.set_ylabel("fgsea NES (synovial fluid vs paired blood)")
    ax.set_xlabel("Sorted population")
    ax.set_title("HSR core points into synovial fluid in Treg and out of it elsewhere")

    core_treg = plot_df[(plot_df["signature"] == "HSR_core")
                        & (plot_df["population"] == "Treg")]
    lede = ""
    if len(core_treg):
        nes = float(core_treg["nes"].iloc[0])
        padj = float(core_treg["padj"].iloc[0])
        lede = (f"HSR core in Treg is {nes:+.2f} at {fmt_fdr(padj)} — a trend, not a\n"
                f"significant result. The selectivity of the sign is the observation.")
    ax.text(0.98, 0.97,
            "Positive NES: enriched toward synovial-fluid-up genes.\n" + lede,
            transform=ax.transAxes, ha="right", va="top", fontsize=ANNOT_SIZE)

    handles = [Patch(facecolor=TERM_COL[t], edgecolor=TERM_COL[t],
                     label=t.replace("_", " ")) for t in TERM_ORDER]
    handles += [
        Patch(facecolor="white", edgecolor="grey", linewidth=1.4,
              label=f"solid: FDR < {FDR}"),
        Patch(facecolor="white", edgecolor="grey", linewidth=2.0, hatch="//",
              label=f"open, hatched: FDR ≥ {FDR}"),
    ]
    ax.legend(handles=handles, loc="lower left", frameon=False, fontsize=LEGEND_SIZE)
    fig.tight_layout()
    return fig


# ===========================================================================
# 2. HSR_core running sum — the same flip, read off the ranked lists
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
    traces, rows = {}, []
    for pop, tag in POP_TAG.items():
        path = tdir / f"runsum_interactive_hsr_gsea_{tag}_HSR_core.csv"
        if not path.exists():
            print(f"[10_hsr_lens_viz] {pop}: no HSR_core trace at {path} — skipping")
            continue
        tr = pd.read_csv(path, usecols=["rank", "running_es", "hit"])
        tr["hit"] = _bool_col(tr["hit"])
        traces[pop] = tr
        r = nes.loc[pop]
        rows.append({
            "population": pop,
            "signature": "HSR_core",
            "contrast": "SF_vs_PB",
            "nes": float(r["nes"]),
            "pvalue": float(r["pvalue"]),
            "padj": float(r["padj"]),
            "set_size": int(r["set_size"]),
            "n_ranked_genes": int(len(tr)),
            "evidence_tier": "secondary_annotation",
        })
    return traces, pd.DataFrame(rows)


def plot_running_sum(traces: dict[str, pd.DataFrame], summary: pd.DataFrame):
    # Two stacked panels in the config's running-sum proportions: the enrichment
    # trace over the gene-hit rug. The config `running_sum_ylim` pins a range ACROSS
    # a family of separate per-population figures; here all three populations share
    # one axis, so comparability is intrinsic and the range is data-driven so the
    # sign flip stays legible.
    fig, (ax, rug) = plt.subplots(
        2, 1, sharex=True, layout="constrained",
        height_ratios=[RS_HEIGHTS[0], RS_HEIGHTS[1]])

    span = max(float(np.abs(t["running_es"]).max()) for t in traces.values())
    pops = [p for p in POP_ORDER if p in traces]
    for pop in pops:
        tr = traces[pop]
        row = summary[summary["population"] == pop].iloc[0]
        ax.plot(tr["rank"], tr["running_es"], color=POP_COL[pop], lw=2.0,
                label=f"{pop}   NES {row['nes']:+.2f}, {fmt_fdr(row['padj'])}")
    ax.axhline(0, color="black", lw=1)
    ax.set_ylim(-span * 1.25, span * 1.25)
    ax.set_ylabel("Running enrichment score")
    ax.set_title("HSR core sits at the synovial-fluid end only in Treg")
    ax.legend(frameon=False, fontsize=LEGEND_SIZE, loc="lower left")

    for i, pop in enumerate(pops):
        tr = traces[pop]
        hits = tr.loc[tr["hit"], "rank"].to_numpy()
        y = len(pops) - 1 - i
        rug.vlines(hits, y, y + 0.86, color=POP_COL[pop], lw=1.0)
    rug.set_ylim(-0.1, len(pops))
    rug.set_yticks([len(pops) - 1 - i + 0.43 for i in range(len(pops))])
    rug.set_yticklabels(pops, fontsize=ANNOT_SIZE)
    rug.set_xlabel("Rank in the synovial-fluid-vs-blood ranked list")
    return fig


# ===========================================================================
# 3. Colocalization — do the empirical and curated lenses mark the same cells?
# ===========================================================================
def hsr_wtheatup_colocalization() -> pd.DataFrame:
    tdir = PATHS.tables(STAGE)
    df = pd.read_csv(tdir / "hsr_colocalization.csv")
    plot_df = df[
        (df["hsr_term"] == "HSR_core")
        & (df["level"] == "cell")
        & (df["method"] == "spearman")
    ].copy()
    plot_df["population"] = pd.Categorical(plot_df["population"], POP_ORDER, ordered=True)
    plot_df = plot_df.sort_values("population").reset_index(drop=True)

    # Carry the gene-overlap context into the plotted table: it is what makes a low
    # r interpretable, so it belongs with the bars rather than in a separate figure.
    ov = pd.read_csv(tdir / "hsr_wtheatup_overlap.csv")
    ov = ov[ov["set_b"] == "HSR_core"]
    if len(ov):
        o = ov.iloc[0]
        plot_df["n_wtheatup_genes"] = int(o["n_a"])
        plot_df["n_hsr_core_genes"] = int(o["n_b"])
        plot_df["n_genes_shared"] = int(o["n_intersect"])
        plot_df["genes_shared"] = str(o["genes_intersect"])
    return plot_df


def plot_colocalization(plot_df: pd.DataFrame):
    x = np.arange(len(plot_df))
    fig, ax = plt.subplots()
    vals = plot_df["r"].to_numpy(dtype=float)
    ax.bar(x, vals, color=COLOC_COL, width=0.55)
    ax.axhline(0, color="black", lw=1)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{p}\nn = {int(n):,} cells"
                        for p, n in zip(plot_df["population"].astype(str), plot_df["n"])])
    ax.set_ylim(-0.05, 1.0)
    ax.set_ylabel("Spearman r, within synovial-fluid cells")
    ax.set_xlabel("Sorted population")
    ax.set_title("The mouse and curated lenses mark largely different cells")
    for xp, r in zip(x, vals):
        va, off = ("bottom", 0.02) if r >= 0 else ("top", -0.02)
        ax.text(xp, r + off, f"{r:.2f}", ha="center", va=va, fontsize=ANNOT_SIZE)

    note = ("Full –0.05 to 1 scale, so bar height is the whole story.\n"
            "This channel is a diagnostic of where two scores sit in the\n"
            "same cells; it is not evidence for a temperature program.")
    if "n_genes_shared" in plot_df.columns and len(plot_df):
        r0 = plot_df.iloc[0]
        note += (f"\nThe lenses share only {int(r0['n_genes_shared'])} genes of "
                 f"{int(r0['n_wtheatup_genes'])} and {int(r0['n_hsr_core_genes'])} "
                 f"({r0['genes_shared'].replace(';', ', ')}), so they are\n"
                 "near-independent by construction rather than by accident.")
    ax.text(0.02, 0.96, note, transform=ax.transAxes, ha="left", va="top",
            fontsize=ANNOT_SIZE)
    fig.tight_layout()
    return fig


# ===========================================================================
def main() -> None:
    set_paper_style(config=FIG_CFG)
    purge_figures(STAGE, "hsr_", overview=True, config=FIG_CFG)

    nes_df = hsr_nes_by_population()
    fig = plot_hsr_nes(nes_df)
    save_overview(
        fig, STAGE, "hsr_nes_by_population",
        table=round_numeric_cols(nes_df),
        finding=("The activation-free HSR lens is selective in sign, not in strength: HSR core "
                 "enriches toward synovial fluid in Treg (+1.50) and away from it in Tcon "
                 "(-1.36) and CD8 (-1.10), with Treg at FDR 0.056 — a trend, not a "
                 "significant result."),
        script=SCRIPT, fn="plot_hsr_nes",
        config_kv=(f"thresholds.gsea_fdr={FDR}; gsea_min_size={PARAMS.gsea_min_size}; "
                   "evidence_tier=secondary_annotation"),
        input="03_results/10_hsr_lens/tables/hsr_lens_nes.csv",
        how_to_read=(
            "Bars are fgsea NES for the two curated HSR terms on each population's "
            "synovial-fluid-vs-paired-blood ranked list; positive means enriched toward the "
            "synovial-fluid-up end, negative toward blood. Every bar is labelled with its NES "
            f"and its FDR. A bar is solid only when it clears FDR < {FDR} and is drawn open "
            "and hatched otherwise, which on this data is every bar — read the sign pattern "
            "across populations, not any single bar's magnitude, and read nothing here as "
            "significant. Deep blue is the 56-gene HSR core, pale blue the 176-gene "
            "sensitivity lens that contains it. Annotation tier, firewalled from the "
            "confirmatory WT_heat effect-size spine."),
        config=FIG_CFG,
    )
    plt.close(fig)

    traces, summary = running_sum_traces()
    if traces:
        fig = plot_running_sum(traces, summary)
        save_overview(
            fig, STAGE, "hsr_core_running_sum",
            table=round_numeric_cols(summary),
            finding=("Walking each population's ranked list, HSR core accumulates a positive "
                     "peak near the synovial-fluid end in Treg while Tcon and CD8 run negative "
                     "throughout, so the sign selectivity is a property of the rankings and not "
                     "an artefact of the summary statistic."),
            script=SCRIPT, fn="plot_running_sum",
            config_kv=(f"figures.running_sum_heights={RS_HEIGHTS[:2]}; "
                       f"thresholds.gsea_fdr={FDR}; evidence_tier=secondary_annotation"),
            input=("03_results/10_hsr_lens/tables/runsum_interactive_hsr_gsea_"
                   "{treg,tcon,cd8}_HSR_core.csv, 03_results/10_hsr_lens/tables/hsr_lens_nes.csv"),
            how_to_read=(
                "Top panel: the weighted running enrichment score as each ranked list is walked "
                "from synovial-fluid-up (left) to blood-up (right); a positive, left-shifted "
                "excursion is synovial-fluid enrichment, a negative trace the opposite. Bottom "
                "panel: where each population's HSR core genes sit in its ranking, in matching "
                "colour. Legend labels carry each NES and FDR, so read the Treg trace as a "
                "trend at FDR 0.056, not a significant enrichment. Ranked-list lengths differ "
                "slightly, so compare shapes rather than x positions; the y range is data-driven "
                "because all three curves share one axis. The source table carries the annotated "
                "numbers, the traces are the cited inputs. Annotation tier."),
            config=FIG_CFG, height=7.0,
        )
        plt.close(fig)

    coloc_df = hsr_wtheatup_colocalization()
    fig = plot_colocalization(coloc_df)
    save_overview(
        fig, STAGE, "hsr_wtheatup_colocalization",
        table=round_numeric_cols(coloc_df),
        finding=("Within synovial-fluid cells the empirical mouse lens and the curated HSR core "
                 "correlate at only 0.11 to 0.19 while sharing just two genes, so they label "
                 "largely different cells and the curated lens is an independent probe rather "
                 "than a restatement."),
        script=SCRIPT, fn="plot_colocalization",
        config_kv="level=cell; method=spearman; tissue=synovial_fluid; evidence_tier=secondary_percell",
        input=("03_results/10_hsr_lens/tables/hsr_colocalization.csv, "
               "03_results/10_hsr_lens/tables/hsr_wtheatup_overlap.csv"),
        how_to_read=(
            "Each bar is the within-synovial-fluid, cell-level Spearman r between "
            "WT_heat_up_AUCell and HSR_core_AUCell, with the cell count under each population. "
            "The axis runs the full -0.05 to 1 range on purpose: read how short the bars are. "
            "Positive r means a heat-high cell tends to be HSR-high. The in-figure gene-overlap "
            "line is what makes a low r interpretable — the two lenses share two genes out of "
            "199 and 56, so they are near-independent by construction. Secondary per-cell tier: "
            "this is a diagnostic of where two scores sit in the same cells and is never read as "
            "evidence for the temperature program, which rests on the donor-pseudobulk "
            "enrichment instead."),
        config=FIG_CFG,
    )
    plt.close(fig)
    print("[10_hsr_lens_viz] wrote 3 overviews (NES by population, HSR core running sum, "
          "colocalization)")


if __name__ == "__main__":
    main()
