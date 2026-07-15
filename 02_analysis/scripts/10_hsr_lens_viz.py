#!/usr/bin/env python
"""
10_hsr_lens_viz.py — VIZ ONLY. Curated HSR second-lens figures.
================================================================
Reads committed stage-10 tables and emits two overview figures plus same-stem
source tables. It computes no statistics.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "02_analysis"))
os.chdir(ROOT)

from config import PATHS  # noqa: E402
from helpers.figure_style import (  # noqa: E402
    FIG_CFG,
    round_numeric_cols,
    save_overview,
    set_paper_style,
)

STAGE = "10_hsr_lens"
POP_ORDER = ["Treg", "Tcon", "CD8"]
TERM_ORDER = ["HSR_core", "HSR_sensitivity"]


def okabe_palette_from_config() -> list[str]:
    colors = ((FIG_CFG.get("colors", {}) or {}).get("okabe_ito", {}) or {})
    if isinstance(colors, dict):
        return list(colors.values())
    return list(colors)


def hsr_nes_by_population() -> pd.DataFrame:
    df = pd.read_csv(PATHS.tables(STAGE) / "hsr_lens_nes.csv")
    plot_df = df[df["signature"].isin(TERM_ORDER)].copy()
    plot_df["population"] = pd.Categorical(plot_df["population"], POP_ORDER, ordered=True)
    plot_df["signature"] = pd.Categorical(plot_df["signature"], TERM_ORDER, ordered=True)
    return plot_df.sort_values(["population", "signature"]).reset_index(drop=True)


def plot_hsr_nes(plot_df: pd.DataFrame):
    set_paper_style(config=FIG_CFG)
    palette = okabe_palette_from_config()
    colors = {TERM_ORDER[0]: palette[0], TERM_ORDER[1]: palette[2]}
    x = np.arange(len(POP_ORDER))
    width = 0.34

    fig, ax = plt.subplots()
    for i, term in enumerate(TERM_ORDER):
        vals = (
            plot_df[plot_df["signature"] == term]
            .set_index("population")
            .reindex(POP_ORDER)["nes"]
            .to_numpy(dtype=float)
        )
        xpos = x + (i - 0.5) * width
        ax.bar(xpos, vals, width=width, label=term.replace("_", " "), color=colors[term])
        for xp, val in zip(xpos, vals):
            if np.isfinite(val):
                va = "bottom" if val >= 0 else "top"
                offset = 0.04 if val >= 0 else -0.04
                ax.text(xp, val + offset, f"{val:.2f}", ha="center", va=va, fontsize=10)

    ax.axhline(0, color="black", linewidth=1)
    ax.set_xticks(x)
    ax.set_xticklabels(POP_ORDER)
    ax.set_ylabel("fgsea NES (SF vs PB)")
    ax.set_xlabel("Sorted population")
    ax.set_title("Curated HSR lens enrichment")
    ax.legend(title="HSR term", frameon=False)
    ax.text(
        0.01,
        0.98,
        "NES > 0: enriched toward synovial-fluid-up genes",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10,
    )
    fig.tight_layout()
    return fig


def hsr_wtheatup_colocalization() -> pd.DataFrame:
    df = pd.read_csv(PATHS.tables(STAGE) / "hsr_colocalization.csv")
    plot_df = df[
        (df["hsr_term"] == "HSR_core")
        & (df["level"] == "cell")
        & (df["method"] == "spearman")
    ].copy()
    plot_df["population"] = pd.Categorical(plot_df["population"], POP_ORDER, ordered=True)
    return plot_df.sort_values("population").reset_index(drop=True)


def plot_colocalization(plot_df: pd.DataFrame):
    set_paper_style(config=FIG_CFG)
    palette = okabe_palette_from_config()
    x = np.arange(len(plot_df))

    fig, ax = plt.subplots()
    vals = plot_df["r"].to_numpy(dtype=float)
    ax.bar(x, vals, color=palette[4])
    ax.axhline(0, color="black", linewidth=1)
    ax.set_xticks(x)
    ax.set_xticklabels(plot_df["population"].astype(str).tolist())
    ax.set_ylabel("Spearman r")
    ax.set_xlabel("SF sorted population")
    ax.set_title("WT_heat_up and HSR_core mark the same SF cells?")
    ax.set_ylim(min(-0.05, np.nanmin(vals) - 0.08), max(0.35, np.nanmax(vals) + 0.08))
    for xp, r, n in zip(x, plot_df["r"], plot_df["n"]):
        va = "bottom" if r >= 0 else "top"
        offset = 0.02 if r >= 0 else -0.02
        ax.text(xp, r + offset, f"r={r:.2f}\nn={int(n)}", ha="center", va=va, fontsize=10)
    ax.text(
        0.01,
        0.98,
        "Low r: empirical WT_heat_up and curated HSR_core label different cells",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10,
    )
    fig.tight_layout()
    return fig


def main() -> None:
    nes_df = hsr_nes_by_population()
    nes_fig = plot_hsr_nes(nes_df)
    save_overview(
        nes_fig,
        STAGE,
        "hsr_nes_by_population",
        table=round_numeric_cols(nes_df),
        finding="The curated HSR lens tests whether the SF-vs-PB WT_heat separation is matched by proteostasis enrichment across the three sorted T-cell populations.",
        script="02_analysis/scripts/10_hsr_lens_viz.py",
        fn="hsr_nes_by_population",
        config_kv="gsea_min_size=5; gsea_max_size=500; evidence_tier=secondary_annotation",
        input="03_results/10_hsr_lens/tables/hsr_lens_nes.csv",
        how_to_read="Bars are fgsea NES values for SF-vs-PB ranked lists. Positive NES means the HSR term is enriched toward synovial-fluid-up genes. This is annotation-tier and is not an effect-size claim.",
        config=FIG_CFG,
    )
    plt.close(nes_fig)

    coloc_df = hsr_wtheatup_colocalization()
    coloc_fig = plot_colocalization(coloc_df)
    save_overview(
        coloc_fig,
        STAGE,
        "hsr_wtheatup_colocalization",
        table=round_numeric_cols(coloc_df),
        finding="Within SF cells, the WT_heat_up-vs-HSR_core Spearman correlation asks whether the empirical mouse lens and curated HSR lens mark the same cells.",
        script="02_analysis/scripts/10_hsr_lens_viz.py",
        fn="hsr_wtheatup_colocalization",
        config_kv="level=cell; method=spearman; tissue=synovial_fluid; evidence_tier=secondary_percell",
        input="03_results/10_hsr_lens/tables/hsr_colocalization.csv",
        how_to_read="Each bar is the within-SF cell-level Spearman r between WT_heat_up_AUCell and HSR_core_AUCell. Low r means WT_heat_up-high and HSR-high cells largely differ; high r means the lenses co-localize. Secondary per-cell tier.",
        config=FIG_CFG,
    )
    plt.close(coloc_fig)
    print("[10_hsr_lens_viz] done")


if __name__ == "__main__":
    main()
