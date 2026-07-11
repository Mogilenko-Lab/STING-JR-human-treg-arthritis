#!/usr/bin/env python
"""
05_score_signatures_viz.py — VIZ (no statistics). The go/no-go money figure.
============================================================================
Reads the fgsea tables + score tables from 05_score_signatures.py and renders:
  - the WT_heat NES forest across Treg / Tcon / CD8 (the Treg-preference check —
    the money panel for the cross-species comparison);
  - per-cell WT_heat_up score violins SF-Treg vs PB-Treg vs SF-Tcon vs SF-CD8
    (running-sum figure moved to 05_score_signatures_viz.R).
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "02_analysis"))
os.chdir(ROOT)

from config import PATHS, PARAMS  # noqa: E402
from helpers.figure_style import set_paper_style, save_overview, FIG_CFG  # noqa: E402

STAGE = "05_scoring"
SCRIPT = "02_analysis/scripts/05_score_signatures_viz.py"
POP_TAG = {"Treg": "treg", "Tcon": "tcon", "CD8": "cd8"}
POP_COL = {"Treg": "#009E73", "Tcon": "#E69F00", "CD8": "#CC79A7"}


def main() -> None:
    set_paper_style(config=FIG_CFG)
    tdir = PATHS.tables(STAGE)

    # ---- gather NES ----
    rows = []
    for pop, tag in POP_TAG.items():
        f = tdir / f"gsea_pseudobulk_{tag}.csv"
        if f.exists():
            g = pd.read_csv(f)
            g["cell_state"] = pop
            rows.append(g)
    gsea = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()

    # ---- 1. NES forest (money) ----
    fig, ax = plt.subplots(figsize=(8.5, 6))
    if not gsea.empty:
        gsea["ypos"] = gsea["direction"].map({"up": 1, "down": 0})
        order = [(p, d) for d in ["up", "down"] for p in ["Treg", "Tcon", "CD8"]]
        ylabels, yv = [], []
        for i, (pop, d) in enumerate(order):
            sub = gsea[(gsea["cell_state"] == pop) & (gsea["direction"] == d)]
            if sub.empty:
                continue
            nes = float(sub["nes"].iloc[0]); padj = float(sub["padj"].iloc[0])
            ax.scatter(nes, i, s=140, color=POP_COL[pop],
                       edgecolors="black", zorder=3,
                       marker="o" if d == "up" else "D")
            star = "*" if padj < float(PARAMS.gsea_fdr) else ""
            ax.text(nes, i + 0.18, f"{nes:.2f}{star}", ha="center", fontsize=9)
            ylabels.append(f"{pop} · {d}"); yv.append(i)
        ax.axvline(0, ls="--", c="grey", lw=1)
        ax.set_yticks(yv); ax.set_yticklabels(ylabels)
        ax.set_xlabel("WT_heat NES (SF synovial-fluid vs PB peripheral-blood)")
        ax.set_title("Mouse 39 °C Treg signature — SF-vs-PB NES across sorted populations")
        from matplotlib.lines import Line2D
        # Legend has two keys: colour = population (matches the point colours), shape = up/down set.
        pop_handles = [
            Line2D([0], [0], marker="o", color="w", markerfacecolor=POP_COL[p],
                   markeredgecolor="k", label=p, markersize=9)
            for p in ["Treg", "Tcon", "CD8"]
        ]
        shape_handles = [
            Line2D([0], [0], marker="o", color="w", markerfacecolor="white",
                   markeredgecolor="k", label="up set (circle)", markersize=9),
            Line2D([0], [0], marker="D", color="w", markerfacecolor="white",
                   markeredgecolor="k", label="down set (diamond)", markersize=9),
            Line2D([0], [0], marker="", linestyle="", label=f"* = FDR < {PARAMS.gsea_fdr}"),
        ]
        ax.legend(handles=pop_handles + shape_handles, frameon=True, loc="lower left",
                  fontsize=9, title="population (colour) · set (shape)",
                  title_fontsize=9)  # lower-left is empty (up-points sit at high NES)
    fig.tight_layout()
    save_overview(fig, STAGE, "wt_heat_nes_forest",
                  table=gsea[["cell_state", "pathway_id", "nes", "pvalue", "padj", "set_size"]]
                  if not gsea.empty else pd.DataFrame(),
                  finding=("The go/no-go readout: whether the mouse 39 °C Treg up-program enriches "
                           "in JIA SF-vs-PB Tregs, and whether that enrichment is Treg-preferential "
                           "over Tcon/CD8."),
                  script=SCRIPT, fn="main",
                  config_kv=f"thresholds.gsea_fdr={PARAMS.gsea_fdr}; gsea_min_size={PARAMS.gsea_min_size}",
                  input="03_results/05_scoring/tables/gsea_pseudobulk_{treg,tcon,cd8}.csv",
                  how_to_read=("Points = fgsea NES for WT_heat up (circle) / down (diamond), colored by "
                               "population; x=0 dashed = no enrichment; * = FDR < threshold. Read "
                               "Treg-preference as: is the Treg-up NES the largest and significant? "
                               "Correlative (consistent-with), not causal."),
                  config=FIG_CFG)

    # ---- 2. per-cell score violins ----
    dm = pd.read_csv(tdir / "donor_label_score_means.csv")
    dm["tissue_s"] = dm["tissue"].map({"synovial_fluid": "SF", "peripheral_blood": "PB"})
    dm["group"] = dm["coarse_label"] + " " + dm["tissue_s"]
    groups = ["Treg SF", "Treg PB", "Tcon SF", "Tcon PB", "CD8 SF", "CD8 PB"]
    fig3, ax = plt.subplots(figsize=(9, 6))
    data, positions, colors = [], [], []
    for i, grp in enumerate(groups):
        vals = dm.loc[dm["group"] == grp, "WT_heat_up"].values
        if len(vals):
            data.append(vals); positions.append(i)
            colors.append(POP_COL[grp.split()[0]])
    parts = ax.violinplot(data, positions=positions, showmeans=True, widths=0.8)
    for b, c in zip(parts["bodies"], colors):
        b.set_facecolor(c); b.set_alpha(0.65)
    for i, grp in zip(positions, [groups[p] for p in positions]):
        vals = dm.loc[dm["group"] == grp, "WT_heat_up"].values
        ax.scatter(np.full(len(vals), i), vals, s=18, c="black", zorder=3)
    ax.axhline(0, ls="--", c="grey", lw=0.8)
    ax.set_xticks(range(len(groups))); ax.set_xticklabels(groups, rotation=30, ha="right")
    ax.set_ylabel("donor-mean WT_heat_up score")
    ax.set_title("Per-cell WT_heat_up score, donor means by state × tissue")
    fig3.tight_layout()
    save_overview(fig3, STAGE, "score_violins", table=dm,
                  finding=("Corroborative per-cell view: donor-mean WT_heat score in SF vs PB across "
                           "Treg/Tcon/CD8 — is the SF-vs-PB shift Treg-preferential?"),
                  script=SCRIPT, fn="main",
                  config_kv="signature = WT_heat_up/down (scanpy score_genes)",
                  input="03_results/05_scoring/tables/donor_label_score_means.csv",
                  how_to_read=("Each dot = one donor's mean WT_heat_up score for that state×tissue; "
                               "violins summarise across donors. Scores sit near/below zero because "
                               "scanpy score_genes centres each cell against a random reference gene set "
                               "— the absolute level is arbitrary; only the RELATIVE SF-vs-PB shift is "
                               "read (e.g. Treg SF > Treg PB). This is a different estimand from the "
                               "forest: that NES is fgsea on the donor-pseudobulk SF-vs-PB signed-Wald "
                               "ranking, asking whether the WT_heat set concentrates at the top of the "
                               "whole ranked list — a normalised statistic on its own scale, positive "
                               "here. Secondary tier (percell) — NEVER pooled with the pseudobulk NES. "
                               "Down arm omitted (up and down co-shift in SF). Correlative."),
                  config=FIG_CFG)
    print("[05_scoring_viz] wrote 2 overviews (NES forest + score violins)")


if __name__ == "__main__":
    main()
