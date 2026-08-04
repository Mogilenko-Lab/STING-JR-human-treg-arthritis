#!/usr/bin/env python
"""
00_build_anndata_viz.py — VIZ (no statistics).
==============================================
Reads the tidy tables written by 00_build_anndata.py and renders the design
overview: cells recovered per GSM (donor x tissue x population). Shows the
paired SF/PB design is intact and which strata are thin (esp. SF-Treg p5).

All figures route through the figure-style contract (save_overview).
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "02_analysis"))
os.chdir(ROOT)

from config import PATHS, COARSE_LABEL, TISSUE_COLORS  # noqa: E402
from helpers.figure_style import set_paper_style, save_overview, FIG_CFG  # noqa: E402

STAGE = "00_build"
SCRIPT = "02_analysis/scripts/00_build_anndata_viz.py"
# The one tissue palette, read from analysis_config.yaml::colors.tissue.
TISSUE_COL = TISSUE_COLORS


def main() -> None:
    set_paper_style(config=FIG_CFG)
    tdir = PATHS.tables(STAGE)
    cpg = pd.read_csv(tdir / "cells_per_gsm.csv")
    cpg["pop"] = cpg["population"].map(COARSE_LABEL)
    cpg["donor_id"] = cpg["donor"].str.replace("JIA_patient_", "p", regex=False)

    pops = ["Treg", "Tcon", "CD8"]
    fig, axes = plt.subplots(1, 3, figsize=(13, 5), sharey=True)
    for ax, pop in zip(axes, pops):
        sub = cpg[cpg["pop"] == pop]
        donors = sorted(sub["donor_id"].unique())
        x = range(len(donors))
        w = 0.4
        for i, (tissue, off) in enumerate([("synovial_fluid", -w / 2), ("peripheral_blood", w / 2)]):
            vals = [int(sub[(sub["donor_id"] == d) & (sub["tissue"] == tissue)]["n_cells"].sum())
                    for d in donors]
            ax.bar([xi + off for xi in x], vals, width=w, color=TISSUE_COL[tissue],
                   label=("SF" if tissue == "synovial_fluid" else "PB"))
        ax.set_xticks(list(x))
        ax.set_xticklabels(donors)
        ax.set_title(pop)
        ax.set_xlabel("donor")
    axes[0].set_ylabel("cells recovered")
    axes[-1].legend(title="tissue", frameon=True)
    fig.suptitle("GSE160097 cells per GSM (sorted population x tissue x donor)")
    fig.tight_layout()

    save_overview(
        fig, STAGE, "cells_per_gsm",
        table=cpg[["gsm", "donor", "tissue", "population", "n_cells"]],
        finding=("At ingest, all 7 donors contribute SF+PB Treg samples; Tcon and CD8 lack "
                 "a PB sample for p3 by design. The near-empty SF-Treg p5 sample is later "
                 "removed by QC, leaving 6 paired donors in each analyzed population."),
        script=SCRIPT, fn="main",
        config_kv="design.populations = [CD4_Treg, CD4_Tcon, CD8]",
        input="03_results/00_build/tables/cells_per_gsm.csv",
        how_to_read=("Grouped bars = cells recovered per donor; orange = synovial fluid (SF), "
                     "blue = peripheral blood (PB); one facet per sorted population. A missing "
                     "PB bar (p3 in Tcon/CD8) is an intentionally-absent sample, not a QC drop. "
                     "These are ingest counts before QC; the donor-level analysis uses 6 paired "
                     "donors per population. Descriptive counts only — no claim tier."),
        config=FIG_CFG, wide=True,
    )
    print("[00_build_viz] wrote cells_per_gsm overview")


if __name__ == "__main__":
    main()
