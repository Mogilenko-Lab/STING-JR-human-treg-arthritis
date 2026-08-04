#!/usr/bin/env python
"""
03_pseudobulk_de_viz.py — VIZ (no statistics).
==============================================
Reads the pseudobulk matrices + DE tables from 03_pseudobulk_de.py and renders:
  - pseudobulk PCA colored by tissue + donor (confirms SF/PB separation; checks
    single-donor dominance red-flag);
  - SF-vs-PB Treg volcano;
  - per-population significant-DE count bar.
No statistics are recomputed here (PCA of log-CPM is a display transform).
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

from config import PATHS, PARAMS, POPULATION_COLORS  # noqa: E402
from helpers.figure_style import set_paper_style, save_overview, FIG_CFG  # noqa: E402

STAGE = "03_pseudobulk"
SCRIPT = "02_analysis/scripts/03_pseudobulk_de_viz.py"
TISSUE_MARK = {"synovial_fluid": "o", "peripheral_blood": "s"}
# The one population palette, read from analysis_config.yaml::colors.populations.
LABEL_COL = POPULATION_COLORS


def main() -> None:
    set_paper_style(config=FIG_CFG)
    tdir = PATHS.tables(STAGE)
    counts = pd.read_csv(tdir / "pseudobulk_counts.csv", index_col=0)
    coldata = pd.read_csv(tdir / "pseudobulk_coldata.csv", index_col=0)
    summary = pd.read_csv(tdir / "de_summary.csv")

    # ---- 1. pseudobulk PCA (log-CPM, HVG) ----
    cpm = counts.div(counts.sum(axis=1) + 1, axis=0) * 1e6
    logcpm = np.log1p(cpm)
    hv = logcpm.var(axis=0).sort_values(ascending=False).head(2000).index
    Xc = logcpm[hv].values
    Xc = Xc - Xc.mean(axis=0)
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    pcs = U[:, :2] * S[:2]
    var_expl = (S ** 2 / (S ** 2).sum())[:2]
    coldata = coldata.loc[counts.index]
    donors = sorted(coldata["donor"].unique())
    dpal = {d: plt.cm.tab10(i % 10) for i, d in enumerate(donors)}

    fig, ax = plt.subplots(figsize=(8.5, 7))
    for i, sid in enumerate(counts.index):
        row = coldata.loc[sid]
        ax.scatter(pcs[i, 0], pcs[i, 1], marker=TISSUE_MARK.get(row["tissue"], "o"),
                   c=[dpal[row["donor"]]], s=90, edgecolors="black", linewidths=0.6)
    ax.set_xlabel(f"PC1 ({var_expl[0]:.0%})"); ax.set_ylabel(f"PC2 ({var_expl[1]:.0%})")
    ax.set_title("Pseudobulk PCA — circle=SF, square=PB, color=donor")
    from matplotlib.lines import Line2D
    shape_handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor="grey",
                            markeredgecolor="k", label="SF (synovial fluid)", markersize=9),
                     Line2D([0], [0], marker="s", color="w", markerfacecolor="grey",
                            markeredgecolor="k", label="PB (peripheral blood)", markersize=9)]
    donor_handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor=dpal[d],
                            markeredgecolor="k", label=d.replace("JIA_patient_", "p"), markersize=8)
                     for d in donors]
    leg1 = ax.legend(handles=shape_handles, frameon=True, loc="upper left", title="tissue")
    ax.add_artist(leg1)
    ax.legend(handles=donor_handles, frameon=True, loc="upper right", title="donor",
              fontsize=8, ncol=2)
    fig.tight_layout()
    pca_tab = coldata.assign(PC1=pcs[:, 0], PC2=pcs[:, 1])[
        ["donor", "tissue", "coarse_label", "n_cells", "PC1", "PC2"]]
    save_overview(fig, STAGE, "pseudobulk_pca", table=pca_tab,
                  finding=("Pseudobulk samples separate by tissue and label without a single donor "
                           "dominating an axis, so donor pseudobulk is well-posed for SF-vs-PB DE."),
                  script=SCRIPT, fn="main",
                  config_kv="thresholds.pseudobulk_min_cells (strata filter)",
                  input="03_results/03_pseudobulk/tables/pseudobulk_counts.csv",
                  how_to_read=("Each point = one donor x tissue x label pseudobulk (log-CPM, top-2000 "
                               "var genes). Circle=SF, square=PB, color=donor. Look for tissue "
                               "separation and NO single-donor axis dominance. Display transform only."),
                  config=FIG_CFG)

    # ---- 2. Treg SF-vs-PB volcano ----
    treg_path = tdir / "de_SFvsPB_treg.csv"
    if treg_path.exists():
        de = pd.read_csv(treg_path, index_col=0)
        de = de.dropna(subset=["padj", "log2FoldChange"])
        sig = (de["padj"] < float(PARAMS.de_fdr)) & (de["log2FoldChange"].abs() >= float(PARAMS.de_logfc))
        fig2, ax = plt.subplots(figsize=(8, 7))
        ax.scatter(de.loc[~sig, "log2FoldChange"], -np.log10(de.loc[~sig, "padj"] + 1e-300),
                   s=6, c="grey", alpha=0.4, linewidths=0)
        ax.scatter(de.loc[sig, "log2FoldChange"], -np.log10(de.loc[sig, "padj"] + 1e-300),
                   s=10, c="#D55E00", linewidths=0)
        n_lab = int(FIG_CFG.get("figures", {}).get("volcano_label_top", 10))
        top = de.loc[sig].reindex(
            de.loc[sig, "log2FoldChange"].abs().sort_values(ascending=False).index).head(n_lab)
        texts = [ax.text(r["log2FoldChange"], -np.log10(r["padj"] + 1e-300),
                         str(r["gene_symbol"]), fontsize=9) for _, r in top.iterrows()]
        from adjustText import adjust_text
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="grey", lw=0.5))
        from matplotlib.lines import Line2D
        ax.legend(handles=[
            Line2D([0], [0], marker="o", color="w", markerfacecolor="#D55E00",
                   label=f"significant (FDR<{PARAMS.de_fdr}, |log2FC|≥{PARAMS.de_logfc})", markersize=8),
            Line2D([0], [0], marker="o", color="w", markerfacecolor="grey", label="ns", markersize=8),
        ], frameon=True, loc="upper left", fontsize=8)
        ax.axhline(-np.log10(float(PARAMS.de_fdr)), ls="--", c="k", lw=0.6)
        ax.axvline(float(PARAMS.de_logfc), ls="--", c="k", lw=0.6)
        ax.axvline(-float(PARAMS.de_logfc), ls="--", c="k", lw=0.6)
        ax.set_xlabel("log2FC (SF / PB)"); ax.set_ylabel("-log10 padj")
        ax.set_title(f"SF-vs-PB Treg pseudobulk DE ({int(sig.sum())} significant)")
        fig2.tight_layout()
        vol_tab = de[["gene_symbol", "log2FoldChange", "stat", "pvalue", "padj"]].head(500)
        save_overview(fig2, STAGE, "treg_volcano", table=vol_tab,
                      finding=("Synovial-fluid Tregs carry a reproducible SF-vs-PB transcriptional "
                               "program (significant up/down genes), the substrate the mouse signature is tested against."),
                      script=SCRIPT, fn="main",
                      config_kv=f"thresholds.de_fdr={PARAMS.de_fdr}; de_logfc={PARAMS.de_logfc}",
                      input="03_results/03_pseudobulk/tables/de_SFvsPB_treg.csv",
                      how_to_read=("x=log2FC SF/PB, y=-log10 padj; orange = significant. Dashed lines = "
                                   "FDR + |log2FC| gates. Correlative DE, top-500 genes tabulated."),
                      config=FIG_CFG)

    # ---- 3. DE-count bar ----
    fig3, ax = plt.subplots(figsize=(6.5, 5))
    s = summary.copy()
    s["pop_color"] = s["population"].map(LABEL_COL)
    ax.bar(s["population"], s["n_sig_de"], color=s["pop_color"])
    for i, r in s.iterrows():
        ax.text(i, r["n_sig_de"], str(int(r["n_sig_de"])), ha="center", va="bottom")
    ax.set_ylabel("significant SF-vs-PB DE genes")
    ax.set_title("Pseudobulk DE burden per population")
    fig3.tight_layout()
    save_overview(fig3, STAGE, "de_count_bar", table=summary,
                  finding=("All three sorted populations yield significant SF-vs-PB DE, so each has a "
                           "ranked list powered for signature enrichment."),
                  script=SCRIPT, fn="main",
                  config_kv=f"thresholds.de_fdr={PARAMS.de_fdr}; de_logfc={PARAMS.de_logfc}",
                  input="03_results/03_pseudobulk/tables/de_summary.csv",
                  how_to_read=("Bar = # significant SF-vs-PB DE genes per population (Treg/Tcon/CD8). "
                               "Confirms each arm has enough signal to rank for fgsea. Diagnostic."),
                  config=FIG_CFG)
    print("[03_pseudobulk_viz] wrote 3 overviews")


if __name__ == "__main__":
    main()
