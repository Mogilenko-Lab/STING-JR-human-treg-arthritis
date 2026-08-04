#!/usr/bin/env python
"""
01_qc_mthi_characterize_viz.py — VIZ (no statistics).
=====================================================
Renders the five defensibility panels for the mt-hi effector-like Treg pocket
from the tables written by 01_qc_mthi_characterize.py. Nothing is computed here.

Panels (each a self-contained _overview figure + same-stem source table + README caption):
  A cluster_mt_etreg      cluster-resolved median %mt vs median eTreg (pocket highlighted) + per-cell hexbin
  B identity_retention    FOXP3/IKZF2/CTLA4/IL2RA/TIGIT expression + fraction-expressing, mt-hi vs normal
  C qc_discrimination     complexity / apoptosis / HSP, mt-hi vs normal (not-debris check)
  D heat_honesty          WT_heat_up / _updown per group (secondary_percell; NOT a heat-signal claim)
  E donor_tissue          donor x tissue composition + Treg embedding with pocket highlighted

Correlative language only: the claim is that these are real, viable, identity-retaining
Tregs, legitimately RETAINED — not a fever/HIF/STING driver, not a rescue of Treg-preference.
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

from config import PATHS, PARAMS, POPULATION_COLORS, TISSUE_COLORS  # noqa: E402
from helpers.figure_style import set_paper_style, save_overview, FIG_CFG  # noqa: E402

STAGE = "01_qc"
SCRIPT = "02_analysis/scripts/01_qc_mthi_characterize_viz.py"

_OKABE = (FIG_CFG.get("colors", {}) or {}).get("okabe_ito", {}) or {}

# The mt-hi grouping is a second categorical axis, and the heat-honesty panel draws it in
# the same axes as the sorted populations, so its hues are taken from the categorical
# palette entries the population palette does not use. `normal_Treg` stays a neutral grey:
# it is the baseline the two pockets are read against, not a fourth category.
GROUP_COL = {"mt_hi_effector": _OKABE["vermillion"],
             "mt_hi_noneffector": _OKABE["blue"],
             "normal_Treg": "#999999"}
GROUP_LAB = {"mt_hi_effector": "mt-hi effector (cl6)", "mt_hi_noneffector": "mt-hi non-eff (cl16)",
             "normal_Treg": "normal Treg"}
# The one population palette, read from analysis_config.yaml::colors.populations.
POP_COL = POPULATION_COLORS


def _read(name: str) -> pd.DataFrame:
    return pd.read_csv(PATHS.tables(STAGE) / f"{name}.csv")


def main() -> None:
    set_paper_style(config=FIG_CFG)
    clu = _read("mthi_cluster_enrichment")
    ident = _read("mthi_identity_retention")
    qc = _read("mthi_qc_discrimination")
    heat = _read("mthi_heat_percell")
    comp = _read("mthi_donor_tissue")
    mem = pd.read_csv(PATHS.tables(STAGE) / "mthi_treg_membership.csv", index_col=0)

    # ==================================================================== #
    # A. cluster-resolved median %mt vs median eTreg + per-cell hexbin
    # ==================================================================== #
    fig, (axl, axr) = plt.subplots(1, 2, figsize=(13, 5.5))
    for _, r in clu.iterrows():
        eff = bool(r["is_mthi_effector"])
        mthi = bool(r["is_mthi_cluster"])
        col = (GROUP_COL["mt_hi_effector"] if eff
               else (GROUP_COL["mt_hi_noneffector"] if mthi else "#BBBBBB"))
        size = 40 + 260 * (r["n_treg"] / clu["n_treg"].max())
        axl.scatter(r["median_pct_mt"], r["median_score_eTreg"], s=size, c=col,
                    edgecolors="black", linewidths=0.6, zorder=3 if mthi else 2)
        if mthi:
            axl.annotate(f"cl{r['leiden']}", (r["median_pct_mt"], r["median_score_eTreg"]),
                         textcoords="offset points", xytext=(6, 6), fontsize=10, fontweight="bold")
    axl.axhline(0, color="grey", lw=0.7, ls="--")
    axl.set_xlabel("cluster median % mito")
    axl.set_ylabel("cluster median score_eTreg (EDA)")
    axl.set_title("Treg clusters: median %mt vs eTreg")

    hb = axr.hexbin(mem["pct_counts_mt"], mem["score_eTreg"], gridsize=45,
                    cmap="Greys", bins="log", mincnt=1)
    fig.colorbar(hb, ax=axr, shrink=0.8, label="log10 cells")
    eff_cells = mem[mem["mthi_group"] == "mt_hi_effector"]
    axr.scatter(eff_cells["pct_counts_mt"], eff_cells["score_eTreg"], s=10,
                c=GROUP_COL["mt_hi_effector"],
                linewidths=0, alpha=0.8, label="mt-hi effector (cl6)")
    axr.axhline(0, color="grey", lw=0.7, ls="--")
    axr.set_xlabel("% mito (per Treg cell)")
    axr.set_ylabel("score_eTreg (EDA)")
    axr.set_title("Per-Treg %mt vs eTreg (pocket = orange)")
    axr.legend(loc="upper right", fontsize=9, framealpha=0.8)
    fig.tight_layout()
    save_overview(
        fig, STAGE, "mthi_cluster_mt_etreg", table=clu,
        finding=("Two Treg leiden clusters carry high %mt (~20% vs ~4% rest): cl6 is the "
                 "effector-like pocket (eTreg-high, SF-restricted) and cl16 is mt-hi but "
                 "eTreg-low. The pocket is a discrete, reproducibly-defined region, not scattered noise."),
        script=SCRIPT, fn="main",
        config_kv=f"thresholds.qc_pct_mt_max = {PARAMS.qc_pct_mt_max}; qc_n_mads_mt = null (ceiling-only)",
        input="03_results/01_qc/tables/mthi_cluster_enrichment.csv + mthi_treg_membership.csv",
        how_to_read=("Left: one dot per Treg leiden cluster (size ~ n cells), x = median %mt, "
                     "y = median score_eTreg; orange = mt-hi effector pocket (cl6), purple = mt-hi "
                     "non-effector (cl16). Right: per-Treg %mt vs eTreg hexbin (log density), pocket "
                     "cells overlaid green. secondary_percell / EDA descriptive — not pseudobulk evidence."),
        config=FIG_CFG, wide=True)

    # ==================================================================== #
    # B. identity retention
    # ==================================================================== #
    genes = ident["gene"].tolist()
    x = np.arange(len(genes))
    w = 0.26
    fig, (axm, axf) = plt.subplots(1, 2, figsize=(13, 5.5))
    axm.bar(x - w, ident["median_mthi_eff"], w, color=GROUP_COL["mt_hi_effector"], label=GROUP_LAB["mt_hi_effector"])
    axm.bar(x, ident["median_mthi_noneff"], w, color=GROUP_COL["mt_hi_noneffector"], label=GROUP_LAB["mt_hi_noneffector"])
    axm.bar(x + w, ident["median_normal"], w, color=GROUP_COL["normal_Treg"], label=GROUP_LAB["normal_Treg"])
    axm.set_xticks(x); axm.set_xticklabels(genes)
    axm.set_ylabel("median log-norm expression")
    axm.set_title("Treg identity markers — median expression")
    axm.legend(fontsize=9)
    for i, r in ident.iterrows():
        star = "n.s." if r["fdr"] >= 0.05 else "*"
        axm.annotate(f"rbc={r['rank_biserial']:.2f}\n{star}", (i - w, r["median_mthi_eff"]),
                     textcoords="offset points", xytext=(0, 3), ha="center", fontsize=7)
    axf.bar(x - w, ident["frac_expr_mthi_eff"], w, color=GROUP_COL["mt_hi_effector"])
    axf.bar(x, ident["frac_expr_mthi_noneff"], w, color=GROUP_COL["mt_hi_noneffector"])
    axf.bar(x + w, ident["frac_expr_normal"], w, color=GROUP_COL["normal_Treg"])
    axf.set_xticks(x); axf.set_xticklabels(genes)
    axf.set_ylabel("fraction of cells expressing (>0)")
    axf.set_ylim(0, 1)
    axf.set_title("Fraction expressing")
    fig.suptitle("mt-hi effector Treg RETAINS Treg identity (IKZF2/CTLA4 up; FOXP3 modestly lower = depth)")
    fig.tight_layout()
    save_overview(
        fig, STAGE, "mthi_identity_retention", table=ident,
        finding=("The mt-hi effector pocket (cl6) retains Treg identity: IKZF2 (rbc +0.59) and "
                 "CTLA4 (+0.38) are UP vs normal Treg, IL2RA/TIGIT comparable; FOXP3 is modestly "
                 "lower (frac 0.55 vs 0.80) consistent with the lower sequencing depth of high-mito "
                 "cells, NOT loss of Treg identity. These are still Tregs."),
        script=SCRIPT, fn="main",
        config_kv="decisions.treg_gate.basis = sorting (refined by FOXP3/IL2RA/CTLA4/IKZF2)",
        input="03_results/01_qc/tables/mthi_identity_retention.csv",
        how_to_read=("Grouped bars over 5 canonical Treg markers: mt-hi effector (orange), mt-hi "
                     "non-effector (purple), normal Treg (grey). Left = median log-norm expression "
                     "(rbc = rank-biserial vs normal Treg, Mann-Whitney BH-FDR; * = FDR<0.05); "
                     "right = fraction expressing. secondary_percell tier."),
        config=FIG_CFG, wide=True)

    # ==================================================================== #
    # C. dying / debris discrimination
    # ==================================================================== #
    fig, axes = plt.subplots(1, 4, figsize=(15, 4.5))
    floor = int(qc["qc_min_genes_floor"].iloc[0])
    for ax, (_, r) in zip(axes, qc.iterrows()):
        vals = [r["median_mthi_eff"], r["median_mthi_noneff"], r["median_normal"]]
        cols = [GROUP_COL["mt_hi_effector"], GROUP_COL["mt_hi_noneffector"], GROUP_COL["normal_Treg"]]
        ax.bar([0, 1, 2], vals, color=cols)
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(["eff", "non", "norm"], fontsize=9)
        arrow = "UP" if r["elevated_in_pocket"] else "DOWN"
        ax.set_title(f"{r['metric']}\nrbc={r['rank_biserial']:.2f} ({arrow})\n[{r['metric_tier']}]", fontsize=9)
        if r["metric"] == "n_genes_by_counts":
            ax.axhline(floor, color="red", ls="--", lw=1.2)
            ax.annotate(f"QC floor {floor}", (2, floor), fontsize=8, color="red",
                        va="bottom", ha="right")
    fig.suptitle("Not dying / not debris: complexity healthy (>>200-gene floor); apoptosis & HSP NOT elevated")
    fig.tight_layout()
    save_overview(
        fig, STAGE, "mthi_qc_discrimination", table=qc,
        finding=("The pocket is not debris/dying: median n_genes 1339 (>>200 QC floor), i.e. "
                 "real cells with lower depth (the expected corollary of high mito fraction). "
                 "score_apoptosis and score_HSP are LOWER, not higher, than normal Treg (rbc -0.09, "
                 "-0.22). Doublets: no cell in the pocket is flagged predicted_doublet, but "
                 "doublet_score was not populated this run (Scrublet gap — see reasoning note)."),
        script=SCRIPT, fn="main",
        config_kv=f"thresholds.qc_min_genes = {int(PARAMS.qc_min_genes)}; qc_pct_mt_max = {PARAMS.qc_pct_mt_max}",
        input="03_results/01_qc/tables/mthi_qc_discrimination.csv",
        how_to_read=("Four bars per metric across groups (eff = mt-hi effector cl6, non = cl16, "
                     "norm = normal Treg). rbc = rank-biserial vs normal Treg (Mann-Whitney, all BH-FDR<0.05). "
                     "Red dashed = 200-gene QC floor. n_genes/total_counts are objective QC; "
                     "score_apoptosis/HSP are Tier-3 hand markers (QC-descriptive, NOT evidence). "
                     "secondary_percell tier."),
        config=FIG_CFG, wide=True)

    # ==================================================================== #
    # D. heat honesty
    # ==================================================================== #
    order = ["mt_hi_effector", "mt_hi_noneffector", "normal_Treg", "Tcon", "CD8"]
    heat = heat.set_index("group").loc[order].reset_index()
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for ax, ch in zip(axes, ["WT_heat_up", "WT_heat_updown"]):
        y = np.arange(len(heat))
        med = heat[f"{ch}_median"]
        lo = med - heat[f"{ch}_q25"]
        hi = heat[f"{ch}_q75"] - med
        cols = [GROUP_COL.get(g, POP_COL.get(g, "#444444")) for g in heat["group"]]
        ax.errorbar(med, y, xerr=[lo, hi], fmt="none", ecolor="grey", capsize=3, zorder=1)
        ax.scatter(med, y, c=cols, s=90, edgecolors="black", linewidths=0.6, zorder=2)
        ax.axvline(0, color="grey", ls="--", lw=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels([GROUP_LAB.get(g, g) for g in heat["group"]])
        ax.invert_yaxis()
        ax.set_xlabel(f"{ch} (median; bars = IQR)")
        ax.set_title(ch)
    fig.suptitle("Heat honesty: WT_heat quiet in the pocket (balanced _updown flat) — NOT a heat-signal claim")
    fig.tight_layout()
    save_overview(
        fig, STAGE, "mthi_heat_honesty", table=heat,
        finding=("WT_heat is quiet in the pocket. The balanced WT_heat_updown is essentially flat "
                 "(mt-hi effector median -0.082 vs normal -0.092); WT_heat_up shifts up but is "
                 "confounded (co-varies with the effector/depth axis). This is secondary_percell and "
                 "must NOT be read as the pocket carrying the mouse 39 °C-derived signature. The "
                 "pseudobulk NES (Treg 2.53 / Tcon 2.59 / CD8 2.07, pan-T) was unchanged when these "
                 "high-mito cells were recovered."),
        script=SCRIPT, fn="main",
        config_kv="decisions.go_no_go.primary_signature = WT_heat (pseudobulk, NOT per-cell)",
        input="03_results/01_qc/tables/mthi_heat_percell.csv",
        how_to_read=("Point = group median, bar = IQR, for WT_heat_up (left) and balanced "
                     "WT_heat_updown (right), across the two mt-hi Treg groups, normal Treg, Tcon, CD8. "
                     "secondary_percell / EDA tier — descriptive only, never pooled with pseudobulk NES."),
        config=FIG_CFG, wide=True)

    # ==================================================================== #
    # E. donor & tissue composition + embedding
    # ==================================================================== #
    fig, (axb, axe) = plt.subplots(1, 2, figsize=(13, 5.5))
    donors = sorted(comp["donor"].unique())
    tissues = ["synovial_fluid", "peripheral_blood"]
    tcol = TISSUE_COLORS
    bottom = np.zeros(len(donors))
    for t in tissues:
        vals = [comp[(comp.donor == d) & (comp.tissue == t)]["n_cells"].sum() for d in donors]
        axb.bar(range(len(donors)), vals, bottom=bottom, color=tcol[t],
                label="SF" if t == "synovial_fluid" else "PB")
        bottom += np.array(vals, dtype=float)
    axb.axhline(int(comp["pseudobulk_min_cells"].iloc[0]), color="black", ls="--", lw=1.0,
                label=f"min_cells={int(comp['pseudobulk_min_cells'].iloc[0])}")
    axb.set_xticks(range(len(donors)))
    axb.set_xticklabels([d.replace("JIA_patient_", "p") for d in donors])
    axb.set_ylabel("cells in mt-hi effector pocket")
    sf_d = int(comp["sf_donors_meeting_floor"].iloc[0])
    pb_d = int(comp["pb_donors_meeting_floor"].iloc[0])
    axb.set_title("Pocket cells by donor x tissue")
    axb.text(0.02, 0.80, f"SF donors>=floor: {sf_d}\nPB donors>=floor: {pb_d}\n"
             f"(min_donors={int(comp['pseudobulk_min_donors'].iloc[0])})\n=> no own DE stratum",
             transform=axb.transAxes, fontsize=9, va="top",
             bbox=dict(boxstyle="round", fc="white", ec="grey", alpha=0.85))
    axb.legend(fontsize=9, loc="upper center")

    # embedding: all Tregs, pocket highlighted
    axe.scatter(mem["x"], mem["y"], s=2, c="#DDDDDD", linewidths=0, label="normal Treg")
    for grp in ["mt_hi_noneffector", "mt_hi_effector"]:
        m = mem["mthi_group"] == grp
        axe.scatter(mem.loc[m, "x"], mem.loc[m, "y"], s=6, c=GROUP_COL[grp],
                    linewidths=0, label=GROUP_LAB[grp])
    axe.set_xticks([]); axe.set_yticks([])
    axe.set_title("Treg UMAP — pocket shown, not its own DE stratum")
    axe.legend(fontsize=9, markerscale=3, loc="upper right")
    fig.suptitle("mt-hi effector pocket is SF-restricted & p6-dominated — retained in the main Treg pseudobulk, shown on the embedding")
    fig.tight_layout()
    save_overview(
        fig, STAGE, "mthi_donor_tissue", table=comp,
        finding=("The pocket is 97% synovial fluid and 69% one donor (p6); it clears min_cells "
                 "(20) in 3 SF donors but has an essentially empty PB arm (p6=5 cells, 0 donors "
                 "at floor), so it cannot support its own paired SF-vs-PB pseudobulk contrast. "
                 "It is therefore retained WITHIN the main SF-vs-PB Treg pseudobulk and shown on "
                 "the embedding, never carved out as its own DE stratum."),
        script=SCRIPT, fn="main",
        config_kv=f"thresholds.pseudobulk_min_cells = {int(PARAMS.pseudobulk_min_cells)}; "
                  f"pseudobulk_min_donors = {int(PARAMS.pseudobulk_min_donors)}",
        input="03_results/01_qc/tables/mthi_donor_tissue.csv + mthi_treg_membership.csv",
        how_to_read=("Left: stacked bars of pocket cells per donor (SF orange / PB blue); black "
                     "dashed = per-stratum min_cells. Right: all Tregs on the unsupervised UMAP with "
                     "the mt-hi effector (orange) and non-effector (purple) clusters highlighted. "
                     "secondary_percell tier; the UMAP is a usability lens, not biology."),
        config=FIG_CFG, wide=True)

    print("[mthi_viz] wrote 5 overview figures + source tables + README captions")


if __name__ == "__main__":
    main()
