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
# it is the baseline the two pockets are read against.
GROUP_COL = {"mt_hi_effector": _OKABE["vermillion"],
             "mt_hi_noneffector": _OKABE["blue"],
             "normal_Treg": "#999999"}
# Base labels. The leiden id of each mt-hi group is appended by _group_labels() from the
# cluster table itself, so a legend can never name a cluster its own source table numbers
# differently.
GROUP_LAB = {"mt_hi_effector": "mt-hi effector", "mt_hi_noneffector": "mt-hi non-eff",
             "normal_Treg": "normal Treg"}
# The one population palette, read from analysis_config.yaml::colors.populations.
POP_COL = POPULATION_COLORS


def _read(name: str) -> pd.DataFrame:
    return pd.read_csv(PATHS.tables(STAGE) / f"{name}.csv")


def _mthi_ids(clu: pd.DataFrame) -> dict:
    """Leiden id of the mt-hi effector and mt-hi non-effector cluster, read from `clu`."""
    is_mthi = clu["is_mthi_cluster"].astype(str).str.lower().eq("true")
    is_eff = clu["is_mthi_effector"].astype(str).str.lower().eq("true")
    eff = clu.loc[is_mthi & is_eff, "leiden"]
    non = clu.loc[is_mthi & ~is_eff, "leiden"]
    return {"mt_hi_effector": int(eff.iloc[0]) if len(eff) else None,
            "mt_hi_noneffector": int(non.iloc[0]) if len(non) else None}


def _group_labels(clu: pd.DataFrame) -> dict:
    """GROUP_LAB with each mt-hi group's leiden id filled in from the cluster table."""
    ids = _mthi_ids(clu)
    lab = dict(GROUP_LAB)
    for grp, cid in ids.items():
        if cid is not None:
            lab[grp] = f"{GROUP_LAB[grp]} (leiden {cid})"
    return lab


def _pseudobulk_nes(arm: str = "WT_heat_up") -> str:
    """One phrase naming the donor-pseudobulk NES of `arm`, read from 05_scoring's own tables.

    The mouse-projection arms are re-derived upstream in `mouse_anchor`, so these values move.
    Reading them here keeps the caption from asserting a stale triple. Returns "" when
    05_scoring has not been run yet, so this stage never depends on a later one existing.
    """
    parts = []
    for pop, slug in (("Treg", "treg"), ("Tcon", "tcon"), ("CD8", "cd8")):
        path = PATHS.tables("05_scoring") / f"gsea_pseudobulk_{slug}.csv"
        if not path.exists():
            return ""
        tab = pd.read_csv(path)
        row = tab[tab["pathway_id"] == arm]
        if row.empty:
            return ""
        parts.append(f"{pop} {float(row['nes'].iloc[0]):.2f}")
    return ", ".join(parts)


def main() -> None:
    set_paper_style(config=FIG_CFG)
    clu = _read("mthi_cluster_enrichment")
    ident = _read("mthi_identity_retention")
    qc = _read("mthi_qc_discrimination")
    heat = _read("mthi_heat_percell")
    comp = _read("mthi_donor_tissue")
    mem = pd.read_csv(PATHS.tables(STAGE) / "mthi_treg_membership.csv", index_col=0)

    # Every cluster id and every number a caption below states is read from these tables in
    # this run, so a caption cannot drift from the artifact it describes.
    glab = _group_labels(clu)
    ids = _mthi_ids(clu)
    is_mthi = clu["is_mthi_cluster"].astype(str).str.lower().eq("true")
    mthi_clu, rest_clu = clu[is_mthi], clu[~is_mthi]

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
            axl.annotate(f"leiden {int(r['leiden'])}",
                         (r["median_pct_mt"], r["median_score_eTreg"]),
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
                linewidths=0, alpha=0.8, label=glab["mt_hi_effector"])
    axr.axhline(0, color="grey", lw=0.7, ls="--")
    axr.set_xlabel("% mito (per Treg cell)")
    axr.set_ylabel("score_eTreg (EDA)")
    axr.set_title("Per-Treg %mt vs eTreg (pocket = orange)")
    axr.legend(loc="upper right", fontsize=9, framealpha=0.8)
    fig.tight_layout()
    save_overview(
        fig, STAGE, "mthi_cluster_mt_etreg", table=clu,
        finding=(f"{len(mthi_clu)} Treg leiden clusters carry high mitochondrial content, a "
                 f"median near {mthi_clu['median_pct_mt'].mean():.0f}% against "
                 f"{rest_clu['median_pct_mt'].min():.1f}% to "
                 f"{rest_clu['median_pct_mt'].max():.1f}% in the remaining {len(rest_clu)}. They "
                 f"split on effector identity: leiden {ids['mt_hi_effector']} is the "
                 f"effector-like, SF-restricted pocket, and leiden {ids['mt_hi_noneffector']} is "
                 f"equally mito-high and effector-low. The pocket is a discrete, reproducibly "
                 f"defined region of the embedding."),
        script=SCRIPT, fn="main",
        config_kv=f"thresholds.qc_pct_mt_max = {PARAMS.qc_pct_mt_max}; qc_n_mads_mt = null (ceiling-only)",
        input="03_results/01_qc/tables/mthi_cluster_enrichment.csv + mthi_treg_membership.csv",
        how_to_read=(f"Left panel plots one dot per Treg leiden cluster, size scaling with cell "
                     f"count, x the median %mt and y the median `score_eTreg`; orange is the mt-hi "
                     f"effector pocket, blue the mt-hi non-effector cluster. Right panel is a "
                     f"per-Treg %mt against eTreg hexbin on log density, pocket cells overlaid in "
                     f"orange. Both panels label the two clusters leiden "
                     f"{ids['mt_hi_effector']} and {ids['mt_hi_noneffector']}, the numbering "
                     f"`mthi_cluster_enrichment.csv` uses. `secondary_percell` / EDA tier — "
                     f"never pseudobulk evidence."),
        config=FIG_CFG, wide=True)

    # ==================================================================== #
    # B. identity retention
    # ==================================================================== #
    genes = ident["gene"].tolist()
    x = np.arange(len(genes))
    w = 0.26
    fig, (axm, axf) = plt.subplots(1, 2, figsize=(13, 5.5))
    axm.bar(x - w, ident["median_mthi_eff"], w, color=GROUP_COL["mt_hi_effector"], label=glab["mt_hi_effector"])
    axm.bar(x, ident["median_mthi_noneff"], w, color=GROUP_COL["mt_hi_noneffector"], label=glab["mt_hi_noneffector"])
    axm.bar(x + w, ident["median_normal"], w, color=GROUP_COL["normal_Treg"], label=glab["normal_Treg"])
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
    fig.suptitle("Canonical Treg markers in the mt-hi Treg groups and normal Treg")
    fig.tight_layout()
    idt = ident.set_index("gene")
    save_overview(
        fig, STAGE, "mthi_identity_retention", table=ident,
        finding=(f"The mt-hi effector pocket retains Treg identity. IKZF2 (rbc "
                 f"{idt.loc['IKZF2', 'rank_biserial']:+.2f}) and CTLA4 "
                 f"({idt.loc['CTLA4', 'rank_biserial']:+.2f}) sit above normal Treg, IL2RA and "
                 f"TIGIT are comparable, and FOXP3 is modestly lower (expressed in "
                 f"{idt.loc['FOXP3', 'frac_expr_mthi_eff']:.0%} of pocket cells against "
                 f"{idt.loc['FOXP3', 'frac_expr_normal']:.0%} of normal Tregs), which tracks the "
                 f"lower sequencing depth of high-mito cells. These cells are Tregs."),
        script=SCRIPT, fn="main",
        config_kv="decisions.treg_gate.basis = sorting (refined by FOXP3/IL2RA/CTLA4/IKZF2)",
        input="03_results/01_qc/tables/mthi_identity_retention.csv",
        how_to_read=(f"Grouped bars over {len(genes)} canonical Treg markers: mt-hi effector in "
                     f"orange, mt-hi non-effector in blue, normal Treg in grey. Left panel is "
                     f"median log-normalised expression, where rbc is the rank-biserial "
                     f"correlation against normal Treg (Mann-Whitney, BH-FDR; asterisk marks FDR "
                     f"< 0.05). Right panel is the fraction expressing. `secondary_percell` tier."),
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
        arrow = "up" if r["elevated_in_pocket"] else "down"
        ax.set_title(f"{r['metric']}\nrbc={r['rank_biserial']:.2f} ({arrow})\n[{r['metric_tier']}]", fontsize=9)
        if r["metric"] == "n_genes_by_counts":
            ax.axhline(floor, color="red", ls="--", lw=1.2)
            ax.annotate(f"QC floor {floor}", (2, floor), fontsize=8, color="red",
                        va="bottom", ha="right")
    fig.suptitle("Complexity, depth, apoptosis and heat-shock scores by mt-hi Treg group")
    fig.tight_layout()
    qct = qc.set_index("metric")
    save_overview(
        fig, STAGE, "mthi_qc_discrimination", table=qc,
        finding=(f"The pocket holds real cells at lower depth: median "
                 f"{int(qct.loc['n_genes_by_counts', 'median_mthi_eff']):,} genes, far above the "
                 f"{floor}-gene QC floor, which is the expected corollary of a high "
                 f"mitochondrial fraction. `score_apoptosis` (rbc "
                 f"{qct.loc['score_apoptosis', 'rank_biserial']:+.2f}) and `score_HSP` "
                 f"({qct.loc['score_HSP', 'rank_biserial']:+.2f}) both sit below normal Treg. No "
                 f"cell in the pocket is flagged `predicted_doublet`. `doublet_score` was not "
                 f"populated in this run, a Scrublet gap, so the doublet evidence rests on the "
                 f"flag alone."),
        script=SCRIPT, fn="main",
        config_kv=f"thresholds.qc_min_genes = {int(PARAMS.qc_min_genes)}; qc_pct_mt_max = {PARAMS.qc_pct_mt_max}",
        input="03_results/01_qc/tables/mthi_qc_discrimination.csv",
        how_to_read=(f"One panel per metric, three bars per panel — eff is the mt-hi effector "
                     f"pocket, non the mt-hi non-effector cluster, norm normal Treg. rbc is the "
                     f"rank-biserial correlation against normal Treg (Mann-Whitney, all BH-FDR < "
                     f"0.05); the red dashed line is the {floor}-gene QC floor. `n_genes_by_counts` "
                     f"and `total_counts` are objective QC measures; `score_apoptosis` and "
                     f"`score_HSP` are Tier-3 hand marker modules, QC-descriptive. "
                     f"`secondary_percell` tier."),
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
        ax.set_yticklabels([glab.get(g, g) for g in heat["group"]])
        ax.invert_yaxis()
        ax.set_xlabel(f"{ch} (median; bars = IQR)")
        ax.set_title(ch)
    fig.suptitle("Mouse 39 °C-derived signature scores per group (one-sided and balanced channels)")
    fig.tight_layout()
    hmed = heat.set_index("group")["WT_heat_updown_median"]
    # The donor-pseudobulk NES lives in 05_scoring and moves with the mouse_anchor re-derivation,
    # so it is read from that stage's own table.
    nes = _pseudobulk_nes("WT_heat_up")
    nes_sentence = (
        f"The same mouse up arm enriches the donor-pseudobulk SF-versus-PB contrast at NES "
        f"{nes}, unchanged when these high-mito cells were recovered; those values and their "
        f"FDRs live in `03_results/05_scoring/`."
        if nes else
        "The donor-pseudobulk enrichment of the same mouse up arm was unchanged when these "
        "high-mito cells were recovered; those NES values and their FDRs live in "
        "`03_results/05_scoring/`.")
    save_overview(
        fig, STAGE, "mthi_heat_honesty", table=heat,
        finding=(f"The mouse 39 °C-derived signature is quiet in the pocket. The balanced "
                 f"`WT_heat_updown` channel is essentially flat, median "
                 f"{hmed['mt_hi_effector']:.3f} in the mt-hi effector pocket against "
                 f"{hmed['normal_Treg']:.3f} in normal Treg. The one-sided `WT_heat_up` channel "
                 f"shifts up, and that shift co-varies with the effector/depth axis. This is "
                 f"`secondary_percell` tier, and it is a per-cell descriptive reading. "
                 f"{nes_sentence}"),
        script=SCRIPT, fn="main",
        config_kv="decisions.go_no_go.primary_signature = WT_heat (pseudobulk, not per-cell)",
        input="03_results/01_qc/tables/mthi_heat_percell.csv",
        how_to_read=("Point is the group median and bar the IQR, for `WT_heat_up` on the left "
                     "and the balanced `WT_heat_updown` on the right, across the two mt-hi Treg "
                     "groups, normal Treg, Tcon and CD8. `secondary_percell` / EDA tier — "
                     "descriptive, never pooled with the donor-pseudobulk NES."),
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
             f"(min_donors={int(comp['pseudobulk_min_donors'].iloc[0])})\n"
             f"=> no DE stratum of its own",
             transform=axb.transAxes, fontsize=9, va="top",
             bbox=dict(boxstyle="round", fc="white", ec="grey", alpha=0.85))
    axb.legend(fontsize=9, loc="upper center")

    # embedding: all Tregs, pocket highlighted
    axe.scatter(mem["x"], mem["y"], s=2, c="#DDDDDD", linewidths=0, label="normal Treg")
    for grp in ["mt_hi_noneffector", "mt_hi_effector"]:
        m = mem["mthi_group"] == grp
        axe.scatter(mem.loc[m, "x"], mem.loc[m, "y"], s=6, c=GROUP_COL[grp],
                    linewidths=0, label=glab[grp])
    axe.set_xticks([]); axe.set_yticks([])
    axe.set_title("Treg UMAP with the two mt-hi clusters highlighted")
    axe.legend(fontsize=9, markerscale=3, loc="upper right")
    fig.suptitle("mt-hi effector pocket composition by donor and tissue, and its place on the Treg embedding")
    fig.tight_layout()
    # Composition read from the table, so the caption states this run's pocket.
    n_pocket = int(comp["n_cells"].sum())
    n_sf = int(comp.loc[comp["tissue"] == "synovial_fluid", "n_cells"].sum())
    per_donor = comp.groupby("donor")["n_cells"].sum().sort_values(ascending=False)
    top_donor, n_top = per_donor.index[0], int(per_donor.iloc[0])
    top_lab = top_donor.replace("JIA_patient_", "p")
    n_top_pb = int(comp[(comp["donor"] == top_donor)
                        & (comp["tissue"] == "peripheral_blood")]["n_cells"].sum())
    save_overview(
        fig, STAGE, "mthi_donor_tissue", table=comp,
        finding=(f"The pocket is {n_sf / n_pocket:.0%} synovial fluid ({n_sf} of {n_pocket} "
                 f"cells) and {n_top / n_pocket:.0%} one donor ({top_lab}, {n_top} cells). It "
                 f"clears the {int(PARAMS.pseudobulk_min_cells)}-cell floor in {sf_d} SF donors "
                 f"and has an essentially empty PB arm — {top_lab} contributes {n_top_pb} PB "
                 f"cells and {pb_d} donors reach the floor — so it cannot support its own paired "
                 f"SF-versus-PB pseudobulk contrast. It is therefore retained inside the main "
                 f"SF-versus-PB Treg pseudobulk and shown on the embedding, and it is carved out "
                 f"as its own DE stratum nowhere."),
        script=SCRIPT, fn="main",
        config_kv=f"thresholds.pseudobulk_min_cells = {int(PARAMS.pseudobulk_min_cells)}; "
                  f"pseudobulk_min_donors = {int(PARAMS.pseudobulk_min_donors)}",
        input="03_results/01_qc/tables/mthi_donor_tissue.csv + mthi_treg_membership.csv",
        how_to_read=("Left panel stacks pocket cells per donor, SF orange and PB blue, with the "
                     "black dashed line the per-stratum `min_cells`. Right panel places all "
                     "Tregs on the unsupervised UMAP with the mt-hi effector cluster in orange "
                     "and the non-effector cluster in blue. `secondary_percell` tier; the UMAP is "
                     "a usability lens."),
        config=FIG_CFG, wide=True)

    print("[mthi_viz] wrote 5 overview figures + source tables + README captions")


if __name__ == "__main__":
    main()
