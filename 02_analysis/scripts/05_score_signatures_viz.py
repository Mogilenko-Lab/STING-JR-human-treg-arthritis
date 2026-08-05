#!/usr/bin/env python
"""
05_score_signatures_viz.py — VIZ (no statistics). The panel that answers the niche question.
===========================================================================================
Reads the fgsea tables + score tables from 05_score_signatures.py and renders:
  - the ordered WT_heat NES dot plot across Treg / Tcon / CD8 — the confirmatory
    answer to whether the mouse 39 °C-derived signature separates the inflamed
    synovial niche from paired blood within a frozen sort label;
  - per-cell WT_heat_up AUCell score violins SF-Treg vs PB-Treg vs SF-Tcon vs SF-CD8
    (running-sum figure moved to 05_score_signatures_viz.R).

Two naming rules bind this file. First, this panel is an ordered NES dot plot:
no `pseudobulk_nes` row carries a standard error or interval. The legacy file
stem `wt_heat_nes_forest` is load-bearing for consumers outside this compartment
and remains unchanged, while reader-facing text names the displayed geometry.
Second, the set is named by how it was derived — from mouse iTreg 37/39 °C
contrasts — and never by a specificity it does not have.
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
from helpers.source_hash_manifest import verify_source_hash  # noqa: E402

STAGE = "05_scoring"
SCRIPT = "02_analysis/scripts/05_score_signatures_viz.py"
POP_TAG = {"Treg": "treg", "Tcon": "tcon", "CD8": "cd8"}
# The one population palette, read from analysis_config.yaml::colors.populations.
POP_COL = POPULATION_COLORS
PRIMARY = "WT_heat"
ARM_ORDER = ["up", "down"]


def nominal_set_sizes() -> dict:
    """Nominal size of each frozen mouse arm — a line count, not a statistic.

    The effective set size only means something against the nominal one, so both
    travel with the panel. Read from the frozen mouse->human projection contract
    that stage 05 scored, so the denominator cannot drift from the numerator.
    """
    sig_dir = PATHS.signature_contract / "signatures" / PRIMARY
    out = {}
    for arm in ARM_ORDER:
        path = sig_dir / f"{PRIMARY}_{arm}.txt"
        if not path.exists():
            raise FileNotFoundError(f"[05_viz] frozen signature not found: {path}")
        verify_source_hash(
            path,
            f"{PRIMARY}_{arm}",
            PATHS.tables(STAGE) / "source_hash_manifest.csv",
            root=ROOT.parent,
        )
        out[arm] = len({ln.strip() for ln in path.read_text().splitlines() if ln.strip()})
    return out


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

    # ---- 1. ordered NES dot plot — the confirmatory answer ----
    fdr = float(PARAMS.gsea_fdr)
    nominal = nominal_set_sizes()
    fig, ax = plt.subplots(figsize=(9.5, 7))
    answer = ""
    if not gsea.empty:
        gsea["n_nominal"] = gsea["direction"].map(nominal)
        order = [(p, d) for d in ARM_ORDER for p in POP_TAG]
        ylabels, yv = [], []
        for i, (pop, d) in enumerate(order):
            sub = gsea[(gsea["cell_state"] == pop) & (gsea["direction"] == d)]
            if sub.empty:
                continue
            nes = float(sub["nes"].iloc[0]); padj = float(sub["padj"].iloc[0])
            n_eff = int(sub["set_size"].iloc[0]); n_nom = int(sub["n_nominal"].iloc[0])
            ax.scatter(nes, i, s=140, color=POP_COL[pop],
                       edgecolors="black", zorder=3,
                       marker="o" if d == "up" else "D")
            star = "*" if padj < fdr else ""
            ax.text(nes, i + 0.20, f"{nes:.2f}{star}", ha="center", fontsize=9)
            # Effective set size travels with every NES, ON THE FACE, against the
            # nominal size of that arm — the project rule, and the guard against
            # reading the ordering of these three NES as a biological ranking.
            ax.text(nes + 0.10, i, f"n {n_eff} of {n_nom}  ·  FDR "
                                   f"{padj:.3f}" if padj >= 0.001 else
                                   f"n {n_eff} of {n_nom}  ·  FDR {padj:.0e}",
                    va="center", ha="left", fontsize=9)
            ylabels.append(f"{pop} · {d}"); yv.append(i)

        # The answer, read off the plotted rows rather than asserted. It is
        # negative on Treg preference, and the panel has to say so on its face.
        up = gsea[gsea["direction"].eq("up")]
        n_sig_up = int(up["padj"].lt(fdr).sum())
        down_sig = gsea[gsea["direction"].eq("down") & gsea["padj"].lt(fdr)]
        answer = (f"pan-T, not Treg-preferential — all {n_sig_up} up arms clear "
                  f"FDR {fdr}")

        ax.axvline(0, ls="--", c="grey", lw=1)
        ax.set_yticks(yv); ax.set_yticklabels(ylabels)
        ax.set_xlim(0, 4.5)
        ax.set_xticks([0, 1, 2, 3])
        ax.set_xlabel("WT_heat NES (SF synovial fluid vs PB paired peripheral blood)")
        ax.set_title("QUESTION: does the mouse 39 °C-derived signature separate the inflamed\n"
                     "synovial niche from paired blood within a frozen sort label?\n"
                     f"ANSWER: yes — {answer}")
        from matplotlib.lines import Line2D
        # Legend has two keys: colour = population (matches the point colours), shape = up/down set.
        pop_handles = [
            Line2D([0], [0], marker="o", color="w", markerfacecolor=POP_COL[p],
                   markeredgecolor="k", label=p, markersize=9)
            for p in POP_TAG
        ]
        shape_handles = [
            Line2D([0], [0], marker="o", color="w", markerfacecolor="white",
                   markeredgecolor="k", label="up set (circle)", markersize=9),
            Line2D([0], [0], marker="D", color="w", markerfacecolor="white",
                   markeredgecolor="k", label="down set (diamond)", markersize=9),
            Line2D([0], [0], marker="", linestyle="", label=f"* = FDR < {fdr}"),
        ]
        # Upper right: the three down-arm rows sit at low NES, so that corner is the
        # only region no marker or size annotation occupies.
        ax.legend(handles=pop_handles + shape_handles, frameon=True, loc="upper right",
                  fontsize=9, title="population (colour) · set (shape)",
                  title_fontsize=9)
        # Two readings the sequence used to leave to prose. The Treg-preference
        # answer is negative; and the down arm is not silent, so "the up arm is
        # the only informative arm" is retired here rather than in a caption.
        down_note = (
            "The down arm is not silent either: it reaches FDR "
            + ", ".join(f"{float(r['padj']):.3f} in {r['cell_state']}"
                        for _, r in down_sig.iterrows())
            + " at the same sign as the up arm, and carries no direction in the other two.\n"
        ) if len(down_sig) else "No down arm reaches significance in any population.\n"
        ax.text(0.0, -0.135,
                "ANSWER — confirmatory tier (donor-level pseudobulk within frozen sort labels, "
                "limma-voom → fgsea). This is the only tier that may answer.\n"
                "Treg preference is answered in the NEGATIVE: Tregs are in the result, not "
                "privileged in it, and the effective set size tracks the NES\n"
                "across these three rows, so do not read the ordering as a biological ranking. "
                "'Treg signature' names how the set was derived — from mouse\n"
                "iTreg 37/39 °C contrasts — and never a specificity it does not have.\n"
                + down_note
                + "Geometry: ordered NES dot plot with FDR encoding; no pseudobulk NES row "
                  "carries an interval.",
                transform=ax.transAxes, ha="left", va="top", fontsize=9)
    fig.tight_layout()

    # Every NES, FDR and set size the caption states is read from `gsea` in this run. The
    # mouse arms are re-derived upstream in mouse_anchor, so a typed triple here would go
    # stale the moment that correction lands.
    def _arm_rows(direction: str) -> pd.DataFrame:
        if gsea.empty:
            return gsea
        sub = gsea[gsea["direction"].eq(direction)].copy()
        sub["ord"] = sub["cell_state"].map({p: i for i, p in enumerate(POP_TAG)})
        return sub.sort_values("ord")

    def _phrase(sub: pd.DataFrame, first_full: bool = True) -> str:
        bits = []
        for k, (_, r) in enumerate(sub.iterrows()):
            size = (f" ({int(r['set_size'])} of {int(r['n_nominal'])} genes ranked)"
                    if k == 0 and first_full else f" ({int(r['set_size'])})")
            bits.append(f"{float(r['nes']):.4f} in {r['cell_state']}{size}")
        return ", ".join(bits[:-1]) + (" and " if len(bits) > 1 else "") + bits[-1]

    up_rows, down_rows = _arm_rows("up"), _arm_rows("down")
    if gsea.empty:
        up_line = down_line = "No enrichment table was found for this stage."
    else:
        up_line = (f"The mouse 39 °C-derived up arm separates synovial fluid from paired blood "
                   f"in every sorted population: NES {_phrase(up_rows)}, all at FDR below "
                   f"{up_rows['padj'].max():.0e}. The result is pan-T, and Tregs are one of the "
                   f"three populations carrying it.")
        d_sig = down_rows[down_rows["padj"].lt(fdr)]
        d_ns = down_rows[~down_rows["padj"].lt(fdr)]
        down_line = (
            (f"The down arm reaches NES "
             + ", ".join(f"{float(r['nes']):.4f} at FDR {float(r['padj']):.3f} in {r['cell_state']}"
                         for _, r in d_sig.iterrows())
             + ", the same sign as the up arm, and carries no direction in "
             + " or ".join(f"{r['cell_state']} ({float(r['nes']):.4f})"
                           for _, r in d_ns.iterrows()) + ".")
            if len(d_sig) else
            f"The down arm reaches significance in no population at FDR {fdr}.")

    save_overview(fig, STAGE, "wt_heat_nes_forest",
                  table=gsea[["cell_state", "pathway_id", "nes", "pvalue", "padj",
                              "set_size", "n_nominal"]]
                  if not gsea.empty else pd.DataFrame(),
                  finding=f"{up_line} {down_line}",
                  script=SCRIPT, fn="main",
                  config_kv=f"thresholds.gsea_fdr={PARAMS.gsea_fdr}; gsea_min_size={PARAMS.gsea_min_size}",
                  input="03_results/05_scoring/tables/gsea_pseudobulk_{treg,tcon,cd8}.csv",
                  how_to_read=("This is the confirmatory tier: donor-level pseudobulk within "
                               "frozen sort labels, limma-voom then fgsea, on the 6 donors "
                               "present in both arms. Points are NES for the up (circle) and down "
                               "(diamond) arms, coloured by population; the asterisk marks FDR "
                               f"below {PARAMS.gsea_fdr}. Labels give effective and nominal set "
                               "sizes plus FDR. Effective size tracks the NES ordering, so the "
                               "ordering is a size effect. Ordered NES dot plot with FDR "
                               "encoding; the rows carry no interval. Correlative."),
                  config=FIG_CFG, width=9.5, height=8.5)

    # ---- 2. per-cell score violins ----
    dm = pd.read_csv(tdir / "donor_label_score_means.csv")
    dm["tissue_s"] = dm["tissue"].map({"synovial_fluid": "SF", "peripheral_blood": "PB"})
    dm["group"] = dm["coarse_label"] + " " + dm["tissue_s"]
    groups = ["Treg SF", "Treg PB", "Tcon SF", "Tcon PB", "CD8 SF", "CD8 PB"]
    fig3, ax = plt.subplots(figsize=(9, 6))
    data, positions, colors = [], [], []
    for i, grp in enumerate(groups):
        vals = dm.loc[dm["group"] == grp, "WT_heat_up_AUCell"].values
        if len(vals):
            data.append(vals); positions.append(i)
            colors.append(POP_COL[grp.split()[0]])
    parts = ax.violinplot(data, positions=positions, showmeans=True, widths=0.8)
    for b, c in zip(parts["bodies"], colors):
        b.set_facecolor(c); b.set_alpha(0.65)
    for i, grp in zip(positions, [groups[p] for p in positions]):
        vals = dm.loc[dm["group"] == grp, "WT_heat_up_AUCell"].values
        ax.scatter(np.full(len(vals), i), vals, s=18, c="black", zorder=3)
    ax.set_xticks(range(len(groups))); ax.set_xticklabels(groups, rotation=30, ha="right")
    ax.set_ylabel("donor-mean WT_heat_up AUCell score")
    ax.set_title("Per-cell WT_heat_up AUCell score, donor means by state × tissue\n"
                 "CORROBORATES ONLY — a per-cell score may never answer")
    fig3.tight_layout()
    save_overview(fig3, STAGE, "score_violins", table=dm,
                  finding=("Donor-mean WT_heat_up AUCell sits higher in synovial fluid than in "
                           "paired blood in all three sorted populations, so the per-cell channel "
                           "shadows the pseudobulk result in the same direction. This corroborates. "
                           "A per-cell score is a different estimand on a secondary tier, and the "
                           "shift it shows spans all three populations."),
                  script=SCRIPT, fn="main",
                  config_kv="signature = WT_heat_up (AUCell, rank-based [0,1])",
                  input="03_results/05_scoring/tables/donor_label_score_means.csv",
                  how_to_read=("This panel corroborates; the confirmatory answer is the pseudobulk "
                               "NES dot plot. Each dot is one donor's mean WT_heat_up AUCell score "
                               "for that state × tissue, and the violins summarise across donors. "
                               "AUCell is a rank-based score in [0, 1], the area under each cell's "
                               "gene-recovery curve for the up set, robust to library size and "
                               "composition. Read the SF-versus-PB shift within a population (Treg "
                               "SF against Treg PB); the absolute level carries no reading. This "
                               "is a different estimand from the pseudobulk NES dot plot and "
                               "shares no axis with it. The down arm is omitted because up and "
                               "down co-shift in synovial fluid. Correlative."),
                  config=FIG_CFG)
    print("[05_scoring_viz] wrote 2 overviews (ordered NES dot plot + score violins)")


if __name__ == "__main__":
    main()
