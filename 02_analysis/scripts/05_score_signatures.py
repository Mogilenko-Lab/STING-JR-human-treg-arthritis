#!/usr/bin/env python
"""
05_score_signatures.py — COMPUTE (no plotting). THE go/no-go readout.
=====================================================================
Primary evidence: pre-ranked fgsea of the frozen mouse `WT_heat` up/down sets
against the donor-pseudobulk SF-vs-PB ranked lists (stage 03), per population ->
NES + FDR. Treg is the money number; Tcon + CD8 are the Treg-specificity control.

Secondary (corroborative) evidence: per-cell AUCell + UCell (rank-based,
composition/depth-robust) of the same sets, summarised to donor x label means, with
an SF-vs-PB standardized mean difference (SMD) at the donor level on the AUCell
`WT_heat_up` score. Per the effect-size contract, the NES (primary_pseudobulk) and
the SMD (secondary_percell) are NEVER pooled.

Deferred (until go = yes): comparators KO_heat/Interaction, Tier-1 MSigDB battery,
eTreg correlation, CoReSh, pathway-explorer.

Outputs:
  03_results/05_scoring/tables/gsea_pseudobulk_{treg,tcon,cd8}.csv   (master_gsea schema NES)
  03_results/05_scoring/tables/gsea_pseudobulk_{treg,tcon,cd8}.rds   (clusterProfiler gseaResult S4)
  03_results/05_scoring/tables/runsum_interactive_{treg,tcon,cd8}_WT_heat_{up,down}.csv
                                                    (report interactive-widget substrate)
  03_results/05_scoring/tables/per_cell_scores.csv
  03_results/05_scoring/tables/donor_label_score_means.csv
  03_results/master/effect_sizes_treg_arthritis.csv   (common_effectsizes schema)
  03_results/master/master_effect_sizes.csv           (mirror)

GSEA engine: clusterProfiler::GSEA(by="fgsea") via helpers/fgsea_prerank.R — a real
gseaResult S4 the RNAseq-toolkit running-sum plotter consumes (migration 2026-07-11).
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "02_analysis"))
os.chdir(ROOT)

from config import (CONFIG, PATHS, PARAMS, TISSUE_KEY, DONOR_KEY, TISSUE_NUM,  # noqa: E402
                    TISSUE_DEN, COARSE_LABEL)
from helpers.geneset_utils import (load_alias_map, load_signature,  # noqa: E402
                                   score_cells_aucell_ucell, _symbol_to_varname)
from helpers.figure_style import append_master_table, FIG_CFG  # noqa: E402
from helpers.source_hash_manifest import verify_source_hashes  # noqa: E402

STAGE = "05_scoring"
PRIMARY = "WT_heat"
POP_TAG = {"Treg": "treg", "Tcon": "tcon", "CD8": "cd8"}
FGSEA_R = "02_analysis/helpers/fgsea_prerank.R"
DATASET = "GSE160097"
# This matrix carries hg19-era HGNC symbols; the frozen mouse arms ship current ones. The
# map is the committed resolution between the two, and both the fgsea call and the per-cell
# scoring read it, so the primary and secondary tiers see the same set membership.
ALIAS_MAP_PATH = CONFIG["symbol_alias"]["map_path"]


def effect_size_signoff_state() -> str:
    """Read the owner-controlled integration gate; an absent/invalid key is an error."""
    state = CONFIG.get("decisions", {}).get("effect_sizes", {}).get("signoff")
    if state not in {"pending", "signed_off"}:
        raise ValueError(
            "decisions.effect_sizes.signoff must be recorded as pending or signed_off"
        )
    return state


def run_fgsea(ranked_path: Path, out_csv: Path, contrast: str, sig_dir: Path) -> pd.DataFrame:
    # The alias map travels with the call, so each set is matched against this matrix's own
    # symbol vintage instead of by exact string against current HGNC. It only ever adds:
    # WT_heat_up gains nothing, WT_heat_down gains CRACR2A -> EFCAB4B and CYSRT1 -> C9orf169.
    cmd = [
        "Rscript", FGSEA_R, str(ranked_path), str(out_csv), contrast,
        str(PARAMS.gsea_min_size), str(PARAMS.gsea_max_size), str(PARAMS.gsea_seed),
        str(PARAMS.gsea_nperm), f"--alias-map={ALIAS_MAP_PATH}",
        f"{PRIMARY}_up:mouse_projection={sig_dir / f'{PRIMARY}_up.txt'}",
        f"{PRIMARY}_down:mouse_projection={sig_dir / f'{PRIMARY}_down.txt'}",
    ]
    subprocess.run(cmd, check=True)
    return pd.read_csv(out_csv)


def smd(sf: np.ndarray, pb: np.ndarray) -> dict:
    """Unpaired standardized mean difference (Cohen's d) on donor-level means + 95% CI."""
    n1, n2 = len(sf), len(pb)
    if n1 < 2 or n2 < 2:
        return dict(estimate=np.nan, se=np.nan, ci_low=np.nan, ci_high=np.nan, pvalue=np.nan)
    s1, s2 = sf.std(ddof=1), pb.std(ddof=1)
    pooled = np.sqrt(((n1 - 1) * s1 ** 2 + (n2 - 1) * s2 ** 2) / (n1 + n2 - 2))
    d = (sf.mean() - pb.mean()) / (pooled + 1e-12)
    se = np.sqrt((n1 + n2) / (n1 * n2) + d ** 2 / (2 * (n1 + n2 - 2)))
    from scipy import stats
    t = d / (se + 1e-12)
    p = 2 * stats.t.sf(abs(t), df=n1 + n2 - 2)
    return dict(estimate=float(d), se=float(se), ci_low=float(d - 1.96 * se),
                ci_high=float(d + 1.96 * se), pvalue=float(p))


def main() -> None:
    adata = sc.read_h5ad(PATHS.object("02_annotation"))
    sig_dir = PATHS.signature_contract / "signatures" / PRIMARY
    verify_source_hashes(
        PATHS.tables(STAGE) / "source_hash_manifest.csv",
        [
            (f"{PRIMARY}_up", sig_dir / f"{PRIMARY}_up.txt"),
            (f"{PRIMARY}_down", sig_dir / f"{PRIMARY}_down.txt"),
            (f"{PRIMARY}_ranked", sig_dir / f"{PRIMARY}_ranked.rnk"),
        ],
        root=ROOT.parent,
    )
    # The arms arrive in current HGNC symbols and this object's var carries the hg19-era
    # vintage, so they are resolved into the object's own vocabulary before anything is
    # scored. The applied pairs are printed rather than silently absorbed: a coverage count
    # that moved needs a reason attached to it.
    alias_map = load_alias_map(ALIAS_MAP_PATH)
    sig = load_signature(PATHS.signature_contract, PRIMARY, alias_map,
                         set(_symbol_to_varname(adata, "gene_symbol")))
    print(f"[05_scoring] {PRIMARY}: {len(sig['up'])} up / {len(sig['down'])} down genes")
    for direction, pairs in sig["alias_applied"].items():
        if pairs:
            print(f"[05_scoring] {PRIMARY}_{direction}: +{len(pairs)} via alias "
                  f"({' '.join(f'{a}->{b}' for a, b in pairs)})")

    # --- per-cell scores (secondary tier): AUCell + UCell on lognorm X ---
    # Rank-based, composition/depth-robust; replaces the mean-centred scanpy score_genes.
    # Up/down sets are scored SEPARATELY (AUCell/UCell are unsigned single-list scorers).
    # AUCell is canonical for the secondary SMD; UCell rides alongside as a cross-check.
    AUC_UP = f"{PRIMARY}_up_AUCell"
    gene_sets = {f"{PRIMARY}_up": list(sig["up"]), f"{PRIMARY}_down": list(sig["down"])}
    scores = score_cells_aucell_ucell(
        adata, gene_sets, layer=None, symbol_col="gene_symbol",
        n_cores=int(PARAMS.get("percell_score_ncores", 4)))
    score_cols = list(scores.columns)   # WT_heat_{up,down}_{AUCell,UCell}
    for c in score_cols:
        adata.obs[c] = scores[c].to_numpy()

    sym_to_var = _symbol_to_varname(adata, "gene_symbol")
    n_up = int(sum(s in sym_to_var for s in sig["up"]))
    n_down = int(sum(s in sym_to_var for s in sig["down"]))
    print(f"[05_scoring] per-cell coverage: up {n_up}/{len(sig['up'])}, "
          f"down {n_down}/{len(sig['down'])}")

    tdir = PATHS.tables(STAGE)
    per_cell = adata.obs[[DONOR_KEY, TISSUE_KEY, "coarse_label"] + score_cols].copy()
    per_cell.to_csv(tdir / "per_cell_scores.csv")

    donor_means = (per_cell.groupby([DONOR_KEY, TISSUE_KEY, "coarse_label"], observed=True)
                   [score_cols].mean().reset_index())
    donor_means["n_cells"] = (per_cell.groupby([DONOR_KEY, TISSUE_KEY, "coarse_label"], observed=True)
                              .size().values)
    donor_means.to_csv(tdir / "donor_label_score_means.csv", index=False)

    ranked_dir = PATHS.tables("03_pseudobulk")
    eff_rows = []
    for pop, tag in POP_TAG.items():
        n_cells = int((adata.obs["coarse_label"].astype(str) == pop).sum())

        # ---- PRIMARY: fgsea NES on the pseudobulk ranked list ----
        ranked_path = ranked_dir / f"ranked_{tag}.tsv"
        if ranked_path.exists():
            gsea = run_fgsea(ranked_path, tdir / f"gsea_pseudobulk_{tag}.csv",
                             f"SF_vs_PB_{pop}", sig_dir)
            # The frame arrives from run_fgsea carrying database="mouse_projection", the
            # gene-set COLLECTION name. append_master_table treats `database` as the key
            # naming the rows THIS CALL owns, and dedupes by comparing that argument
            # against the column. While the two disagreed the filter matched nothing, so
            # every run appended a fresh copy instead of replacing the previous one and
            # the accumulator grew without bound. The column is therefore set to the
            # per-call owner key before appending, which is the convention 10_hsr_lens.py
            # already uses on this same accumulator. The collection name is kept in this
            # stage's own gsea_pseudobulk_<tag>.csv, and `pathway_id` plus `stage`
            # identify these rows in the master table.
            # The key must be UNIQUE PER CALL: this loop appends once per population, so a
            # key shared across the three calls (either "mouse_projection" for all, or the
            # `stage` value) would have each call delete the previous population's rows.
            gsea_master = gsea.copy()
            gsea_master["database"] = f"WT_heat_{tag}"
            append_master_table(gsea_master, database=f"WT_heat_{tag}", stage=STAGE,
                                name="master_gsea_pseudobulk", config=FIG_CFG)
            # donors contributing to this pop's pseudobulk contrast
            dm = donor_means[donor_means["coarse_label"] == pop]
            n_don_sf = dm[dm[TISSUE_KEY] == TISSUE_NUM][DONOR_KEY].nunique()
            n_don_pb = dm[dm[TISSUE_KEY] == TISSUE_DEN][DONOR_KEY].nunique()
            for _, r in gsea.iterrows():
                nes = r["nes"]
                eff_rows.append(dict(
                    dataset=DATASET, signature=r["pathway_id"], cell_state=pop,
                    contrast="SF_vs_PB", effect_metric="pseudobulk_nes",
                    evidence_tier="primary_pseudobulk", estimate=nes,
                    se=np.nan, ci_low=np.nan, ci_high=np.nan,
                    direction=("up" if (pd.notna(nes) and nes > 0) else "down" if pd.notna(nes) else "ns"),
                    pvalue=r["pvalue"], padj=r["padj"],
                    n_donors=min(n_don_sf, n_don_pb), n_cells=n_cells))
        else:
            print(f"[05_scoring] {pop}: no ranked list (DE skipped) — no primary NES")

        # ---- SECONDARY: per-cell AUCell SMD (SF vs PB) on WT_heat_up ----
        dm = donor_means[donor_means["coarse_label"] == pop]
        sf = dm[dm[TISSUE_KEY] == TISSUE_NUM][AUC_UP].values
        pb = dm[dm[TISSUE_KEY] == TISSUE_DEN][AUC_UP].values
        st = smd(sf, pb)
        eff_rows.append(dict(
            dataset=DATASET, signature=f"{PRIMARY}_up", cell_state=pop,
            contrast="SF_vs_PB", effect_metric="percell_auc_smd",
            evidence_tier="secondary_percell", estimate=st["estimate"], se=st["se"],
            ci_low=st["ci_low"], ci_high=st["ci_high"],
            direction=("up" if (pd.notna(st["estimate"]) and st["estimate"] > 0)
                       else "down" if pd.notna(st["estimate"]) else "ns"),
            pvalue=st["pvalue"], padj=np.nan,
            n_donors=min(len(sf), len(pb)), n_cells=n_cells))

    eff = pd.DataFrame(eff_rows)
    # BH-FDR within each evidence tier for the secondary rows (primary padj comes from fgsea).
    from scipy.stats import false_discovery_control
    sec = eff["evidence_tier"] == "secondary_percell"
    if sec.any() and eff.loc[sec, "pvalue"].notna().any():
        p = eff.loc[sec, "pvalue"].fillna(1.0).values
        eff.loc[sec, "padj"] = false_discovery_control(p)

    cols = ["dataset", "signature", "cell_state", "contrast", "effect_metric", "evidence_tier",
            "estimate", "se", "ci_low", "ci_high", "direction", "pvalue", "padj",
            "n_donors", "n_cells", "signoff_state"]
    eff["signoff_state"] = effect_size_signoff_state()
    eff = eff[cols]
    eff.to_csv(PATHS.master_file("effect_sizes_treg_arthritis.csv"), index=False)
    append_master_table(eff, database=DATASET, stage=STAGE, name="master_effect_sizes",
                        config=FIG_CFG)

    print("[05_scoring] effect sizes:\n",
          eff[["signature", "cell_state", "effect_metric", "estimate", "pvalue", "padj"]]
          .to_string(index=False))
    print("[05_scoring] done")


if __name__ == "__main__":
    main()
