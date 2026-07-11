#!/usr/bin/env python
"""
01_qc_mthi_characterize.py — COMPUTE (no plotting).
===================================================
Turns the interactively-found "mito-high (mt-hi), effector-like Treg pocket"
(01_qc_explore.qmd, saved lasso `eda/selection_sf_ctla4_pocket.csv`, 335 cells)
into a REPRODUCIBLE population rule + honest statistics, so the go/no-go can
defensibly RETAIN all Tregs including this high-mito subset (not QC it out as
dying/doublets).

Substrate: the frozen explorer parquet `03_results/interactive/01_qc_explore.parquet`
(per-cell; NOT the h5ad). Everything is computed among SORTED Treg cells.

Population rule (two independent constructions, cross-checked against each other
AND the saved 335 barcodes):
  A. CLUSTER rule — among Treg, leiden_unsupervised clusters enriched for the
     saved lasso by one-sided Fisher (BH-FDR<0.05). Two clusters clear: {6, 16},
     both median %mt ~20 (vs ~4 rest). The EFFECTOR-LIKE pocket is the one that is
     also effector (median score_eTreg > global Treg median) = cluster 6. Cluster
     16 is mt-hi but eTreg-LOW, PB-dominated, low-complexity — a DISTINCT cluster,
     NOT the effector pocket.
  B. THRESHOLD rule — among Treg, pct_counts_mt >= 97.5th percentile (=10.0%,
     the %mt that cleanly separates clusters {6,16}: ~93% of them vs ~0.8% of other
     Treg exceed it) AND score_eTreg >= median. Recovers cluster 6 (Jaccard 0.64).

Primary comparison for the retention claim: mt-hi effector Treg (cluster 6) vs
normal Treg (Treg not in {6,16}). Cluster 16 carried as a labelled reference group.

Evidence tier: ALL per-cell comparisons here are `secondary_percell` — never pooled
with pseudobulk NES. score_eTreg is an EDA-derived signature (descriptive), and
score_HSP / score_apoptosis are Tier-3 hand markers (QC-descriptive, NOT evidence);
the decisive not-debris datum is n_genes_by_counts (objective complexity).

Outputs (03_results/01_qc/tables/):
  mthi_population_rule.csv       — rule definitions + A/B/saved concordance
  mthi_cluster_enrichment.csv    — per Treg cluster: Fisher enrichment + %mt/eTreg (Panel A)
  mthi_identity_retention.csv    — FOXP3/IKZF2/CTLA4/IL2RA/TIGIT, mt-hi vs rest (Panel B)
  mthi_qc_discrimination.csv     — complexity/apoptosis/HSP + doublet availability (Panel C)
  mthi_heat_percell.csv          — WT_heat_up/_updown per group (Panel D)
  mthi_donor_tissue.csv          — donor x tissue composition + DE-floor verdict (Panel E)
  mthi_treg_membership.csv       — per-Treg reproducible label + coords (viz input)
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "02_analysis"))
os.chdir(ROOT)

from config import PATHS, PARAMS  # noqa: E402

STAGE = "01_qc"
SCRIPT = "02_analysis/scripts/01_qc_mthi_characterize.py"

# Frozen inputs.
PARQUET = PATHS.interactive / "01_qc_explore.parquet"
SAVED_POCKET = (ROOT / "02_analysis/notebooks/01_qc_explore/eda/"
                "selection_sf_ctla4_pocket.csv")

MT_PCTILE = 97.5            # threshold-rule %mt percentile (within Treg)
IDENTITY_GENES = ["FOXP3", "IKZF2", "CTLA4", "IL2RA", "TIGIT"]
QC_METRICS = ["n_genes_by_counts", "total_counts", "score_apoptosis", "score_HSP"]
TIER = "secondary_percell"


# --------------------------------------------------------------------------- #
# stat helpers                                                                #
# --------------------------------------------------------------------------- #
def _bh_fdr(pvals: list[float]) -> np.ndarray:
    """Benjamini-Hochberg FDR (stdlib-safe; no statsmodels dependency)."""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def _mann_whitney(x: pd.Series, y: pd.Series) -> tuple[float, float, float]:
    """Two-sided Mann-Whitney U + rank-biserial correlation effect size."""
    x = x.astype(float).values
    y = y.astype(float).values
    U, p = stats.mannwhitneyu(x, y, alternative="two-sided")
    rbc = 2.0 * U / (len(x) * len(y)) - 1.0  # rank-biserial in [-1, 1]
    return float(U), float(p), float(rbc)


def _jaccard(a: set, b: set) -> float:
    return len(a & b) / len(a | b) if (a | b) else 0.0


# --------------------------------------------------------------------------- #
def main() -> None:
    df = pd.read_parquet(PARQUET)
    saved = pd.read_csv(SAVED_POCKET, index_col=0)
    tdir = PATHS.tables(STAGE)

    treg = df[df["population_short"] == "Treg"].copy()
    treg["in_saved335"] = treg.index.isin(saved.index)
    n_treg = len(treg)
    saved_treg = set(treg.index[treg["in_saved335"]])

    # --- Method A: per-cluster Fisher enrichment for the saved lasso ------- #
    n_pocket = int(treg["in_saved335"].sum())
    global_etreg_med = float(treg["score_eTreg"].median())
    rows = []
    for cl, g in treg.groupby("leiden_unsupervised", observed=True):
        a = int(g["in_saved335"].sum())
        b = len(g) - a
        c = n_pocket - a
        d = (n_treg - len(g)) - c
        OR, p = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
        rows.append({
            "leiden": str(cl),
            "n_treg": len(g),
            "n_saved_in_cluster": a,
            "frac_saved_captured": a / n_pocket if n_pocket else 0.0,
            "median_pct_mt": float(g["pct_counts_mt"].median()),
            "median_score_eTreg": float(g["score_eTreg"].median()),
            "median_n_genes": float(g["n_genes_by_counts"].median()),
            "fisher_or": float(OR) if np.isfinite(OR) else np.nan,
            "fisher_p": float(p),
        })
    clu = pd.DataFrame(rows)
    clu["fisher_fdr"] = _bh_fdr(clu["fisher_p"].tolist())
    # A cluster is mt-hi if Fisher-enriched (FDR<0.05) AND median %mt > 2x global Treg median.
    global_mt_med = float(treg["pct_counts_mt"].median())
    clu["is_mthi_cluster"] = (clu["fisher_fdr"] < 0.05) & (clu["median_pct_mt"] > 2 * global_mt_med)
    # The effector-like pocket = mt-hi cluster(s) that are ALSO effector (eTreg above global median).
    clu["is_mthi_effector"] = clu["is_mthi_cluster"] & (clu["median_score_eTreg"] > global_etreg_med)
    clu = clu.sort_values("n_saved_in_cluster", ascending=False).reset_index(drop=True)
    clu["evidence_tier"] = TIER
    clu.to_csv(tdir / "mthi_cluster_enrichment.csv", index=False)

    mthi_clusters = clu.loc[clu["is_mthi_cluster"], "leiden"].tolist()          # {6,16}
    effector_clusters = clu.loc[clu["is_mthi_effector"], "leiden"].tolist()     # {6}

    # --- per-Treg group label --------------------------------------------- #
    def _group(cl: str) -> str:
        if cl in effector_clusters:
            return "mt_hi_effector"
        if cl in mthi_clusters:
            return "mt_hi_noneffector"
        return "normal_Treg"

    treg["mthi_group"] = treg["leiden_unsupervised"].astype(str).map(_group)

    # --- Method B: threshold rule ----------------------------------------- #
    mt_thr = float(np.percentile(treg["pct_counts_mt"], MT_PCTILE))
    treg["rule_B_threshold"] = (treg["pct_counts_mt"] >= mt_thr) & (treg["score_eTreg"] >= global_etreg_med)
    treg["rule_A_effector"] = treg["mthi_group"] == "mt_hi_effector"

    set_A = set(treg.index[treg["rule_A_effector"]])
    set_B = set(treg.index[treg["rule_B_threshold"]])

    # --- population-rule + concordance table ------------------------------ #
    rule_rows = [
        {"rule": "A_cluster_effector",
         "definition": f"Treg in leiden cluster(s) {effector_clusters} "
                       "(Fisher-enriched for saved lasso, FDR<0.05, median %mt>2x global, "
                       "median score_eTreg>global median)",
         "n_cells": len(set_A)},
        {"rule": "A_cluster_all_mthi",
         "definition": f"Treg in leiden cluster(s) {mthi_clusters} (all mt-hi enriched clusters; "
                       "includes the eTreg-LOW non-effector cluster)",
         "n_cells": int(treg["leiden_unsupervised"].astype(str).isin(mthi_clusters).sum())},
        {"rule": "B_threshold",
         "definition": f"Treg with pct_counts_mt>=P{MT_PCTILE}({mt_thr:.2f}) AND "
                       f"score_eTreg>=median({global_etreg_med:.4f})",
         "n_cells": len(set_B)},
        {"rule": "saved_lasso_335_Treg",
         "definition": "Treg subset of the interactively saved 335-cell lasso "
                       "(selection_sf_ctla4_pocket.csv)",
         "n_cells": len(saved_treg)},
    ]
    rule = pd.DataFrame(rule_rows)
    rule["jaccard_vs_A_effector"] = [
        np.nan,
        _jaccard(set_A, set(treg.index[treg["leiden_unsupervised"].astype(str).isin(mthi_clusters)])),
        _jaccard(set_A, set_B),
        _jaccard(set_A, saved_treg),
    ]
    rule["frac_saved_recovered"] = [
        len(set_A & saved_treg) / len(saved_treg),
        len(set(treg.index[treg["leiden_unsupervised"].astype(str).isin(mthi_clusters)]) & saved_treg) / len(saved_treg),
        len(set_B & saved_treg) / len(saved_treg),
        1.0,
    ]
    rule["mt_pctile"] = MT_PCTILE
    rule["mt_threshold"] = mt_thr
    rule["global_etreg_median"] = global_etreg_med
    rule["evidence_tier"] = TIER
    rule.to_csv(tdir / "mthi_population_rule.csv", index=False)

    # --- groups for the comparison stats ---------------------------------- #
    eff = treg[treg["mthi_group"] == "mt_hi_effector"]      # cluster 6
    non = treg[treg["mthi_group"] == "mt_hi_noneffector"]   # cluster 16
    rest = treg[treg["mthi_group"] == "normal_Treg"]

    # --- Panel B: identity retention (mt-hi effector vs normal Treg) ------- #
    rows = []
    for g in IDENTITY_GENES:
        U, p, rbc = _mann_whitney(eff[g], rest[g])
        rows.append({
            "gene": g,
            "median_mthi_eff": float(eff[g].median()),
            "median_mthi_noneff": float(non[g].median()),
            "median_normal": float(rest[g].median()),
            "frac_expr_mthi_eff": float((eff[g] > 0).mean()),
            "frac_expr_mthi_noneff": float((non[g] > 0).mean()),
            "frac_expr_normal": float((rest[g] > 0).mean()),
            "mannwhitney_u": U,
            "rank_biserial": rbc,
            "pvalue": p,
        })
    ident = pd.DataFrame(rows)
    ident["fdr"] = _bh_fdr(ident["pvalue"].tolist())
    ident["comparison"] = "mt_hi_effector_vs_normal_Treg"
    ident["evidence_tier"] = TIER
    ident.to_csv(tdir / "mthi_identity_retention.csv", index=False)

    # --- Panel C: dying/debris discrimination ----------------------------- #
    # Tier tags: n_genes_by_counts/total_counts = objective QC; apoptosis/HSP = Tier-3 hand marker.
    tier_map = {"n_genes_by_counts": "objective_qc", "total_counts": "objective_qc",
                "score_apoptosis": "tier3_handmarker_descriptive",
                "score_HSP": "tier3_handmarker_descriptive"}
    rows = []
    for m in QC_METRICS:
        U, p, rbc = _mann_whitney(eff[m], rest[m])
        rows.append({
            "metric": m,
            "metric_tier": tier_map[m],
            "median_mthi_eff": float(eff[m].median()),
            "median_mthi_noneff": float(non[m].median()),
            "median_normal": float(rest[m].median()),
            "mannwhitney_u": U,
            "rank_biserial": rbc,
            "pvalue": p,
            "elevated_in_pocket": bool(rbc > 0),
        })
    qc = pd.DataFrame(rows)
    qc["fdr"] = _bh_fdr(qc["pvalue"].tolist())
    # Doublet availability: Scrublet did not populate doublet_score in this run.
    n_ds = int(treg["doublet_score"].notna().sum()) if "doublet_score" in treg else 0
    qc.attrs["min_genes_floor"] = int(PARAMS.qc_min_genes)
    qc["comparison"] = "mt_hi_effector_vs_normal_Treg"
    qc["evidence_tier"] = TIER
    qc["doublet_score_n_nonnull"] = n_ds
    qc["qc_min_genes_floor"] = int(PARAMS.qc_min_genes)
    qc.to_csv(tdir / "mthi_qc_discrimination.csv", index=False)

    # --- Panel D: heat honesty (secondary_percell; must NOT read as heat signal) --- #
    groups = {
        "mt_hi_effector": eff, "mt_hi_noneffector": non, "normal_Treg": rest,
        "Tcon": df[df["population_short"] == "Tcon"],
        "CD8": df[df["population_short"] == "CD8"],
    }
    rows = []
    for name, gd in groups.items():
        rec = {"group": name, "n_cells": len(gd)}
        for h in ["WT_heat_up", "WT_heat_updown"]:
            rec[f"{h}_median"] = float(gd[h].median())
            rec[f"{h}_q25"] = float(gd[h].quantile(0.25))
            rec[f"{h}_q75"] = float(gd[h].quantile(0.75))
        rows.append(rec)
    heat = pd.DataFrame(rows)
    # mt-hi effector vs normal Treg tests (both channels).
    for h in ["WT_heat_up", "WT_heat_updown"]:
        U, p, rbc = _mann_whitney(eff[h], rest[h])
        heat.loc[heat["group"] == "mt_hi_effector", f"{h}_rbc_vs_normal"] = rbc
        heat.loc[heat["group"] == "mt_hi_effector", f"{h}_p_vs_normal"] = p
    heat["evidence_tier"] = TIER
    heat["caveat"] = ("secondary_percell; WT_heat_up shift is confounded (correlates with "
                      "effector/depth); balanced WT_heat_updown is flat; pseudobulk go/no-go NES "
                      "unchanged when these cells were recovered. NOT evidence the pocket carries heat.")
    heat.to_csv(tdir / "mthi_heat_percell.csv", index=False)

    # --- Panel E: donor x tissue composition + DE-floor verdict ------------ #
    min_cells = int(PARAMS.pseudobulk_min_cells)
    min_donors = int(PARAMS.pseudobulk_min_donors)
    rows = []
    for (donor, tissue), g in eff.groupby(["donor", "tissue"], observed=True):
        rows.append({"donor": donor, "tissue": tissue, "n_cells": len(g)})
    comp = pd.DataFrame(rows)
    comp["frac_of_pocket"] = comp["n_cells"] / len(eff)
    comp["meets_min_cells"] = comp["n_cells"] >= min_cells
    # Per-tissue donor floor: how many donors clear min_cells in each arm?
    sf = comp[(comp["tissue"] == "synovial_fluid") & comp["meets_min_cells"]]["donor"].nunique()
    pb = comp[(comp["tissue"] == "peripheral_blood") & comp["meets_min_cells"]]["donor"].nunique()
    comp["evidence_tier"] = TIER
    comp["pseudobulk_min_cells"] = min_cells
    comp["pseudobulk_min_donors"] = min_donors
    comp["sf_donors_meeting_floor"] = sf
    comp["pb_donors_meeting_floor"] = pb
    # A paired SF-vs-PB contrast needs >=min_donors in BOTH arms.
    comp["supports_own_de_stratum"] = (sf >= min_donors) and (pb >= min_donors)
    comp = comp.sort_values(["tissue", "n_cells"], ascending=[True, False])
    comp.to_csv(tdir / "mthi_donor_tissue.csv", index=False)

    # --- per-Treg membership (reproducible label + viz coords) ------------- #
    membership = treg[[
        "leiden_unsupervised", "mthi_group", "rule_A_effector", "rule_B_threshold",
        "in_saved335", "x", "y", "pct_counts_mt", "score_eTreg",
        "n_genes_by_counts", "donor", "tissue",
    ]].copy()
    membership.to_csv(tdir / "mthi_treg_membership.csv")

    # --- console summary --------------------------------------------------- #
    print(f"[mthi] Treg n={n_treg}; mt-hi clusters={mthi_clusters}; effector pocket={effector_clusters}")
    print(f"[mthi] rule A (cluster {effector_clusters}) n={len(set_A)}; "
          f"rule B (mt>=P{MT_PCTILE}={mt_thr:.1f} & eTreg>=med) n={len(set_B)}; "
          f"A-B Jaccard={_jaccard(set_A, set_B):.3f}; saved recovered by A={len(set_A & saved_treg)/len(saved_treg):.3f}")
    print(f"[mthi] Panel B FOXP3 frac {ident.loc[ident.gene=='FOXP3','frac_expr_mthi_eff'].iat[0]:.2f} "
          f"vs {ident.loc[ident.gene=='FOXP3','frac_expr_normal'].iat[0]:.2f}; "
          f"IKZF2 rbc {ident.loc[ident.gene=='IKZF2','rank_biserial'].iat[0]:.2f}")
    print(f"[mthi] Panel C n_genes {qc.loc[qc.metric=='n_genes_by_counts','median_mthi_eff'].iat[0]:.0f} "
          f"(floor {int(PARAMS.qc_min_genes)}); apoptosis rbc "
          f"{qc.loc[qc.metric=='score_apoptosis','rank_biserial'].iat[0]:.2f}; "
          f"doublet_score nonnull={n_ds}")
    print(f"[mthi] Panel E SF donors>=floor={sf}, PB donors>=floor={pb}, "
          f"own DE stratum? {comp['supports_own_de_stratum'].iat[0]}")
    print("[mthi] wrote 7 tables to", tdir)


if __name__ == "__main__":
    main()
