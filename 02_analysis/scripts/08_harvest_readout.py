#!/usr/bin/env python
"""
08_harvest_readout.py — COMPUTE (no plotting). Exploratory READOUT harvest.
===========================================================================
Adds two MSigDB Hallmark per-cell readouts — HALLMARK_HYPOXIA and
HALLMARK_UNFOLDED_PROTEIN_RESPONSE — to the frozen JIA SF/PB T-cell annotation, aligned
column-for-column with the reference compartment's Phase-1 readout so the two are directly
comparable. These are READOUTS: a hypoxia score is "consistent with a low-O2 /
metabolically-stressed state", and HIF causality is a further claim. This whole stage is
SECONDARY / annotation tier — per-cell, rank-based, composition-robust — and stays out of the
confirmatory donor-pseudobulk NES spine (stage 05) and the OR-gated POI harvest (stage 07). It
feeds the reactive review and the cross-species harvest question.

Engine: `score_cells_aucell_ucell` (helpers/geneset_utils.py → percell_score.R), the same
AUCell+UCell rank-based scorer stage 05 and the reference compartment use. AUCell is the
canonical readout; UCell rides alongside as a cross-check. The already-derived per-cell
readouts of THIS compartment — score_HSP, score_eTreg, and the mouse-anchor WT_heat up / down
/ up-minus-down channels — are carried in verbatim from the frozen explorer parquets, giving
the review one tidy table on one scorer run.

The mouse anchor enters as THREE readouts. `WT_heat_updown` is the balanced up-minus-down
composite, so a coordinated rise in both arms cancels inside it and the composite reads flat
whatever the arms do. Carrying `WT_heat_up` and `WT_heat_down` as their own readouts alongside
it makes each arm visible, which is what the pseudobulk ranked-list enrichment scores — the up
and down sets are scored separately there too, since AUCell/UCell are unsigned single-list
scorers.

Two levels of SF-vs-PB contrast are emitted, and they answer different questions:
  * per-cell (`harvest_readout_summary.csv`) — pooled over all cells, so it is
    pseudoreplicated: the unit of replication is the cell.
  * per-donor (`harvest_readout_donor_contrast.csv`) — each donor contributes one mean per
    (label × tissue), and the contrast is taken ACROSS donors, paired within donor. This is
    the level the donor-pseudobulk spine works at, so it is the one to read next to it. It
    stays annotation tier and produces no effect-size row.

Outputs:
  03_results/interactive/08_harvest_readout.parquet          (per-cell; gitignored/regenerable)
  03_results/08_harvest_readout/tables/harvest_readout_summary.csv
  03_results/08_harvest_readout/tables/harvest_readout_donor_means.csv
  03_results/08_harvest_readout/tables/harvest_readout_donor_contrast.csv
  03_results/08_harvest_readout/README.md                    (caption; hand-written)

Run in-container from the compartment root (or anywhere, via `from config import`):
    python 02_analysis/scripts/08_harvest_readout.py
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "02_analysis"))
os.chdir(ROOT)

from config import (PATHS, PARAMS, TISSUE_KEY, DONOR_KEY, TISSUE_NUM,  # noqa: E402
                    TISSUE_DEN)
from helpers.geneset_utils import score_cells_aucell_ucell, _symbol_to_varname  # noqa: E402

STAGE = "08_harvest_readout"
COARSE = "coarse_label"

# Frozen Hallmark sets (freeze_hallmark_sets.R). Scored per-cell → <set>_{AUCell,UCell}.
HALLMARK_DIR = ROOT / "00_data" / "references" / "msigdb_hallmark"
HALLMARK_SETS = ["HALLMARK_HYPOXIA", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"]

# Readout name → source column in the assembled per-cell table. AUCell is canonical
# for the two new Hallmark sets (aligned with the sting compartment).
READOUTS = {
    "WT_heat_up": "WT_heat_up",                               # mouse 39C anchor, up arm (annotation only)
    "WT_heat_down": "WT_heat_down",                           # mouse 39C anchor, down arm (annotation only)
    "WT_heat_updown": "WT_heat_updown",                       # balanced composite (arms can cancel)
    "score_eTreg": "score_eTreg",                             # effector-Treg (GSE161426)
    "score_HSP": "score_HSP",                                 # heat-shock / proteostasis
    "HALLMARK_HYPOXIA": "HALLMARK_HYPOXIA_AUCell",
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE": "HALLMARK_UNFOLDED_PROTEIN_RESPONSE_AUCell",
}


def _read_gene_list(path: Path) -> list[str]:
    with open(path) as fh:
        return [ln.strip() for ln in fh if ln.strip()]


def load_hallmark_sets() -> dict[str, list[str]]:
    missing = [s for s in HALLMARK_SETS if not (HALLMARK_DIR / f"{s}.txt").exists()]
    if missing:
        raise FileNotFoundError(
            f"frozen Hallmark set(s) absent: {missing} under {HALLMARK_DIR}. "
            "Regenerate with: Rscript 02_analysis/scripts/freeze_hallmark_sets.R")
    return {s: _read_gene_list(HALLMARK_DIR / f"{s}.txt") for s in HALLMARK_SETS}


def carry_frozen_readouts(index: pd.Index) -> pd.DataFrame:
    """Join already-derived per-cell readouts (score_HSP, score_eTreg, and the
    WT_heat up / down / up-minus-down channels) from the frozen explorer parquets —
    reindexed to the annotation object's barcodes. These are carried verbatim,
    NEVER re-derived here (guardrail)."""
    idir = PATHS.interactive_dir()
    out = pd.DataFrame(index=index)

    gonogo = pd.read_parquet(idir / "05_gonogo_explore.parquet")
    gonogo.index = gonogo.index.astype(str)
    for col in ("WT_heat_updown", "WT_heat_up", "WT_heat_down", "score_eTreg"):
        out[col] = gonogo[col].reindex(index).to_numpy() if col in gonogo else np.nan

    qc = pd.read_parquet(idir / "01_qc_explore.parquet")
    qc.index = qc.index.astype(str)
    out["score_HSP"] = qc["score_HSP"].reindex(index).to_numpy() if "score_HSP" in qc else np.nan
    return out


def build_parquet(adata, scored: pd.DataFrame, carried: pd.DataFrame) -> pd.DataFrame:
    """One row per cell: barcode + UMAP + coarse_label/tissue/donor + %mt + the two
    new Hallmark AUCell/UCell scores + the carried readouts."""
    idx = adata.obs_names.astype(str)
    df = pd.DataFrame(index=idx)
    df["barcode"] = idx
    um = np.asarray(adata.obsm["X_umap_unsupervised"])
    df["x"] = um[:, 0]
    df["y"] = um[:, 1]
    df[COARSE] = adata.obs[COARSE].astype(str).to_numpy()
    df[TISSUE_KEY] = adata.obs[TISSUE_KEY].astype(str).to_numpy()
    df[DONOR_KEY] = adata.obs[DONOR_KEY].astype(str).to_numpy()
    df["pct_counts_mt"] = adata.obs["pct_counts_mt"].to_numpy()
    for col in scored.columns:                       # HALLMARK_*_{AUCell,UCell}
        df[col] = scored[col].reindex(idx).to_numpy()
    for col in carried.columns:                      # WT_heat_{up,down,updown}, score_eTreg, score_HSP
        df[col] = carried[col].reindex(idx).to_numpy()
    return df.reset_index(drop=True)


def summarise(df: pd.DataFrame) -> pd.DataFrame:
    """Per (coarse_label × tissue × readout): mean/spread + mean %mt + n, plus a compact
    verdict — fraction above the readout's global P90 (high-pocket signal), and the
    per-cell SF-vs-PB shift (mean difference + standardized d), attached per cell-state.
    Long/tidy, mirroring the sting Phase-1 summary. Descriptive, correlative only.

    The `sf_minus_pb_*` columns here pool cells across donors, so their unit of
    replication is the cell and they are pseudoreplicated relative to the donor-level
    spine. `donor_contrast()` is the donor-level companion; read that one alongside
    donor-pseudobulk results."""
    # Global P90 per readout (across all cells) → the high-pocket threshold.
    p90 = {name: np.nanpercentile(df[col].to_numpy(dtype=float), 90)
           for name, col in READOUTS.items()}

    rows = []
    for label, gl in df.groupby(COARSE, observed=True):
        # per cell-state SF-vs-PB per-cell shift (mean diff + Cohen's d), per readout.
        sf_all = gl[gl[TISSUE_KEY] == TISSUE_NUM]
        pb_all = gl[gl[TISSUE_KEY] == TISSUE_DEN]
        shift = {}
        for name, col in READOUTS.items():
            sf = sf_all[col].to_numpy(dtype=float)
            pb = pb_all[col].to_numpy(dtype=float)
            sf = sf[np.isfinite(sf)]
            pb = pb[np.isfinite(pb)]
            mean_shift = (float(np.mean(sf)) - float(np.mean(pb))
                          if len(sf) and len(pb) else np.nan)
            d = np.nan
            if len(sf) >= 2 and len(pb) >= 2:
                n1, n2 = len(sf), len(pb)
                pooled = np.sqrt(((n1 - 1) * sf.std(ddof=1) ** 2
                                  + (n2 - 1) * pb.std(ddof=1) ** 2) / (n1 + n2 - 2))
                d = (sf.mean() - pb.mean()) / (pooled + 1e-12)
            shift[name] = (mean_shift, d)

        for tissue, gt in gl.groupby(TISSUE_KEY, observed=True):
            for name, col in READOUTS.items():
                v = gt[col].to_numpy(dtype=float)
                v = v[np.isfinite(v)]
                rows.append(dict(
                    coarse_label=label, tissue=tissue, readout=name,
                    n_cells=int(len(v)),
                    mean=float(np.mean(v)) if len(v) else np.nan,
                    sd=float(np.std(v, ddof=1)) if len(v) > 1 else np.nan,
                    median=float(np.median(v)) if len(v) else np.nan,
                    frac_above_p90=(float(np.mean(v > p90[name])) if len(v) else np.nan),
                    mean_pct_counts_mt=float(np.nanmean(gt["pct_counts_mt"].to_numpy(dtype=float)))
                    if len(gt) else np.nan,
                    sf_minus_pb_mean=shift[name][0],
                    sf_minus_pb_smd=shift[name][1],
                    global_p90=float(p90[name])))
    out = pd.DataFrame(rows)
    order = {"Treg": 0, "Tcon": 1, "CD8": 2}
    out = out.sort_values(
        by=["coarse_label", "readout", "tissue"],
        key=lambda s: s.map(order) if s.name == "coarse_label" else s
    ).reset_index(drop=True)
    return out


def donor_means(df: pd.DataFrame, min_cells: int) -> pd.DataFrame:
    """One row per (donor × tissue × coarse_label): the mean of every readout over that
    donor's cells, plus `n_cells`. Strata thinner than `min_cells` are dropped, reusing
    the same donor-stratum floor the donor-pseudobulk aggregation applies, so the two
    levels see the same donors."""
    keys = [DONOR_KEY, TISSUE_KEY, COARSE]
    cols = list(READOUTS.values())
    g = df.groupby(keys, observed=True)
    out = g[cols].mean()
    out["n_cells"] = g.size()
    out = out.reset_index().rename(columns={v: k for k, v in READOUTS.items() if v != k})
    dropped = out[out["n_cells"] < min_cells]
    if len(dropped):
        print(f"[08_harvest] donor strata below min_cells={min_cells}, dropped: "
              + ", ".join(f"{r[DONOR_KEY]}/{r[TISSUE_KEY]}/{r[COARSE]}(n={r['n_cells']})"
                          for _, r in dropped.iterrows()))
    out = out[out["n_cells"] >= min_cells].reset_index(drop=True)
    order = {"Treg": 0, "Tcon": 1, "CD8": 2}
    return out.sort_values(
        by=[COARSE, DONOR_KEY, TISSUE_KEY],
        key=lambda s: s.map(order) if s.name == COARSE else s).reset_index(drop=True)


def paired_effect(diff: np.ndarray) -> dict:
    """Paired standardized mean difference (Cohen's dz) on the within-donor SF-minus-PB
    differences, with a 95% CI and a paired-t p-value. dz = mean(diff)/sd(diff); its
    approximate SE is sqrt(1/n + dz^2/(2n)). Mirrors the unpaired donor-level SMD used
    for the secondary per-cell effect sizes, but keeps the donor pairing."""
    diff = diff[np.isfinite(diff)]
    n = len(diff)
    nan = dict(estimate=np.nan, se=np.nan, ci_low=np.nan, ci_high=np.nan, pvalue=np.nan)
    if n < 2:
        return nan
    sd = float(diff.std(ddof=1))
    if not np.isfinite(sd) or sd == 0.0:
        return nan
    dz = float(diff.mean()) / sd
    se = float(np.sqrt(1.0 / n + dz ** 2 / (2.0 * n)))
    t = dz * np.sqrt(n)
    p = float(2 * stats.t.sf(abs(t), df=n - 1))
    return dict(estimate=dz, se=se, ci_low=dz - 1.96 * se, ci_high=dz + 1.96 * se,
                pvalue=p)


def unpaired_smd(sf: np.ndarray, pb: np.ndarray) -> float:
    """Donor-level Cohen's d ignoring the pairing — the same formula the secondary
    per-cell effect sizes use, kept here only so the paired estimate can be read against
    an unpaired one on the same donor means."""
    sf, pb = sf[np.isfinite(sf)], pb[np.isfinite(pb)]
    n1, n2 = len(sf), len(pb)
    if n1 < 2 or n2 < 2:
        return np.nan
    pooled = np.sqrt(((n1 - 1) * sf.std(ddof=1) ** 2
                      + (n2 - 1) * pb.std(ddof=1) ** 2) / (n1 + n2 - 2))
    return float((sf.mean() - pb.mean()) / (pooled + 1e-12))


def donor_contrast(dm: pd.DataFrame) -> pd.DataFrame:
    """Per (coarse_label × readout): the SF-minus-PB contrast taken ACROSS DONORS, on
    the per-donor means from `donor_means()`. Donors carrying both tissues are paired
    within donor (GSE160097 is a paired SF/PB design); `n_donors_paired` records how
    many actually pair, alongside the per-arm donor counts. Positive = higher in
    synovial fluid. Annotation tier — never an effect-size row."""
    rows = []
    order = {"Treg": 0, "Tcon": 1, "CD8": 2}
    for label, gl in dm.groupby(COARSE, observed=True):
        sf_arm = gl[gl[TISSUE_KEY] == TISSUE_NUM].set_index(DONOR_KEY)
        pb_arm = gl[gl[TISSUE_KEY] == TISSUE_DEN].set_index(DONOR_KEY)
        both = sorted(set(sf_arm.index) & set(pb_arm.index))
        for name in READOUTS:
            sf = sf_arm[name].to_numpy(dtype=float)
            pb = pb_arm[name].to_numpy(dtype=float)
            diff = (sf_arm.loc[both, name].to_numpy(dtype=float)
                    - pb_arm.loc[both, name].to_numpy(dtype=float)) if both else np.array([])
            eff = paired_effect(diff)
            finite = diff[np.isfinite(diff)]
            rows.append(dict(
                coarse_label=label, readout=name,
                n_donors_paired=int(len(finite)),
                n_donors_sf=int(np.isfinite(sf).sum()),
                n_donors_pb=int(np.isfinite(pb).sum()),
                donor_mean_sf=float(np.nanmean(sf)) if len(sf) else np.nan,
                donor_mean_pb=float(np.nanmean(pb)) if len(pb) else np.nan,
                sf_minus_pb_mean=float(finite.mean()) if len(finite) else np.nan,
                sf_minus_pb_sd=float(finite.std(ddof=1)) if len(finite) > 1 else np.nan,
                sf_minus_pb_dz=eff["estimate"], dz_se=eff["se"],
                dz_ci_low=eff["ci_low"], dz_ci_high=eff["ci_high"],
                dz_pvalue=eff["pvalue"],
                sf_minus_pb_smd_unpaired=unpaired_smd(sf, pb),
                min_cells_per_donor_stratum=int(gl["n_cells"].min())))
    out = pd.DataFrame(rows)
    return out.sort_values(
        by=["coarse_label", "readout"],
        key=lambda s: s.map(order) if s.name == "coarse_label" else s
    ).reset_index(drop=True)


def main() -> None:
    n_cores = int(PARAMS.get("percell_score_ncores", 4))
    hallmark = load_hallmark_sets()
    print(f"[08_harvest] Hallmark sets: "
          + ", ".join(f"{k}={len(v)}" for k, v in hallmark.items()))

    adata = sc.read_h5ad(PATHS.object("02_annotation"))
    print(f"[08_harvest] annotation object: {adata.n_obs} cells x {adata.n_vars} genes")

    # per-cell coverage (how many set genes are present in the object)
    sym_to_var = _symbol_to_varname(adata, "gene_symbol")
    for k, v in hallmark.items():
        print(f"[08_harvest] {k} coverage: "
              f"{sum(g in sym_to_var for g in v)}/{len(v)} genes present")

    # --- score the two Hallmark sets per cell (AUCell + UCell on lognorm X) ---
    scored = score_cells_aucell_ucell(adata, hallmark, layer=None,
                                      symbol_col="gene_symbol", n_cores=n_cores)
    auc_cols = [c for c in scored.columns if c.endswith("_AUCell")]
    for c in auc_cols:
        vals = scored[c].to_numpy(dtype=float)
        vals = vals[np.isfinite(vals)]
        assert vals.min() >= 0.0 and vals.max() <= 1.0, f"{c} out of [0,1]: {vals.min()}..{vals.max()}"
        print(f"[08_harvest] {c} in [{vals.min():.4f}, {vals.max():.4f}] (n={len(vals)}) OK")

    # --- carry already-derived readouts from the frozen explorers (never re-derived) ---
    carried = carry_frozen_readouts(adata.obs_names.astype(str))
    for col in carried.columns:
        print(f"[08_harvest] carried {col}: NaN frac {carried[col].isna().mean():.4f}")

    # --- per-cell parquet (gitignored/regenerable) ---
    df = build_parquet(adata, scored, carried)
    pq_path = PATHS.interactive_dir() / "08_harvest_readout.parquet"
    df.to_parquet(pq_path, index=False)
    print(f"[08_harvest] wrote {pq_path.name} ({df.shape[0]} cells, {df.shape[1]} cols)")

    # --- summary (per coarse_label × tissue × readout) ---
    summary = summarise(df)
    sdir = PATHS.tables(STAGE)
    summary.to_csv(sdir / "harvest_readout_summary.csv", index=False)
    print(f"[08_harvest] wrote harvest_readout_summary.csv ({summary.shape[0]} rows)")
    print("\n[08_harvest] summary (mean, frac>P90, SF-PB shift):\n",
          summary[["coarse_label", "tissue", "readout", "n_cells", "mean",
                   "frac_above_p90", "sf_minus_pb_mean", "sf_minus_pb_smd"]]
          .to_string(index=False))

    # --- donor-level companion (per-donor means, then SF-vs-PB across donors) ---
    # The per-cell shift above pools cells and is pseudoreplicated; this is the level
    # the donor-pseudobulk spine works at, so it is the comparable one. Still
    # annotation tier: nothing here is written to master/ or to an effect-size row.
    dm = donor_means(df, min_cells=int(PARAMS.get("pseudobulk_min_cells", 20)))
    dm.to_csv(sdir / "harvest_readout_donor_means.csv", index=False)
    print(f"[08_harvest] wrote harvest_readout_donor_means.csv ({dm.shape[0]} rows, "
          f"{dm[DONOR_KEY].nunique()} donors)")

    contrast = donor_contrast(dm)
    contrast.to_csv(sdir / "harvest_readout_donor_contrast.csv", index=False)
    print(f"[08_harvest] wrote harvest_readout_donor_contrast.csv ({contrast.shape[0]} rows)")
    print("\n[08_harvest] donor-level SF-vs-PB (paired within donor):\n",
          contrast[["coarse_label", "readout", "n_donors_paired", "sf_minus_pb_mean",
                    "sf_minus_pb_dz", "dz_ci_low", "dz_ci_high", "dz_pvalue"]]
          .to_string(index=False))

    print("\n[08_harvest] DONE. Exploratory / annotation tier — never pooled with the "
          "pseudobulk NES spine.")


if __name__ == "__main__":
    main()
