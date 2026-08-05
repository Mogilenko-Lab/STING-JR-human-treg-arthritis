#!/usr/bin/env python
"""
16_narrative_scoring.py — COMPUTE ONLY (no plotting).
=============================================================================
Builds ONE per-cell substrate that can colour the frozen JIA embedding by an
empirical mouse-derived up arm and by a curated program lens, so the two
colourings can be read side by side on the same map.

SCOPE AND TIER. This is SECONDARY / annotation tier. Nothing here is pooled with
the donor-pseudobulk NES spine and no row is written to
`03_results/master/effect_sizes_treg_arthritis.csv`. A per-cell score localises a
program on a map; it does not test it. Confirmatory claims stay with donor-level
pseudobulk differential expression.

ONE METRIC: AUCell. The project computes both AUCell and UCell elsewhere
(`helpers/geneset_utils.score_cells_aucell_ucell`). This stage ships **AUCell
only** — a single metric so that two colourings of the same map differ by gene
set and by nothing else. UCell is computed as a by-product of the shared scorer
and is dropped before the substrate is written.

THIRTEEN GENE SETS, UP ARMS ONLY. Four mouse-derived human-projected up arms and
nine curated, anchor-independent program lenses (see `MOUSE_ARMS` / `CURATED`).
Down arms are deliberately absent: they are not scored, not carried, and not
named. A per-cell colouring answers "where on this map is this program high",
and a down arm inverted onto that question reads as an absence, which a
continuous colour scale cannot show honestly.

THE EMBEDDING IS NOT RECOMPUTED. `barcode, x, y, coarse_label, tissue, donor,
pct_counts_mt` are taken verbatim from the published per-cell readout
`03_results/interactive/08_harvest_readout.parquet` so this substrate lands on
exactly the map the published figures use. Three already-published score columns
ride along under a `published_` prefix so the notebook needs one file, not two.

PROVENANCE SEAM GUARD — an in-run assertion, not a published artifact. Three of
the sets scored here already have a per-cell column in the published readout, so
re-deriving them is a scale check on this whole substrate. Recomputing AUCell over
the same matrix with the same gene set is deterministic, so the guard demands
r = 1.000000 to within SEAM_R_EXACT on every same-metric row and halts the run
otherwise. It carries no biological reading, so it is asserted and printed here
rather than written into `03_results/` as a table or a figure.

  * `HALLMARK_HYPOXIA` and `HALLMARK_UNFOLDED_PROTEIN_RESPONSE` reproduce the
    published AUCell columns at r = 1.000000. The scorer is bit-for-bit stable.
  * `published_WT_heat_up` reproduces at only r ~ 0.755, because that column is
    **not AUCell**. It is a stale mean-centred scanpy `score_genes` module score,
    carried verbatim into the published readout from `05_gonogo_explore.parquet`,
    which predates the migration of stage 05 to AUCell. Its values run negative
    (min -0.154, mean -0.052); AUCell is bounded in [0, 1]. That row is the one
    `cross_metric` comparison and is exempt from the guard by construction.

A fourth row compares the newly scored `WT_heat_up_AUCell` against stage 05's
`per_cell_scores.csv` — the apples-to-apples anchor for the mouse arm, which does
reproduce at r = 1.000000, and which the guard therefore holds to. **Colour the
mouse up arm with `WT_heat_up_AUCell`, never with `published_WT_heat_up`.** The
stale column is carried only so the discrepancy stays visible rather than being
quietly dropped.

Outputs:
  03_results/interactive/16_narrative_embedding.parquet             (per-cell substrate)
  03_results/16_narrative_scoring/tables/narrative_scoring_manifest.csv
  03_results/16_narrative_scoring/tables/narrative_score_summary.csv

Run in-container from the compartment root:
    python 02_analysis/scripts/16_narrative_scoring.py
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

COMPARTMENT_ROOT = Path(__file__).resolve().parents[2]
REPO_ROOT = COMPARTMENT_ROOT.parent

sys.path.insert(0, str(COMPARTMENT_ROOT))
sys.path.insert(0, str(COMPARTMENT_ROOT / "02_analysis"))
os.chdir(COMPARTMENT_ROOT)

from config import CONFIG, PATHS, PARAMS, TISSUE_KEY, DONOR_KEY  # noqa: E402
from helpers.geneset_utils import (load_alias_map, resolve_symbols,  # noqa: E402
                                   score_cells_aucell_ucell, _symbol_to_varname)

STAGE = "16_narrative_scoring"
COARSE = "coarse_label"
METRIC = "AUCell"

# The published per-cell readout supplying the embedding + the carried-forward columns.
PUBLISHED_PARQUET = "08_harvest_readout.parquet"
SUBSTRATE_PARQUET = "16_narrative_embedding.parquet"

# Columns taken verbatim from the published readout — the map is NOT recomputed here.
EMBEDDING_COLS = ["barcode", "x", "y", COARSE, TISSUE_KEY, DONOR_KEY, "pct_counts_mt"]

# Already-published per-cell score columns carried forward under a `published_` prefix,
# so the provenance of every column is legible from its own name.
CARRY_FORWARD = [
    "WT_heat_up",                                   # NOT AUCell — see the seam check
    "HALLMARK_HYPOXIA_AUCell",
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE_AUCell",
]
CARRY_PREFIX = "published_"

_MOUSE_SIG_DIR = REPO_ROOT / "mouse_anchor/03_results/human_projection/signatures"
_HALLMARK_DIR = COMPARTMENT_ROOT / "00_data/references/msigdb_hallmark"
_HSR_DIR = COMPARTMENT_ROOT / "00_data/references/temp_hsr_lens"
_STING_SIG_DIR = REPO_ROOT / "sting_positive_control/03_results/06_reference_axis/signatures"

# --- the four mouse-derived, human-projected UP arms (empirical; anchor-dependent) ---
MOUSE_ARMS = {
    "WT_heat_up": _MOUSE_SIG_DIR / "WT_heat/WT_heat_up.txt",
    "KO_heat_up": _MOUSE_SIG_DIR / "KO_heat/KO_heat_up.txt",
    "Interaction_up": _MOUSE_SIG_DIR / "Interaction/Interaction_up.txt",
    "Interaction_fdrOnly_up": _MOUSE_SIG_DIR / "Interaction/Interaction_fdrOnly_up.txt",
}

# --- the nine curated program lenses (versioned/published; anchor-independent) ---
CURATED = {
    "HALLMARK_HYPOXIA": _HALLMARK_DIR / "HALLMARK_HYPOXIA.txt",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB": _HALLMARK_DIR / "HALLMARK_TNFA_SIGNALING_VIA_NFKB.txt",
    "HALLMARK_INFLAMMATORY_RESPONSE": _HALLMARK_DIR / "HALLMARK_INFLAMMATORY_RESPONSE.txt",
    "HALLMARK_IL2_STAT5_SIGNALING": _HALLMARK_DIR / "HALLMARK_IL2_STAT5_SIGNALING.txt",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE": _HALLMARK_DIR / "HALLMARK_INTERFERON_ALPHA_RESPONSE.txt",
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE": _HALLMARK_DIR / "HALLMARK_UNFOLDED_PROTEIN_RESPONSE.txt",
    "HSR_core": _HSR_DIR / "HSR_core.txt",
    "sting_specific_published": _STING_SIG_DIR / "sting_specific_up.txt",
    "ifn_generic_axis": _STING_SIG_DIR / "ifn_only_up.txt",
}

# Expected nominal set sizes — reproduction checks, NOT targets. A disagreement is
# reported and raised, never reconciled to.
EXPECTED_SIZES = {
    "WT_heat_up": 202, "KO_heat_up": 221, "Interaction_up": 7,
    "Interaction_fdrOnly_up": 19, "HSR_core": 56,
    "sting_specific_published": 21, "ifn_generic_axis": 200,
}
EXPECTED_N_CELLS = 99_915

# A set this thin cannot support a per-cell reading; the same power bands stage 12 uses.
GATE_TESTABLE = 15
GATE_UNDERPOWERED = 5

# Seam-guard thresholds. SEAM_R_FLOOR is the loose question: is this substrate on the
# published AUCell scale at all. Below it, every colouring built on the substrate would be
# misleading. SEAM_R_EXACT is the strict one: recomputing AUCell over the same matrix with
# the same gene set is deterministic, and both same-metric references land at r = 1.0
# exactly, so a same-metric row that misses 1.0 by more than this tolerance is a real change
# in the scorer rather than float noise, and the run must stop on it.
SEAM_R_FLOOR = 0.98
SEAM_R_EXACT = 1e-6


def read_gene_list(path: Path) -> list[str]:
    if not path.exists():
        raise FileNotFoundError(f"gene set absent at {path}")
    return [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]


def assign_gate(n_found: int) -> str:
    """Power band on the number of set genes actually present in the object — the
    size that a per-cell score is really computed over, not the nominal size."""
    if n_found >= GATE_TESTABLE:
        return "testable"
    if n_found >= GATE_UNDERPOWERED:
        return "underpowered_reported"
    return "untestable"


def load_gene_sets() -> tuple[dict[str, list[str]], dict[str, str], dict[str, Path]]:
    """Load all 13 sets, returning {name: genes}, {name: kind}, {name: source_path}."""
    gene_sets: dict[str, list[str]] = {}
    kinds: dict[str, str] = {}
    sources: dict[str, Path] = {}
    for kind, spec in (("mouse_derived_arm", MOUSE_ARMS), ("curated_lens", CURATED)):
        for name, path in spec.items():
            gene_sets[name] = read_gene_list(path)
            kinds[name] = kind
            sources[name] = path
    mismatched = {n: (len(gene_sets[n]), e) for n, e in EXPECTED_SIZES.items()
                  if len(gene_sets.get(n, [])) != e}
    if mismatched:
        raise ValueError(
            "nominal gene-set size disagrees with the expected reproduction check "
            f"(observed, expected): {mismatched}. Reported, not reconciled — inspect "
            "the source files before re-running.")
    return gene_sets, kinds, sources


def load_published_embedding() -> pd.DataFrame:
    """The published per-cell readout, indexed by barcode. The embedding coordinates
    and the carried-forward score columns both come from here; no UMAP is recomputed."""
    path = PATHS.interactive_dir() / PUBLISHED_PARQUET
    if not path.exists():
        raise FileNotFoundError(
            f"published per-cell readout absent at {path}; run 08_harvest_readout.py first")
    pub = pd.read_parquet(path)
    missing = [c for c in EMBEDDING_COLS + CARRY_FORWARD if c not in pub.columns]
    if missing:
        raise ValueError(f"{PUBLISHED_PARQUET} is missing expected column(s): {missing}")
    pub["barcode"] = pub["barcode"].astype(str)
    if pub["barcode"].duplicated().any():
        raise ValueError(f"{PUBLISHED_PARQUET} carries duplicate barcodes")
    if len(pub) != EXPECTED_N_CELLS:
        raise ValueError(
            f"{PUBLISHED_PARQUET} has {len(pub)} rows, expected {EXPECTED_N_CELLS}. "
            "Reported, not reconciled.")
    return pub.set_index("barcode", drop=False)


def check_barcode_coverage(adata_barcodes: pd.Index, pub_barcodes: pd.Index) -> None:
    """The h5ad must cover every barcode in the published readout. A missing barcode is
    reported loudly and fatally rather than silently dropping a row from the substrate."""
    missing = pub_barcodes.difference(adata_barcodes)
    extra = adata_barcodes.difference(pub_barcodes)
    n_pub, n_cov = len(pub_barcodes), len(pub_barcodes) - len(missing)
    print(f"[{STAGE}] barcode coverage: {n_cov}/{n_pub} published barcodes present in the "
          f"annotation object ({100.0 * n_cov / n_pub:.4f}%); "
          f"{len(extra)} object barcodes absent from the published readout")
    if len(missing):
        raise ValueError(
            f"{len(missing)} published barcode(s) absent from 02_annotation.h5ad, e.g. "
            f"{list(missing[:5])}. Reported, not silently dropped — the substrate must "
            "cover the published embedding row-for-row.")


def build_manifest(gene_sets, kinds, sources, sym_to_var) -> pd.DataFrame:
    """One row per scored set. `n_genes_found_in_object` is where a reader checks that a
    set is thick enough to carry a reading at all — a set that barely intersects the
    object cannot, whatever its nominal size."""
    rows = []
    for name, genes in gene_sets.items():
        n_set = len(genes)
        n_found = int(sum(g in sym_to_var for g in genes))
        rows.append(dict(
            set_name=name,
            kind=kinds[name],
            source_path=str(sources[name].resolve().relative_to(REPO_ROOT.resolve())),
            n_genes_in_set=n_set,
            n_genes_found_in_object=n_found,
            frac_found=(n_found / n_set if n_set else np.nan),
            metric=METRIC,
            gate=assign_gate(n_found),
        ))
    order = {"mouse_derived_arm": 0, "curated_lens": 1}
    return (pd.DataFrame(rows)
            .sort_values(["kind", "set_name"], key=lambda s: s.map(order) if s.name == "kind" else s)
            .reset_index(drop=True))


def build_substrate(pub: pd.DataFrame, auc: pd.DataFrame) -> pd.DataFrame:
    """The one per-cell substrate: the published embedding + the 13 new AUCell columns +
    the three carried-forward published columns. Row order follows the published readout
    so the substrate and the published figures index the same cells identically."""
    df = pub[EMBEDDING_COLS].copy()
    for col in auc.columns:                                  # <set>_AUCell
        df[col] = auc[col].reindex(df.index).to_numpy()
    for col in CARRY_FORWARD:
        df[f"{CARRY_PREFIX}{col}"] = pub[col].to_numpy()
    return df.reset_index(drop=True)


def summarise(df: pd.DataFrame, set_names: list[str]) -> pd.DataFrame:
    """Mean / median / sd AUCell per (set_name x coarse_label x tissue) with cell counts,
    so the substrate's content is readable without loading the parquet. Descriptive and
    correlative: pooling cells across donors, so the unit here is the cell, not the donor."""
    rows = []
    order = {"Treg": 0, "Tcon": 1, "CD8": 2}
    for (label, tissue), g in df.groupby([COARSE, TISSUE_KEY], observed=True):
        for name in set_names:
            v = g[f"{name}_{METRIC}"].to_numpy(dtype=float)
            v = v[np.isfinite(v)]
            rows.append(dict(
                set_name=name, coarse_label=label, tissue=tissue,
                n_cells=int(len(v)),
                n_donors=int(g[DONOR_KEY].nunique()),
                mean=float(np.mean(v)) if len(v) else np.nan,
                median=float(np.median(v)) if len(v) else np.nan,
                sd=float(np.std(v, ddof=1)) if len(v) > 1 else np.nan,
                metric=METRIC,
                evidence_tier="secondary_percell",
            ))
    out = pd.DataFrame(rows)
    return out.sort_values(
        ["set_name", COARSE, TISSUE_KEY],
        key=lambda s: s.map(order) if s.name == COARSE else s).reset_index(drop=True)


def _corr(a: np.ndarray, b: np.ndarray) -> tuple[float, float, int]:
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3:
        return np.nan, np.nan, int(ok.sum())
    return (float(stats.pearsonr(a[ok], b[ok])[0]),
            float(stats.spearmanr(a[ok], b[ok])[0]),
            int(ok.sum()))


def seam_check(df: pd.DataFrame, stage05: pd.DataFrame) -> pd.DataFrame:
    """Re-derivation check: does this substrate sit on the same scale as the already-
    published per-cell columns? Rows 1-2 compare against genuine published AUCell and are
    expected at r = 1. Row 3 compares against `published_WT_heat_up`, which is a stale
    scanpy `score_genes` module score rather than AUCell, so it is expected to FAIL the
    floor; row 4 is the apples-to-apples anchor for the same mouse arm, against stage 05's
    canonical AUCell. Row 3 is kept and printed on purpose: the log must state why the mouse
    arm has to be coloured with the new column. The frame is consumed by
    assert_seam_reproduces() and is not written to disk."""
    comparisons = [
        ("HALLMARK_HYPOXIA", "HALLMARK_HYPOXIA_AUCell",
         f"{CARRY_PREFIX}HALLMARK_HYPOXIA_AUCell", f"interactive/{PUBLISHED_PARQUET}",
         "AUCell", "same_metric"),
        ("HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE_AUCell",
         f"{CARRY_PREFIX}HALLMARK_UNFOLDED_PROTEIN_RESPONSE_AUCell",
         f"interactive/{PUBLISHED_PARQUET}", "AUCell", "same_metric"),
        ("WT_heat_up", "WT_heat_up_AUCell", f"{CARRY_PREFIX}WT_heat_up",
         f"interactive/{PUBLISHED_PARQUET}", "scanpy_score_genes", "cross_metric"),
        ("WT_heat_up", "WT_heat_up_AUCell", "stage05_WT_heat_up_AUCell",
         "05_scoring/tables/per_cell_scores.csv", "AUCell", "same_metric"),
    ]
    rows = []
    for set_name, new_col, ref_col, ref_source, ref_metric, kind in comparisons:
        new = df[new_col].to_numpy(dtype=float)
        if ref_col.startswith("stage05_"):
            ref = stage05["WT_heat_up_AUCell"].reindex(df["barcode"]).to_numpy(dtype=float)
        else:
            ref = df[ref_col].to_numpy(dtype=float)
        r_p, r_s, n = _corr(new, ref)
        rows.append(dict(
            set_name=set_name, new_column=new_col, reference_column=ref_col,
            reference_source=ref_source, new_metric=METRIC, reference_metric=ref_metric,
            comparison_kind=kind, n_shared_cells=n,
            pearson_r=r_p, spearman_r=r_s, r_floor=SEAM_R_FLOOR,
            passes_floor=bool(np.isfinite(r_p) and r_p >= SEAM_R_FLOOR)))
    return pd.DataFrame(rows)


def assert_seam_reproduces(seam: pd.DataFrame) -> None:
    """Halt the run unless every same-metric row reproduces its reference exactly.

    This is the guard, not a result: it carries no biological reading and writes nothing to
    `03_results/`. A cross-metric row is not evidence that this scorer drifted, so only the
    same-metric rows are held to a threshold. Both thresholds are checked, the loose scale
    floor first so its message is the one a genuinely rescaled substrate raises.
    """
    print(f"\n[{STAGE}] PROVENANCE SEAM GUARD (new vs already-published per-cell columns):")
    print(seam[["set_name", "reference_column", "reference_metric", "comparison_kind",
                "n_shared_cells", "pearson_r", "spearman_r", "passes_floor"]]
          .to_string(index=False))
    same = seam[seam["comparison_kind"] == "same_metric"]
    if not len(same):
        raise ValueError(
            "seam guard found no same_metric row to check — the guard cannot pass vacuously; "
            "inspect seam_check() before trusting this substrate.")
    failed_floor = same[~same["passes_floor"]]
    if len(failed_floor):
        raise ValueError(
            f"SAME-METRIC seam guard FAILED the r >= {SEAM_R_FLOOR} scale floor — this "
            "substrate is not on the published AUCell scale and nothing downstream may be "
            f"built on it:\n{failed_floor.to_string(index=False)}")
    drifted = same[(1.0 - same["pearson_r"] > SEAM_R_EXACT)
                   | (1.0 - same["spearman_r"] > SEAM_R_EXACT)]
    if len(drifted):
        raise ValueError(
            "SAME-METRIC seam guard FAILED: recomputed AUCell no longer reproduces its "
            f"reference to within {SEAM_R_EXACT:g} of r = 1. Scoring the same matrix with the "
            "same gene set is deterministic, so this is a change in the scorer or in its "
            f"input, reported and not reconciled:\n{drifted.to_string(index=False)}")
    print(f"[{STAGE}] all {len(same)} same-metric seam rows reproduce at Pearson and Spearman "
          f"r = 1 to within {SEAM_R_EXACT:g}: the AUCell scorer is stable against the published "
          "readout.")
    for _, r in seam[seam["comparison_kind"] == "cross_metric"].iterrows():
        print(f"[{STAGE}] NOTE — `{r['reference_column']}` is a {r['reference_metric']} score, "
              f"NOT AUCell, so its r = {r['pearson_r']:.6f} against `{r['new_column']}` is a "
              "metric difference and not a scoring drift. Colour the mouse up arm with "
              f"`{r['new_column']}`; `{r['reference_column']}` is carried only to keep the "
              "discrepancy visible.")


def main() -> None:
    n_cores = int(PARAMS.get("percell_score_ncores", 4))

    gene_sets, kinds, sources = load_gene_sets()
    print(f"[{STAGE}] {len(gene_sets)} gene sets (up arms only; no down arm is scored):")
    for name, genes in gene_sets.items():
        print(f"[{STAGE}]   {kinds[name]:17s} {name:38s} {len(genes):4d} genes")

    pub = load_published_embedding()
    print(f"[{STAGE}] published readout: {len(pub)} cells, embedding carried verbatim "
          "(no UMAP recomputed)")

    adata = sc.read_h5ad(PATHS.object("02_annotation"))
    print(f"[{STAGE}] annotation object: {adata.n_obs} cells x {adata.n_vars} genes")
    check_barcode_coverage(pd.Index(adata.obs_names.astype(str)), pd.Index(pub.index))

    # Resolved into this object's symbol vintage here rather than in load_gene_sets(),
    # whose EXPECTED_SIZES check is a fact about the source files and must stay a check on
    # them. Where the sets meet the matrix is where the vintage matters, and the gate band
    # below is read off `n_genes_found_in_object`, so an unresolved count would understate
    # a set's power for a vocabulary reason.
    sym_to_var = _symbol_to_varname(adata, "gene_symbol")
    alias_map = load_alias_map(CONFIG["symbol_alias"]["map_path"])
    for name, genes in list(gene_sets.items()):
        genes_r, applied = resolve_symbols(genes, alias_map, set(sym_to_var))
        gene_sets[name] = genes_r
        if applied:
            print(f"[16_narrative_scoring] {name}: +{len(applied)} via alias "
                  f"({' '.join(f'{a}->{b}' for a, b in applied)})")
    manifest = build_manifest(gene_sets, kinds, sources, sym_to_var)
    thin = manifest[manifest["gate"] != "testable"]
    if len(thin):
        print(f"[{STAGE}] sets below the testable floor ({GATE_TESTABLE} genes found) — "
              "reported WITH their size, not dropped:")
        print(thin[["set_name", "n_genes_in_set", "n_genes_found_in_object", "gate"]]
              .to_string(index=False))

    # --- score all 13 sets in one scorer call, then keep AUCell only ---
    print(f"[{STAGE}] scoring {len(gene_sets)} sets on log-normalised X "
          f"(n_cores={n_cores}) ...")
    scores = score_cells_aucell_ucell(adata, gene_sets, layer=None,
                                      symbol_col="gene_symbol", n_cores=n_cores)
    auc = scores[[f"{n}_{METRIC}" for n in gene_sets]].copy()
    auc.index = auc.index.astype(str)
    for c in auc.columns:
        v = auc[c].to_numpy(dtype=float)
        v = v[np.isfinite(v)]
        if not (v.min() >= 0.0 and v.max() <= 1.0):
            raise ValueError(f"{c} outside [0,1]: {v.min()}..{v.max()} — not an AUCell score")
        print(f"[{STAGE}] {c:48s} [{v.min():.6f}, {v.max():.6f}] mean {v.mean():.6f} OK")
    del adata

    # --- the one substrate ---
    df = build_substrate(pub, auc)
    pq_path = PATHS.interactive_dir() / SUBSTRATE_PARQUET
    df.to_parquet(pq_path, index=False)
    print(f"[{STAGE}] wrote {pq_path.name}: {df.shape[0]} cells x {df.shape[1]} cols")
    print(f"[{STAGE}] columns: {list(df.columns)}")
    if df.shape[0] != EXPECTED_N_CELLS:
        raise ValueError(f"substrate has {df.shape[0]} rows, expected {EXPECTED_N_CELLS}")

    tdir = PATHS.tables(STAGE)
    manifest.to_csv(tdir / "narrative_scoring_manifest.csv", index=False)
    print(f"\n[{STAGE}] narrative_scoring_manifest.csv:")
    print(manifest.to_string(index=False))

    summary = summarise(df, list(gene_sets))
    summary.to_csv(tdir / "narrative_score_summary.csv", index=False)
    print(f"\n[{STAGE}] wrote narrative_score_summary.csv ({summary.shape[0]} rows)")

    # --- seam guard: asserted in-run, published nowhere ---
    stage05 = pd.read_csv(PATHS.tables("05_scoring") / "per_cell_scores.csv", index_col=0)
    stage05.index = stage05.index.astype(str)
    assert_seam_reproduces(seam_check(df, stage05))

    print(f"\n[{STAGE}] COMPUTE DONE. Secondary / annotation tier — never pooled with the "
          "donor-pseudobulk NES spine; no effect-size row written.")


if __name__ == "__main__":
    main()
