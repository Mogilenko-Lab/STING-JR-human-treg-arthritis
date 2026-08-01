#!/usr/bin/env python
"""
17_treg_reembedding.py — COMPUTE ONLY (no plotting).
=============================================================================
A second map of the sorted JIA T-cell compartment, computed on the Treg gate
alone. The full-object embedding lays out three sorted lineages in one geometry,
so the variance separating Treg from Tcon from CD8 sets the axes and Treg
substructure occupies what room is left. Here the subset gets the whole canvas.

SCOPE AND TIER. Annotation / visualisation only. Nothing here is pooled with the
donor-pseudobulk NES spine and no row is written to
`03_results/master/effect_sizes_treg_arthritis.csv`. An embedding places cells; it
tests nothing. Confirmatory claims stay with donor-level pseudobulk differential
expression.

THE RECIPE IS COPIED FROM THE FULL-OBJECT RUN. `01_qc_filter.py:169-174` is:

    highly_variable_genes(n_top_genes=hvg_n_top)
      -> subset to HVGs -> scale(max_value=10) -> pca(n_comps=n_pcs)
      -> neighbors(n_pcs=n_pcs) -> umap

The same six calls run here on the subset. `hvg_n_top` and `n_pcs` are read from
`thresholds:`, the same config block the full-object run read, so the two maps
cannot drift apart on a parameter. `n_neighbors` stays at the scanpy default
because the full-object call left it there. `random_state` is set explicitly from
`treg_reembedding.random_seed`, which is the value scanpy would otherwise have
used, so fixing it changes no coordinate and makes the fixing visible. The stored
log-normalised `X` is used as it stands: normalisation is per-cell, so it does not
depend on which other cells are in the object, and re-deriving it on the subset
would reproduce the same matrix.

TWO COORDINATE PAIRS, ONE PCA. That recipe applies no batch correction, and on this
subset the uncorrected map resolves donor-by-tissue sample of origin (see the
mixing table: 66% of a cell's 30 nearest neighbours share its donor, against 14.6%
expected). So the neighbours-and-UMAP tail is run twice off the SAME PCA:

  x_uncorrected / y_uncorrected   neighbors(X_pca)         -> umap
  x / y                           neighbors(X_pca_harmony) -> umap

`X_pca_harmony` is `scanpy.external.pp.harmony_integrate` over
`harmony_batch_key`, at harmonypy standard settings, seeded. Harmony is sanctioned
for annotation and visualisation by the umbrella embeddings guardrail; it makes no
claim and its output never reaches a statistical test. Both pairs are carried in
one parquet so they are comparable cell for cell.

SCORES ARE JOINED, NEVER RECOMPUTED. Every score column of the narrative substrate
`03_results/interactive/16_narrative_embedding.parquet` is joined onto the new
coordinates on `barcode`, and the join is asserted complete: row count equals the
subset size, no barcode is missing on either side, and no carried column gains a
NaN. Re-scoring here would put the maps on separate scorer runs and make a colour
difference ambiguous between scoring and layout.

HOW MUCH DOES DONOR STRUCTURE EACH MAP. For each cell, the fraction of its k
nearest neighbours carrying the same `donor` (then the same `tissue`) is averaged
over cells and compared with the fraction expected from the group proportions;
`excess_over_chance` rescales that difference to read 0 at dataset composition and
1 when every neighbour shares the group. The statistic is computed in the 2D map
coordinates for all three maps, because the full-object embedding has no stored
latent space in `02_annotation.h5ad` (`obsm` holds `X_umap_unsupervised` only) and
because 2D is the space a widget draws. Both latent spaces of the subset are
measured as additional rows, which separates "the layout crowds donors together"
from "the donors sit apart in the representation the layout came from".

Outputs:
  03_results/interactive/17_treg_reembedding.parquet          (per-cell widget substrate)
  03_results/objects/17_treg_reembedding.parquet              (coordinate checkpoint)
  03_results/17_treg_reembedding/tables/treg_reembedding_manifest.csv
  03_results/17_treg_reembedding/tables/treg_reembedding_mixing.csv

Run in-container from the compartment root:
    python 02_analysis/scripts/17_treg_reembedding.py
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import harmonypy
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.neighbors import NearestNeighbors

COMPARTMENT_ROOT = Path(__file__).resolve().parents[2]

sys.path.insert(0, str(COMPARTMENT_ROOT))
sys.path.insert(0, str(COMPARTMENT_ROOT / "02_analysis"))
os.chdir(COMPARTMENT_ROOT)

from config import CONFIG, PATHS, PARAMS, DONOR_KEY, TISSUE_KEY  # noqa: E402

STAGE = "17_treg_reembedding"
SOURCE_STAGE_OBJECT = "02_annotation"

# Coordinate pairs carried by the substrate. `x`/`y` is the batch-corrected map, which
# is the pair a widget should draw; the uncorrected pair rides along for comparison.
COORD_CORRECTED = ["x", "y"]
COORD_UNCORRECTED = ["x_uncorrected", "y_uncorrected"]
ALL_COORDS = COORD_CORRECTED + COORD_UNCORRECTED

# Latent-representation column prefixes in the checkpoint.
PC_PREFIX = "PC"
HARMONY_PREFIX = "HPC"

# The columns the substrate carries beyond the coordinates and the scores.
META_COLS = ["coarse_label", TISSUE_KEY, DONOR_KEY, "pct_counts_mt"]

# Embedding labels used in the mixing table.
E_UNCORR, E_HARMONY, E_REF = "treg_only", "treg_only_harmony", "full_object_restricted"

_CFG = CONFIG.get("treg_reembedding", {}) or {}


def cfg(key: str):
    """A stage parameter, or a loud failure. No stage parameter is inlined here."""
    if key not in _CFG:
        raise KeyError(
            f"analysis_config.yaml::treg_reembedding is missing '{key}'; add it there "
            "rather than inlining a value in this script")
    return _CFG[key]


def _pc_names(prefix: str, n: int) -> list[str]:
    return [f"{prefix}{i + 1}" for i in range(n)]


def run_harmony(rep: np.ndarray, obs: pd.DataFrame, batch_key: str, max_iter: int,
                seed: int) -> np.ndarray:
    """Harmony over `batch_key` on a cells x dims representation, returning the corrected
    representation in the same orientation.

    `harmonypy` is called directly rather than through
    `scanpy.external.pp.harmony_integrate`, because that wrapper ends in
    `harmony_out.Z_corr.T`: harmonypy 1.x returned `Z_corr` as dims x cells, harmonypy
    2.0.0 returns it as cells x dims, and the unconditional transpose makes the wrapper
    raise a shape error against 2.0.0. Orientation is resolved here from the shape and
    asserted, so this stage keeps working across either harmonypy generation."""
    n_cells, n_dims = rep.shape
    if n_cells == n_dims:
        raise ValueError(
            f"representation is square ({n_cells}x{n_dims}); the orientation of the "
            "Harmony output cannot be resolved from its shape")
    out = harmonypy.run_harmony(rep.astype(np.float64), obs, batch_key,
                                max_iter_harmony=max_iter, random_state=seed)
    z = np.asarray(out.Z_corr)
    if z.shape == (n_dims, n_cells):
        z = z.T
    if z.shape != (n_cells, n_dims):
        raise ValueError(
            f"Harmony returned Z_corr with shape {np.asarray(out.Z_corr).shape}, which is "
            f"neither {(n_cells, n_dims)} nor its transpose")
    print(f"[{STAGE}] harmonypy {harmonypy.__version__} returned "
          f"{np.asarray(out.Z_corr).shape}, used as {z.shape} (cells x dims)")
    return np.ascontiguousarray(z, dtype=np.float32)


# --------------------------------------------------------------------------- #
# the embedding
# --------------------------------------------------------------------------- #
def recompute_embedding(adata, hvg_n_top: int, n_pcs: int, seed: int, batch_key: str,
                        basis: str, adjusted_basis: str, max_iter: int) -> pd.DataFrame:
    """The `01_qc_filter.py` recipe on whatever cells `adata` holds, with the
    neighbours-and-UMAP tail run twice off the same PCA: once on the uncorrected
    representation and once on the Harmony-corrected one.

    Returns a barcode-indexed frame of both `n_pcs`-dimensional representations and
    both coordinate pairs."""
    emb = adata.copy()
    sc.pp.highly_variable_genes(emb, n_top_genes=hvg_n_top)
    n_hvg = int(emb.var["highly_variable"].sum())
    print(f"[{STAGE}] highly variable genes selected on the subset: {n_hvg} "
          f"(requested {hvg_n_top})")

    emb = emb[:, emb.var["highly_variable"]].copy()
    sc.pp.scale(emb, max_value=10)
    sc.tl.pca(emb, n_comps=n_pcs, random_state=seed)

    # --- variant 1: the recipe as written, no batch correction ---
    sc.pp.neighbors(emb, n_pcs=n_pcs, random_state=seed)
    sc.tl.umap(emb, random_state=seed)
    umap_uncorrected = np.asarray(emb.obsm["X_umap"], dtype=np.float32).copy()

    # --- variant 2: the same PCA, corrected over the batch key, then the same tail ---
    n_batches = int(emb.obs[batch_key].astype(str).nunique())
    print(f"[{STAGE}] Harmony over '{batch_key}' ({n_batches} levels) on {basis}, "
          f"standard settings (max_iter_harmony={max_iter}, random_state={seed}) ...")
    emb.obsm[adjusted_basis] = run_harmony(np.asarray(emb.obsm[basis]), emb.obs,
                                           batch_key, max_iter, seed)
    if emb.obsm[adjusted_basis].shape[1] != n_pcs:
        raise ValueError(
            f"{adjusted_basis} has {emb.obsm[adjusted_basis].shape[1]} dimensions, "
            f"expected {n_pcs}")
    sc.pp.neighbors(emb, use_rep=adjusted_basis, n_pcs=n_pcs, random_state=seed)
    sc.tl.umap(emb, random_state=seed)
    umap_corrected = np.asarray(emb.obsm["X_umap"], dtype=np.float32).copy()

    pcs = np.asarray(emb.obsm[basis], dtype=np.float32)
    hpcs = np.asarray(emb.obsm[adjusted_basis], dtype=np.float32)
    idx = pd.Index(emb.obs_names.astype(str), name="barcode")
    out = pd.DataFrame(pcs, index=idx, columns=_pc_names(PC_PREFIX, n_pcs))
    for j, name in enumerate(_pc_names(HARMONY_PREFIX, n_pcs)):
        out[name] = hpcs[:, j]
    out[COORD_CORRECTED[0]], out[COORD_CORRECTED[1]] = umap_corrected[:, 0], umap_corrected[:, 1]
    out[COORD_UNCORRECTED[0]], out[COORD_UNCORRECTED[1]] = (umap_uncorrected[:, 0],
                                                            umap_uncorrected[:, 1])
    return out


def load_or_build_coords(adata, hvg_n_top: int, n_pcs: int, seed: int, batch_key: str,
                         basis: str, adjusted_basis: str,
                         max_iter: int) -> tuple[pd.DataFrame, str]:
    """Coordinate checkpoint. HVG -> scale -> PCA -> Harmony -> two neighbour graphs ->
    two UMAPs on 27k cells runs well past a minute, so the coordinates are cached. The
    cache is accepted only when it covers exactly the barcodes of the current subset and
    carries both representations at the configured dimensionality and both coordinate
    pairs, so an older cache from a single-variant run is recomputed."""
    ckpt = PATHS.object(str(cfg("checkpoint")))
    want = pd.Index(adata.obs_names.astype(str))
    needed = _pc_names(PC_PREFIX, n_pcs) + _pc_names(HARMONY_PREFIX, n_pcs) + ALL_COORDS
    if ckpt.exists():
        cached = pd.read_parquet(ckpt).set_index("barcode")
        same_cells = len(cached) == len(want) and cached.index.sort_values().equals(
            want.sort_values())
        missing = [c for c in needed if c not in cached.columns]
        if same_cells and not missing:
            print(f"[{STAGE}] reusing coordinate checkpoint {ckpt.name} "
                  f"({len(cached)} cells, both variants present)")
            return cached.loc[want], "checkpoint_reused"
        why = (f"{len(cached)} cells vs {len(want)}" if not same_cells
               else f"missing column(s) {missing[:4]}")
        print(f"[{STAGE}] coordinate checkpoint {ckpt.name} does not match the current "
              f"request ({why}); recomputing")
    coords = recompute_embedding(adata, hvg_n_top, n_pcs, seed, batch_key, basis,
                                 adjusted_basis, max_iter)
    coords.reset_index().to_parquet(ckpt, index=False)
    print(f"[{STAGE}] wrote coordinate checkpoint {ckpt} ({len(coords)} cells)")
    return coords.loc[want], "recomputed"


# --------------------------------------------------------------------------- #
# neighbourhood composition
# --------------------------------------------------------------------------- #
def neighbour_composition(coords: np.ndarray, groups: pd.Series, k: int) -> pd.DataFrame:
    """Per-cell fraction of the k nearest neighbours sharing the cell's group value,
    against the fraction expected from the group proportions.

    `expected` is the same-group fraction a random relabelling would give,
    sum_g p_g * (n_g - 1) / (N - 1), so unequal group sizes are already absorbed.
    `excess_over_chance` = (observed - expected) / (1 - expected) reads 0 when
    neighbourhoods look like the dataset composition and 1 when every neighbour
    shares the group. One row per group value plus an `_all_` row."""
    g = pd.Series(np.asarray(groups, dtype=object)).astype(str).to_numpy()
    n = len(g)
    if n <= k:
        raise ValueError(f"k={k} is not smaller than the {n} cells being measured")

    nn = NearestNeighbors(n_neighbors=k + 1).fit(coords)
    _, idx = nn.kneighbors(coords)
    idx = idx[:, 1:]                                   # drop self
    same = (g[idx] == g[:, None])
    per_cell = same.mean(axis=1)

    counts = pd.Series(g).value_counts()
    p = counts / n
    exp_by_group = (counts - 1) / (n - 1)              # expected for a cell OF that group

    rows = []
    for grp in counts.index:
        m = g == grp
        e = float(exp_by_group[grp])
        obs = float(per_cell[m].mean())
        rows.append(dict(group=str(grp), n_cells=int(m.sum()),
                         group_frac=float(p[grp]), observed_same_frac=obs,
                         expected_same_frac=e,
                         excess_over_chance=(obs - e) / (1.0 - e) if e < 1.0 else np.nan))
    e_all = float((p * exp_by_group).sum())
    obs_all = float(per_cell.mean())
    rows.insert(0, dict(group="_all_", n_cells=int(n), group_frac=1.0,
                        observed_same_frac=obs_all, expected_same_frac=e_all,
                        excess_over_chance=(obs_all - e_all) / (1.0 - e_all)))
    return pd.DataFrame(rows)


def mixing_table(spaces: dict[str, tuple[str, np.ndarray]], obs: pd.DataFrame,
                 keys: list[str], k: int) -> pd.DataFrame:
    """`neighbour_composition` over every (embedding x space x key) combination."""
    rows = []
    for label, (space, coords) in spaces.items():
        for key in keys:
            part = neighbour_composition(coords, obs[key], k)
            part.insert(0, "grouping_key", key)
            part.insert(0, "space", space)
            part.insert(0, "embedding", label)
            part["k"] = k
            rows.append(part)
    out = pd.concat(rows, ignore_index=True)
    out["evidence_tier"] = "annotation_embedding"
    return out


def report_mixing(mix: pd.DataFrame) -> None:
    """Print the three-way comparison in the 2D space every map shares, then state for
    each factor how far the corrected map moved and where it landed relative to chance."""
    k = int(mix["k"].iloc[0])
    two_d = mix[(mix["group"] == "_all_") & (mix["space"] == "umap_2d")]
    print(f"\n[{STAGE}] NEIGHBOURHOOD COMPOSITION at k = {k}, in the 2D coordinates all "
          "three maps share:")
    print(two_d[["embedding", "grouping_key", "n_cells", "observed_same_frac",
                 "expected_same_frac", "excess_over_chance"]].to_string(index=False))
    print(f"\n[{STAGE}] latent-space rows (the representations the layouts came from):")
    lat = mix[(mix["group"] == "_all_") & (mix["space"] != "umap_2d")]
    print(lat[["embedding", "space", "grouping_key", "observed_same_frac",
               "expected_same_frac", "excess_over_chance"]].to_string(index=False))

    def _row(label: str, key: str):
        r = two_d[(two_d["embedding"] == label) & (two_d["grouping_key"] == key)]
        return None if not len(r) else r.iloc[0]

    print()
    for key in two_d["grouping_key"].unique():
        unc, har, ref = _row(E_UNCORR, key), _row(E_HARMONY, key), _row(E_REF, key)
        if unc is None or har is None or ref is None:
            continue
        exp = float(unc["expected_same_frac"])
        print(f"[{STAGE}] {key}: same-{key} neighbour fraction {float(unc['observed_same_frac']):.3f} "
              f"uncorrected, {float(har['observed_same_frac']):.3f} Harmony-corrected, "
              f"{float(ref['observed_same_frac']):.3f} full-object restricted; "
              f"{exp:.3f} expected from the {key} proportions "
              f"(excess over chance {float(unc['excess_over_chance']):.3f} / "
              f"{float(har['excess_over_chance']):.3f} / {float(ref['excess_over_chance']):.3f}).")
        # Where the corrected map lands relative to chance and to the two references.
        share = float(har["observed_same_frac"]) / exp if exp > 0 else np.nan
        print(f"[{STAGE}]   corrected map sits at {share:.2f}x the chance fraction for "
              f"{key}; {float(har['observed_same_frac']) - float(unc['observed_same_frac']):+.3f} "
              f"against the uncorrected map and "
              f"{float(har['observed_same_frac']) - float(ref['observed_same_frac']):+.3f} "
              "against the full-object map.")


# --------------------------------------------------------------------------- #
# the substrate
# --------------------------------------------------------------------------- #
def load_source_substrate(subset_key: str, subset_value: str) -> tuple[pd.DataFrame, list[str]]:
    """The narrative substrate, restricted to the subset. Returns the frame indexed by
    barcode plus the list of its score columns (everything that is not a coordinate or
    a carried metadata column), so the widget can colour by any set the main map can."""
    path = PATHS.interactive_dir() / str(cfg("source_parquet"))
    if not path.exists():
        raise FileNotFoundError(
            f"per-cell score substrate absent at {path}; run 16_narrative_scoring.py first")
    src = pd.read_parquet(path)
    needed = ["barcode"] + COORD_CORRECTED + META_COLS
    missing = [c for c in needed if c not in src.columns]
    if missing:
        raise ValueError(f"{path.name} is missing expected column(s): {missing}")
    src["barcode"] = src["barcode"].astype(str)
    if src["barcode"].duplicated().any():
        raise ValueError(f"{path.name} carries duplicate barcodes")
    score_cols = [c for c in src.columns if c not in needed]
    print(f"[{STAGE}] score substrate {path.name}: {len(src)} cells, "
          f"{len(score_cols)} score columns carried forward")
    sub = src[src[subset_key].astype(str) == subset_value].copy()
    print(f"[{STAGE}] restricted to {subset_key} == '{subset_value}': {len(sub)} cells")
    return sub.set_index("barcode", drop=False), score_cols


def build_substrate(coords: pd.DataFrame, src: pd.DataFrame,
                    score_cols: list[str]) -> pd.DataFrame:
    """Both coordinate pairs + the carried metadata and scores, joined on barcode. The
    join is asserted complete on both sides and checked for NaNs it could introduce."""
    only_coords = coords.index.difference(src.index)
    only_src = src.index.difference(coords.index)
    if len(only_coords) or len(only_src):
        raise ValueError(
            f"barcode join is incomplete: {len(only_coords)} embedded cell(s) absent from "
            f"the score substrate (e.g. {list(only_coords[:5])}), {len(only_src)} substrate "
            f"cell(s) absent from the embedding (e.g. {list(only_src[:5])}). Reported, not "
            "silently dropped.")

    carried = META_COLS + score_cols
    df = pd.DataFrame({"barcode": coords.index.astype(str)})
    for c in ALL_COORDS:
        df[c] = coords[c].to_numpy(dtype=np.float32)
    aligned = src.reindex(coords.index)
    for c in carried:
        df[c] = aligned[c].to_numpy()

    nan_before = src[carried].isna().sum()
    nan_after = df[carried].isna().sum()
    gained = {c: int(nan_after[c] - nan_before[c]) for c in carried
              if nan_after[c] > nan_before[c]}
    if gained:
        raise ValueError(f"the join introduced NaNs in {gained}; reported, not filled")
    print(f"[{STAGE}] join complete: {len(df)} rows, no NaN introduced "
          f"({int(nan_before.sum())} NaN already present upstream)")
    return df[["barcode"] + ALL_COORDS + carried]


def build_manifest(n_cells: int, n_expected: int, hvg_n_top: int, n_pcs: int, seed: int,
                   k: int, coord_state: str, score_cols: list[str], subset_key: str,
                   subset_value: str, batch_key: str, max_iter: int,
                   adjusted_basis: str) -> pd.DataFrame:
    """One row per recorded parameter, so the substrate's provenance is readable without
    the config file open beside it."""
    src = str(cfg("source_parquet"))
    rows = [
        ("subset_key", subset_key, "treg_reembedding.subset_key"),
        ("subset_value", subset_value, "treg_reembedding.subset_value"),
        ("n_cells_embedded", n_cells, "derived"),
        ("n_cells_expected", n_expected, "treg_reembedding.expected_n_cells"),
        ("hvg_n_top", hvg_n_top, "thresholds.hvg_n_top"),
        ("n_pcs", n_pcs, "thresholds.n_pcs"),
        ("n_neighbors", "scanpy default (unset, as in 01_qc_filter.py)", "recipe"),
        ("scale_max_value", 10, "recipe (01_qc_filter.py)"),
        ("random_seed", seed, "treg_reembedding.random_seed"),
        ("harmony_batch_key", batch_key, "treg_reembedding.harmony_batch_key"),
        ("harmony_max_iter", max_iter, "treg_reembedding.harmony_max_iter"),
        ("harmony_adjusted_basis", adjusted_basis,
         "treg_reembedding.harmony_adjusted_basis"),
        ("corrected_coord_columns", "/".join(COORD_CORRECTED), "derived"),
        ("uncorrected_coord_columns", "/".join(COORD_UNCORRECTED), "derived"),
        ("mixing_k", k, "treg_reembedding.mixing_k"),
        ("reference_embedding_key", str(cfg("reference_embedding_key")),
         "treg_reembedding.reference_embedding_key"),
        ("source_object", f"03_results/objects/{SOURCE_STAGE_OBJECT}.h5ad", "PATHS.object"),
        ("score_source_parquet", f"03_results/interactive/{src}",
         "treg_reembedding.source_parquet"),
        ("n_score_columns_carried", len(score_cols), "derived"),
        ("coordinate_state", coord_state, "derived"),
        ("evidence_tier", "annotation_embedding", "fixed"),
    ]
    return pd.DataFrame(rows, columns=["parameter", "value", "source"])


# --------------------------------------------------------------------------- #
def main() -> None:
    subset_key = str(cfg("subset_key"))
    subset_value = str(cfg("subset_value"))
    n_expected = int(cfg("expected_n_cells"))
    seed = int(cfg("random_seed"))
    k = int(cfg("mixing_k"))
    mixing_keys = [str(x) for x in cfg("mixing_keys")]
    ref_key = str(cfg("reference_embedding_key"))
    batch_key = str(cfg("harmony_batch_key"))
    basis = str(cfg("harmony_basis"))
    adjusted_basis = str(cfg("harmony_adjusted_basis"))
    max_iter = int(cfg("harmony_max_iter"))
    hvg_n_top = int(PARAMS.hvg_n_top)
    n_pcs = int(PARAMS.n_pcs)

    print(f"[{STAGE}] recipe parameters: hvg_n_top={hvg_n_top} (thresholds), "
          f"n_pcs={n_pcs} (thresholds), random_seed={seed}, scale max_value=10, "
          "n_neighbors at the scanpy default — the same values the full-object "
          "embedding used")
    print(f"[{STAGE}] two coordinate pairs off one PCA: {'/'.join(COORD_UNCORRECTED)} "
          f"uncorrected, {'/'.join(COORD_CORRECTED)} after Harmony over '{batch_key}'")

    src, score_cols = load_source_substrate(subset_key, subset_value)

    adata = sc.read_h5ad(PATHS.object(SOURCE_STAGE_OBJECT))
    print(f"[{STAGE}] annotation object: {adata.n_obs} cells x {adata.n_vars} genes")
    if ref_key not in adata.obsm:
        raise ValueError(f"{ref_key} absent from obsm; available: {list(adata.obsm)}")
    if batch_key not in adata.obs:
        raise ValueError(f"batch key '{batch_key}' absent from obs")
    keep = adata.obs[subset_key].astype(str) == subset_value
    adata = adata[keep.to_numpy()].copy()
    print(f"[{STAGE}] subset {subset_key} == '{subset_value}': {adata.n_obs} cells")
    if adata.n_obs != n_expected:
        raise ValueError(
            f"subset holds {adata.n_obs} cells, expected {n_expected}. Reported, not "
            "reconciled — check the frozen labels before re-running.")

    # --- the two new maps (checkpointed) + the full-object map on the same cells ---
    coords, coord_state = load_or_build_coords(adata, hvg_n_top, n_pcs, seed, batch_key,
                                               basis, adjusted_basis, max_iter)
    ref_2d = pd.DataFrame(
        np.asarray(adata.obsm[ref_key], dtype=np.float32)[:, :2],
        index=pd.Index(adata.obs_names.astype(str), name="barcode"),
        columns=COORD_CORRECTED)
    obs = adata.obs.set_index(pd.Index(adata.obs_names.astype(str))).loc[coords.index]
    ref_2d = ref_2d.loc[coords.index]
    del adata

    # --- the widget substrate ---
    df = build_substrate(coords, src, score_cols)
    pq_path = PATHS.interactive_dir() / str(cfg("output_parquet"))
    df.to_parquet(pq_path, index=False)
    print(f"[{STAGE}] wrote {pq_path.name}: {df.shape[0]} cells x {df.shape[1]} cols")
    print(f"[{STAGE}] columns: {list(df.columns)}")
    if df.shape[0] != n_expected:
        raise ValueError(f"substrate has {df.shape[0]} rows, expected {n_expected}")

    # --- how much donor structures each map ---
    spaces = {
        E_HARMONY: ("umap_2d", df[COORD_CORRECTED].to_numpy(dtype=np.float32)),
        E_UNCORR: ("umap_2d", df[COORD_UNCORRECTED].to_numpy(dtype=np.float32)),
        E_REF: ("umap_2d", ref_2d.to_numpy(dtype=np.float32)),
        f"{E_UNCORR}_latent": (f"pca_{n_pcs}",
                               coords[_pc_names(PC_PREFIX, n_pcs)].to_numpy(dtype=np.float32)),
        f"{E_HARMONY}_latent": (f"harmony_{n_pcs}",
                                coords[_pc_names(HARMONY_PREFIX, n_pcs)].to_numpy(dtype=np.float32)),
    }
    mix = mixing_table(spaces, obs, mixing_keys, k)

    tdir = PATHS.tables(STAGE)
    mix.to_csv(tdir / "treg_reembedding_mixing.csv", index=False)
    print(f"\n[{STAGE}] wrote treg_reembedding_mixing.csv ({len(mix)} rows)")
    report_mixing(mix)

    manifest = build_manifest(len(df), n_expected, hvg_n_top, n_pcs, seed, k, coord_state,
                              score_cols, subset_key, subset_value, batch_key, max_iter,
                              adjusted_basis)
    manifest.to_csv(tdir / "treg_reembedding_manifest.csv", index=False)
    print(f"\n[{STAGE}] treg_reembedding_manifest.csv:")
    print(manifest.to_string(index=False))

    print(f"\n[{STAGE}] COMPUTE DONE. Annotation / visualisation tier — never pooled with "
          "the donor-pseudobulk NES spine; no effect-size row written.")


if __name__ == "__main__":
    main()
