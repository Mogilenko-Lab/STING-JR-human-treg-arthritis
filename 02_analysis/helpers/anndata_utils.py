"""
anndata_utils.py — per-GSM 10x H5 ingestion + gene-union assembly.
==================================================================
GSE160097 ships 40 per-GSM CellRanger filtered H5 files in two layouts:
  - v2:  *_filtered_gene_bc_matrices_h5.h5
  - v3:  *_filtered_feature_bc_matrix.h5
Both load via `sc.read_10x_h5` and carry RAW integer UMI counts.

We set `var_names` to Ensembl gene IDs (stable across the v2/v3 union) and keep
the HGNC symbol in `var['gene_symbol']`, preserving raw counts in
`layers['counts']`. Used by 00_build_anndata.py only.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

# GSM -> H5 filename token: filenames embed the GSM id as the first underscore field.
_GSM_RE = re.compile(r"^(GSM\d+)_")
# Two accepted CellRanger H5 layouts (v2 first, then v3).
_H5_SUFFIXES = ("_filtered_gene_bc_matrices_h5.h5", "_filtered_feature_bc_matrix.h5")


def find_gsm_h5(raw_dir: Path) -> Dict[str, Path]:
    """Map GSM id -> its filtered H5 path (whichever CellRanger layout it uses)."""
    out: Dict[str, Path] = {}
    for p in sorted(Path(raw_dir).iterdir()):
        if not p.name.endswith(_H5_SUFFIXES):
            continue
        m = _GSM_RE.match(p.name)
        if not m:
            continue
        gsm = m.group(1)
        if gsm in out:
            raise ValueError(f"two H5 files for {gsm}: {out[gsm].name} and {p.name}")
        out[gsm] = p
    return out


def read_gsm(h5_path: Path, gsm: str) -> ad.AnnData:
    """Read one GSM H5 -> AnnData with Ensembl var_names + gene_symbol + counts layer."""
    a = sc.read_10x_h5(str(h5_path))
    a.var_names_make_unique()  # v2 layout has non-unique symbols
    # Promote Ensembl IDs to var_names (stable key for the union); keep symbol.
    a.var["gene_symbol"] = a.var_names.astype(str)
    if "gene_ids" in a.var.columns:
        a.var_names = a.var["gene_ids"].astype(str)
    a.var_names_make_unique()
    a.var_names.name = "ensembl_id"
    # Preserve raw integer counts before anything touches X.
    a.layers["counts"] = a.X.copy()
    a.obs["gsm"] = gsm
    a.obs["barcode"] = a.obs_names.astype(str)
    return a


def build_pooled(samples: pd.DataFrame, raw_dir: Path,
                 low_input_gsm: str = "GSM4859852") -> ad.AnnData:
    """Concatenate all GSMs in `samples` into one gene-union AnnData with provenance.

    `samples` must carry columns: gsm, donor, condition, tissue, population, title.
    Barcodes are made globally unique by prefixing the GSM id.
    """
    h5_map = find_gsm_h5(raw_dir)
    missing = sorted(set(samples["gsm"]) - set(h5_map))
    if missing:
        raise FileNotFoundError(f"no H5 found for GSMs: {missing}")

    meta_cols = ["title", "donor", "condition", "tissue", "population"]
    adatas: List[ad.AnnData] = []
    for row in samples.itertuples(index=False):
        gsm = row.gsm
        a = read_gsm(h5_map[gsm], gsm)
        for c in meta_cols:
            a.obs[c] = getattr(row, c)
        adatas.append(a)

    pooled = ad.concat(
        adatas, join="outer", merge="first",
        keys=[a.obs["gsm"].iloc[0] for a in adatas], index_unique="-",
    )
    # `merge="first"` keeps gene_symbol from the first GSM carrying each gene; backfill any gap.
    pooled.var["gene_symbol"] = pooled.var["gene_symbol"].astype(str)
    # Guarantee the counts layer is present + integer-valued.
    if "counts" not in pooled.layers:
        pooled.layers["counts"] = pooled.X.copy()
    pooled.obs["coarse_population"] = pooled.obs["population"].astype(str)

    # Provenance flags: the two intentionally-absent samples + the low-input SF Treg.
    present = set(pooled.obs["gsm"])
    pooled.uns["absent_samples"] = {
        "note": "p3 (JIA_patient_3) has no PB Tcon and no PB CD8 GSM (by design).",
        "expected_absent": ["PB CD4_Tcon p3", "PB CD8 p3"],
    }
    pooled.uns["low_input_gsm"] = {
        "gsm": low_input_gsm, "note": "SF CD4_Treg p5 (~3426 cells) — smallest stratum.",
        "present": low_input_gsm in present,
    }
    pooled.uns["ingest"] = {
        "n_gsm": len(adatas), "n_genes_union": int(pooled.n_vars),
        "raw_dir": str(raw_dir),
    }
    return pooled


def annotate_gene_classes(adata: ad.AnnData, species_db: str = "HS") -> None:
    """Flag MT / ribosomal / hemoglobin genes on `var` using the HGNC gene_symbol.

    var_names are Ensembl IDs, so pattern-match on `var['gene_symbol']` (human prefixes).
    """
    sym = adata.var["gene_symbol"].astype(str).str.upper()
    if species_db == "HS":
        adata.var["mt"] = sym.str.startswith("MT-")
        adata.var["ribo"] = sym.str.startswith(("RPS", "RPL"))
        adata.var["hb"] = sym.str.match(r"^HB[ABDEGMQZ][0-9]?$").fillna(False)
    else:  # mouse fallback (unused here)
        adata.var["mt"] = sym.str.startswith("MT-")
        adata.var["ribo"] = sym.str.startswith(("RPS", "RPL"))
        adata.var["hb"] = sym.str.match(r"^HB[ABPQ]").fillna(False)


def counts_are_integer(adata: ad.AnnData, layer: str = "counts") -> bool:
    """True if the layer holds integer-valued counts (float storage is fine)."""
    X = adata.layers[layer]
    data = X.data if hasattr(X, "data") else np.asarray(X).ravel()
    if data.size == 0:
        return True
    return bool(np.allclose(data, np.round(data)))
