"""
interactive_utils.py — loaders for the jscatter breakpoint explorers.
=====================================================================
The interactive notebooks (`02_analysis/notebooks/*_explore/`) are live-kernel
steering tools: open in VS Code (Quarto + Jupyter) or JupyterLab and run the cells
to explore the embedding with jupyter-scatter (brush, color-by, linked views).

They load COMPACT per-inflection-point tables (2D coords + a few obs columns + a
handful of gene/score columns) written by `scripts/export_explorers.py` to
`03_results/interactive/`. That keeps the interactive kernel light — it never
touches the multi-GB `.h5ad` checkpoints. Nothing here computes biology; it is a
read-only projection for visualization.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd


def find_root(sentinel: str = "02_analysis/config/analysis_config.yaml") -> Path:
    """Walk up from CWD to the compartment root (works from any notebook location)."""
    d = Path.cwd().resolve()
    while not (d / sentinel).exists():
        if d.parent == d:
            raise FileNotFoundError(f"compartment root ({sentinel}) not found above {Path.cwd()}")
        d = d.parent
    return d


def load_explorer(name: str) -> pd.DataFrame:
    """Load a compact explorer table, e.g. load_explorer('01_qc').

    Reads `03_results/interactive/<name>_explore.parquet` (built by
    scripts/export_explorers.py). Index = cell barcode.
    """
    path = find_root() / "03_results" / "interactive" / f"{name}_explore.parquet"
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found — run `python 02_analysis/scripts/export_explorers.py` first.")
    return pd.read_parquet(path)


def first_selection(scatters) -> "np.ndarray":
    """Return the brushed row indices from the FIRST scatter that has a non-empty selection.

    Robust to the common trap where several cells reuse a variable name: register every
    panel in a list and pass it here, so `df.iloc[first_selection(PANELS)]` works no matter
    which grid you lassoed in. Returns an empty int array if nothing is selected.
    """
    import numpy as np
    for sc in scatters:
        try:
            s = np.asarray(sc.selection())
        except Exception:
            continue
        if s.size:
            return s.astype(int)
    return np.array([], dtype=int)


def color_key(df: pd.DataFrame, col: str) -> dict:
    """A jscatter-friendly categorical palette dict for `col` (glasbey if many levels)."""
    import jscatter
    levels = sorted(df[col].astype(str).unique())
    pal = list(jscatter.glasbey_dark) if len(levels) > 8 else list(jscatter.okabe_ito)
    return {lvl: pal[i % len(pal)] for i, lvl in enumerate(levels)}


def _eda_dir(subdir: str = "01_qc_explore") -> Path:
    d = find_root() / "02_analysis" / "notebooks" / subdir / "eda"
    d.mkdir(parents=True, exist_ok=True)
    return d


def save_selection(df: pd.DataFrame, indices, label: str,
                   subdir: str = "01_qc_explore", cols=None) -> Path:
    """Persist a lasso selection's BARCODES + key columns for durable reference.

    Writes `02_analysis/notebooks/<subdir>/eda/selection_<label>.csv` (index = cell
    barcode) and appends a one-row summary to `.../eda/selections_index.csv`. `indices`
    is what `first_selection(PANELS)` / `scatter.selection()` returns.
    """
    import numpy as np
    idx = np.asarray(indices, dtype=int)
    if idx.size == 0:
        raise ValueError("empty selection — lasso a cluster first")
    default = ["x", "y", "population_short", "coarse_label", "tissue", "donor",
               "pct_counts_mt", "n_genes_by_counts", "FOXP3", "CTLA4", "IL2RA",
               "score_HSP", "score_apoptosis", "score_eTreg", "WT_heat_updown"]
    cols = [c for c in (cols or default) if c in df.columns]
    sub = df.iloc[idx]
    out = _eda_dir(subdir) / f"selection_{label}.csv"
    sub[cols].to_csv(out)  # index (barcode) is written

    # append a summary row (composition + means) for quick cross-selection comparison
    summ = {"label": label, "n_cells": len(sub),
            "pct_SF": round((sub.get("tissue") == "synovial_fluid").mean(), 3)
            if "tissue" in sub else None,
            "top_donor": sub["donor"].value_counts().idxmax() if "donor" in sub else None,
            "top_donor_frac": round(sub["donor"].value_counts(normalize=True).max(), 3)
            if "donor" in sub else None}
    for c in ["FOXP3", "CTLA4", "IL2RA", "pct_counts_mt", "score_HSP",
              "score_apoptosis", "score_eTreg", "WT_heat_updown"]:
        if c in sub:
            summ[c] = round(float(sub[c].mean()), 4)
    index_csv = _eda_dir(subdir) / "selections_index.csv"
    prev = pd.read_csv(index_csv) if index_csv.exists() else pd.DataFrame()
    prev = prev[prev.get("label") != label] if "label" in prev else prev
    pd.concat([prev, pd.DataFrame([summ])], ignore_index=True).to_csv(index_csv, index=False)
    print(f"saved {len(sub)} barcodes -> {out.relative_to(find_root())}")
    return out


def snapshot(df: pd.DataFrame, color: str, name: str, indices=None,
             subdir: str = "01_qc_explore", cat: bool | None = None,
             cmap: str = "magma", title: str | None = None):
    """Render a LABELLED matplotlib snapshot of the embedding (legend/axes/title that
    jscatter's viewport export cannot capture) to `.../<subdir>/eda/<name>.{png,pdf}`.

    Colors by `color`; if `indices` is given, non-selected cells are greyed and the
    selection outlined. EDA-tier (not a 03_results deliverable), so it saves directly.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    if cat is None:
        cat = not pd.api.types.is_numeric_dtype(df[color])
    fig, ax = plt.subplots(figsize=(8, 8))
    xy = df[["x", "y"]].to_numpy()
    if indices is not None and np.asarray(indices).size:
        mask = np.zeros(len(df), bool); mask[np.asarray(indices, int)] = True
        ax.scatter(xy[~mask, 0], xy[~mask, 1], s=2, c="lightgrey", linewidths=0)
        base = df.iloc[np.asarray(indices, int)]; pts = xy[mask]
    else:
        base = df; pts = xy
    if cat:
        pal = color_key(df, color)
        for lvl, col in pal.items():
            m = (base[color].astype(str).to_numpy() == lvl)
            ax.scatter(pts[m, 0], pts[m, 1], s=3, c=col, label=str(lvl), linewidths=0)
        ax.legend(markerscale=4, fontsize=8, loc="best", frameon=True, title=color)
    else:
        sctr = ax.scatter(pts[:, 0], pts[:, 1], s=3, c=base[color].to_numpy(), cmap=cmap, linewidths=0)
        fig.colorbar(sctr, ax=ax, shrink=0.7, label=color)
    ax.set_title(title or f"{color}" + (f"  (n={len(base)} selected)" if indices is not None else ""))
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("UMAP-1"); ax.set_ylabel("UMAP-2")
    fig.tight_layout()
    d = _eda_dir(subdir)
    for ext in ("png", "pdf"):
        fig.savefig(d / f"{name}.{ext}", dpi=200, bbox_inches="tight")
    print(f"saved snapshot -> {(d / f'{name}.png').relative_to(find_root())}")
    return fig
