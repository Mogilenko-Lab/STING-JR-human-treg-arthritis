"""
geneset_utils.py — load the frozen Phase-0 signature contract + per-cell scoring.
================================================================================
Phase-0 exports the mouse->human signature at
`mouse_anchor/03_results/human_projection/` (human HGNC symbols). We LOAD it,
never re-derive it. Primary axis = `WT_heat`.

- `load_signature(contract_dir, contrast)` -> {up, down, ranked} human-symbol sets.
- `score_cells(...)` -> per-cell up/down module scores via scanpy score_genes.
Used by 05_score_signatures.py.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List

import pandas as pd


def load_signature(contract_dir: Path, contrast: str = "WT_heat") -> Dict[str, object]:
    """Load one contrast's up/down gene lists + ranked list from the frozen contract."""
    sig_dir = Path(contract_dir) / "signatures" / contrast
    up = _read_gene_list(sig_dir / f"{contrast}_up.txt")
    down = _read_gene_list(sig_dir / f"{contrast}_down.txt")
    ranked_path = sig_dir / f"{contrast}_ranked.rnk"
    ranked = None
    if ranked_path.exists():
        ranked = pd.read_csv(ranked_path, sep="\t", header=None, names=["symbol", "t"])
    return {"contrast": contrast, "up": up, "down": down, "ranked": ranked}


def _read_gene_list(path: Path) -> List[str]:
    if not Path(path).exists():
        return []
    with open(path) as fh:
        return [ln.strip() for ln in fh if ln.strip()]


def load_manifest(contract_dir: Path) -> pd.DataFrame:
    return pd.read_csv(Path(contract_dir) / "manifest.csv")


def score_cells(adata, sig: Dict[str, object], layer: str | None = None,
                symbol_col: str = "gene_symbol", seed: int = 0,
                prefix: str | None = None) -> List[str]:
    """Add per-cell up/down module scores (scanpy score_genes) on human symbols.

    Scores on `adata.X` (assumed log-normalized) unless a `layer` is given. Genes
    are matched against `adata.var[symbol_col]`; only present genes are scored.
    Adds obs columns `<prefix>_up`, `<prefix>_down`, `<prefix>_updown`.
    Returns the list of obs columns written.
    """
    import numpy as np
    import scanpy as sc

    prefix = prefix or str(sig["contrast"])
    sym_to_var = _symbol_to_varname(adata, symbol_col)
    written: List[str] = []
    for direction in ("up", "down"):
        genes = [sym_to_var[s] for s in sig[direction] if s in sym_to_var]  # type: ignore[operator]
        col = f"{prefix}_{direction}"
        if len(genes) >= 2:
            kw = {"layer": layer} if layer is not None else {}
            sc.tl.score_genes(adata, gene_list=genes, score_name=col,
                              use_raw=False, random_state=seed, **kw)
        else:
            adata.obs[col] = np.nan
        written.append(col)
    adata.obs[f"{prefix}_updown"] = adata.obs[f"{prefix}_up"] - adata.obs[f"{prefix}_down"]
    written.append(f"{prefix}_updown")
    adata.uns[f"{prefix}_score_coverage"] = {
        "n_up_in_data": int(sum(s in sym_to_var for s in sig["up"])),   # type: ignore[operator]
        "n_up_total": len(sig["up"]),                                    # type: ignore[arg-type]
        "n_down_in_data": int(sum(s in sym_to_var for s in sig["down"])),# type: ignore[operator]
        "n_down_total": len(sig["down"]),                                # type: ignore[arg-type]
    }
    return written


def derive_etreg_from_xlsx(xlsx_path, min_l2fc: float = 1.0, max_p: float = 0.05,
                           n_cap: int = 200) -> Dict[str, object]:
    """Derive the human effector-Treg (eTreg) signature from the GSE161426 log2 matrix.

    Contrast = SF Treg vs PB Treg (Mijnheer/Lutter 2021, Nat Commun 12:2710 — "conserved
    human effector Treg core signature in arthritic joint inflammation"). This is a
    lightweight derivation on the provided log2-normalised table (mean-difference + Welch t)
    intended for EDA scoring; a rigorous version (limma/DESeq2 on the GEO counts, or the
    paper's Supp core-signature list) is the post-go/no-go deliverable.

    Returns {up, down, ranked(gene,t), n_sf, n_pb}. up/down are HGNC symbols.
    """
    import numpy as np
    import pandas as pd
    from scipy import stats

    x = pd.read_excel(xlsx_path)
    sample_cols = list(x.columns[8:])  # first 8 cols are gene + annotation
    sf = [c for c in sample_cols if "SFTreg" in c]
    pb = [c for c in sample_cols if "PBTreg" in c]
    genes = x["gene"].astype(str)
    SF = x[sf].to_numpy(float)
    PB = x[pb].to_numpy(float)
    l2fc = SF.mean(axis=1) - PB.mean(axis=1)
    t, p = stats.ttest_ind(SF, PB, axis=1, equal_var=False)

    d = pd.DataFrame({"gene": genes, "l2fc": l2fc, "t": t, "p": p}).dropna()
    d = d[d["gene"] != "nan"].drop_duplicates("gene")
    sig = d[(d["p"] < max_p) & (d["l2fc"].abs() >= min_l2fc)]
    up = sig[sig["l2fc"] > 0].sort_values("l2fc", ascending=False).head(n_cap)["gene"].tolist()
    down = sig[sig["l2fc"] < 0].sort_values("l2fc").head(n_cap)["gene"].tolist()
    ranked = d.sort_values("t", ascending=False)[["gene", "t"]].reset_index(drop=True)
    return {"up": up, "down": down, "ranked": ranked, "n_sf": len(sf), "n_pb": len(pb)}


def _symbol_to_varname(adata, symbol_col: str) -> Dict[str, str]:
    """Map HGNC symbol -> var_name (Ensembl). First occurrence wins on duplicate symbols."""
    out: Dict[str, str] = {}
    for vn, sym in zip(adata.var_names.astype(str), adata.var[symbol_col].astype(str)):
        if sym and sym != "nan" and sym not in out:
            out[sym] = vn
    return out
