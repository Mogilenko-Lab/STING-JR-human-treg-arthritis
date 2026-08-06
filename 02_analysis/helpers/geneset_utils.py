"""
geneset_utils.py — load the frozen Phase-0 signature contract + per-cell scoring.
================================================================================
Phase-0 exports the mouse->human signature at
`mouse_anchor/03_results/human_projection/` (human HGNC symbols). We LOAD it,
never re-derive it. Primary axis = `WT_heat`.

- `load_signature(contract_dir, contrast)` -> {up, down, ranked} human-symbol sets.
- `load_alias_map(path)` / `resolve_symbols(...)` -> the symbol-vintage fix (below).
- `score_cells(...)` -> per-cell up/down module scores via scanpy score_genes.
Used by 05_score_signatures.py.

THE SYMBOL-VINTAGE SEAM. GSE160097 was quantified against a CellRanger hg19
reference, so this compartment's matrix carries that build's HGNC vintage: cGAS is
`MB21D1`, STING is `TMEM173`, MARCHF5 is `MARCH5`, MRE11 is `MRE11A`. Reference sets
ship current symbols, and every match here is an exact string match, so a renamed
gene leaves a set silently and the loss reads as biological absence. Python cannot
reach org.Hs.eg.db, so the resolution is precomputed by
`02_analysis/scripts/00_symbol_alias_map.R` into a committed CSV that both languages
read; pass it through `alias_map=`/`vocabulary=` wherever a reference set meets this
matrix. Only `accepted` pairs are applied — `flagged_for_review` is withheld by
construction, so no consumer can apply a pair a human has not signed off on.
"""
from __future__ import annotations

import os
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Sequence, Set, Tuple

import pandas as pd


def load_signature(contract_dir: Path, contrast: str = "WT_heat",
                   alias_map: Dict[str, str] | None = None,
                   vocabulary: Set[str] | None = None) -> Dict[str, object]:
    """Load one contrast's up/down gene lists + ranked list from the frozen contract.

    `alias_map` + `vocabulary` resolve the up/down lists into this matrix's symbol
    vintage; the applied pairs are returned under `alias_applied` so a caller can
    report the recovery alongside the benefit. The ranked list is left verbatim, being
    the mouse contract's own ordering rather than a set matched against
    this matrix.
    """
    sig_dir = Path(contract_dir) / "signatures" / contrast
    up, up_pairs = resolve_symbols(
        _read_gene_list(sig_dir / f"{contrast}_up.txt"), alias_map, vocabulary)
    down, down_pairs = resolve_symbols(
        _read_gene_list(sig_dir / f"{contrast}_down.txt"), alias_map, vocabulary)
    ranked_path = sig_dir / f"{contrast}_ranked.rnk"
    ranked = None
    if ranked_path.exists():
        ranked = pd.read_csv(ranked_path, sep="\t", header=None, names=["symbol", "t"])
    return {"contrast": contrast, "up": up, "down": down, "ranked": ranked,
            "alias_applied": {"up": up_pairs, "down": down_pairs}}


def load_alias_map(path: str | Path) -> Dict[str, str]:
    """Read the committed reference_symbol -> matrix_symbol map, accepted pairs only.

    Rows carrying any other `resolution` are withheld at this point:
    `flagged_for_review` is a human decision to exclude, and the rejection classes are
    candidates the ownership guard refused (a retired symbol that now names a
    DIFFERENT gene, e.g. ACOD1 -> CAD or IL17F -> IL17A). Applying one would attach
    one gene's expression to another gene's set membership.
    """
    df = pd.read_csv(path)
    for col in ("reference_symbol", "matrix_symbol", "resolution"):
        if col not in df.columns:
            raise ValueError(f"{path} is not a symbol alias map (no `{col}` column)")
    acc = df[df["resolution"] == "accepted"]
    return dict(zip(acc["reference_symbol"].astype(str), acc["matrix_symbol"].astype(str)))


def resolve_symbols(genes: Sequence[str], alias_map: Dict[str, str] | None,
                    vocabulary: Set[str] | None) -> Tuple[List[str], List[Tuple[str, str]]]:
    """Add each gene's matrix-vintage name where the reference name is not in `vocabulary`.

    Only ever ADDS: the reference symbol is kept alongside its resolved twin, so a
    caller's exact-match count cannot move. Order is preserved and the result is
    de-duplicated, because a set carrying both vintages of one gene would otherwise
    weight it twice.

    Returns (resolved genes, [(reference_symbol, matrix_symbol) applied]).
    """
    out: List[str] = []
    seen: Set[str] = set()
    applied: List[Tuple[str, str]] = []
    for g in genes:
        if g in seen:
            continue
        out.append(g)
        seen.add(g)
        if not alias_map or vocabulary is None or g in vocabulary:
            continue
        tgt = alias_map.get(g)
        if tgt is None or tgt not in vocabulary:
            continue
        applied.append((g, tgt))
        if tgt not in seen:
            out.append(tgt)
            seen.add(tgt)
    return out, applied


def _read_gene_list(path: Path) -> List[str]:
    """Read a one-symbol-per-line gene list verbatim; resolution is a separate step.

    Kept unresolved on purpose: a caller has to name the vocabulary it is resolving
    into, and this function cannot know one. Pass the result through
    `resolve_symbols()` at the point the list meets a matrix.
    """
    if not Path(path).exists():
        return []
    with open(path) as fh:
        return [ln.strip() for ln in fh if ln.strip()]


def load_manifest(contract_dir: Path) -> pd.DataFrame:
    return pd.read_csv(Path(contract_dir) / "manifest.csv")


def score_cells(adata, sig: Dict[str, object], layer: str | None = None,
                symbol_col: str = "gene_symbol", seed: int = 0,
                prefix: str | None = None) -> List[str]:
    """DEPRECATED per-cell scorer (scanpy score_genes module score).

    Superseded by `score_cells_aucell_ucell` (rank-based AUCell + UCell, the
    rigorous composition-robust secondary lens). Kept intact so nothing breaks
    until the pipeline scripts switch over; do not use in new code.

    Adds per-cell up/down module scores on human symbols. Scores on `adata.X`
    (assumed log-normalized) unless a `layer` is given. Genes are matched against
    `adata.var[symbol_col]`; only present genes are scored. Adds obs columns
    `<prefix>_up`, `<prefix>_down`, `<prefix>_updown`. Returns the columns written.
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


def score_cells_aucell_ucell(adata, gene_sets: Dict[str, List[str]],
                             layer: str | None = None, symbol_col: str = "gene_symbol",
                             n_cores: int = 4, tmp_dir: str | Path | None = None,
                             keep_tmp: bool = False) -> pd.DataFrame:
    """Per-cell AUCell + UCell scores — the rigorous, composition-robust secondary lens.

    Rank-based scoring (AUCell AUC over the per-cell gene ranking; UCell Mann-Whitney
    U statistic). Both are invariant to per-cell library size and to the exact
    normalization, unlike the mean-centred scanpy `score_genes` module score this
    replaces. This is the SECONDARY tier; donor-pseudobulk fgsea NES stays primary.

    Interop mirrors `fgsea_prerank.R`: the expression matrix and gene sets are
    written to a scratch dir, `percell_score.R` (same helpers/ dir) is called via
    subprocess, and the per-cell CSV is read back.

    Parameters
    ----------
    adata : AnnData with HGNC symbols in `var[symbol_col]`.
    gene_sets : {set_name: [HGNC symbols]}. Set names are preserved verbatim as
        column prefixes. Score up/down sets separately (e.g. "WT_heat_up",
        "WT_heat_down") — AUCell/UCell are unsigned, single-list scorers.
    layer : expression layer to score; None -> `adata.X` (must be LOG-NORMALIZED,
        NOT raw counts). Rank-based scores are monotone-invariant, but stage-05
        scores on log-normalized X for consistency.
    n_cores : parallel workers (AUCell block BPPARAM + UCell ncores).

    Returns
    -------
    DataFrame indexed by `adata.obs_names`, with columns `<set>_AUCell` and
    `<set>_UCell` for each gene set (float in [0, 1]). AUCell is the canonical
    source for `effect_metric=percell_auc_smd`; UCell rides alongside as a
    cross-check. Coverage (n genes matched per set) is logged to stderr.
    """
    import numpy as np
    import scipy.io as sio
    import scipy.sparse as sp

    helper_dir = Path(__file__).resolve().parent
    r_script = helper_dir / "percell_score.R"
    if not r_script.exists():
        raise FileNotFoundError(f"percell_score.R not found at {r_script}")

    sym_to_var = _symbol_to_varname(adata, symbol_col)
    # Restrict to genes carrying a unique HGNC symbol (first occurrence wins),
    # relabelling matrix rows by symbol so the R side matches gene sets directly.
    symbols = list(sym_to_var.keys())
    var_order = [sym_to_var[s] for s in symbols]
    sub = adata[:, var_order]
    X = sub.layers[layer] if layer is not None else sub.X
    # genes x cells (AUCell/UCell convention), sparse CSC
    Xt = sp.csc_matrix(X.T if sp.issparse(X) else np.asarray(X).T)

    ctx = tempfile.TemporaryDirectory(dir=str(tmp_dir) if tmp_dir else None,
                                      prefix="percell_score_")
    work = Path(ctx.name)
    try:
        sio.mmwrite(str(work / "expr.mtx"), Xt, field="real")
        (work / "genes.txt").write_text("\n".join(symbols) + "\n")
        (work / "barcodes.txt").write_text("\n".join(map(str, adata.obs_names)) + "\n")
        set_specs = []
        for name, genes in gene_sets.items():
            present = [g for g in genes if g in sym_to_var]
            sp_path = work / f"set__{name}.txt"
            sp_path.write_text("\n".join(present) + ("\n" if present else ""))
            set_specs.append(f"{name}={sp_path}")
        out_csv = work / "percell_scores.csv"
        rscript = os.environ.get("RSCRIPT", "Rscript")
        cmd = [rscript, str(r_script), str(work / "expr.mtx"), str(work / "genes.txt"),
               str(work / "barcodes.txt"), str(out_csv), str(int(n_cores)), *set_specs]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(
                f"percell_score.R failed (rc={proc.returncode}):\n{proc.stdout}\n{proc.stderr}")
        if proc.stderr.strip():
            print(proc.stderr.strip())
        df = pd.read_csv(out_csv).set_index("cell")
        df.index = df.index.astype(str)
        df = df.reindex(adata.obs_names.astype(str))
        df.index.name = adata.obs_names.name
        return df
    finally:
        if keep_tmp:
            print(f"[score_cells_aucell_ucell] kept scratch at {work}")
        else:
            ctx.cleanup()


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
