#!/usr/bin/env python
"""
11_heat_decomposition.py — COMPUTE ONLY. Where does each part of the mouse
39 C-derived signature sit in each population's ranking?
=========================================================================
The JIA SF-vs-PB enrichment of the mouse 39 °C `WT_heat_up` set survives the hypoxia purge and
is carried by 199 genes doing many different things. This stage splits those 199 (and the 94
down genes) into subcomponents defined by CURATED, VERSIONED, ANCHOR-INDEPENDENT gene sets,
then scores every subcomponent against the same donor-pseudobulk ranked lists the whole
signature was scored on. The question is which parts of the mouse program carry the
synovial-fluid shift.

Two design decisions are load-bearing:

1. The partition uses curated MSigDB sets plus the frozen HSR core, keeping the
   `WT_heat_up` leading-edge taxonomy out of it. That taxonomy covers the 66 genes that are
   the union of the three populations' leading edges, so testing subsets of it would test
   enrichment of genes selected because they enriched. It is a post-hoc annotation of a
   result.
2. Subcomponents OVERLAP — a gene may sit in two curated programs — and the overlap is kept.
   Any priority-ordered disjoint partition is an extra assumption that silently decides which
   program gets credit for a shared gene. `decomposition_gene_assignment.csv` records every
   gene's full membership so the sharing is auditable, and the subset of each mouse arm that
   no curated set claims is reported as its own `unassigned` subcomponent.

Inputs:
  - 03_results/09_heat_hypoxia/tables/_signatures_full/WT_heat_{up,down}.txt
  - 00_data/references/msigdb_hallmark/HALLMARK_*.txt
  - 00_data/references/temp_hsr_lens/HSR_core.txt
  - 03_results/03_pseudobulk/tables/ranked_{treg,tcon,cd8}.tsv
  - ../sting_positive_control/03_results/06_reference_axis/signatures/sting_specific_up.txt

Outputs (all under 03_results/11_heat_decomposition/tables/):
  - _signatures_decomp/<subcomponent>_<arm>.txt
  - decomposition_overlap.csv
  - decomposition_gene_assignment.csv
  - decomposition_assignment_multiplicity.csv
  - decomp_gsea_{treg,tcon,cd8}.csv (+ .rds + runsum_interactive_*.csv)
  - decomposition_nes.csv
  - sting_axis_overlap.csv

Tier note: SECONDARY / annotation tier, firewalled from the confirmatory `WT_heat` claim
spine. No row is written to effect_sizes_treg_arthritis.csv and no row is appended to any
03_results/master/ accumulator.

Size floor: a subcomponent whose testable size (its intersection with a ranked list) falls
below `gsea_min_size` gets NO NES. It is reported as untestable WITH its size and the reason,
since silent truncation would read as full coverage.

Sign convention: NES > 0 means the subcomponent is enriched toward genes up in synovial fluid
versus paired peripheral blood.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "02_analysis"))
os.chdir(ROOT)

from config import PARAMS, PATHS  # noqa: E402
from helpers.source_hash_manifest import verify_source_hash  # noqa: E402

STAGE = "11_heat_decomposition"
PRIMARY = "WT_heat"
POP_TAG = {"Treg": "treg", "Tcon": "tcon", "CD8": "cd8"}
ARMS = ("up", "down")
FGSEA_R = "02_analysis/helpers/fgsea_prerank.R"

# --- reference gene-set locations (frozen, reproducible, committed reproducers) ---
HALLMARK_DIR = Path("00_data/references/msigdb_hallmark")
HSR_DIR = Path("00_data/references/temp_hsr_lens")
# The mouse arms exactly as handed to fgsea for the whole-signature run, so this
# decomposition and the whole-signature NES share one gene universe.
MOUSE_ARM_DIR = PATHS.tables("09_heat_hypoxia") / "_signatures_full"
MOUSE_ARM_SIZE = {"up": 202, "down": 96}

# The published de Cevins Table S6 IFN-independent STING signature, from the STING
# positive-control compartment. Carried as a gene-overlap TALLY only: its overlap with
# the mouse arms is far below the size floor, so fabricating a GSEA arm for it would be
# a statistic with no support.
STING_SIG_PATH = Path("../sting_positive_control/03_results/06_reference_axis/"
                      "signatures/sting_specific_up.txt")
STING_SIG_LABEL = "de_Cevins_sting_specific_up"

# The curated presumptions. Each is a versioned, anchor-independent public gene set;
# the mouse signature knows nothing about any of them, so an intersection is a genuine
# partition of the signature, independent of any result.
PROGRAMS: list[tuple[str, str, Path, str]] = [
    ("hsr_curated", "HSR_core", HSR_DIR / "HSR_core.txt",
     "curated heat-shock-response core, Reactome/GO-derived"),
    ("upr_er", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
     HALLMARK_DIR / "HALLMARK_UNFOLDED_PROTEIN_RESPONSE.txt",
     "ER-side proteostasis"),
    ("hypoxia", "HALLMARK_HYPOXIA", HALLMARK_DIR / "HALLMARK_HYPOXIA.txt",
     "MSigDB Hallmark hypoxia response"),
    ("nfkb_tnfa", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
     HALLMARK_DIR / "HALLMARK_TNFA_SIGNALING_VIA_NFKB.txt",
     "TNFA / NF-kB inflammatory signalling"),
    ("ifn_type_i", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
     HALLMARK_DIR / "HALLMARK_INTERFERON_ALPHA_RESPONSE.txt",
     "type-I interferon, the cGAS/STING-adjacent arm"),
    ("inflammatory", "HALLMARK_INFLAMMATORY_RESPONSE",
     HALLMARK_DIR / "HALLMARK_INFLAMMATORY_RESPONSE.txt",
     "broad inflammatory response"),
    ("t_activation", "HALLMARK_IL2_STAT5_SIGNALING",
     HALLMARK_DIR / "HALLMARK_IL2_STAT5_SIGNALING.txt",
     "T-cell activation / IL2-STAT5 growth signalling"),
]
RESIDUAL = "unassigned"
RESIDUAL_DESC = "the part of the mouse arm no curated presumption claims"
SUBCOMPONENTS = [p[0] for p in PROGRAMS] + [RESIDUAL]
DATABASE_BY_SUBCOMPONENT = {
    "hsr_curated": "curated_hsr_reactome_go",
    "upr_er": "msigdb_hallmark",
    "hypoxia": "msigdb_hallmark",
    "nfkb_tnfa": "msigdb_hallmark",
    "ifn_type_i": "msigdb_hallmark",
    "inflammatory": "msigdb_hallmark",
    "t_activation": "msigdb_hallmark",
    RESIDUAL: "mouse_projection",
}


def read_gene_list(path: Path) -> list[str]:
    if not path.exists():
        raise FileNotFoundError(f"gene list absent: {path}")
    return [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]


def write_gene_list(path: Path, genes: Iterable[str]) -> None:
    """One symbol per line. An empty part writes a genuinely empty file, zero bytes."""
    path.parent.mkdir(parents=True, exist_ok=True)
    listed = list(genes)
    path.write_text("\n".join(listed) + "\n" if listed else "")


def load_mouse_arms() -> dict[str, list[str]]:
    """The full mouse arms, with a size guard so a contract change is loud."""
    arms = {}
    for arm in ARMS:
        genes = read_gene_list(MOUSE_ARM_DIR / f"{PRIMARY}_{arm}.txt")
        want = MOUSE_ARM_SIZE[arm]
        if len(genes) != want:
            raise ValueError(
                f"{PRIMARY}_{arm} size drift: expected {want} genes, found {len(genes)} "
                f"in {MOUSE_ARM_DIR}. Re-run 09_heat_hypoxia.py before decomposing."
            )
        arms[arm] = genes
    return arms


def load_curated_sets() -> dict[str, set[str]]:
    return {key: set(read_gene_list(path)) for key, _, path, _ in PROGRAMS}


# ===========================================================================
# 1. Partition each mouse arm by intersection with the curated presumptions
# ===========================================================================
def build_subsignatures(arms: dict[str, list[str]],
                        curated: dict[str, set[str]]) -> dict[tuple[str, str], list[str]]:
    """Return {(subcomponent, arm): genes}. Subcomponents overlap; `unassigned` does not."""
    subs: dict[tuple[str, str], list[str]] = {}
    for arm, arm_genes in arms.items():
        claimed: set[str] = set()
        for key, _, _, _ in PROGRAMS:
            hit = [g for g in arm_genes if g in curated[key]]
            subs[(key, arm)] = sorted(hit)
            claimed.update(hit)
        subs[(RESIDUAL, arm)] = sorted(g for g in arm_genes if g not in claimed)
    return subs


def prepare_signature_dirs(tables_dir: Path,
                           subs: dict[tuple[str, str], list[str]]) -> Path:
    """Materialise every sub-signature as its own one-symbol-per-line list."""
    sig_dir = tables_dir / "_signatures_decomp"
    if sig_dir.exists():
        shutil.rmtree(sig_dir)
    sig_dir.mkdir(parents=True, exist_ok=True)
    for (key, arm), genes in subs.items():
        write_gene_list(sig_dir / f"{key}_{arm}.txt", genes)
    return sig_dir


def overlap_tallies(tables_dir: Path, arms: dict[str, list[str]],
                    curated: dict[str, set[str]],
                    subs: dict[tuple[str, str], list[str]]) -> pd.DataFrame:
    """Plain, non-exclusive tallies: how much of each mouse arm each presumption claims."""
    desc = {key: d for key, _, _, d in PROGRAMS}
    label = {key: lab for key, lab, _, _ in PROGRAMS}
    rows = []
    for arm in ARMS:
        n_arm = len(arms[arm])
        for key, _, _, _ in PROGRAMS:
            genes = subs[(key, arm)]
            rows.append({
                "subcomponent": key,
                "curated_set": label[key],
                "role": desc[key],
                "mouse_arm": f"{PRIMARY}_{arm}",
                "n_curated_set": len(curated[key]),
                "n_mouse_arm": n_arm,
                "n_intersect": len(genes),
                "frac_of_mouse_arm": len(genes) / n_arm if n_arm else np.nan,
                "frac_of_curated_set": len(genes) / len(curated[key]) if curated[key] else np.nan,
                "genes": ";".join(genes),
            })
        res = subs[(RESIDUAL, arm)]
        rows.append({
            "subcomponent": RESIDUAL,
            "curated_set": "(none)",
            "role": RESIDUAL_DESC,
            "mouse_arm": f"{PRIMARY}_{arm}",
            "n_curated_set": np.nan,
            "n_mouse_arm": n_arm,
            "n_intersect": len(res),
            "frac_of_mouse_arm": len(res) / n_arm if n_arm else np.nan,
            "frac_of_curated_set": np.nan,
            "genes": ";".join(res),
        })
    out = pd.DataFrame(rows)
    out.to_csv(tables_dir / "decomposition_overlap.csv", index=False)
    return out


def gene_assignment(tables_dir: Path, arms: dict[str, list[str]],
                    curated: dict[str, set[str]]) -> pd.DataFrame:
    """One row per mouse-signature gene, listing every subcomponent that claims it.

    This is what makes the overlapping partition auditable: a gene in two curated
    programs appears in both arms of the decomposition, and the reader can see it.
    """
    rows = []
    for arm in ARMS:
        for gene in arms[arm]:
            hits = [key for key, _, _, _ in PROGRAMS if gene in curated[key]]
            rows.append({
                "gene": gene,
                "mouse_arm": f"{PRIMARY}_{arm}",
                "n_subcomponents": len(hits),
                "subcomponents": ";".join(hits) if hits else RESIDUAL,
            })
    out = pd.DataFrame(rows).sort_values(["mouse_arm", "gene"]).reset_index(drop=True)
    out.to_csv(tables_dir / "decomposition_gene_assignment.csv", index=False)
    return out


def assignment_multiplicity(tables_dir: Path, assignment: pd.DataFrame) -> pd.DataFrame:
    """One row per arm: how far the overlapping assignment is from a partition.

    The coverage figure draws one bar per curated set, and overlapping bars read
    as a partition unless the panel says otherwise. The quantity that shows they
    are not is the multiplicity: how many of the claimed genes carry more than one
    assignment, and therefore how many more claims exist than claimed genes. Summing
    the bars double-counts exactly that excess.

    Computed here, so the number the figure prints is read from a committed table.
    """
    rows = []
    for arm in ARMS:
        sub = assignment[assignment["mouse_arm"] == f"{PRIMARY}_{arm}"]
        mult = sub["n_subcomponents"].astype(int)
        claimed = mult > 0
        rows.append({
            "mouse_arm": f"{PRIMARY}_{arm}",
            "arm": arm,
            "n_arm": int(len(sub)),
            "n_unassigned": int((~claimed).sum()),
            "n_claimed": int(claimed.sum()),
            "n_claimed_once": int((mult == 1).sum()),
            "n_claimed_multiply": int((mult >= 2).sum()),
            "max_subcomponents_per_gene": int(mult.max()) if len(mult) else 0,
            "n_claims_total": int(mult.sum()),
            "n_excess_claims": int(mult.sum() - claimed.sum()),
            "is_partition": bool(int((mult >= 2).sum()) == 0),
            "evidence_tier": "secondary_annotation",
        })
    out = pd.DataFrame(rows)
    out.to_csv(tables_dir / "decomposition_assignment_multiplicity.csv", index=False)
    return out


# ===========================================================================
# 2. Score every subcomponent against the same donor-pseudobulk ranked lists
# ===========================================================================
def run_fgsea(ranked_path: Path, out_csv: Path, contrast: str, sig_dir: Path,
              set_names: list[str]) -> pd.DataFrame:
    cmd = [
        "Rscript",
        FGSEA_R,
        str(ranked_path),
        str(out_csv),
        contrast,
        str(PARAMS.gsea_min_size),
        str(PARAMS.gsea_max_size),
        str(PARAMS.gsea_seed),
        str(PARAMS.gsea_nperm),
    ] + [
        f"{nm}:{DATABASE_BY_SUBCOMPONENT[nm.rsplit('_', 1)[0]]}={sig_dir / f'{nm}.txt'}"
        for nm in set_names
    ]
    subprocess.run(cmd, check=True)
    return pd.read_csv(out_csv)


def leading_edge_from_row(row: pd.Series) -> str:
    for col in ("leading_edge", "core_enrichment"):
        if col in row and pd.notna(row[col]):
            return str(row[col]).replace("/", ";")
    return ""


def decomposition_nes(tables_dir: Path, sig_dir: Path,
                      subs: dict[tuple[str, str], list[str]]) -> pd.DataFrame:
    """Run one fgsea per population over every non-empty sub-signature, then report.

    Every requested subcomponent gets a row, whether or not fgsea could score it.
    `testable` is False when the subcomponent has no genes in the mouse arm at all or
    when its intersection with the ranked list falls under `gsea_min_size`; the reason
    and the size are carried in the row.
    """
    min_size = int(PARAMS.gsea_min_size)
    label = {key: lab for key, lab, _, _ in PROGRAMS}
    label[RESIDUAL] = "(none)"

    # Empty sub-signatures cannot be handed to fgsea at all.
    requested = [(key, arm) for arm in ARMS for key in SUBCOMPONENTS if subs[(key, arm)]]
    empty = [(key, arm) for arm in ARMS for key in SUBCOMPONENTS if not subs[(key, arm)]]
    for key, arm in empty:
        print(f"[11_heat_decomposition] {key}_{arm}: 0 genes in {PRIMARY}_{arm} — "
              "not passed to fgsea, reported as untestable")
    set_names = [f"{key}_{arm}" for key, arm in requested]

    ranked_dir = PATHS.tables("03_pseudobulk")
    rows = []
    for pop, tag in POP_TAG.items():
        ranked_path = ranked_dir / f"ranked_{tag}.tsv"
        if not ranked_path.exists():
            print(f"[11_heat_decomposition] {pop}: no ranked list at {ranked_path} — skipping")
            continue
        gsea = run_fgsea(
            ranked_path,
            tables_dir / f"decomp_gsea_{tag}.csv",
            f"SF_vs_PB_{pop}_heat_decomposition",
            sig_dir,
            set_names,
        ).set_index("pathway_id")

        for arm in ARMS:
            for key in SUBCOMPONENTS:
                genes = subs[(key, arm)]
                sid = f"{key}_{arm}"
                if sid in gsea.index:
                    r = gsea.loc[sid]
                    in_ranked = int(r["set_size"])
                    nes, pval, padj = r["nes"], r["pvalue"], r["padj"]
                    le = leading_edge_from_row(r)
                else:
                    in_ranked, nes, pval, padj, le = 0, np.nan, np.nan, np.nan, ""
                testable = bool(pd.notna(nes))
                if testable:
                    reason = ""
                elif not genes:
                    reason = f"no {PRIMARY}_{arm} gene belongs to this curated set"
                else:
                    reason = (f"{in_ranked} of {len(genes)} gene(s) present in the ranked "
                              f"list, below the gsea_min_size floor of {min_size}")
                    print(f"[11_heat_decomposition] {pop} {sid}: untestable — {reason}")
                rows.append({
                    "population": pop,
                    "mouse_arm": f"{PRIMARY}_{arm}",
                    "subcomponent": key,
                    "curated_set": label[key],
                    "contrast": "SF_vs_PB",
                    "n_genes": len(genes),
                    "set_size_in_ranked": in_ranked,
                    "gsea_min_size": min_size,
                    "testable": testable,
                    "untestable_reason": reason,
                    "nes": nes,
                    "pvalue": pval,
                    "padj": padj,
                    "leading_edge": le,
                    "evidence_tier": "secondary_annotation",
                })
    out = pd.DataFrame(rows)
    out.to_csv(tables_dir / "decomposition_nes.csv", index=False)
    return out


# ===========================================================================
# 3. The cGAS/STING axis — carried as a tally
# ===========================================================================
def sting_axis_overlap(tables_dir: Path, arms: dict[str, list[str]],
                       curated: dict[str, set[str]]) -> pd.DataFrame:
    """Tally the published IFN-independent STING signature against the mouse arms.

    Its overlap is a handful of genes, far under the size floor, so this stays a gene
    tally. The row records how many of the shared genes the hypoxia purge also removes,
    which is what makes the tally interpretable.
    """
    if not STING_SIG_PATH.exists():
        print(f"[11_heat_decomposition] STING signature absent at {STING_SIG_PATH} — "
              "tally skipped (the positive-control compartment is not checked out)")
        return pd.DataFrame()
    verify_source_hash(
        STING_SIG_PATH,
        "savi_sting_specific_up",
        tables_dir / "source_hash_manifest.csv",
        root=ROOT.parent,
    )
    sting = set(read_gene_list(STING_SIG_PATH))
    hypoxia = curated["hypoxia"]
    rows = []
    for arm in ARMS:
        shared = sorted(set(arms[arm]) & sting)
        also_hypoxia = sorted(g for g in shared if g in hypoxia)
        rows.append({
            "signature_a": f"{PRIMARY}_{arm}",
            "signature_b": STING_SIG_LABEL,
            "n_a": len(arms[arm]),
            "n_b": len(sting),
            "n_intersect": len(shared),
            "genes_intersect": ";".join(shared),
            "n_intersect_also_in_hypoxia": len(also_hypoxia),
            "genes_intersect_also_in_hypoxia": ";".join(also_hypoxia),
            "gsea_min_size": int(PARAMS.gsea_min_size),
            "testable_as_gsea_arm": len(shared) >= int(PARAMS.gsea_min_size),
            "evidence_tier": "secondary_annotation",
        })
    out = pd.DataFrame(rows)
    out.to_csv(tables_dir / "sting_axis_overlap.csv", index=False)
    return out


# ===========================================================================
def main() -> None:
    tables_dir = PATHS.tables(STAGE)
    arms = load_mouse_arms()
    curated = load_curated_sets()
    print("[11_heat_decomposition] mouse arms: "
          + ", ".join(f"{PRIMARY}_{a}={len(g)}" for a, g in arms.items()))
    print("[11_heat_decomposition] curated sets: "
          + ", ".join(f"{k}={len(v)}" for k, v in curated.items()))

    subs = build_subsignatures(arms, curated)
    sig_dir = prepare_signature_dirs(tables_dir, subs)
    overlap = overlap_tallies(tables_dir, arms, curated, subs)
    assignment = gene_assignment(tables_dir, arms, curated)
    multiplicity = assignment_multiplicity(tables_dir, assignment)
    nes = decomposition_nes(tables_dir, sig_dir, subs)
    sting = sting_axis_overlap(tables_dir, arms, curated)

    print("[11_heat_decomposition] how much of each arm each presumption claims:")
    print(overlap[["mouse_arm", "subcomponent", "n_intersect", "frac_of_mouse_arm"]]
          .to_string(index=False))
    print("[11_heat_decomposition] the assignment is NOT a partition:")
    print(multiplicity[["mouse_arm", "n_arm", "n_unassigned", "n_claimed",
                        "n_claimed_multiply", "n_claims_total", "is_partition"]]
          .to_string(index=False))
    if len(nes):
        up = nes[nes["mouse_arm"] == f"{PRIMARY}_up"]
        print("[11_heat_decomposition] up-arm NES:")
        print(up[["population", "subcomponent", "n_genes", "set_size_in_ranked",
                  "testable", "nes", "padj"]].to_string(index=False))
        n_untestable = int((~nes["testable"]).sum())
        print(f"[11_heat_decomposition] {n_untestable} of {len(nes)} population x "
              "subcomponent cells are untestable and carry their reason in the table")
    if len(sting):
        print("[11_heat_decomposition] STING axis tally:")
        print(sting[["signature_a", "n_intersect", "genes_intersect",
                     "testable_as_gsea_arm"]].to_string(index=False))
    print("[11_heat_decomposition] done — annotation tier only; "
          "effect_sizes_treg_arthritis.csv and 03_results/master/ untouched")


if __name__ == "__main__":
    main()
