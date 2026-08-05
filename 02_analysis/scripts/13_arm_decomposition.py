#!/usr/bin/env python
"""
13_arm_decomposition.py — COMPUTE ONLY. What are the mouse-derived UP arms MADE OF?
==================================================================================
The `11_heat_decomposition` results establish the curated-set composition of one arm,
`WT_heat_up`. This script asks the same MEMBERSHIP question of all three mouse-derived up
arms in human projection space, against a wider lens panel, and emits the result in a tidy
(arm, program, gene) shape an alluvial or stacked-bar display can draw without recomputing
anything.

WHAT THIS MEASURES, AND WHAT IT DOES NOT
----------------------------------------
Membership is not enrichment. A curated lens CONTAINING a gene of an arm is arithmetic over
two committed text files; a lens ENRICHING in a ranked list is a separate measurement with a
different failure mode. This script performs only the first. Nothing here carries an NES, a
p-value, an effect size or a direction of shift, and no row reaches
`effect_sizes_treg_arthritis.csv` or any `03_results/master/` accumulator. The enrichment
question for `WT_heat_up` is answered by the `11_heat_decomposition` results and by the
donor-pseudobulk panels; it is deliberately not re-answered here.

FOUR THINGS THE OUTPUT IS NOT
-----------------------------
1. NOT a partition of genes. The lenses overlap, so a gene can be claimed by several
   programs. Every such gene appears once per claiming program in `arm_program_gene.csv`,
   with `n_programs_for_gene` recording the multiplicity and `weight_fractional =
   1 / n_programs_for_gene` so a consumer can choose duplicated accounting (sum the rows)
   or fractional accounting (sum the weights, which totals the arm exactly). Forcing a
   priority-ordered disjoint assignment would silently decide which program gets credit for
   a shared gene, so the sharing is published instead.
2. NOT four independent arms. The three mouse contrasts are linearly dependent by
   construction — WT_heat = KO_heat + Interaction — so agreement between `WT_heat_up` and
   `KO_heat_up` is expected arithmetic, not corroboration. `Interaction_up_fdrOnly` is the
   same Interaction contrast read at a relaxed gate, so it is not independent of
   `Interaction_up` either; both are carried, and the `gate` column marks which is which.
3. NOT a statement about the down arms. Only up arms are in scope. No down arm is read,
   tallied or mentioned in the outputs.
4. NOT a re-curation of anything. Every lens is read verbatim from a frozen committed file.

Sizes are contracts, not observations: each arm's gene count is checked against the mouse
anchor's own `manifest.csv`, and a drift is a hard stop rather than a quietly different
denominator.

Inputs (all read-only):
  - ../mouse_anchor/03_results/human_projection/manifest.csv          (arm sizes + gates)
  - ../mouse_anchor/03_results/human_projection/signatures/WT_heat/WT_heat_up.txt
  - ../mouse_anchor/03_results/human_projection/signatures/KO_heat/KO_heat_up.txt
  - ../mouse_anchor/03_results/human_projection/signatures/Interaction/Interaction_up.txt
  - ../mouse_anchor/03_results/human_projection/signatures/Interaction/Interaction_fdrOnly_up.txt
  - 00_data/references/msigdb_hallmark/HALLMARK_*.txt                 (6 frozen Hallmark sets)
  - 00_data/references/temp_hsr_lens/HSR_core.txt                     (frozen curated HSR core)
  - ../sting_positive_control/03_results/06_reference_axis/signatures/sting_specific_up.txt
  - ../sting_positive_control/03_results/06_reference_axis/signatures/ifn_only_up.txt

Outputs (all under 03_results/13_arm_decomposition/tables/):
  - arm_program_gene.csv          one row per (arm, program, gene) — the alluvial substrate
  - arm_program_summary.csv       one row per (arm, program) — the tallies
  - arm_program_multiplicity.csv  one row per (arm, gene) — the double-counting, inspectable

Run from the compartment root, before 13_arm_decomposition_viz.py:
  python 02_analysis/scripts/13_arm_decomposition.py
"""
from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Iterable

import pandas as pd

COMPARTMENT_ROOT = Path(__file__).resolve().parents[2]
REPO_ROOT = COMPARTMENT_ROOT.parent

sys.path.insert(0, str(COMPARTMENT_ROOT))
sys.path.insert(0, str(COMPARTMENT_ROOT / "02_analysis"))
os.chdir(COMPARTMENT_ROOT)

from config import PATHS  # noqa: E402
from helpers.source_hash_manifest import verify_source_hash  # noqa: E402

STAGE = "13_arm_decomposition"
SCRIPT = "02_analysis/scripts/13_arm_decomposition.py"
TIER = "secondary_annotation"
RESIDUAL = "unassigned"
RESIDUAL_SET = "(none)"

# --- the mouse-derived UP arms, in human projection space --------------------------------
# `arm` is the label carried in every output; `contrast` and `gate` are the mouse anchor's own
# manifest keys, so the size guard reads the contract rather than a number retyped here.
# The `Interaction_up_fdrOnly` label names the arm as "the Interaction up arm at the relaxed
# gate"; its frozen file is `Interaction_fdrOnly_up.txt`.
PROJECTION_DIR = REPO_ROOT / "mouse_anchor/03_results/human_projection"
PROJECTION_MANIFEST = PROJECTION_DIR / "manifest.csv"

ARMS: list[tuple[str, str, str, Path]] = [
    ("WT_heat_up", "WT_heat", "fdr_logfc",
     PROJECTION_DIR / "signatures/WT_heat/WT_heat_up.txt"),
    ("KO_heat_up", "KO_heat", "fdr_logfc",
     PROJECTION_DIR / "signatures/KO_heat/KO_heat_up.txt"),
    ("Interaction_up", "Interaction", "fdr_logfc",
     PROJECTION_DIR / "signatures/Interaction/Interaction_up.txt"),
    ("Interaction_up_fdrOnly", "Interaction", "fdr_only",
     PROJECTION_DIR / "signatures/Interaction/Interaction_fdrOnly_up.txt"),
]
ARM_ORDER = [a[0] for a in ARMS]

# --- the lens panel: nine frozen, versioned, anchor-independent gene sets ----------------
# Read verbatim. The first seven keys are the ones the 11_heat_decomposition results already
# use, kept byte-identical so the two sets of tables join on `program`. The last two are the
# STING positive-control compartment's two reference axes.
HALLMARK_DIR = COMPARTMENT_ROOT / "00_data/references/msigdb_hallmark"
HSR_DIR = COMPARTMENT_ROOT / "00_data/references/temp_hsr_lens"
AXIS_DIR = REPO_ROOT / "sting_positive_control/03_results/06_reference_axis/signatures"

PROGRAMS: list[tuple[str, str, Path]] = [
    ("hypoxia", "HALLMARK_HYPOXIA", HALLMARK_DIR / "HALLMARK_HYPOXIA.txt"),
    ("nfkb_tnfa", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
     HALLMARK_DIR / "HALLMARK_TNFA_SIGNALING_VIA_NFKB.txt"),
    ("inflammatory", "HALLMARK_INFLAMMATORY_RESPONSE",
     HALLMARK_DIR / "HALLMARK_INFLAMMATORY_RESPONSE.txt"),
    ("t_activation", "HALLMARK_IL2_STAT5_SIGNALING",
     HALLMARK_DIR / "HALLMARK_IL2_STAT5_SIGNALING.txt"),
    ("ifn_type_i", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
     HALLMARK_DIR / "HALLMARK_INTERFERON_ALPHA_RESPONSE.txt"),
    ("upr_er", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
     HALLMARK_DIR / "HALLMARK_UNFOLDED_PROTEIN_RESPONSE.txt"),
    ("hsr_curated", "HSR_core", HSR_DIR / "HSR_core.txt"),
    ("sting_specific_published", "de_Cevins_sting_specific_up",
     AXIS_DIR / "sting_specific_up.txt"),
    ("ifn_generic_axis", "ifn_only_up", AXIS_DIR / "ifn_only_up.txt"),
]
PROGRAM_ORDER = [p[0] for p in PROGRAMS]

# The already-committed pin for the published STING gene set, owned by the
# 11_heat_decomposition tables. Re-using that pin rather than minting a fresh one means a
# change to the positive-control axis stops BOTH readers, instead of letting a new stage
# quietly re-pin whatever it happens to find on disk.
STING_PIN_MANIFEST = PATHS.tables("11_heat_decomposition") / "source_hash_manifest.csv"
STING_PIN_LABEL = "savi_sting_specific_up"


def read_gene_list(path: Path) -> list[str]:
    """Newline-delimited HGNC symbols, order preserved, blank lines dropped."""
    if not path.exists():
        raise FileNotFoundError(f"gene list absent: {path}")
    return [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]


def _rel(path: Path) -> str:
    """Repo-root-relative path string, for provenance columns."""
    try:
        return str(Path(path).resolve().relative_to(REPO_ROOT.resolve()))
    except ValueError:
        return str(path)


# ===========================================================================
# 1. Load the arms against the mouse anchor's own size contract
# ===========================================================================
def expected_arm_sizes() -> dict[tuple[str, str, str], int]:
    """`{(contrast, direction, gate): n_human}` read from the projection manifest.

    The arm sizes are a frozen contract of the mouse anchor, so they are read from the
    anchor's manifest rather than retyped here. A retyped number silently becomes a
    different denominator the day the contract changes.
    """
    if not PROJECTION_MANIFEST.exists():
        raise FileNotFoundError(
            f"projection manifest absent: {PROJECTION_MANIFEST}. The mouse anchor "
            "sub-project is not checked out, so the arm sizes cannot be verified.")
    man = pd.read_csv(PROJECTION_MANIFEST)
    return {
        (str(r["contrast"]), str(r["direction"]), str(r["gate"])): int(r["n_human"])
        for _, r in man.iterrows()
    }


def load_arms() -> dict[str, list[str]]:
    """Return `{arm: genes}`, hard-stopping on any size drift from the frozen contract."""
    want = expected_arm_sizes()
    arms: dict[str, list[str]] = {}
    for arm, contrast, gate, path in ARMS:
        genes = read_gene_list(path)
        key = (contrast, "up", gate)
        if key not in want:
            raise ValueError(
                f"{arm}: no manifest row for contrast={contrast}, direction=up, gate={gate} "
                f"in {PROJECTION_MANIFEST}")
        if len(genes) != want[key]:
            raise ValueError(
                f"{arm} size drift: the projection manifest declares {want[key]} human "
                f"genes at gate {gate}, the file holds {len(genes)} ({path}). Reconcile the "
                "mouse anchor contract before decomposing.")
        if len(set(genes)) != len(genes):
            raise ValueError(f"{arm}: duplicate symbols in {path}; the arm is not a set.")
        arms[arm] = genes
        print(f"[{STAGE}] arm {arm}: {len(genes)} genes, gate {gate}, from {_rel(path)}")
    return arms


def load_programs() -> dict[str, list[str]]:
    """Return `{program: genes}` for the nine frozen lenses, verbatim from disk."""
    out: dict[str, list[str]] = {}
    for key, label, path in PROGRAMS:
        genes = read_gene_list(path)
        out[key] = genes
        print(f"[{STAGE}] lens {key} ({label}): {len(genes)} genes from {_rel(path)}")
    return out


def gate_of(arm: str) -> str:
    return next(g for a, _, g, _ in ARMS if a == arm)


def curated_set_of(program: str) -> str:
    if program == RESIDUAL:
        return RESIDUAL_SET
    return next(lab for k, lab, _ in PROGRAMS if k == program)


# ===========================================================================
# 2. Membership — which lenses contain each gene of each arm
# ===========================================================================
def claims(arms: dict[str, list[str]],
           programs: dict[str, list[str]]) -> dict[str, dict[str, list[str]]]:
    """`{arm: {gene: [programs containing it]}}`, in the declared program order."""
    sets = {k: set(v) for k, v in programs.items()}
    return {
        arm: {gene: [p for p in PROGRAM_ORDER if gene in sets[p]] for gene in genes}
        for arm, genes in arms.items()
    }


def arm_program_gene(tables_dir: Path, arms: dict[str, list[str]],
                     claimed: dict[str, dict[str, list[str]]]) -> pd.DataFrame:
    """The alluvial substrate: one row per (arm, program, gene).

    A gene claimed by k lenses yields k rows, each carrying `n_programs_for_gene = k` and
    `weight_fractional = 1/k`. A gene no lens claims yields ONE row on the `unassigned`
    program with `n_programs_for_gene = 0` and `weight_fractional = 1`. So summing
    `weight_fractional` over an arm returns the arm size exactly, whichever way a
    downstream display chooses to account for sharing.
    """
    rows = []
    for arm in ARM_ORDER:
        gate = gate_of(arm)
        for gene in arms[arm]:
            hits = claimed[arm][gene]
            if hits:
                for program in hits:
                    rows.append({
                        "arm": arm,
                        "gate": gate,
                        "program": program,
                        "curated_set": curated_set_of(program),
                        "gene": gene,
                        "n_programs_for_gene": len(hits),
                        "weight_fractional": 1.0 / len(hits),
                    })
            else:
                rows.append({
                    "arm": arm,
                    "gate": gate,
                    "program": RESIDUAL,
                    "curated_set": RESIDUAL_SET,
                    "gene": gene,
                    "n_programs_for_gene": 0,
                    "weight_fractional": 1.0,
                })
    out = pd.DataFrame(rows)
    out["_arm"] = out["arm"].map(ARM_ORDER.index)
    out["_prog"] = out["program"].map(lambda p: len(PROGRAM_ORDER) if p == RESIDUAL
                                     else PROGRAM_ORDER.index(p))
    out = (out.sort_values(["_arm", "_prog", "gene"])
              .drop(columns=["_arm", "_prog"]).reset_index(drop=True))
    out.to_csv(tables_dir / "arm_program_gene.csv", index=False)
    return out


def arm_program_summary(tables_dir: Path, arms: dict[str, list[str]],
                        programs: dict[str, list[str]],
                        claimed: dict[str, dict[str, list[str]]]) -> pd.DataFrame:
    """One row per (arm, program), plus the `unassigned` remainder row per arm.

    `n_intersect` is a plain containment count and `frac_of_arm` its share of the arm. Rows
    within an arm overlap, so `n_intersect` does not sum to `n_arm` — that is what
    `arm_program_multiplicity.csv` and the fractional weights are for.
    """
    rows = []
    for arm in ARM_ORDER:
        gate, n_arm = gate_of(arm), len(arms[arm])
        for program in PROGRAM_ORDER:
            genes = sorted(g for g in arms[arm] if program in claimed[arm][g])
            rows.append({
                "arm": arm,
                "gate": gate,
                "n_arm": n_arm,
                "program": program,
                "curated_set": curated_set_of(program),
                "n_curated_set": len(programs[program]),
                "n_intersect": len(genes),
                "frac_of_arm": len(genes) / n_arm,
                "genes": ";".join(genes),
            })
        residual = sorted(g for g in arms[arm] if not claimed[arm][g])
        rows.append({
            "arm": arm,
            "gate": gate,
            "n_arm": n_arm,
            "program": RESIDUAL,
            "curated_set": RESIDUAL_SET,
            "n_curated_set": 0,
            "n_intersect": len(residual),
            "frac_of_arm": len(residual) / n_arm,
            "genes": ";".join(residual),
        })
    out = pd.DataFrame(rows)
    out.to_csv(tables_dir / "arm_program_summary.csv", index=False)
    return out


def arm_program_multiplicity(tables_dir: Path, arms: dict[str, list[str]],
                             claimed: dict[str, dict[str, list[str]]]) -> pd.DataFrame:
    """One row per (arm, gene): how many lenses claim it, and which.

    This is the audit trail that makes the double-counting inspectable. Genes with
    `n_programs = 0` read `programs = unassigned`, so the remainder is visible here too
    rather than only as an absence.
    """
    rows = []
    for arm in ARM_ORDER:
        for gene in sorted(arms[arm]):
            hits = claimed[arm][gene]
            rows.append({
                "arm": arm,
                "gene": gene,
                "n_programs": len(hits),
                "programs": ";".join(hits) if hits else RESIDUAL,
            })
    out = pd.DataFrame(rows)
    out.to_csv(tables_dir / "arm_program_multiplicity.csv", index=False)
    return out


# ===========================================================================
# 3. Provenance — the one source that is pinned
# ===========================================================================
def verify_pinned_sources() -> None:
    """Hard-check the one input this stage pins, against a hash already committed.

    Only the published STING lens is pinned, and it is checked against the hash
    `11_heat_decomposition` carries, so a changed source stops the run. The other
    inputs are read as-is. Recording their hashes beside this one produced a table
    of provenance that nothing enforced, which reads as a guarantee it does not
    give; the enforced check is kept and the record is not.
    """
    for key, _, path in PROGRAMS:
        if key == "sting_specific_published":
            verify_source_hash(path, STING_PIN_LABEL, STING_PIN_MANIFEST, root=REPO_ROOT)
            print(f"[{STAGE}] {key}: SHA-256 matches the committed pin in "
                  f"{_rel(STING_PIN_MANIFEST)}")


# ===========================================================================
# 4. Self-checks — the invariants that make the tables readable
# ===========================================================================
def check_invariants(gene: pd.DataFrame, summary: pd.DataFrame,
                     mult: pd.DataFrame, arms: dict[str, list[str]]) -> None:
    """Assert the arithmetic the three tables promise each other."""
    for arm in ARM_ORDER:
        n_arm = len(arms[arm])
        g = gene[gene["arm"] == arm]
        # Fractional weights partition the arm exactly; duplicated rows do not.
        total = float(g["weight_fractional"].sum())
        if abs(total - n_arm) > 1e-9:
            raise AssertionError(f"{arm}: fractional weights sum to {total}, not {n_arm}")
        # Every gene of the arm is represented, exactly once per claiming program.
        if set(g["gene"]) != set(arms[arm]):
            raise AssertionError(f"{arm}: arm_program_gene.csv does not cover the arm")
        if len(mult[mult["arm"] == arm]) != n_arm:
            raise AssertionError(f"{arm}: multiplicity table is not one row per gene")
        # The summary's per-program counts are the gene table's per-program row counts.
        want = g.groupby("program").size().to_dict()
        got = (summary[(summary["arm"] == arm) & (summary["n_intersect"] > 0)]
               .set_index("program")["n_intersect"].to_dict())
        if want != got:
            raise AssertionError(f"{arm}: summary counts disagree with the gene table")
        # The remainder is exactly the genes no lens claims.
        n_res = int(summary[(summary["arm"] == arm)
                            & (summary["program"] == RESIDUAL)]["n_intersect"].iloc[0])
        n_zero = int((mult[mult["arm"] == arm]["n_programs"] == 0).sum())
        if n_res != n_zero:
            raise AssertionError(f"{arm}: remainder {n_res} != unclaimed genes {n_zero}")
    print(f"[{STAGE}] invariants hold: fractional weights total each arm exactly, "
          "duplicated counts and the remainder agree across all three tables")


def check_against_11(summary: pd.DataFrame, mult: pd.DataFrame) -> None:
    """Reproduce the `11_heat_decomposition` WT_heat_up tallies from this stage's tables.

    The seven lenses that stage uses are a SUBSET of the nine here, so restricting to them
    must return its published numbers gene for gene. This reads the committed
    `decomposition_overlap.csv` rather than retyping the counts, so the check cannot drift
    into agreement with a remembered number.

    The one quantity that legitimately differs is the remainder. `unassigned` here means
    "claimed by none of NINE lenses", which is a smaller set than "claimed by none of
    SEVEN". That gap is printed with the genes responsible, because a remainder quoted
    against the wrong panel is the easiest number in these tables to misread.
    """
    seven = set(PROGRAM_ORDER) - {"sting_specific_published", "ifn_generic_axis"}
    ref_path = PATHS.tables("11_heat_decomposition") / "decomposition_overlap.csv"
    wt = mult[mult["arm"] == "WT_heat_up"].copy()
    wt["_n7"] = wt["programs"].map(
        lambda s: 0 if s == RESIDUAL else len(seven & set(str(s).split(";"))))

    if not ref_path.exists():
        print(f"[{STAGE}] {_rel(ref_path)} absent — the seven-lens cross-check is skipped, "
              "so this run's WT_heat_up tallies are unconfirmed against the published ones")
        return
    ref = pd.read_csv(ref_path)
    ref = ref[ref["mouse_arm"] == "WT_heat_up"].set_index("subcomponent")["n_intersect"]

    bad = []
    for program in sorted(seven):
        here = int(summary[(summary["arm"] == "WT_heat_up")
                           & (summary["program"] == program)]["n_intersect"].iloc[0])
        there = int(ref.loc[program])
        if here != there:
            bad.append(f"{program}: this stage {here}, 11_heat_decomposition {there}")
    res_here = int((wt["_n7"] == 0).sum())
    res_there = int(ref.loc[RESIDUAL])
    if res_here != res_there:
        bad.append(f"{RESIDUAL} over the seven shared lenses: this stage {res_here}, "
                   f"11_heat_decomposition {res_there}")
    if bad:
        raise AssertionError(
            "WT_heat_up membership disagrees with the published seven-lens tallies:\n  "
            + "\n  ".join(bad)
            + "\nOne of the two stages is reading a different arm or a different frozen "
              "lens. Diagnose the input files before touching either script.")

    res_nine = int(summary[(summary["arm"] == "WT_heat_up")
                           & (summary["program"] == RESIDUAL)]["n_intersect"].iloc[0])
    print(f"[{STAGE}] WT_heat_up reproduces every published seven-lens tally, including a "
          f"{res_there}-gene remainder over those seven")
    if res_nine != res_there:
        newly = sorted(wt[(wt["_n7"] == 0) & (wt["programs"] != RESIDUAL)]["gene"])
        print(f"[{STAGE}] over all NINE lenses that remainder is {res_nine}, because the two "
              f"reference axes claim {len(newly)} gene(s) the seven leave unclaimed: "
              + ", ".join(newly))


def report_arm_overlap(arms: dict[str, list[str]]) -> None:
    """Print how much gene content the four arms actually share.

    The mouse anchor's contract states the three contrasts are linearly dependent —
    WT_heat = KO_heat + Interaction. That is a statement about the MODEL COEFFICIENTS, and
    it does not carry to the thresholded gene lists: a gene can pass one contrast's gate and
    fail another's. So the dependence has to be read off the actual overlaps rather than
    assumed from the algebra, and printing them keeps a reader from doing either — treating
    the arms as independent, or treating them as literal set sums.
    """
    print(f"[{STAGE}] gene content shared between arms (the arms are neither independent "
          "nor literal set sums of one another):")
    for i, a in enumerate(ARM_ORDER):
        for b in ARM_ORDER[i + 1:]:
            shared = set(arms[a]) & set(arms[b])
            print(f"  {a} & {b}: {len(shared)} shared "
                  f"({len(shared)}/{len(arms[a])} of {a}, {len(shared)}/{len(arms[b])} of {b})")


def report(summary: pd.DataFrame, mult: pd.DataFrame) -> None:
    """Print the membership tallies and the not-a-partition counts."""
    print(f"[{STAGE}] containment of each arm by each lens (NOT a partition — "
          "these counts overlap and must not be summed):")
    wide = (summary.pivot(index="program", columns="arm", values="n_intersect")
            .reindex(PROGRAM_ORDER + [RESIDUAL])[ARM_ORDER])
    print(wide.to_string())
    print(f"[{STAGE}] how far the assignment is from a partition:")
    for arm in ARM_ORDER:
        m = mult[mult["arm"] == arm]["n_programs"].astype(int)
        print(f"  {arm}: n_arm={len(m)}  unclaimed={int((m == 0).sum())}  "
              f"claimed={int((m > 0).sum())}  claimed_by_2plus={int((m >= 2).sum())}  "
              f"max_lenses_per_gene={int(m.max())}  total_claims={int(m.sum())}  "
              f"excess_claims={int(m.sum() - (m > 0).sum())}")


def main() -> None:
    tables_dir = PATHS.tables(STAGE)
    arms = load_arms()
    programs = load_programs()
    claimed = claims(arms, programs)

    gene = arm_program_gene(tables_dir, arms, claimed)
    summary = arm_program_summary(tables_dir, arms, programs, claimed)
    mult = arm_program_multiplicity(tables_dir, arms, claimed)
    verify_pinned_sources()

    check_invariants(gene, summary, mult, arms)
    check_against_11(summary, mult)
    report_arm_overlap(arms)
    report(summary, mult)
    print(f"[{STAGE}] done — membership only; no NES, no p-value, no effect size, "
          "and effect_sizes_treg_arthritis.csv plus 03_results/master/ untouched")


if __name__ == "__main__":
    main()
