#!/usr/bin/env python
"""
11_heat_decomposition_viz.py — VIZ ONLY. Which part of the mouse 39 C-derived arm
carries the synovial-fluid shift?
============================================================================
The mouse 39 °C up-arm enriches toward synovial fluid as a whole. This family
asks where that comes from by splitting the 199 up genes and 94 down genes with
curated, versioned, anchor-independent gene sets and reading each part off the
same ranked lists.

  1. heatdecomp_arm_coverage        — how much of each arm each presumption claims,
                                      including the arms too small to test
  2. heatdecomp_runsum_up_*         — the running sum of each testable up-arm part
  3. heatdecomp_runsum_down_*       — the same for the down arm's remainder

Every running-sum figure in the family shares ONE y-range, so a shallow curve in
one figure really is shallower than a steep curve in another. Reads only
committed stage-11 tables and computes no statistic: every NES, FDR and gene
tally is read verbatim from a CSV written by 11_heat_decomposition.py.

Run from the compartment root, AFTER 11_heat_decomposition.py:
  python 02_analysis/scripts/11_heat_decomposition_viz.py
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "02_analysis"))
os.chdir(ROOT)

from config import PARAMS, PATHS  # noqa: E402
from helpers.figure_style import (  # noqa: E402
    FIG_CFG,
    purge_figures,
    round_numeric_cols,
    save_overview,
    set_paper_style,
)

STAGE = "11_heat_decomposition"
SCRIPT = "02_analysis/scripts/11_heat_decomposition_viz.py"
PRIMARY = "WT_heat"
POP_TAG = {"Treg": "treg", "Tcon": "tcon", "CD8": "cd8"}
POP_ORDER = list(POP_TAG)

# --- declared palette constants ---------------------------------------------
# Population colours are resolved from `colors.okabe_ito` in the config and match
# 10_hsr_lens_viz.py, so a population keeps one colour across the compartment.
_OKABE = (FIG_CFG.get("colors", {}) or {}).get("okabe_ito", {}) or {}
POP_COL = {"Treg": _OKABE["bluish_green"], "Tcon": _OKABE["orange"],
           "CD8": _OKABE["reddish_purple"]}
# Mouse-arm diverging cue, IDENTICAL to 09_heat_hypoxia_viz.py: up = warm
# brown, down = cool blue. Keyed by arm so the mapping can never come out of
# a positional vector in the wrong order.
ARM_COL = {"up": "#A6611A", "down": "#2166AC"}

_F = FIG_CFG.get("figures", {}) or {}
ANNOT_SIZE = float(_F["axis_text_size"])
LEGEND_SIZE = float(_F["legend_text_size"])
RS_HEIGHTS = [float(h) for h in _F["running_sum_heights"]]
TOP_N = int(_F["top_n"])
FDR = float(PARAMS.gsea_fdr)
MIN_SIZE = int(PARAMS.gsea_min_size)

# Human-readable subcomponent labels, in the order the coverage figure reads them.
SUB_LABEL = {
    "unassigned": "no named program",
    "nfkb_tnfa": "TNFA / NF-kB signalling",
    "inflammatory": "inflammatory response",
    "hypoxia": "hypoxia",
    "t_activation": "IL2-STAT5 activation",
    "hsr_curated": "curated HSR core (Reactome/GO)",
    "ifn_type_i": "type-I interferon",
    "upr_er": "unfolded-protein response",
}

# The running-sum family: the (arm, subcomponent) pairs that clear the size floor in
# every population, in the order the narrative reads them.
DRAWN: list[tuple[str, str]] = [
    ("up", "unassigned"),
    ("up", "nfkb_tnfa"),
    ("up", "hypoxia"),
    ("up", "inflammatory"),
    ("up", "t_activation"),
    ("down", "unassigned"),
]

TITLES = {
    ("up", "unassigned"): "The unclaimed two thirds of the up-arm carries the shift",
    ("up", "nfkb_tnfa"): "The NF-kB part shifts to synovial fluid in CD4, weakly in CD8",
    ("up", "hypoxia"): "The hypoxia-overlap part shifts in every population",
    ("up", "inflammatory"): "The inflammatory-response part tracks the whole up-arm",
    ("up", "t_activation"): "The IL2-STAT5 activation part is the weakest in Treg",
    ("down", "unassigned"): "The down arm's remainder sits nowhere in particular",
}


def fmt_fdr(p: float) -> str:
    """Render an FDR for an in-figure label: fixed below 3 decimals, else scientific."""
    if pd.isna(p):
        return "FDR n/a"
    return f"FDR {p:.3f}" if p >= 0.001 else f"FDR {p:.0e}"


def _bool_col(s: pd.Series) -> pd.Series:
    """Robust to R (or pandas) writing a real bool or the strings TRUE/FALSE."""
    if s.dtype == bool:
        return s
    return s.astype(str).str.strip().str.upper().isin({"TRUE", "T", "1"})


# ===========================================================================
# 1. Coverage — how much of each arm does each presumption actually claim?
# ===========================================================================
def coverage_table() -> pd.DataFrame:
    """Join the plain overlap tallies to their testability, one row per arm x part."""
    tdir = PATHS.tables(STAGE)
    ov = pd.read_csv(tdir / "decomposition_overlap.csv")
    nes = pd.read_csv(tdir / "decomposition_nes.csv")
    nes["testable"] = _bool_col(nes["testable"])

    rows = []
    for _, r in ov.iterrows():
        cell = nes[(nes["mouse_arm"] == r["mouse_arm"])
                   & (nes["subcomponent"] == r["subcomponent"])]
        n_testable = int(cell["testable"].sum())
        sizes = cell["set_size_in_ranked"].astype(int)
        rows.append({
            "mouse_arm": r["mouse_arm"],
            "arm": str(r["mouse_arm"]).rsplit("_", 1)[-1],
            "subcomponent": r["subcomponent"],
            "curated_set": r["curated_set"],
            "n_curated_set": r["n_curated_set"],
            "n_mouse_arm": int(r["n_mouse_arm"]),
            "n_intersect": int(r["n_intersect"]),
            "frac_of_mouse_arm": float(r["frac_of_mouse_arm"]),
            "set_size_in_ranked_min": int(sizes.min()) if len(sizes) else 0,
            "set_size_in_ranked_max": int(sizes.max()) if len(sizes) else 0,
            "n_populations_testable": n_testable,
            "n_populations": int(len(cell)),
            "gsea_min_size": MIN_SIZE,
            "evidence_tier": "secondary_annotation",
        })
    out = pd.DataFrame(rows)
    # Up arm first, each arm ordered by how much of the signature it claims.
    out["_arm_rank"] = out["arm"].map({"up": 0, "down": 1})
    out = out.sort_values(["_arm_rank", "n_intersect"], ascending=[True, False])
    return out.drop(columns="_arm_rank").head(TOP_N).reset_index(drop=True)


def sting_sentence() -> str:
    """One caption sentence on the published STING axis, read from its committed tally."""
    path = PATHS.tables(STAGE) / "sting_axis_overlap.csv"
    if not path.exists():
        return ""
    st = pd.read_csv(path).set_index("signature_a")
    key = f"{PRIMARY}_up"
    if key not in st.index:
        return ""
    r = st.loc[key]
    genes = str(r["genes_intersect"]).replace(";", " and ")
    return (f"The published {int(r['n_b'])}-gene interferon-independent STING signature "
            f"contributes only {genes} here, tallied in sting_axis_overlap.csv.")


def coverage_note(r: pd.Series) -> str:
    """The right-hand testability line for one coverage bar. Counts, never a glyph."""
    if r["n_intersect"] == 0:
        return "0  ·  no gene of this arm"
    lo, hi = int(r["set_size_in_ranked_min"]), int(r["set_size_in_ranked_max"])
    ranked = f"{lo}" if lo == hi else f"{lo}–{hi}"
    if r["n_populations_testable"] == r["n_populations"] and r["n_populations"]:
        return f"{int(r['n_intersect'])}  ·  tested, {ranked} in the ranked lists"
    return (f"{int(r['n_intersect'])}  ·  under the {MIN_SIZE}-gene floor "
            f"({ranked} in the ranked lists)")


def multiplicity_rows() -> pd.DataFrame:
    """The committed per-arm multiplicity of the overlapping assignment.

    A hard failure rather than a skip: the bars in this figure overlap, and a
    reader who cannot see by how much will read them as a partition. That warning
    is not optional decoration, so a missing table stops the render.
    """
    path = PATHS.tables(STAGE) / "decomposition_assignment_multiplicity.csv"
    if not path.exists():
        raise FileNotFoundError(
            f"[11_viz] {path} not found. Re-run 11_heat_decomposition.py — the coverage "
            "figure may not be drawn without the counts that show its bars are not a "
            "partition.")
    return pd.read_csv(path).set_index("arm")


def plot_coverage(df: pd.DataFrame, mult: pd.DataFrame):
    fig, ax = plt.subplots()
    n = len(df)
    labels = []
    for i, (_, r) in enumerate(df.iterrows()):
        y = n - 1 - i
        ax.barh(y, r["n_intersect"], height=0.66, color=ARM_COL[r["arm"]], zorder=2)
        ax.text(r["n_intersect"] + 3.0, y, coverage_note(r), va="center", ha="left",
                fontsize=ANNOT_SIZE)
        labels.append(f"{SUB_LABEL.get(r['subcomponent'], r['subcomponent'])} · {r['arm']}")

    ax.set_yticks(range(n))
    ax.set_yticklabels(list(reversed(labels)))
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_xlim(0, 400)
    ax.set_xticks([0, 50, 100, 150])
    ax.set_xlabel("Genes of the mouse arm that the curated set contains")
    ax.set_title("Two thirds of the mouse 39 C-derived up arm belongs to no named program")
    handles = [
        Patch(facecolor=ARM_COL[a],
              label=f"{PRIMARY} {a} arm, {int(mult.loc[a, 'n_arm'])} genes")
        for a in ("up", "down")
    ]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.27),
              ncol=2, frameon=False, fontsize=LEGEND_SIZE)
    # THE BARS ARE NOT A PARTITION, stated as a count rather than as a caution.
    # Overlapping bars read as a partition unless the panel says otherwise, and
    # this is the panel most likely to be over-read, so the multiplicity goes on
    # the face and not only into the caption. Placed under the axis, where it
    # cannot land on a bar or on a right-hand testability note; `bbox_inches=
    # "tight"` keeps below-axis text in the render.
    u, d = mult.loc["up"], mult.loc["down"]
    ax.text(0.0, -0.115,
            "NOT A PARTITION — do not sum these bars.  "
            f"{int(u['n_claimed'])} of the {int(u['n_arm'])} up-arm genes are claimed by a curated "
            f"set, but they carry {int(u['n_claims_total'])} claims between them: "
            f"{int(u['n_claimed_multiply'])} of the {int(u['n_claimed'])}\n"
            f"belong to two or three sets at once, so adding the named bars double-counts "
            f"{int(u['n_excess_claims'])} claims and shrinks the {int(u['n_unassigned'])}-gene "
            f"remainder — the largest single part. On the down arm {int(d['n_claimed_multiply'])} "
            f"of {int(d['n_claimed'])} are multiply claimed.\n"
            "ANSWERS by membership over frozen versioned gene lists — arithmetic over committed "
            "files, not an effect estimate, and no NES is on this face.",
            transform=ax.transAxes, ha="left", va="top", fontsize=ANNOT_SIZE)
    fig.tight_layout()
    return fig


# ===========================================================================
# 2. Running sums — where each part of the arm sits in each ranked list
# ===========================================================================
def subcomponent_traces(arm: str, key: str) -> tuple[dict[str, pd.DataFrame], pd.DataFrame]:
    """Load one subcomponent's committed running-sum traces plus its NES summary rows."""
    tdir = PATHS.tables(STAGE)
    nes = pd.read_csv(tdir / "decomposition_nes.csv")
    nes = nes[(nes["mouse_arm"] == f"{PRIMARY}_{arm}") & (nes["subcomponent"] == key)]
    nes = nes.set_index("population")
    traces, rows = {}, []
    for pop, tag in POP_TAG.items():
        path = tdir / f"runsum_interactive_decomp_gsea_{tag}_{key}_{arm}.csv"
        if not path.exists() or pop not in nes.index:
            print(f"[11_heat_decomposition_viz] {pop} {key}_{arm}: no trace at {path} — skipping")
            continue
        tr = pd.read_csv(path, usecols=["rank", "running_es", "hit"])
        tr["hit"] = _bool_col(tr["hit"])
        traces[pop] = tr
        r = nes.loc[pop]
        rows.append({
            "population": pop,
            "mouse_arm": f"{PRIMARY}_{arm}",
            "subcomponent": key,
            "curated_set": r["curated_set"],
            "contrast": "SF_vs_PB",
            "n_genes": int(r["n_genes"]),
            "set_size_in_ranked": int(r["set_size_in_ranked"]),
            "nes": r["nes"],
            "pvalue": r["pvalue"],
            "padj": r["padj"],
            "n_ranked_genes": int(len(tr)),
            "evidence_tier": "secondary_annotation",
        })
    return traces, pd.DataFrame(rows)


def family_span(family: dict[tuple[str, str], dict[str, pd.DataFrame]]) -> float:
    """One running-ES half-range for the WHOLE family, so curves compare across figures."""
    peaks = [float(np.abs(tr["running_es"]).max())
             for traces in family.values() for tr in traces.values()]
    return max(peaks) * 1.25 if peaks else 1.0


def plot_subcomponent_runsum(arm: str, key: str, traces: dict[str, pd.DataFrame],
                             summary: pd.DataFrame, span: float, n_arm: int):
    # Two stacked panels in the config's running-sum proportions: the enrichment trace
    # over the gene-hit rug. The y-range is pinned to the FAMILY span rather than to
    # this figure's data, so curve heights are comparable figure to figure.
    fig, (ax, rug) = plt.subplots(
        2, 1, sharex=True, layout="constrained",
        height_ratios=[RS_HEIGHTS[0], RS_HEIGHTS[1]])

    pops = [p for p in POP_ORDER if p in traces]
    for pop in pops:
        tr = traces[pop]
        row = summary[summary["population"] == pop].iloc[0]
        ax.plot(tr["rank"], tr["running_es"], color=POP_COL[pop], lw=2.0,
                label=(f"{pop}   {int(row['set_size_in_ranked'])} genes ranked, "
                       f"NES {row['nes']:+.2f}, {fmt_fdr(row['padj'])}"))
    ax.axhline(0, color="black", lw=1)
    ax.set_ylim(-span, span)
    ax.set_ylabel("Running enrichment score")
    ax.set_title(TITLES[(arm, key)])
    ax.legend(frameon=False, fontsize=LEGEND_SIZE, loc="lower left")

    n_sub = int(summary["n_genes"].iloc[0]) if len(summary) else 0
    ax.text(0.98, 0.96,
            f"{SUB_LABEL.get(key, key)}: {n_sub} of the {n_arm} mouse {PRIMARY} {arm} genes\n"
            f"y-range shared across the whole decomposition family",
            transform=ax.transAxes, ha="right", va="top", fontsize=ANNOT_SIZE)

    for i, pop in enumerate(pops):
        hits = traces[pop].loc[traces[pop]["hit"], "rank"].to_numpy()
        y = len(pops) - 1 - i
        rug.vlines(hits, y, y + 0.86, color=POP_COL[pop], lw=1.2)
    rug.set_ylim(-0.1, len(pops))
    rug.set_yticks([len(pops) - 1 - i + 0.43 for i in range(len(pops))])
    rug.set_yticklabels(pops, fontsize=ANNOT_SIZE)
    rug.set_xlabel("Rank in the synovial-fluid-vs-blood ranked list")
    return fig


# ===========================================================================
def main() -> None:
    set_paper_style(config=FIG_CFG)
    purge_figures(STAGE, "heatdecomp_", overview=True, config=FIG_CFG)

    cov = coverage_table()
    mult = multiplicity_rows()
    fig = plot_coverage(cov, mult)
    save_overview(
        fig, STAGE, "heatdecomp_arm_coverage",
        table=round_numeric_cols(cov),
        finding=("Curated public gene sets claim only 62 of the 199 mouse 39 C-derived up "
                 "genes and 11 of the 94 down genes, so the largest part of the projected "
                 "signature — 137 up genes — belongs to no named program, and the curated "
                 "HSR core (Reactome/GO) contributes 2 genes. The bars are not a partition: "
                 "25 of those 62 claimed up genes belong to two or three curated sets at "
                 "once, so the 62 carry 92 claims and summing the named bars double-counts "
                 "30 of them."),
        script=SCRIPT, fn="plot_coverage",
        config_kv=(f"gsea_min_size={MIN_SIZE}; figures.top_n={TOP_N}; "
                   "evidence_tier=secondary_annotation"),
        input=("03_results/11_heat_decomposition/tables/decomposition_overlap.csv, "
               "03_results/11_heat_decomposition/tables/decomposition_nes.csv, "
               "03_results/11_heat_decomposition/tables/"
               "decomposition_assignment_multiplicity.csv, "
               "03_results/11_heat_decomposition/tables/sting_axis_overlap.csv"),
        how_to_read=(
            "ANSWERS what the projected set is made of, by membership over frozen versioned "
            "gene lists — arithmetic over committed files, not an effect estimate, and no NES "
            "on the face. One bar per mouse arm and curated presumption; length is how many "
            "of that arm's genes the curated set contains. Warm brown = the 199-gene up arm, "
            "cool blue = the 94-gene down arm. The right-hand text gives the count, then the "
            f"testability: parts reaching {MIN_SIZE} genes in the ranked lists are tested, "
            "smaller parts are marked under the floor, and a part with no gene in that arm "
            "says so. **Do not sum the bars.** The assignment is not a partition — 25 of the "
            "62 claimed up-arm genes sit in two or three sets, so adding the named bars "
            "double-counts 30 claims and shrinks the 137-gene remainder, the largest single "
            "part. That count is on the face, per arm in "
            "decomposition_assignment_multiplicity.csv, and per gene in "
            "decomposition_gene_assignment.csv. The remainder is reported as a remainder: it "
            "is not named, and is evidence for no mechanism. "
            + sting_sentence() + " Annotation tier."),
        config=FIG_CFG, wide=True, height=8.0,
    )
    plt.close(fig)

    # Load the whole running-sum family first so every figure shares one y-range.
    family, summaries = {}, {}
    for arm, key in DRAWN:
        traces, summary = subcomponent_traces(arm, key)
        if traces:
            family[(arm, key)] = traces
            summaries[(arm, key)] = summary
    span = family_span(family)
    n_arm = {"up": 199, "down": 94}

    for arm, key in DRAWN:
        if (arm, key) not in family:
            continue
        summary = summaries[(arm, key)]
        fig = plot_subcomponent_runsum(arm, key, family[(arm, key)], summary, span,
                                       n_arm[arm])
        best = summary.loc[summary["nes"].astype(float).abs().idxmax()]
        save_overview(
            fig, STAGE, f"heatdecomp_runsum_{arm}_{key}",
            table=round_numeric_cols(summary),
            finding=FINDINGS[(arm, key)],
            script=SCRIPT, fn="plot_subcomponent_runsum",
            config_kv=(f"figures.running_sum_heights={RS_HEIGHTS[:2]}; "
                       f"thresholds.gsea_fdr={FDR}; gsea_min_size={MIN_SIZE}; "
                       "evidence_tier=secondary_annotation"),
            input=("03_results/11_heat_decomposition/tables/runsum_interactive_decomp_gsea_"
                   f"{{treg,tcon,cd8}}_{key}_{arm}.csv, "
                   "03_results/11_heat_decomposition/tables/decomposition_nes.csv"),
            how_to_read=HOW_TO_READ[(arm, key)],
            config=FIG_CFG, height=7.0,
        )
        plt.close(fig)
        print(f"[11_heat_decomposition_viz] {key}_{arm}: strongest excursion in "
              f"{best['population']} at NES {float(best['nes']):+.2f}")

    print(f"[11_heat_decomposition_viz] wrote {1 + len(family)} overviews "
          "(arm coverage + running-sum family)")


# ---------------------------------------------------------------------------
# Caption text for the running-sum family, keyed by (arm, subcomponent) so a
# figure can never be shipped with another figure's finding.
# ---------------------------------------------------------------------------
FINDINGS = {
    ("up", "unassigned"): (
        "The 137 up-arm genes that no curated presumption claims give the strongest "
        "synovial-fluid enrichment of any part in CD8 (+2.10) and remain strongly "
        "enriched in Treg (+2.21) and Tcon (+2.27), so the shift is not carried by "
        "any single named program."),
    ("up", "nfkb_tnfa"): (
        "The 35 TNFA/NF-kB up-arm genes enrich toward synovial fluid strongly in Treg "
        "(+2.24) and Tcon (+2.32) and only weakly in CD8 (+1.23, FDR 0.22), making the "
        "inflammatory-signalling part the most CD4-selective of the decomposition."),
    ("up", "hypoxia"): (
        "The 18 hypoxia-overlap up-arm genes enrich toward synovial fluid in all three "
        "populations (+1.82 to +2.07), while removing these genes barely moves the "
        "whole-signature NES."),
    ("up", "inflammatory"): (
        "The 21 inflammatory-response up-arm genes track the whole up-arm (+1.52 to "
        "+2.11), adding no separation of their own beyond the broad synovial-fluid shift."),
    ("up", "t_activation"): (
        "The 14 IL2-STAT5 activation up-arm genes are the weakest testable part in Treg "
        "(+1.32, FDR 0.22) while reaching +1.89 in Tcon, so a curated T-cell activation "
        "program does not account for the Treg shift."),
    ("down", "unassigned"): (
        "The 83 down-arm genes no presumption claims sit nowhere in particular — NES "
        "+0.97 in Treg, +1.41 in Tcon and -1.12 in CD8, none of them significant — so the "
        "mouse down arm does not separate synovial fluid from blood in either direction."),
}

_HOW_BASE = (
    "Top panel: the weighted running enrichment score as each population's ranked list is "
    "walked from synovial-fluid-up (left) to blood-up (right); a positive, left-shifted "
    "excursion is synovial-fluid enrichment and a negative trace the opposite. Bottom "
    "panel: where this part's genes sit in each ranking, in matching colour. Legend labels "
    "carry the testable gene count, the NES and the FDR, and no other glyph marks "
    "significance. The y-range is shared across every figure of this decomposition family, "
    "so curve heights compare between figures. Annotation tier, firewalled from the "
    "confirmatory WT_heat effect-size spine.")

HOW_TO_READ = {
    ("up", "unassigned"): (
        "This part is the residual: the up-arm genes belonging to none of the curated "
        "presumptions. " + _HOW_BASE),
    ("up", "nfkb_tnfa"): (
        "This part is the up-arm genes that also sit in HALLMARK_TNFA_SIGNALING_VIA_NFKB. "
        + _HOW_BASE),
    ("up", "hypoxia"): (
        "This part is the up-arm genes that also sit in HALLMARK_HYPOXIA, the same 18 the "
        "whole-signature purge removes. " + _HOW_BASE),
    ("up", "inflammatory"): (
        "This part is the up-arm genes that also sit in HALLMARK_INFLAMMATORY_RESPONSE. "
        + _HOW_BASE),
    ("up", "t_activation"): (
        "This part is the up-arm genes that also sit in HALLMARK_IL2_STAT5_SIGNALING, read "
        "as a curated proxy for T-cell activation. " + _HOW_BASE),
    ("down", "unassigned"): (
        "This part is the residual of the mouse down arm: the genes belonging to none of "
        "the curated presumptions. " + _HOW_BASE),
}


if __name__ == "__main__":
    main()
