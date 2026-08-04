#!/usr/bin/env python
"""
14d_runsum_overview_viz.py — VIZ ONLY. One running sum per set, three populations.
==================================================================================
The per-database battery in 14c answers "within one population, which sets move?".
This script asks the transposed question: for ONE set, where along the
synovial-fluid-versus-paired-blood ranking do its genes sit in Treg, in Tcon and in
CD8 — on one axis, so the comparison is made by the eye rather than by holding three
figures in mind. That transposition is the whole point: the constraint this
compartment already reports for the mouse-derived arm is that its synovial
enrichment is NOT cell-subset-selective, and a figure that puts the three
populations on separate pages cannot show that either way.

WHY THE X AXIS IS A FRACTION AND NOT A RANK
-------------------------------------------
The three populations have different ranked lists — different lengths, because
filterByExpr keeps a different gene set in each — so a raw rank axis would end at a
different x for each curve and the shortest list would read as a truncated curve
rather than a shorter one. x is therefore rank / n_ranked, and each curve's own
n_ranked travels in the legend and in the same-stem table.

WHICH SETS GET A PANEL, AND WHY THE RULE MATTERS
------------------------------------------------
Only sets whose running-sum trace exists in ALL THREE populations. A trace is written
by 14_unbiased_enrichment.R for every mouse-derived arm, every set named in
`unbiased_enrichment.runsum_always`, and the top-N curated per population — so a set
that ranked top-N in Treg alone HAS a Treg curve and no other. Drawing that as a
one-curve panel in this family would read as "the set does not enrich in Tcon or
CD8", when it means "no trace was written there". The partial-coverage sets are
listed in the run summary and in the family caption, named, with that reason.

Reads committed stage-14 tables and computes no statistic. Every NES and adjusted p
is read from gsea_all.csv; the curve is read from the trace table the compute stage
wrote. The only arithmetic here is rank / n_ranked and the peak of a column.

Input (03_results/14_unbiased_enrichment/tables/):
  runsum_interactive_index.csv                 which trace tables exist, and why
  runsum_interactive_<pop>_<set>.csv           the traces themselves
  gsea_all.csv                                 NES, pooled FDR, set and leading-edge size

Output (03_results/14_unbiased_enrichment/):
  figures/_overview/runsum_<set>.{pdf,png}     one per complete-coverage set
  tables/_overview/runsum_<set>.csv            its three rows
  README.md                                    one caption each, plus a family caption

Run from the compartment root, AFTER 14_unbiased_enrichment.R:
  python 02_analysis/scripts/14d_runsum_overview_viz.py
"""
from __future__ import annotations

import os
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "02_analysis"))
os.chdir(ROOT)

from config import PARAMS, PATHS, POPULATION_COLORS  # noqa: E402
from helpers.figure_style import (  # noqa: E402
    FIG_CFG,
    purge_figures,
    round_numeric_cols,
    save_overview,
    set_paper_style,
    write_caption,
)

STAGE = "14_unbiased_enrichment"
SCRIPT = "02_analysis/scripts/14d_runsum_overview_viz.py"

POP_TAG = {"Treg": "treg", "Tcon": "tcon", "CD8": "cd8"}
POP_ORDER = list(POP_TAG)
POP_COL = POPULATION_COLORS

_F = FIG_CFG.get("figures", {}) or {}
ANNOT_SIZE = float(_F["axis_text_size"])
LEGEND_SIZE = float(_F["legend_text_size"])
RS_HEIGHTS = [float(h) for h in _F["running_sum_heights"]]
FDR = float(PARAMS.gsea_fdr)

# Sets this compartment owns, named by HOW THEY WERE DERIVED rather than by the
# mechanism they are hoped to represent, so a panel title cannot smuggle the
# conclusion. Anything absent here keeps its identifier verbatim.
DERIVATION = {
    "WT_heat_up": "mouse 39 °C-derived up arm, wild type",
    "KO_heat_up": "mouse 39 °C-derived up arm, cGAS knockout",
    "Interaction_up": "mouse 39 °C x genotype interaction up arm",
    "HSR_core": "curated heat-shock-response core, independent of the mouse anchor",
    "sting_specific_up": "SAVI-derived, STING-attributable",
    "ifn_only_up": "SAVI-derived, generic type-I interferon",
}

# MSigDB collection prefixes, stripped for display only; the full identifier stays in
# the table, the caption and the file stem.
_PREFIX = re.compile(r"^(HALLMARK|REACTOME|KEGG|KEGG_LEGACY|WP|GOBP|GOCC|GOMF|MITOPATHWAYS)_")


def fmt_fdr(p: float) -> str:
    if pd.isna(p):
        return "FDR n/a"
    return f"FDR {p:.3f}" if p >= 0.001 else f"FDR {p:.0e}"


def display_name(set_id: str) -> str:
    """A legible panel title for a set id, without ever truncating it."""
    if set_id in DERIVATION:
        return set_id
    stripped = _PREFIX.sub("", set_id)
    if stripped == set_id:
        return set_id
    return stripped.replace("_", " ").capitalize()


def _bool_col(s: pd.Series) -> pd.Series:
    """Robust to R writing a real bool or the strings TRUE/FALSE."""
    if s.dtype == bool:
        return s
    return s.astype(str).str.strip().str.upper().isin({"TRUE", "T", "1"})


# ===========================================================================
# Load: the trace index, the traces, and the published statistics
# ===========================================================================
def load_family() -> tuple[dict[str, dict[str, pd.DataFrame]], pd.DataFrame, pd.DataFrame]:
    """Return {set_id: {population: trace}}, the summary rows, and the skipped sets.

    A set is admitted only with a trace in every population, for the reason in the
    module docstring. The skipped frame is returned rather than discarded so the
    family caption can name what is missing and why.
    """
    tdir = PATHS.tables(STAGE)
    index_path = tdir / "runsum_interactive_index.csv"
    if not index_path.exists():
        raise FileNotFoundError(
            f"[14d] {index_path} not found. Run 02_analysis/scripts/14_unbiased_enrichment.R first."
        )
    idx = pd.read_csv(index_path)
    sweep = pd.read_csv(tdir / "gsea_all.csv")

    coverage = idx.groupby("gene_set")["population"].apply(lambda s: set(s))
    complete = sorted(k for k, v in coverage.items() if set(POP_ORDER) <= v)
    partial = sorted(k for k, v in coverage.items() if not set(POP_ORDER) <= v)

    skipped = (
        idx[idx["gene_set"].isin(partial)]
        .groupby(["gene_set", "database"], as_index=False)
        .agg(populations_with_a_trace=("population", lambda s: " ".join(sorted(s))))
        .assign(reason="no trace in every population: not named in runsum_always and "
                       "outside the top-N curated quota there")
        .sort_values("gene_set")
    )

    traces: dict[str, dict[str, pd.DataFrame]] = {}
    rows = []
    for set_id in complete:
        per_pop = {}
        for pop in POP_ORDER:
            rec = idx[(idx["gene_set"] == set_id) & (idx["population"] == pop)].iloc[0]
            path = tdir / str(rec["file"])
            if not path.exists():
                raise FileNotFoundError(
                    f"[14d] {set_id} / {pop}: the index lists {path.name} but it is not on disk. "
                    "The index and the trace tables are from different runs."
                )
            tr = pd.read_csv(path, usecols=["rank", "stat", "running_es", "hit", "leading_edge"])
            tr["hit"] = _bool_col(tr["hit"])
            tr["leading_edge"] = _bool_col(tr["leading_edge"])
            # The one piece of arithmetic in this script: a fraction of the list, so
            # three rankings of different length share one axis.
            tr["rank_fraction"] = tr["rank"] / len(tr)
            per_pop[pop] = tr

            s = sweep[(sweep["population"] == pop)
                      & (sweep["pathway_id"] == set_id)
                      & (sweep["database"] == rec["database"])]
            if s.empty:
                raise ValueError(
                    f"[14d] {set_id} / {pop}: no row in gsea_all.csv under database "
                    f"{rec['database']}. Every number on a panel must be readable from the "
                    "published table."
                )
            s = s.iloc[0]
            peak_i = int(tr["running_es"].abs().idxmax())
            rows.append({
                "gene_set": set_id,
                "display_name": display_name(set_id),
                "database": rec["database"],
                "population": pop,
                "contrast": rec["contrast"],
                "nes": float(s["nes"]),
                "pvalue": float(s["pvalue"]),
                "padj_in_database": float(s["padj"]),
                "padj_pooled": float(s["padj_pooled"]),
                "direction": str(s["direction"]),
                "set_size": int(s["set_size"]),
                "leading_edge_size": int(s["leading_edge_size"]),
                "n_ranked": int(len(tr)),
                "peak_running_es": float(tr["running_es"].iloc[peak_i]),
                "peak_rank": int(tr["rank"].iloc[peak_i]),
                "peak_rank_fraction": float(tr["rank_fraction"].iloc[peak_i]),
                "emitted_because_named": bool(_bool_col(pd.Series([rec["always_emitted"]])).iloc[0]),
            })
        traces[set_id] = per_pop
    return traces, pd.DataFrame(rows), skipped


# ===========================================================================
# The panel
# ===========================================================================
def plot_set(set_id: str, per_pop: dict[str, pd.DataFrame], summary: pd.DataFrame):
    """Three populations' running sums for one set, on a shared fractional-rank axis.

    Three stacked panels in the config's running-sum proportions: the enrichment
    traces, the gene-hit rug with one NAMED row per population, and the three ranked
    metrics. The metric panel is kept even though it repeats between panels of this
    family, because it is what shows the three rankings crossing zero at comparable
    fractions — the assumption the shared x axis rests on.
    """
    fig, (ax, rug, met) = plt.subplots(
        3, 1, sharex=True, layout="constrained", height_ratios=RS_HEIGHTS)

    sub = summary[summary["gene_set"] == set_id].set_index("population")

    for pop in POP_ORDER:
        tr, r = per_pop[pop], sub.loc[pop]
        ax.plot(tr["rank_fraction"], tr["running_es"], color=POP_COL[pop], lw=2.0,
                label=(f"{pop}   NES {r['nes']:+.2f}, pooled {fmt_fdr(r['padj_pooled'])}, "
                       f"{int(r['set_size'])} genes ranked of "
                       f"{int(r['n_ranked']):,}, {int(r['leading_edge_size'])} in the leading edge"))
    ax.axhline(0, color="black", lw=1)
    ax.set_ylabel("Running enrichment score")
    note = DERIVATION.get(set_id)
    ax.set_title(
        f"{display_name(set_id)} along the synovial-fluid-versus-paired-blood ranking"
        + (f"\n{note}" if note else f"\n{set_id}"))

    # The y range is DATA-DRIVEN, and asymmetric on purpose. The config's
    # running_sum_ylim pins one fixed range across a family of separate per-population
    # figures so their shapes stay comparable; here all three populations already share
    # one axis, so comparability is intrinsic and a fixed range would spend most of the
    # panel on white space. Zero is always included — a running-sum panel that cropped
    # it would hide the sign — and the padding is split unevenly so the legend gets the
    # room it needs on the side the curves leave empty, decided from the sign of the
    # largest excursion rather than by eye. Without that, a three-line legend either
    # overlaps the traces or forces the symmetric range back.
    hi = max(0.0, max(float(per_pop[p]["running_es"].max()) for p in POP_ORDER))
    lo = min(0.0, min(float(per_pop[p]["running_es"].min()) for p in POP_ORDER))
    rng = max(hi - lo, 1e-9)
    peak_positive = sub.loc[sub["peak_running_es"].abs().idxmax(), "peak_running_es"] > 0
    if peak_positive:
        ax.set_ylim(lo - 0.45 * rng, hi + 0.12 * rng)
    else:
        ax.set_ylim(lo - 0.12 * rng, hi + 0.45 * rng)
    ax.legend(frameon=False, fontsize=LEGEND_SIZE,
              loc="lower left" if peak_positive else "upper left")

    for i, pop in enumerate(POP_ORDER):
        tr = per_pop[pop]
        hits = tr.loc[tr["hit"], "rank_fraction"].to_numpy()
        y = len(POP_ORDER) - 1 - i
        rug.vlines(hits, y, y + 0.86, color=POP_COL[pop], lw=1.0)
    rug.set_ylim(-0.1, len(POP_ORDER))
    rug.set_yticks([len(POP_ORDER) - 1 - i + 0.43 for i in range(len(POP_ORDER))])
    rug.set_yticklabels(POP_ORDER, fontsize=ANNOT_SIZE)

    for pop in POP_ORDER:
        tr = per_pop[pop]
        met.plot(tr["rank_fraction"], tr["stat"], color=POP_COL[pop], lw=1.0, alpha=0.75)
    met.axhline(0, color="black", lw=0.8)
    met.set_ylabel("Moderated t")
    met.set_xlabel("Position in the ranked list, as a fraction of its length "
                   "(0 = most synovial-fluid-up, 1 = most blood-up)")
    met.set_xlim(0, 1)
    return fig


# ===========================================================================
def main() -> None:
    set_paper_style(config=FIG_CFG)
    purge_figures(STAGE, "runsum_", overview=True, config=FIG_CFG)

    traces, summary, skipped = load_family()
    if not traces:
        print("[14d] no set has a trace in all three populations — nothing drawn")
        return

    CONFIG_KV = (f"figures.running_sum_heights={RS_HEIGHTS}; thresholds.gsea_fdr={FDR}; "
                 f"x_axis=rank/n_ranked; y_range=data-driven symmetric; "
                 f"unbiased_enrichment.runsum_always="
                 f"{list((FIG_CFG.get('unbiased_enrichment', {}) or {}).get('runsum_always', []))}")

    for set_id, per_pop in traces.items():
        sub = summary[summary["gene_set"] == set_id]
        stem = "runsum_" + re.sub(r"[^A-Za-z0-9._-]", "_", set_id)
        fig = plot_set(set_id, per_pop, summary)

        by_pop = {r["population"]: r for _, r in sub.iterrows()}
        parts = ", ".join(
            f"{p} NES {by_pop[p]['nes']:+.2f} at pooled {fmt_fdr(by_pop[p]['padj_pooled'])}"
            for p in POP_ORDER)
        signs = {np.sign(by_pop[p]["nes"]) for p in POP_ORDER}
        sig_pops = [p for p in POP_ORDER if by_pop[p]["padj_pooled"] < FDR]
        shape = ("all three populations carry the same sign" if len(signs) == 1
                 else "the sign is not shared across the three populations")
        finding = (
            f"On one fractional-rank axis, {display_name(set_id)} gives {parts}, so {shape}"
            + (f" and {'it' if len(sig_pops) == 1 else 'they'} reach pooled FDR < {FDR:g} in "
               f"{', '.join(sig_pops)}." if sig_pops
               else f" and no population reaches pooled FDR < {FDR:g}.")
            + " Where a curve peaks is where the set's genes concentrate in that population's "
              "ranking; a difference in peak height between populations is a difference in "
              "concentration, not a measured difference in the program the set is named for.")

        save_overview(
            fig, STAGE, stem, table=round_numeric_cols(sub),
            finding=finding, script=SCRIPT, fn="plot_set",
            config_kv=CONFIG_KV,
            input=("03_results/14_unbiased_enrichment/tables/runsum_interactive_index.csv, "
                   "runsum_interactive_{treg,tcon,cd8}_" + set_id + ".csv, gsea_all.csv"),
            how_to_read=(
                "One set, three populations, one axis. TOP: the running enrichment score as "
                "each ranked list is walked from its most synovial-fluid-up gene (left) to "
                "its most blood-up gene (right), so a positive left-shifted excursion is "
                "concentration on the synovial-fluid side and a negative trace the opposite. "
                "MIDDLE: where that population's members of the set sit in its own ranking, "
                "one named row each. BOTTOM: the three moderated-t rankings, which show them "
                "crossing zero at comparable fractions and so justify the shared axis. "
                "X IS A FRACTION of each list's length, not a rank, because the lists differ "
                "in length. Y IS DATA-DRIVEN, not the fixed range the by_contrast panels use, "
                "so do not compare curve heights across the two families. The legend carries "
                "each NES, pooled adjusted p, genes reaching the ranked list and leading-edge "
                "count, so a tall curve resting on few genes shows as one. Correlative: this "
                "says where gene content sits in a ranking, not that the program the set is "
                "named for is present. Claim tier: L3; no row reaches an effect-size "
                "accumulator."),
            config=FIG_CFG, width=11.0, height=8.5,
        )
        plt.close(fig)

    # ---- family caption: what the family is, and what it deliberately omits ----
    skipped_txt = (
        "Every set with a trace is drawn." if skipped.empty else
        "Not drawn, because a trace exists in only some populations and a one-curve panel in "
        "this family would read as absence of enrichment rather than absence of a trace: "
        + "; ".join(f"{r['gene_set']} [{r['database']}] ({r['populations_with_a_trace']})"
                    for _, r in skipped.iterrows())
        + ". Their statistics are all in gsea_all.csv, and any of them can be added to the "
          "family by naming it in unbiased_enrichment.runsum_always.")
    n_named = int(summary.groupby("gene_set")["emitted_because_named"].max().sum())
    write_caption(
        stage=STAGE,
        filename="figures/_overview/runsum_&lt;set&gt;.png (cross-population running-sum family)",
        finding=(
            f"The transposed view of the sweep: {summary['gene_set'].nunique()} sets, each drawn "
            f"once with all three sorted populations on one fractional-rank axis, so whether a "
            f"set's enrichment is shared across Treg, Tcon and CD8 or specific to one of them is "
            f"read off a single panel. {n_named} of them are in the family because they are named "
            f"in unbiased_enrichment.runsum_always — the comparators this compartment's niche "
            f"question is read against, including the mouse-derived arms, HALLMARK_HYPOXIA, the "
            f"curated HSR core and the two SAVI-derived axes — and the rest because they ranked "
            f"into the top-N curated quota in every population. {skipped_txt}"),
        script=SCRIPT, fn="load_family",
        config_kv=CONFIG_KV,
        input=("03_results/14_unbiased_enrichment/tables/{runsum_interactive_index,gsea_all}.csv "
               "+ runsum_interactive_&lt;population&gt;_&lt;set&gt;.csv"),
        how_to_read=(
            "LAYOUT. One figure per set, figures/_overview/runsum_&lt;set&gt;.{pdf,png}, with its "
            "three rows under tables/_overview/. This family is the transpose of "
            "figures/by_contrast/, where a panel holds one population and many sets. "
            "MEMBERSHIP IS A REPORTING CHOICE, NOT A RESULT: a set is here because it was named "
            "in the config or ranked top-N in all three populations, so the family is not a "
            "ranking and a set's presence in it privileges nothing. READ IT FOR SHARED VERSUS "
            "SELECTIVE. Three curves of similar shape and height mean the set's genes sit in a "
            "comparable part of all three rankings; one curve standing apart means it does not. "
            "Neither reading is evidence about the program a set is named for, and neither "
            "reaches the confirmatory spine, which stays the donor-pseudobulk effect sizes. "
            "Correlative. Claim tier: L3 (enrichment statistics)."),
        config=FIG_CFG)

    print(f"[14d] wrote {summary['gene_set'].nunique()} cross-population running sums "
          f"({n_named} named, {summary['gene_set'].nunique() - n_named} top-N), "
          f"{len(skipped)} set(s) skipped for partial coverage")
    if not skipped.empty:
        for _, r in skipped.iterrows():
            print(f"    skipped {r['gene_set']:<45s} trace only in {r['populations_with_a_trace']}")


if __name__ == "__main__":
    main()
