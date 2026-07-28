#!/usr/bin/env python
"""
09_heat_hypoxia_viz.py — VIZ ONLY. The heat-versus-hypoxia narrative in four beats.
==================================================================================
The inflamed synovial niche is hypoxic as well as inflamed, so one bounded
question can be asked of the mouse 39 °C-derived enrichment: is it reducible to
the set's own HALLMARK_HYPOXIA-overlap gene content? That is a membership
question and it is answered by deleting those genes from the mouse sets and
re-running the same fgsea engine. It is not a question about temperature, and it
is not a question about whether hypoxia is a confound or a co-exposure — those
are not separable in cross-sectional human data, and nothing here licenses a
statement about either. These figures walk the bounded answer, then put the
mouse signature back into the plain Treg differential-expression view so its
size and its specificity are both visible:

  1. heat_purge_nes_paired         — how much of the enrichment survives the purge
  2. heat_hypoxia_colocalization   — do heat-high and hypoxia-high mark the same cells
  3. heat_treg_volcano_signature   — where the mouse arms sit in the Treg SF-vs-PB volcano
  4. heat_treg_volcano_programs    — which programs those genes are, and how little
                                     they share with the published STING axis

The model-assigned leading-edge composition table remains an exploratory compute
resource, including for the reactive review notebook, but its visualization is
withdrawn from the published overview. Whole-arm composition is carried only by
11_heat_decomposition/heatdecomp_arm_coverage.

Computes no statistic. Every NES, FDR, correlation, logFC and p is read verbatim
from a committed CSV; the only derived quantities are set membership over frozen
gene lists, plain counts and differences of those memberships, and the usual
-log10(FDR) volcano display transform.

Run from the compartment root, AFTER 09_heat_hypoxia.py:
  python 02_analysis/scripts/09_heat_hypoxia_viz.py
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
from adjustText import adjust_text  # noqa: E402
from matplotlib.colors import to_rgba  # noqa: E402
from matplotlib.legend_handler import HandlerTuple  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
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
from helpers.source_hash_manifest import verify_source_hash  # noqa: E402

STAGE = "09_heat_hypoxia"
SCRIPT = "02_analysis/scripts/09_heat_hypoxia_viz.py"
POP_TAG = {"Treg": "treg", "Tcon": "tcon", "CD8": "cd8"}
ARM_ORDER = ["up", "down"]
ARM_SET = {"up": "WT_heat_up", "down": "WT_heat_down"}

# --- committed inputs, declared once ---------------------------------------
# Same module-constant idiom 09_heat_hypoxia.py uses for its reference inputs.
DE_TREG_PATH = Path("03_results/03_pseudobulk/tables/de_SFvsPB_treg.csv")
TAXONOMY_PATH = Path("00_data/references/heat_leadingedge_taxonomy/leadingedge_gene_taxonomy.csv")
# Frozen reference axes published by the STING positive-control compartment: the
# 21-gene IFN-independent STING-activation signature of de Cevins et al. 2023
# (Cell Rep Med, PMID 38118407) and the 200-gene generic type-I IFN program from
# an IFN-beta donor pseudobulk. Read-only, scored nowhere here — membership only.
AXES_DIR = Path("../sting_positive_control/03_results/06_reference_axis/signatures")
AXIS_LABEL = {"sting_specific": "STING-specific axis (21 published genes)",
              "ifn_only": "generic type-I IFN axis (200 genes)"}

# Declared palette constants -------------------------------------------------
# Mouse-arm diverging cue, IDENTICAL to 05_score_signatures_viz.R: heat-up =
# warm brown, heat-down = cool blue. Keyed by arm so the mapping can never come
# out of a positional vector in the wrong order.
ARM_COL = {"up": "#A6611A", "down": "#2166AC"}

# Leading-edge program colours come from the config Okabe-Ito palette, so the
# only literal here is the config key, not a hex string.
_OKABE = (FIG_CFG.get("colors", {}) or {}).get("okabe_ito", {}) or {}
LE_PROGRAMS = [
    ("heat_shock_proteostasis", "heat shock / proteostasis", _OKABE["vermillion"]),
    ("hypoxia_HIF", "hypoxia / HIF overlap", _OKABE["sky_blue"]),
    ("immediate_early_stress", "immediate-early stress", _OKABE["orange"]),
    ("effector_activation", "effector / activation", _OKABE["bluish_green"]),
    ("other", "other", _OKABE["reddish_purple"]),
]

# Reference-axis cues use the two Okabe-Ito hues LE_PROGRAMS leaves free, so a
# program colour can never be mistaken for an axis colour inside one figure.
AXIS_COL = {"sting_specific": _OKABE["black"], "ifn_only": _OKABE["blue"]}
# Named greys (not hex) for the genes that carry no membership.
BG_COL = "lightgrey"
EDGE_HALO = "white"
EDGE_CONTRAST = _OKABE["black"]

_F = FIG_CFG.get("figures", {}) or {}
ANNOT_SIZE = float(_F["axis_text_size"])
LEGEND_SIZE = float(_F["legend_text_size"])
LABEL_TOP = int(_F["volcano_label_top"])
FDR = float(PARAMS.gsea_fdr)
DE_FDR = float(PARAMS.de_fdr)
DE_LFC = float(PARAMS.de_logfc)
# One translucency for every paired/overlapping marker in this module, so
# coincident markers read through each other instead of hiding one another.
MARKER_ALPHA = 0.65


def fmt_fdr(p: float) -> str:
    """Render an FDR for an in-figure label: fixed below 3 decimals, else scientific."""
    if pd.isna(p):
        return "FDR n/a"
    return f"FDR {p:.3f}" if p >= 0.001 else f"FDR {p:.0e}"


# ===========================================================================
# 1. Paired full-versus-purged NES — does the enrichment survive the purge?
# ===========================================================================
def _count_symbols(path: Path) -> int:
    """Plain line count of a frozen newline-delimited symbol file. No statistic."""
    if not path.exists():
        raise FileNotFoundError(f"[09_viz] frozen signature file missing: {path}")
    return len({ln.strip() for ln in path.read_text().splitlines() if ln.strip()})


def purge_paired_table() -> pd.DataFrame:
    """Marshal the per-population full/purged fgsea rows into one row per arm.

    `gene_purge_nes_comparison.csv` carries only the up arm, so the down arm is read
    from the same `gsea_{full,purged}_*.csv` files that table was built from. No
    number is recomputed — NES, p, FDR and set_size are all read as written.

    Two membership counts travel with every row so the caption can quote the
    nominal-versus-effective distinction from its own same-stem table rather than
    from a neighbour:

      `n_removed_nominal`  — genes the purge deletes from the frozen set FILE
                             (a line-count difference over two committed files);
      `n_removed_testable` — how many of those were in that population's ranked
                             list at all, i.e. the number the purge could actually
                             move (`set_size_full - set_size_purged`).

    They differ, and quoting the nominal count as "genes removed" overstates the
    purge. `delta_nes` is likewise carried so the caption never has to subtract
    two columns in prose.
    """
    tdir = PATHS.tables(STAGE)
    n_nominal = {arm: _count_symbols(tdir / "_signatures_full" / f"{ARM_SET[arm]}.txt")
                 for arm in ARM_ORDER}
    nominal = {
        arm: n_nominal[arm]
        - _count_symbols(tdir / "_signatures_purged" / f"{ARM_SET[arm]}.txt")
        for arm in ARM_ORDER
    }
    rows = []
    for pop, tag in POP_TAG.items():
        full = pd.read_csv(tdir / f"gsea_full_{tag}.csv").set_index("pathway_id")
        purged = pd.read_csv(tdir / f"gsea_purged_{tag}.csv").set_index("pathway_id")
        for arm in ARM_ORDER:
            sid = ARM_SET[arm]
            if sid not in full.index or sid not in purged.index:
                continue
            size_full = int(full.loc[sid, "set_size"])
            size_purged = int(purged.loc[sid, "set_size"])
            rows.append({
                "population": pop,
                "arm": arm,
                "signature": sid,
                "contrast": "SF_vs_PB",
                "nes_full": float(full.loc[sid, "nes"]),
                "padj_full": float(full.loc[sid, "padj"]),
                "set_size_full": size_full,
                "nes_purged": float(purged.loc[sid, "nes"]),
                "padj_purged": float(purged.loc[sid, "padj"]),
                "set_size_purged": size_purged,
                "delta_nes": float(purged.loc[sid, "nes"]) - float(full.loc[sid, "nes"]),
                "n_nominal": n_nominal[arm],
                "n_removed_nominal": nominal[arm],
                "n_removed_testable": size_full - size_purged,
                "evidence_tier": "primary_pseudobulk",
            })
    return pd.DataFrame(rows)


def _pair_handles(colour: str) -> tuple:
    """Legend handle for one arm: the full-set diamond beside the purged circle.

    Both are drawn at the plotted translucency and with no contrasting edge, so
    the key shows the default (FDR-not-passing) marker state and the outline key
    below it is the only thing that carries significance.
    """
    face = to_rgba(colour, MARKER_ALPHA)
    return (Line2D([0], [0], linestyle="none", marker="D", markerfacecolor=face,
                   markeredgecolor=face, markersize=15),
            Line2D([0], [0], linestyle="none", marker="o", markerfacecolor=face,
                   markeredgecolor=face, markersize=10))


def plot_purge_paired(df: pd.DataFrame):
    fig, ax = plt.subplots()
    n = len(df)
    ylabels = []
    for i, (_, r) in enumerate(df.iterrows()):
        y = n - 1 - i
        col = ARM_COL[r["arm"]]
        ax.plot([r["nes_full"], r["nes_purged"]], [y, y], color=col, lw=2.6, zorder=2,
                solid_capstyle="round")
        # Every marker is filled and translucent, larger full-set diamond first and
        # smaller purged circle on top: where the purge removes nothing and both land
        # on the same NES (the whole down arm) the pair reads as a darker circle
        # inside a lighter diamond rather than one marker hiding the other.
        for nes, padj, marker, size in ((r["nes_full"], r["padj_full"], "D", 300),
                                        (r["nes_purged"], r["padj_purged"], "o", 120)):
            passes = padj < FDR
            ax.scatter(nes, y, marker=marker, s=size, zorder=3,
                       facecolors=to_rgba(col, MARKER_ALPHA),
                       edgecolors=EDGE_CONTRAST if passes else to_rgba(col, MARKER_ALPHA),
                       linewidths=2.2 if passes else 0.0)
        # Effective set size against this arm's OWN nominal size, then what the
        # purge cost. Where the purge takes nothing there is no NES to quote, so
        # the row says that instead of printing a rounded zero.
        if int(r["n_removed_testable"]) == 0:
            cost = "purge removes no gene"
        else:
            cost = f"ΔNES {r['delta_nes']:+.3f}"
        ax.text(2.95, y,
                f"n {r['set_size_full']}→{r['set_size_purged']} of {r['n_nominal']}  ·  "
                f"{cost}  ·  {fmt_fdr(r['padj_purged'])}",
                va="center", ha="left", fontsize=ANNOT_SIZE)
        ylabels.append(f"{r['population']} · {r['arm']}")

    ax.axvline(0, color="black", lw=1)
    ax.set_yticks(range(n))
    ax.set_yticklabels(list(reversed(ylabels)))
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_xlim(0, 4.6)
    ax.set_xticks([0, 1, 2, 3])
    ax.set_xlabel("fgsea NES (synovial fluid vs paired blood)")
    # The title states the bounded answer in membership terms. It may not say
    # anything about temperature, and it may not call hypoxia a confound or a
    # co-exposure: the purge tests gene content and licenses nothing else.
    ax.set_title("The synovial-fluid enrichment is not reducible to the set's\n"
                 "HALLMARK_HYPOXIA-overlap gene content")
    # Q -> A marker, placed under the axis so it cannot collide with a marker:
    # this panel answers, and it says on its face what it answers, at what tier,
    # and where the answer stops. `bbox_inches="tight"` keeps below-axis text.
    up = df[df["arm"].eq("up")]
    ax.text(0.0, -0.135,
            "ANSWER — confirmatory tier (donor-level pseudobulk, limma-voom → fgsea).\n"
            f"Deleting the {int(up['n_removed_nominal'].iloc[0])} HALLMARK_HYPOXIA-overlap genes "
            f"takes {int(up['n_removed_testable'].min())}–{int(up['n_removed_testable'].max())} "
            "testable genes out of the up arm and costs "
            f"{abs(float(up['delta_nes'].max())):.3f}–{abs(float(up['delta_nes'].min())):.3f} NES;\n"
            "all three up arms stay significant, so the enrichment is not reducible to that gene "
            "content. That is the whole claim. It says nothing\n"
            "about temperature, and nothing about whether hypoxia is a confound or a co-exposure — "
            "those are not separable in cross-sectional human data.",
            transform=ax.transAxes, ha="left", va="top", fontsize=ANNOT_SIZE)
    # Three keys: each arm key carries its colour AND its two shapes, and the one
    # remaining key explains the outline. No open/filled convention anywhere.
    handles = [_pair_handles(ARM_COL["up"]), _pair_handles(ARM_COL["down"]),
               Line2D([0], [0], linestyle="none", marker="o",
                      markerfacecolor=to_rgba("grey", MARKER_ALPHA),
                      markeredgecolor=EDGE_CONTRAST, markeredgewidth=2.2, markersize=11)]
    # The down-arm key names where that arm is significant, read from the table
    # rather than asserted, so the panel cannot outlive "the up arm is the only
    # informative arm" by carrying it in a legend.
    down_sig = df[df["arm"].eq("down") & df["padj_full"].lt(FDR)]["population"].tolist()
    down_where = (", ".join(down_sig) + " only") if down_sig else "no population"
    labels = ["WT_heat up arm — diamond = full mouse set, circle = hypoxia-purged set",
              "WT_heat down arm — same pair, and the purge removes no gene, so they coincide; "
              f"this arm is not silent, reaching FDR below {FDR} in {down_where}",
              f"dark outline = FDR below {FDR} (every FDR is printed at right)"]
    ax.legend(handles, labels, handler_map={tuple: HandlerTuple(ndivide=None, pad=0.7)},
              loc="upper center", bbox_to_anchor=(0.5, -0.30), ncol=1, frameon=False,
              fontsize=LEGEND_SIZE, handlelength=3.2, handletextpad=1.0, labelspacing=0.7)
    fig.tight_layout()
    return fig


# ===========================================================================
# 2. Per-cell co-localization — are these the same cells?
# ===========================================================================
def colocalization_table() -> pd.DataFrame:
    df = pd.read_csv(PATHS.tables(STAGE) / "heat_hypoxia_colocalization.csv")
    cell = df[df["level"] == "cell"].copy()
    cell["population"] = pd.Categorical(cell["population"], list(POP_TAG), ordered=True)
    return cell.sort_values(["population", "method"]).reset_index(drop=True)


def plot_colocalization(df: pd.DataFrame):
    fig, ax = plt.subplots()
    methods = ["spearman", "pearson"]
    m_col = {"spearman": _OKABE["blue"], "pearson": _OKABE["sky_blue"]}
    width = 0.34
    pops = list(POP_TAG)
    for i, method in enumerate(methods):
        sub = df[df["method"] == method].set_index("population").reindex(pops)
        xpos = [x + (i - 0.5) * width for x in range(len(pops))]
        ax.bar(xpos, sub["r"].astype(float), width=width, color=m_col[method],
               label=method.capitalize())
        for xp, r in zip(xpos, sub["r"].astype(float)):
            ax.text(xp, r + 0.025, f"{r:.2f}", ha="center", va="bottom", fontsize=ANNOT_SIZE)

    ns = df[df["method"] == "spearman"].set_index("population").reindex(pops)["n"]
    ax.set_xticks(range(len(pops)))
    ax.set_xticklabels([f"{p}\nn = {int(v):,} cells" for p, v in zip(pops, ns)])
    ax.axhline(0, color="black", lw=1)
    ax.set_ylim(-0.05, 1.0)
    ax.set_ylabel("Correlation of WT_heat_up and HALLMARK_HYPOXIA score, within SF cells")
    ax.set_title("WT_heat_up-high and HALLMARK_HYPOXIA-high are largely different cells")
    ax.legend(frameon=False, fontsize=LEGEND_SIZE, loc="upper right")
    ax.text(0.02, 0.94,
            "CORROBORATES ONLY — per-cell scores never answer. Full –0.05 to 1\n"
            "scale, so the height of a bar is the whole story: the two scores\n"
            "co-localize weakly, i.e. they mark largely different cells.\n"
            "Donor-level SF means (n = 6–7) are unpowered and reported in the\n"
            "stage table only.",
            transform=ax.transAxes, ha="left", va="top", fontsize=ANNOT_SIZE)
    fig.tight_layout()
    return fig


# ===========================================================================
# 3/4. The mouse signature inside the plain Treg SF-vs-PB volcano
# ===========================================================================
def _read_symbols(path: Path) -> set:
    """Read a newline-delimited, header-free gene-symbol list."""
    return {ln.strip() for ln in Path(path).read_text().splitlines() if ln.strip()}


def volcano_annotated_de() -> pd.DataFrame:
    """Annotate the committed Treg SF-vs-PB pseudobulk DE table by set membership.

    `log2FoldChange`, `pvalue` and `padj` are read verbatim from the pseudobulk DE
    table and `-log10(padj)` is the usual volcano display transform. Everything
    else added here is plain membership over frozen gene lists — the two mouse
    arms, the 18 purged `HALLMARK_HYPOXIA` genes, the leading-edge taxonomy, and
    the two published reference axes. No model is refitted, nothing is re-ranked.
    """
    de = pd.read_csv(DE_TREG_PATH).dropna(subset=["padj", "log2FoldChange"]).copy()
    sig_dir = PATHS.tables(STAGE) / "_signatures_full"
    arms = {arm: _read_symbols(sig_dir / f"{ARM_SET[arm]}.txt") for arm in ARM_ORDER}
    cmp_tab = pd.read_csv(PATHS.tables(STAGE) / "gene_purge_nes_comparison.csv")
    purged = set(str(cmp_tab.loc[cmp_tab["population"] == "Treg",
                                 "genes_removed"].iloc[0]).split(";"))
    tax = pd.read_csv(TAXONOMY_PATH).set_index("gene")["category"]
    axes = {}
    for k in AXIS_LABEL:
        source_path = AXES_DIR / f"{k}_up.txt"
        verify_source_hash(
            source_path,
            f"savi_{k}_up",
            PATHS.tables(STAGE) / "source_hash_manifest.csv",
            root=ROOT.parent,
        )
        axes[k] = _read_symbols(source_path)

    sym = de["gene_symbol"].astype(str)
    de["arm"] = np.where(sym.isin(arms["up"]), "up",
                         np.where(sym.isin(arms["down"]), "down", "not_in_signature"))
    de["hypoxia_purged"] = sym.isin(purged) & de["arm"].eq("up")
    de["le_program"] = sym.map(tax).fillna("not_annotated")
    de["reference_axis"] = np.where(sym.isin(axes["sting_specific"]), "sting_specific",
                                    np.where(sym.isin(axes["ifn_only"]), "ifn_only", "none"))
    de["neg_log10_padj"] = -np.log10(de["padj"] + 1e-300)
    de["sf_direction"] = np.where(de["log2FoldChange"] > 0, "SF_high", "blood_high")
    de["passes_de_gates"] = de["padj"].lt(DE_FDR) & de["log2FoldChange"].abs().ge(DE_LFC)
    de.attrs["sets"] = {"arms": arms, "axes": axes, "taxonomy": set(tax.index),
                        "purged": purged}
    return de


def volcano_tallies(de: pd.DataFrame) -> dict:
    """Plain counts over the annotated table — set sizes and gate tallies only."""
    sets = de.attrs["sets"]
    gated = de[de["passes_de_gates"]]
    t = {"n_tested": len(de),
         "n_gates_sf": int(gated["sf_direction"].eq("SF_high").sum()),
         "n_gates_blood": int(gated["sf_direction"].eq("blood_high").sum())}
    for arm in ARM_ORDER:
        sub, g = de[de["arm"].eq(arm)], gated[gated["arm"].eq(arm)]
        t[f"n_{arm}_set"] = len(sets["arms"][arm])
        t[f"n_{arm}_tested"] = len(sub)
        t[f"n_{arm}_gates_sf"] = int(g["sf_direction"].eq("SF_high").sum())
        t[f"n_{arm}_gates_blood"] = int(g["sf_direction"].eq("blood_high").sum())
    t["pct_sf_response_from_up"] = 100.0 * t["n_up_gates_sf"] / max(t["n_gates_sf"], 1)
    t["n_taxonomy"] = len(sets["taxonomy"])
    t["n_taxonomy_in_up"] = len(sets["taxonomy"] & sets["arms"]["up"])
    t["n_taxonomy_in_down"] = len(sets["taxonomy"] & sets["arms"]["down"])
    t["n_taxonomy_tested"] = int(de["le_program"].ne("not_annotated").sum())
    for key, members in sets["axes"].items():
        sub, g = de[de["reference_axis"].eq(key)], gated[gated["reference_axis"].eq(key)]
        t[f"n_{key}_set"] = len(members)
        t[f"n_{key}_tested"] = len(sub)
        t[f"n_{key}_gates_sf"] = int(g["sf_direction"].eq("SF_high").sum())
        t[f"n_{key}_in_up"] = len(members & sets["arms"]["up"])
        t[f"n_{key}_in_down"] = len(members & sets["arms"]["down"])
        t[f"genes_{key}_in_up"] = sorted(members & sets["arms"]["up"])
    return t


def _label_points(ax, rows: pd.DataFrame, avoid: pd.DataFrame | None = None):
    """Draw gene-symbol labels for `rows` and return the symbols drawn (never truncated).

    `avoid` are the foreground points the labels must also clear, so a name never
    lands on top of a coloured marker.
    """
    texts = [ax.text(r["log2FoldChange"], r["neg_log10_padj"], str(r["gene_symbol"]),
                     fontsize=ANNOT_SIZE) for _, r in rows.iterrows()]
    if texts:
        static = {} if avoid is None or avoid.empty else {
            "x": avoid["log2FoldChange"].to_numpy(),
            "y": avoid["neg_log10_padj"].to_numpy()}
        # Generous repulsion from other labels AND from the foreground points, upward
        # only: a gene name is never abbreviated, so a crowded label is lifted into
        # the headroom above the cloud and leadered back to its point rather than
        # shortened or pushed down into the dense band along the FDR line.
        adjust_text(texts, ax=ax, expand=(1.5, 2.1), force_text=(0.6, 1.1),
                    force_static=(1.0, 1.8), force_pull=(0.002, 0.002),
                    only_move={"text": "xy+", "static": "xy+", "explode": "xy+",
                               "pull": "xy"},
                    min_arrow_len=6, time_lim=15, **static,
                    arrowprops=dict(arrowstyle="-", color="grey", lw=0.6))
    return [str(g) for g in rows["gene_symbol"]]


def _panel_note(ax, text: str) -> None:
    """Printed tallies inside a panel, on an opaque plate so the gate lines never
    cut through the digits."""
    ax.text(0.02, 0.98, text, transform=ax.transAxes, ha="left", va="top",
            fontsize=ANNOT_SIZE, zorder=6,
            bbox=dict(boxstyle="square,pad=0.4", facecolor="white", edgecolor="none",
                      alpha=0.9))


def _panel_title(ax, text: str) -> None:
    """Panel label at the config strip size, weighted like every other theme title."""
    ax.set_title(text, fontsize=float(_F["strip_size"]),
                 fontweight=plt.rcParams["axes.titleweight"])


def _outside_legend(ax, handles, ncol: int = 2):
    """Legend below the panel. Call AFTER fig.tight_layout() so the solver sizes the
    panel on its own content and the legend simply hangs below it — the saved tight
    bounding box still includes it."""
    return ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.11),
                     ncol=ncol, frameon=False, fontsize=LEGEND_SIZE)


# Vertical headroom above the most significant gene, so the printed tallies never
# sit on top of a data point, and the one canvas height both volcanoes are laid
# out AND saved at (widths come from the config: one-column and two-column).
VOLCANO_HEADROOM = 1.30
VOLCANO_HEIGHT = 7.6


def _volcano_frame(ax, de: pd.DataFrame) -> None:
    """Shared volcano furniture: background cloud, config gates, axis labels."""
    ax.scatter(de["log2FoldChange"], de["neg_log10_padj"], s=6, c=BG_COL,
               linewidths=0, alpha=0.55, zorder=1, rasterized=True)
    ax.axhline(-np.log10(DE_FDR), ls="--", color="black", lw=0.8, zorder=2)
    for x in (-DE_LFC, DE_LFC):
        ax.axvline(x, ls="--", color="black", lw=0.8, zorder=2)
    ax.set_xlabel("log2 fold change (synovial fluid / paired blood)")
    ax.set_ylim(-4, de["neg_log10_padj"].max() * VOLCANO_HEADROOM)


def plot_signature_volcano(de: pd.DataFrame, tal: dict):
    """Beat 4: the two mouse arms drawn on the plain Treg SF-vs-PB volcano."""
    fig, ax = plt.subplots(figsize=(float(_F["width"]), VOLCANO_HEIGHT))
    _volcano_frame(ax, de)
    ax.set_ylabel("-log10 FDR")
    for arm, z in (("down", 3), ("up", 4)):
        sub = de[de["arm"].eq(arm)]
        for purged, marker, size in ((False, "o", 62), (True, "^", 110)):
            s2 = sub[sub["hypoxia_purged"].eq(purged)]
            if s2.empty:
                continue
            ax.scatter(s2["log2FoldChange"], s2["neg_log10_padj"], marker=marker, s=size,
                       facecolors=to_rgba(ARM_COL[arm], 0.85), edgecolors=EDGE_HALO,
                       linewidths=0.5, zorder=z)
    _panel_note(ax,
                f"Gates FDR < {DE_FDR:g} and |log2FC| ≥ {DE_LFC:g}: "
                f"{tal['n_gates_sf']} SF-high, {tal['n_gates_blood']} blood-high\n"
                f"Up arm: {tal['n_up_gates_sf']} of {tal['n_up_tested']} tested clear them "
                f"SF-high, {tal['n_up_gates_blood']} blood-high\n"
                f"Down arm: {tal['n_down_gates_sf']} of {tal['n_down_tested']} SF-high, "
                f"{tal['n_down_gates_blood']} blood-high — no direction\n"
                f"Up-arm genes are {tal['n_up_gates_sf']} of the {tal['n_gates_sf']} "
                f"SF-high genes ({tal['pct_sf_response_from_up']:.1f}%)")
    ax.set_title("The mouse 39 °C up-arm sits in the SF-high corner of the Treg volcano")
    handles = [
        Line2D([0], [0], linestyle="none", marker="o", markersize=9,
               markerfacecolor=to_rgba(ARM_COL["up"], 0.85), markeredgecolor=EDGE_HALO,
               label=f"WT_heat up arm — {tal['n_up_tested']} of {tal['n_up_set']} tested here"),
        Line2D([0], [0], linestyle="none", marker="^", markersize=11,
               markerfacecolor=to_rgba(ARM_COL["up"], 0.85), markeredgecolor=EDGE_HALO,
               label="up-arm gene the hypoxia purge removes"),
        Line2D([0], [0], linestyle="none", marker="o", markersize=9,
               markerfacecolor=to_rgba(ARM_COL["down"], 0.85), markeredgecolor=EDGE_HALO,
               label=f"WT_heat down arm — {tal['n_down_tested']} of {tal['n_down_set']} tested here"),
        Line2D([0], [0], linestyle="none", marker="o", markersize=8, markerfacecolor=BG_COL,
               markeredgecolor=BG_COL,
               label=f"every other tested gene ({tal['n_tested']:,} in all)"),
    ]
    # Lay out first, then label and legend: adjust_text needs the final panel
    # geometry to decide what actually overlaps.
    fig.tight_layout()
    # Cap the labels at the config volcano_label_top, taken from the up arm's
    # gate-passing SF-high genes by significance — the claim's own genes.
    cand = de[de["arm"].eq("up") & de["passes_de_gates"] & de["sf_direction"].eq("SF_high")]
    labelled = _label_points(ax, cand.nlargest(LABEL_TOP, "neg_log10_padj"),
                             avoid=de[de["arm"].ne("not_in_signature")])
    _outside_legend(ax, handles, ncol=2)
    return fig, labelled


def plot_programs_axes_volcano(de: pd.DataFrame, tal: dict):
    """Beat 5: the same volcano twice — by leading-edge program, then by reference axis."""
    fig, (ax_p, ax_a) = plt.subplots(1, 2, sharex=True, sharey=True,
                                     figsize=(float(_F["width_wide"]), VOLCANO_HEIGHT))

    # --- left: which programs the up-arm genes belong to ---------------------
    _volcano_frame(ax_p, de)
    ax_p.set_ylabel("-log10 FDR")
    unann = de[de["arm"].eq("up") & de["le_program"].eq("not_annotated")]
    ax_p.scatter(unann["log2FoldChange"], unann["neg_log10_padj"], s=55,
                 facecolors=to_rgba(ARM_COL["up"], 0.35), edgecolors=EDGE_HALO,
                 linewidths=0.4, zorder=3)
    for key, _lab, colour in LE_PROGRAMS:
        sub = de[de["le_program"].eq(key)]
        if sub.empty:
            continue
        ax_p.scatter(sub["log2FoldChange"], sub["neg_log10_padj"], s=95,
                     facecolors=to_rgba(colour, 0.9), edgecolors=EDGE_HALO,
                     linewidths=0.6, zorder=4)
    # The heat-shock minority IS this panel's point, and its genes sit in the densest
    # part of the cloud, so they are named in the note rather than leadered from three
    # markers that a label would then cover. Every gene is in the source table.
    hsp = de[de["le_program"].eq("heat_shock_proteostasis")].nlargest(
        LABEL_TOP, "neg_log10_padj")
    named = [str(g) for g in hsp["gene_symbol"]]
    _panel_note(ax_p,
                f"The leading-edge annotation covers {tal['n_taxonomy_in_up']} of the "
                f"{tal['n_up_set']} up-arm genes\n"
                f"({tal['n_taxonomy_tested']} of them tested here) and "
                f"{tal['n_taxonomy_in_down']} down-arm genes, so it\n"
                "annotates the leading edge only and decomposes no whole set\n"
                f"Heat shock and proteostasis is {len(named)} genes: {', '.join(named)}")
    _panel_title(ax_p, "Effector and activation genes carry the surviving signal")
    handles_p = [Patch(facecolor=c, label=lab) for _, lab, c in LE_PROGRAMS]
    handles_p += [Patch(facecolor=to_rgba(ARM_COL["up"], 0.35),
                        label="up-arm gene the annotation does not cover"),
                  Patch(facecolor=BG_COL, label="every other tested gene")]

    # --- right: the two published reference axes ----------------------------
    _volcano_frame(ax_a, de)
    for key, marker, size in (("ifn_only", "o", 70), ("sting_specific", "s", 150)):
        sub = de[de["reference_axis"].eq(key)]
        in_arm = sub["arm"].ne("not_in_signature")
        for mouse, edge, lw in ((False, EDGE_HALO, 0.5), (True, ARM_COL["up"], 2.4)):
            s2 = sub[in_arm.eq(mouse)]
            if s2.empty:
                continue
            ax_a.scatter(s2["log2FoldChange"], s2["neg_log10_padj"], marker=marker, s=size,
                         facecolors=to_rgba(AXIS_COL[key], 0.8), edgecolors=edge,
                         linewidths=lw, zorder=5 if key == "sting_specific" else 4)
    _panel_note(ax_a,
                f"{tal['n_sting_specific_in_up']} of the {tal['n_sting_specific_set']} "
                f"STING-specific genes are mouse up-arm members\n"
                f"({', '.join(tal['genes_sting_specific_in_up'])}), against "
                f"{tal['n_ifn_only_in_up']} of {tal['n_ifn_only_set']} generic IFN genes\n"
                f"{tal['n_sting_specific_gates_sf']} of the "
                f"{tal['n_sting_specific_tested']} tested STING genes are SF-high, and "
                "they are\nalmost all genes the mouse program never contained")
    _panel_title(ax_a, "Almost none of those genes is a published STING gene")
    handles_a = [
        Line2D([0], [0], linestyle="none", marker="s", markersize=11,
               markerfacecolor=to_rgba(AXIS_COL["sting_specific"], 0.8),
               markeredgecolor=EDGE_HALO,
               label=f"{AXIS_LABEL['sting_specific']} — {tal['n_sting_specific_tested']} tested"),
        Line2D([0], [0], linestyle="none", marker="o", markersize=9,
               markerfacecolor=to_rgba(AXIS_COL["ifn_only"], 0.8), markeredgecolor=EDGE_HALO,
               label=f"{AXIS_LABEL['ifn_only']} — {tal['n_ifn_only_tested']} tested"),
        Line2D([0], [0], linestyle="none", marker="o", markersize=9,
               markerfacecolor=to_rgba(AXIS_COL["ifn_only"], 0.8),
               markeredgecolor=ARM_COL["up"], markeredgewidth=2.4,
               label="brown outline: also a mouse WT_heat gene"),
        Line2D([0], [0], linestyle="none", marker="o", markersize=8, markerfacecolor=BG_COL,
               markeredgecolor=BG_COL, label="every other tested gene"),
    ]
    fig.suptitle("The mouse 39 °C program that survives the hypoxia purge is an "
                 "activation program, and it is not a STING program",
                 fontsize=plt.rcParams["axes.titlesize"],
                 fontweight=plt.rcParams["axes.titleweight"])
    # Lay out first, then label and legend: adjust_text needs the final panel
    # geometry to decide what actually overlaps.
    fig.tight_layout()
    # The testable STING axis is 11 genes, so every one is named: a complete
    # enumeration rather than a top-N that hides members.
    labelled_a = _label_points(ax_a, de[de["reference_axis"].eq("sting_specific")],
                               avoid=de[de["reference_axis"].ne("none")])
    _outside_legend(ax_p, handles_p, ncol=2)
    _outside_legend(ax_a, handles_a, ncol=1)
    return fig, named + labelled_a


def volcano_source_table(de: pd.DataFrame, tal: dict, subset: pd.Series,
                         labelled: list) -> pd.DataFrame:
    """The per-gene rows behind a volcano, plus the printed denominators per row."""
    cols = ["gene_symbol", "log2FoldChange", "pvalue", "padj", "neg_log10_padj",
            "sf_direction", "passes_de_gates", "arm", "hypoxia_purged", "le_program",
            "reference_axis"]
    out = de.loc[subset, cols].copy()
    out["label_drawn"] = out["gene_symbol"].astype(str).isin(set(labelled))
    # Constant per row, so the counts printed in the figure are checkable here
    # without reopening the pseudobulk DE table.
    out["n_genes_tested"] = tal["n_tested"]
    out["n_gates_sf_high_all_genes"] = tal["n_gates_sf"]
    out["n_gates_blood_high_all_genes"] = tal["n_gates_blood"]
    out["contrast"] = "SF_vs_PB"
    out["evidence_tier"] = "primary_pseudobulk"
    return out.sort_values(["arm", "neg_log10_padj"], ascending=[True, False])


# ===========================================================================
def main() -> None:
    set_paper_style(config=FIG_CFG)
    purge_figures(STAGE, "heat_", overview=True, config=FIG_CFG)

    paired = purge_paired_table()
    fig = plot_purge_paired(paired)
    save_overview(
        fig, STAGE, "heat_purge_nes_paired",
        table=round_numeric_cols(paired),
        finding=("Deleting the 18 HALLMARK_HYPOXIA-overlap genes from the mouse 39 °C-derived "
                 "up-set takes 12 to 15 testable genes out of the arm and costs 0.129 to 0.165 "
                 "NES — 2.5915 to 2.4268 in Treg, 2.6809 to 2.5516 in Tcon, 2.0710 to 1.9261 in "
                 "CD8 — leaving all three significant, so the synovial-fluid enrichment is not "
                 "reducible to its HALLMARK_HYPOXIA-overlap gene content. That is a statement "
                 "about gene content and nothing else: it says nothing about temperature, and "
                 "nothing about whether hypoxia is a confound or a co-exposure, which are not "
                 "separable in cross-sectional human data."),
        script=SCRIPT, fn="plot_purge_paired",
        config_kv=(f"thresholds.gsea_fdr={FDR}; gsea_min_size={PARAMS.gsea_min_size}; "
                   f"gsea_nperm={PARAMS.gsea_nperm}"),
        input="03_results/09_heat_hypoxia/tables/gsea_{full,purged}_{treg,tcon,cd8}.csv",
        how_to_read=(
            "ANSWERS at confirmatory tier: donor-level pseudobulk within frozen sort labels, "
            "limma-voom then fgsea. Positive NES points toward synovial fluid. Each row pairs "
            "the full set (large diamond) with its purged form (small circle); the connecting "
            "bar is the NES cost. Warm brown is the up arm and cool blue the down arm. A dark "
            f"outline marks FDR below {FDR}. Right-hand text reports effective and nominal set "
            "sizes, the NES cost, and purged FDR. Distinguish the 18 genes removed from the "
            "frozen set from the 12 to 15 present in a ranked list. The Tcon down arm remains "
            "significant at the up arm's sign. This licenses a membership statement only. "
            "Correlative."),
        config=FIG_CFG, height=7.6,
    )
    plt.close(fig)

    coloc = colocalization_table()
    fig = plot_colocalization(coloc)
    save_overview(
        fig, STAGE, "heat_hypoxia_colocalization",
        table=round_numeric_cols(coloc),
        finding=("Within synovial-fluid cells the per-cell WT_heat_up score and the "
                 "HALLMARK_HYPOXIA score correlate only weakly (Spearman 0.08 to 0.20), so the "
                 "two scores are carried by largely different cells rather than reading out one "
                 "shared cell state. Per-cell tier: this corroborates the membership result and "
                 "cannot answer anything on its own."),
        script=SCRIPT, fn="plot_colocalization",
        config_kv="level=cell; tissue=synovial_fluid; evidence_tier=secondary_percell",
        input="03_results/09_heat_hypoxia/tables/heat_hypoxia_colocalization.csv",
        how_to_read=(
            "CORROBORATES and never answers: this per-cell tier cannot support a claim. Bars "
            "show within-SF cell-level correlation between WT_heat_up and HALLMARK_HYPOXIA "
            "AUCell scores, with Spearman dark and Pearson light; cell counts sit below each "
            "population. The y-axis spans -0.05 to 1, so read bar height rather than rank. "
            "Positive r means the scores tend to coincide. Donor-level SF means rest on 6 to "
            "7 donors and remain in the stage table. Never pool this diagnostic with "
            "pseudobulk NES or read it as directional evidence."),
        config=FIG_CFG,
    )
    plt.close(fig)

    # --- 3/4. the same mouse sets, seen in the plain Treg DE volcano ---------
    de = volcano_annotated_de()
    tal = volcano_tallies(de)

    fig, labelled = plot_signature_volcano(de, tal)
    save_overview(
        fig, STAGE, "heat_treg_volcano_signature",
        table=volcano_source_table(de, tal, de["arm"].ne("not_in_signature"), labelled),
        finding=(f"{tal['n_up_gates_sf']} of the {tal['n_up_tested']} testable mouse 39 °C "
                 f"up-arm genes clear the Treg SF-vs-PB significance gates on the "
                 f"synovial-fluid side against {tal['n_up_gates_blood']} on the blood side, so "
                 "the enrichment is visible in the plain differential-expression view, while "
                 f"accounting for {tal['pct_sf_response_from_up']:.1f}% of the SF-high "
                 "response."),
        script=SCRIPT, fn="plot_signature_volcano",
        config_kv=(f"thresholds.de_fdr={DE_FDR}; de_logfc={DE_LFC}; "
                   f"figures.volcano_label_top={LABEL_TOP}"),
        input="03_results/03_pseudobulk/tables/de_SFvsPB_treg.csv",
        how_to_read=(
            "Every tested Treg gene, x = log2 fold change synovial fluid over paired blood, "
            "y = -log10 FDR, dashed lines = the config gates. Warm brown = mouse up-arm "
            "member, cool blue = down-arm member, grey = every other gene, triangle = an "
            "up-arm gene the hypoxia purge removes. Membership is the only thing added to the "
            "committed DE table, and the printed tallies are counts of it. Gene names are "
            f"capped at the top {LABEL_TOP} up-arm genes by FDR and the rest are in the "
            "source table. The down arm scattering both ways is the honest caveat. Primary "
            "donor-pseudobulk tier, correlative."),
        config=FIG_CFG, height=VOLCANO_HEIGHT,
    )
    plt.close(fig)

    fig, labelled = plot_programs_axes_volcano(de, tal)
    keep = (de["arm"].eq("up") | de["le_program"].ne("not_annotated")
            | de["reference_axis"].ne("none"))
    save_overview(
        fig, STAGE, "heat_treg_volcano_programs",
        table=volcano_source_table(de, tal, keep, labelled),
        finding=(f"Only {tal['n_sting_specific_in_up']} of the "
                 f"{tal['n_sting_specific_set']} published IFN-independent STING-activation "
                 f"genes and {tal['n_ifn_only_in_up']} of {tal['n_ifn_only_set']} generic "
                 "type-I IFN genes are in the mouse 39 °C up-arm, so the SF-high program the "
                 "purge leaves standing is an effector and activation program that shares "
                 "almost nothing with the STING reference axis."),
        script=SCRIPT, fn="plot_programs_axes_volcano",
        config_kv=(f"thresholds.de_fdr={DE_FDR}; de_logfc={DE_LFC}; "
                   "taxonomy=00_data/references/heat_leadingedge_taxonomy"),
        input=("03_results/03_pseudobulk/tables/de_SFvsPB_treg.csv, "
               "00_data/references/heat_leadingedge_taxonomy/leadingedge_gene_taxonomy.csv, "
               "../sting_positive_control/03_results/06_reference_axis/signatures/"),
        how_to_read=(
            "Two views of one volcano, same axes as the signature volcano. Left colours the "
            f"up-arm genes by leading-edge program, with the {tal['n_taxonomy_in_up']}-gene "
            "annotation covering the leading edge only and pale brown marking the up-arm "
            "genes it leaves unlabelled, so read it as an annotation and not a decomposition "
            "of the 199-gene set. Right draws the two frozen reference axes: black squares "
            "are the published IFN-independent STING-activation genes, all named, blue "
            "circles the generic type-I IFN program, and a brown outline means the gene is "
            "also a mouse signature member. The heat-shock trio is named in the left panel "
            "note rather than at its markers. Primary donor-pseudobulk tier, correlative."),
        config=FIG_CFG, height=VOLCANO_HEIGHT, wide=True,
    )
    plt.close(fig)
    print("[09_heat_hypoxia_viz] wrote 4 overviews (purge pairing, co-localization, "
          "Treg volcano, programs vs reference axes)")


if __name__ == "__main__":
    main()
