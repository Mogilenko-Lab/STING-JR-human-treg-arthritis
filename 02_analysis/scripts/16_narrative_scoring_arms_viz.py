#!/usr/bin/env python
"""
16_narrative_scoring_arms_viz.py: VIZ ONLY (no statistics).
=============================================================================
The distributions behind `umap_full_arms`, in that map's panel order: the three mouse
39 °C-derived up arms, then the three lenses ruled off beside them there. Six violins
per panel (three frozen sort labels x two tissues) over a Cliff's delta row.

The lenses share the figure because they share the map: all six lean synovial, which is
what makes the arms' colouring undistinctive, and the lenses supply what an arm's offset
is read against.

eTreg_up reads as a ceiling: another cohort's synovial-versus-blood contrast, asking this
contrast's own question, and scoring near the top of the whole enrichment sweep.

TIER. Annotation, one vote per cell over 7 donors of unequal yield. Ranking the sort
labels belongs to
03_results/14_unbiased_enrichment/figures/_overview/arm_nes_by_cell_state.png.

Panel order comes from `analysis_config.yaml::percell_map_panels.arm_strip`, which the
map, its Treg-only counterpart, this figure and the sweep-coverage audit all read.

Input:
  03_results/interactive/16_narrative_embedding.parquet     99,915 cells x score columns
  03_results/16_narrative_scoring/tables/narrative_scoring_manifest.csv   genes scored

Output (03_results/16_narrative_scoring/):
  figures/_overview/arm_score_violins.{pdf,png}
  tables/_overview/arm_score_violins.csv    36 rows: set x sort label x tissue summary
  README.md                                 caption (via save_overview)

Run in-container from the compartment root, AFTER 16_narrative_scoring.py:
  python 02_analysis/scripts/16_narrative_scoring_arms_viz.py
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import pandas as pd  # noqa: E402

COMPARTMENT_ROOT = Path(__file__).resolve().parents[2]

sys.path.insert(0, str(COMPARTMENT_ROOT))
sys.path.insert(0, str(COMPARTMENT_ROOT / "02_analysis"))
os.chdir(COMPARTMENT_ROOT)

from config import (  # noqa: E402
    ARM_STRIP_SETS,
    PATHS,
    POPULATION_COLORS,
    TISSUE_COLORS,
    TISSUE_DEN,
    TISSUE_NUM,
)
from helpers import score_violins as sv  # noqa: E402
from helpers.figure_style import (  # noqa: E402
    FIG_CFG,
    purge_figures,
    save_overview,
    set_paper_style,
)

STAGE = "16_narrative_scoring"
SCRIPT = "02_analysis/scripts/16_narrative_scoring_arms_viz.py"
STEM = "arm_score_violins"
SUBSTRATE_PARQUET = "16_narrative_embedding.parquet"
MANIFEST_CSV = "narrative_scoring_manifest.csv"

MAP_FIGURE = "figures/_overview/umap_full_arms.png in this directory"
CONFIRMATORY_FIGURE = ("03_results/14_unbiased_enrichment/figures/_overview/"
                       "arm_nes_by_cell_state.png")

# Display label per tissue; the hue comes from the one tissue palette in the config.
TISSUES = [(TISSUE_NUM, "synovial fluid"), (TISSUE_DEN, "paired blood")]


def panels() -> list:
    """[(column, set name, genes scored)] in the map's panel order.

    The size is `n_genes_found_in_object` from the scoring manifest — the number of genes
    the score was really computed over, and the same number the map's panel titles carry,
    so a title reads identically in both figures.
    """
    manifest = pd.read_csv(PATHS.tables(STAGE) / MANIFEST_CSV)
    sizes = dict(zip(manifest["set_name"], manifest["n_genes_found_in_object"].astype(int)))
    missing = [s for s in ARM_STRIP_SETS if s not in sizes]
    if missing:
        raise KeyError(
            f"[16_arms_viz] {MANIFEST_CSV} carries no row for {missing}, so the arm strip "
            "names a set this compartment never scored per cell. Run "
            "02_analysis/scripts/16_narrative_scoring.py, or fix "
            "analysis_config.yaml::percell_map_panels.arm_strip.")
    return [(f"{s}_AUCell", s, sizes[s]) for s in ARM_STRIP_SETS]


def main() -> None:
    set_paper_style(config=FIG_CFG)
    purge_figures(STAGE, STEM, overview=True, config=FIG_CFG)

    pq_path = PATHS.interactive_dir() / SUBSTRATE_PARQUET
    if not pq_path.exists():
        raise FileNotFoundError(
            f"per-cell substrate missing at {pq_path}; "
            "run 02_analysis/scripts/16_narrative_scoring.py first")

    spec = panels()
    need = ["coarse_label", "tissue", "donor"] + [c for c, _s, _n in spec]
    df = pd.read_parquet(pq_path, columns=need)
    missing = sorted(set(need) - set(df.columns))
    if missing:
        raise KeyError(f"[16_arms_viz] substrate is missing columns: {missing}")

    summary = sv.summary_table(df, spec, TISSUES)
    width, height = 16.0, 11.4
    fig, dlim = sv.build_figure(
        df, spec, TISSUES, TISSUE_COLORS, POPULATION_COLORS, summary, FIG_CFG,
        title=("Per-cell scores by sort label and tissue, "
               "the mouse up arms beside three lenses"),
        subtitle=("Annotation tier. AUCell runs 0 to 1, one value per cell, each panel on its "
                  "own y axis. The grey row gives the synovial-against-blood separation "
                  "inside each sort label."),
        width=width, height=height)

    def cell(set_name: str, state: str) -> pd.Series:
        hit = summary[summary["gene_set"].eq(set_name) & summary["cell_state"].eq(state)]
        return hit.iloc[0]

    wt = {s: float(cell("WT_heat_up", s)["delta_cliffs_sf_vs_pb"]) for s in sv.CELL_STATES}
    wt_treg_top = wt["Treg"] >= max(wt.values())
    lens_treg = {s: float(cell(s, "Treg")["delta_cliffs_sf_vs_pb"])
                 for _c, s, _n in spec[3:]}
    best_lens = max(lens_treg, key=lens_treg.get)
    inter = cell("Interaction_up", "Treg")
    zero_pct = 100.0 * float(summary[summary["gene_set"].eq("Interaction_up")]
                             ["frac_at_zero"].mean())
    agree = summary[summary["gene_set"].eq("WT_heat_up")]
    all_donors_agree = bool((agree["delta_donors_agreeing_in_sign"]
                             == agree["delta_n_donors_paired"]).all())

    save_overview(
        fig, STAGE, STEM,
        table=summary,
        finding=(
            f"Every one of the six sets drawn on the score map scores higher in synovial "
            f"fluid than in paired blood in all three sort labels, and the mouse "
            f"39 °C-derived up arm's separation is Cliff's δ "
            f"{wt['Treg']:.3f} in Treg, {wt['Tcon']:.3f} in Tcon and {wt['CD8']:.3f} in "
            f"CD8 — "
            + ("largest in Treg of the three, and the three lenses beside it separate "
               if wt_treg_top else
               "so the separation is shared across the sort labels, and the lenses beside "
               "it separate ")
            + f"the same tissues at least as far ({best_lens} reaches "
            f"{lens_treg[best_lens]:.3f} in Treg). "
            + ("Every donor that carries both tissues agrees in sign with the pooled value "
               "on that arm, in all three sort labels. " if all_donors_agree else
               "The donor markers show where that agreement breaks, panel by panel. ")
            + "This per-cell ordering of the sort labels stands on its own and the donor-level "
              "panel is the one that ranks them, where the same arm reaches NES 2.68 in Tcon "
              "against 2.59 in Treg. "
            + f"The {int(inter['genes_scored'])}-gene Interaction_up arm leaves "
            f"{zero_pct:.0f}% of cells at exactly zero, which is what a set that thin does "
            "when none of its genes reaches a cell's top-ranked genes."
        ),
        script=SCRIPT, fn="build_figure",
        config_kv=("metric = AUCell (rank-based, 0 to 1); panels = "
                   "percell_map_panels.arm_strip = "
                   + ", ".join(f"{s} ({n} genes scored)" for _c, s, n in spec)
                   + f"; y_pad_frac = {sv.Y_PAD_FRAC}; delta axis = ±{dlim:g} shared"),
        input=("03_results/interactive/16_narrative_embedding.parquet, "
               f"03_results/{STAGE}/tables/{MANIFEST_CSV}"),
        how_to_read=sv.how_to_read(df, summary, spec, dlim, MAP_FIGURE,
                                   CONFIRMATORY_FIGURE),
        config=FIG_CFG, width=width, height=height)

    print(f"[16_narrative_scoring_arms_viz] wrote {STEM}: {len(spec)} panels, "
          f"{len(summary)} summarised violins, shared delta axis ±{dlim:g}")


if __name__ == "__main__":
    main()
