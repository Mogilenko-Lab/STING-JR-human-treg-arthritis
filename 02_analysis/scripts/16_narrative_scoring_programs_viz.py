#!/usr/bin/env python
"""
16_narrative_scoring_programs_viz.py: VIZ ONLY (no statistics).
=============================================================================
The distributions behind `umap_full_programs`, on the geometry of `arm_score_violins` so
the two violin figures pair as their two maps do. Six panels in the map's panel order:
the published interferon-independent STING panel, the generic type-I interferon axis and
Hallmark interferon alpha, ruled off from three inflammation and activation programs.

What it adds over the map: the size of each tissue difference, and whether that size
changes between sort labels. The delta row carries both on one shared axis, so a
Treg-selective effect shows as a marker standing away from its neighbours in one panel and
sitting level with them in the rest.

It also shows the zero-inflated baseline the map cannot: sting_specific_published is 18
genes scored and leaves at least half of Treg blood cells at exactly zero, which draws as
a body pinned to the axis. Zero fractions are in the source table.

TIER. Annotation, one vote per cell over 7 donors of unequal yield. Testing belongs to
03_results/14_unbiased_enrichment/figures/_overview/program_nes_by_cell_state.png.

Panel order comes from `analysis_config.yaml::percell_map_panels.program_strip`, which
the map, its Treg-only counterpart, this figure and the sweep-coverage audit all read.

Input:
  03_results/interactive/16_narrative_embedding.parquet     99,915 cells x score columns
  03_results/16_narrative_scoring/tables/narrative_scoring_manifest.csv   genes scored

Output (03_results/16_narrative_scoring/):
  figures/_overview/program_score_violins.{pdf,png}
  tables/_overview/program_score_violins.csv   36 rows: set x sort label x tissue summary
  README.md                                    caption (via save_overview)

Run in-container from the compartment root, AFTER 16_narrative_scoring.py:
  python 02_analysis/scripts/16_narrative_scoring_programs_viz.py
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
    PATHS,
    POPULATION_COLORS,
    PROGRAM_STRIP_SETS,
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
SCRIPT = "02_analysis/scripts/16_narrative_scoring_programs_viz.py"
STEM = "program_score_violins"
SUBSTRATE_PARQUET = "16_narrative_embedding.parquet"
MANIFEST_CSV = "narrative_scoring_manifest.csv"

MAP_FIGURE = "figures/_overview/umap_full_programs.png in this directory"
CONFIRMATORY_FIGURE = ("03_results/14_unbiased_enrichment/figures/_overview/"
                       "program_nes_by_cell_state.png")

TISSUES = [(TISSUE_NUM, "synovial fluid"), (TISSUE_DEN, "paired blood")]


def panels() -> list:
    """[(column, set name, genes scored)] in the map's panel order.

    The size is `n_genes_found_in_object` from the scoring manifest — the number of genes
    the score was really computed over, and the same number the map's panel titles carry.
    """
    manifest = pd.read_csv(PATHS.tables(STAGE) / MANIFEST_CSV)
    sizes = dict(zip(manifest["set_name"], manifest["n_genes_found_in_object"].astype(int)))
    missing = [s for s in PROGRAM_STRIP_SETS if s not in sizes]
    if missing:
        raise KeyError(
            f"[16_programs_viz] {MANIFEST_CSV} carries no row for {missing}, so the program "
            "strip names a set this compartment never scored per cell. Run "
            "02_analysis/scripts/16_narrative_scoring.py, or fix "
            "analysis_config.yaml::percell_map_panels.program_strip.")
    return [(f"{s}_AUCell", s, sizes[s]) for s in PROGRAM_STRIP_SETS]


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
        raise KeyError(f"[16_programs_viz] substrate is missing columns: {missing}")

    summary = sv.summary_table(df, spec, TISSUES)
    width, height = 16.0, 11.4
    fig, dlim = sv.build_figure(
        df, spec, TISSUES, TISSUE_COLORS, POPULATION_COLORS, summary, FIG_CFG,
        title=("Per-cell scores by sort label and tissue, "
               "cGAS-STING and interferon against inflammation"),
        subtitle=("Annotation tier. AUCell runs 0 to 1, one value per cell, each panel on its "
                  "own y axis. The grey row gives the synovial-against-blood separation "
                  "inside each sort label."),
        width=width, height=height)

    def cell(set_name: str, state: str) -> pd.Series:
        return summary[summary["gene_set"].eq(set_name)
                       & summary["cell_state"].eq(state)].iloc[0]

    sting = "sting_specific_published"
    sting_zero = cell(sting, "Treg")
    sting_zero_pb = float(summary[summary["gene_set"].eq(sting)
                                  & summary["cell_state"].eq("Treg")
                                  & summary["tissue"].eq("paired blood")]["frac_at_zero"].iloc[0])
    deltas = {(s, st): float(cell(s, st)["delta_cliffs_sf_vs_pb"])
              for _c, s, _n in spec for st in sv.CELL_STATES}
    strongest = max(deltas, key=deltas.get)
    weakest = min(deltas, key=deltas.get)
    treg_gap = {s: deltas[(s, "Treg")] - max(deltas[(s, "Tcon")], deltas[(s, "CD8")])
                for _c, s, _n in spec}
    treg_lead = max(treg_gap, key=treg_gap.get)
    means = {s: (float(cell(s, "Treg")["mean"]), float(cell(s, "Tcon")["mean"]),
                 float(cell(s, "CD8")["mean"])) for _c, s, _n in spec}

    save_overview(
        fig, STAGE, STEM,
        table=summary,
        finding=(
            f"All six lenses score higher in synovial fluid than in paired blood in all "
            f"three sort labels, so a synovial-side colouring is shared by the cGAS-STING "
            f"family and by the inflammation programs it has to be told apart from, and "
            f"the separations span Cliff's δ {deltas[weakest]:.3f} "
            f"({weakest[0]} in {weakest[1]}) to {deltas[strongest]:.3f} "
            f"({strongest[0]} in {strongest[1]}). Reading the tissue separation by sort "
            f"label, the widest Treg lead over both other labels belongs to {treg_lead} "
            f"(δ {deltas[(treg_lead, 'Treg')]:.3f} in Treg against "
            f"{deltas[(treg_lead, 'Tcon')]:.3f} and {deltas[(treg_lead, 'CD8')]:.3f}), a "
            f"gap of {treg_gap[treg_lead]:.3f}, which is the scale of Treg-selectivity this "
            f"per-cell channel carries at all. The {int(sting_zero['genes_scored'])}-gene "
            f"published STING panel sits lowest of the six in the Treg gate "
            f"(per-cell mean {means[sting][0]:.4f} against {means[sting][1]:.4f} in Tcon "
            f"and {means[sting][2]:.4f} in CD8) and leaves "
            f"{100 * sting_zero_pb:.0f}% of Treg blood cells at exactly zero, so its Treg "
            "tissue difference rests on a zero-inflated baseline."
        ),
        script=SCRIPT, fn="build_figure",
        config_kv=("metric = AUCell (rank-based, 0 to 1); panels = "
                   "percell_map_panels.program_strip = "
                   + ", ".join(f"{s} ({n} genes scored)" for _c, s, n in spec)
                   + f"; y_pad_frac = {sv.Y_PAD_FRAC}; delta axis = ±{dlim:g} shared"),
        input=("03_results/interactive/16_narrative_embedding.parquet, "
               f"03_results/{STAGE}/tables/{MANIFEST_CSV}"),
        how_to_read=(
            sv.how_to_read(df, summary, spec, dlim, MAP_FIGURE, CONFIRMATORY_FIGURE)
            + " One reading this figure supports and its map does not: the zero fraction. "
            "A thin set leaves many cells at exactly zero, which draws as a violin body "
            "pinned to the axis and is given per violin in the source table, and a tissue "
            "difference resting on that baseline tracks how many cells score anything at "
            "all."
        ),
        config=FIG_CFG, width=width, height=height)

    print(f"[16_narrative_scoring_programs_viz] wrote {STEM}: {len(spec)} panels, "
          f"{len(summary)} summarised violins, shared delta axis ±{dlim:g}")


if __name__ == "__main__":
    main()
