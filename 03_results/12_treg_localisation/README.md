# 12_treg_localisation -- Treg-only niche localisation

Secondary and corroborative tier evaluation of per-cell AUCell score distributions across synovial fluid (SF) versus peripheral blood (PB) niches in sorted JIA CD4+ Tregs.

## figures/_overview/treg_localisation.png

Across sorted JIA Tregs, per-cell AUCell scores for WT_heat_up,
HALLMARK_HYPOXIA, and WT_heat_up_purged_hypoxia are consistently
higher in synovial fluid than in peripheral blood; the cGAS-dependent
Interaction_fdrOnly_up (18-gene gate) shows modest SF elevation,
whereas the 7-gene Interaction_up gate is underpowered and dominated
by detection noise at per-cell resolution.

**How to read:** Box plots show per-cell AUCell score distributions for sorted CD4_Treg
cells in peripheral blood (PB, blue) versus synovial fluid (SF,
vermillion). Solid lines indicate medians; dashed red lines indicate
means. Set sizes (nominal and in-dataset effective) and power-band
classifications are declared for each panel. This panel is hypothesis-
generating tier; primary statistical claims are carried by donor-level
pseudobulk DE.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_treg_localisation_viz.py` | `build_figure` | `gsea_fdr = 0.05, percell_score_ncores = 4` | `03_results/12_treg_localisation/tables/treg_per_cell_scores.csv` |

## tables/treg_localisation_summary.csv

Summary statistics (cell counts, donor counts, mean, median, IQR, set sizes, and power bands) per signature and tissue arm for sorted CD4+ Tregs.

**How to read:** Each row gives summary metrics for one signature and tissue niche (`synovial_fluid` vs `peripheral_blood`) in sorted CD4+ Tregs (`coarse_label == "Treg"`). Set sizes reflect nominal gene list counts and effective in-dataset matches. Power bands follow standard thresholds: `testable` (≥15 genes), `underpowered_reported` (5–14 genes), and `untestable` (<5 genes).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_treg_localisation.py` | `main` | `percell_score_ncores = 4` | `03_results/objects/02_annotation.h5ad` |

## tables/treg_per_cell_scores.csv

Per-cell AUCell score table for all cells in the substrate across five evaluated signatures.

**How to read:** Per-cell rank-based AUCell scores for each cell (indexed by barcode), including donor, tissue, coarse label, and coordinate metadata.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_treg_localisation.py` | `main` | `percell_score_ncores = 4` | `03_results/objects/02_annotation.h5ad` |
