# 12_treg_localisation -- Treg-only niche localisation

Secondary and corroborative tier evaluation of per-cell AUCell score distributions across synovial fluid (SF) versus peripheral blood (PB) niches in sorted JIA CD4+ Tregs.

This stage is a compute resource and publishes no figure. Its two tables carry the per-cell scores and their per-signature summaries; the statistical claims about the same contrast rest on donor-level pseudobulk DE, not on per-cell distributions.

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
