# 03_pseudobulk — artifact captions

_**Abbreviations:** SF = synovial fluid (inflamed joint); PB = peripheral blood. The SF-vs-PB contrast is paired within each of the 7 JIA donors. Treg = CD4⁺CD127ˡᵒCD25⁺ regulatory; Tcon = CD4⁺CD25⁻ conventional; CD8 = CD8⁺CD45RO⁺ memory._

## figures/_overview/pseudobulk_pca.png

Pseudobulk samples separate by tissue and label without a single donor
dominating an axis, so donor pseudobulk is well-posed for SF-vs-PB DE.

**How to read:** Each point = one donor x tissue x label pseudobulk (log-CPM, top-2000
var genes). Circle=SF, square=PB, color=donor. Look for tissue
separation and NO single-donor axis dominance. Display transform only.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03_pseudobulk_de_viz.py` | `main` | `thresholds.pseudobulk_min_cells (strata filter)` | `03_results/03_pseudobulk/tables/pseudobulk_counts.csv` |

## figures/_overview/treg_volcano.png

Synovial-fluid Tregs carry a reproducible SF-vs-PB transcriptional
program (significant up/down genes), the substrate the mouse signature
is tested against.

**How to read:** x=log2FC SF/PB, y=-log10 padj; orange = significant. Dashed lines =
FDR + |log2FC| gates. Correlative DE, top-500 genes tabulated.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03_pseudobulk_de_viz.py` | `main` | `thresholds.de_fdr=0.05; de_logfc=1.0` | `03_results/03_pseudobulk/tables/de_SFvsPB_treg.csv` |

## figures/_overview/de_count_bar.png

All three sorted populations yield significant SF-vs-PB DE, so each
has a ranked list powered for signature enrichment.

**How to read:** Bar = # significant SF-vs-PB DE genes per population (Treg/Tcon/CD8).
Confirms each arm has enough signal to rank for fgsea. Diagnostic.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03_pseudobulk_de_viz.py` | `main` | `thresholds.de_fdr=0.05; de_logfc=1.0` | `03_results/03_pseudobulk/tables/de_summary.csv` |
