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

## tables/pseudobulk_coldata.csv

39 of the 42 possible donor x tissue x label strata survive aggregation, each
carrying 266–4,454 cells with none falling under the 20-cell floor; the three
absences (SF Treg in patient 5, PB Tcon and PB CD8 in patient 3) are why every
population's contrast runs 6 donors in one tissue arm against 7 in the other
rather than a clean 7-vs-7.

**How to read:** One row per column of `pseudobulk_counts.csv`, keyed
`donor|tissue|coarse_label`. `donor` is the blocking factor, `tissue` the contrast
factor (`synovial_fluid` = numerator, `peripheral_blood` = denominator),
`coarse_label` the frozen population, `n_cells` the cells summed into that column.
Cell counts are aggregation weights, not a result — read them to judge whether a
stratum is thin. This is the design table behind the primary donor-pseudobulk tier:
a stratum missing here removes that donor from the population's contrast.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03_pseudobulk_de.py` | `main` | `thresholds.pseudobulk_min_cells=20` | `03_results/objects/02_annotation.h5ad` |

## tables/de_summary.csv

Every sorted population clears the bar for a powered enrichment test: the paired
`~ donor + tissue` model returns 1,584 significant SF-vs-PB genes in Treg out of
14,714 ranked, 1,560 of 15,229 in Tcon, and 1,795 of 15,070 in CD8, with no arm
dropped for lack of donors.

**How to read:** One row per sorted population. `n_sf`/`n_pb` count donor strata in
each tissue arm (strata, not cells); `model` is the fitted design, or `skipped` if
an arm fell under the donor floor; `n_sig_de` counts genes passing both the FDR and
the |log2FC| gate; `n_ranked` is the length of the sign-preserving ranked list
handed to pre-ranked GSEA. Diagnostic/powering table — the confirmatory statistics
are the donor-pseudobulk NES computed downstream from `ranked_{treg,tcon,cd8}.tsv`,
not the gene counts here.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03_pseudobulk_de.py` | `main` | `thresholds.de_fdr=0.05; de_logfc=1.0; pseudobulk_min_donors=3` | `03_results/objects/02_annotation.h5ad` |
