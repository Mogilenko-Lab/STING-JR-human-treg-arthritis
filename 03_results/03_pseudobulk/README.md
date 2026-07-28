# 03_pseudobulk — artifact captions

_**Abbreviations:** SF = synovial fluid (inflamed joint); PB = peripheral blood. The cohort contains 7 JIA donors, of whom 6 span both arms in each analyzed population after QC. Treg = CD4⁺CD127ˡᵒCD25⁺ regulatory; Tcon = CD4⁺CD25⁻ conventional; CD8 = CD8⁺CD45RO⁺ memory._

Differential expression runs in R. Aggregation to donor x tissue x label pseudobulk is done in
Python (`03a_pseudobulk_export.py`, counts only), and the model is fitted in R
(`03b_pseudobulk_de.R`: `filterByExpr` -> TMM -> voom -> `lmFit` -> `eBayes(robust=TRUE)`).
The ranking metric handed to every downstream enrichment stage is the signed moderated
t-statistic. The DE tables carry an engine-agnostic column schema (`log2FoldChange`, `stat`,
`pvalue`, `padj`) plus a `de_engine` column recording which engine produced them, so a future
engine change does not ripple into the consumers.

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
`~ donor + tissue` model returns 1,797 significant SF-vs-PB genes in Treg out of
13,999 ranked, 1,949 of 14,411 in Tcon, and 1,695 of 14,014 in CD8, with no arm
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

## tables/pseudobulk_counts.csv

Summing raw UMI counts within donor x tissue x label gives 39 pseudobulk strata across 21,740
genes, with no stratum falling below the per-stratum cell floor.

**How to read:** Rows are strata (`<donor>_<tissue>_<label>`), columns are Ensembl gene ids, values
are summed raw integer UMIs — not normalised, because the DE engine expects raw counts and does its
own library-size normalisation. Pair with `pseudobulk_coldata.csv` for the design factors. This is
the matrix handed across the Python-to-R seam.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03a_pseudobulk_export.py` | `main` | `pseudobulk_min_cells` | `03_results/objects/02_annotation.h5ad` (`counts` layer) |

## tables/gene_symbols.csv

All 21,740 exported Ensembl ids carry a gene symbol and every symbol is distinct, so the
Ensembl-to-symbol rename loses no gene and creates no collision in this dataset.

**How to read:** Two columns, `ensembl_id` and `gene_symbol`, in the same order as the columns of
`pseudobulk_counts.csv`. This file exists because the counts matrix is keyed by Ensembl id while
every reference gene set the ranked lists are tested against — the mouse-projection signatures,
MSigDB Hallmark, the curated HSR lens — matches on HGNC symbol. The R stage asserts this file is
present before ranking; without it the ranked lists would carry Ensembl ids and intersect the
reference sets at approximately zero, which enrichment tools report as an empty result rather than
an error.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03a_pseudobulk_export.py` | `main` | — | `03_results/objects/02_annotation.h5ad` (`var['gene_symbol']`) |

## tables/de_SFvsPB_{treg,tcon,cd8}.csv

Paired SF-vs-PB differential expression recovers 1,797 / 1,949 / 1,695 genes at FDR < 0.05 and
|log2FC| >= 1 in Treg / Tcon / CD8, all three fitted on the donor-paired model.

**How to read:** One file per sorted population, one row per gene, sorted by p-value. `stat` is the
signed moderated t-statistic and `log2FoldChange` the log2 fold change, both positive when the gene
is higher in synovial fluid; `padj` is BH across genes within that population. `avg_expr` is
limma's average log2-CPM — deliberately not named `baseMean`, which denotes a differently scaled
quantity. `model` records the fitted design and `n_paired_donors` how many donors appeared in both
tissues; `de_engine` records the engine. Primary donor-pseudobulk tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03b_pseudobulk_de.R` | `(top-level)` | `design.tissue_key`, `design.donor_key`, `thresholds.de_fdr` | `03_results/03_pseudobulk/tables/pseudobulk_counts.csv`, `pseudobulk_coldata.csv`, `gene_symbols.csv` |

## tables/ranked_{treg,tcon,cd8}.tsv

The signed ranked lists carry 13,999 / 14,411 / 14,014 genes for Treg / Tcon / CD8 after
expression filtering, and are the single input every downstream enrichment stage reads.

**How to read:** Two columns, no header: HGNC symbol and the signed moderated t-statistic, sorted
descending, so the top of the file is most synovial-fluid-up and the bottom most blood-up. One row
per symbol — where several Ensembl ids share a symbol the most extreme |t| is kept, because a
duplicated symbol would corrupt a pre-ranked enrichment run. Gene count differs per population
because `filterByExpr` is applied within each population's own design.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03b_pseudobulk_de.R` | `(top-level)` | `design.tissue_key`, `design.donor_key` | `03_results/03_pseudobulk/tables/pseudobulk_counts.csv`, `pseudobulk_coldata.csv`, `gene_symbols.csv` |

## tables/de_engine_migration.csv

The limma-voom migration validation reported in commit `799e4fe` cannot be recomputed from tracked
ranked-list artifacts in this repository.

**How to read:** This is a provenance table, not a new validation. `not_reproducible` means the
reported value exists only in the commit message and the required pre-migration ranked lists are not
tracked. Current `ranked_*.tsv` files exist on disk only, so they are insufficient to reconstruct the
before/after comparison.

| Script | Function | Config | Input |
|---|---|---|---|
| repository audit | `git show`, `git ls-files`, `find` | `commit=799e4fe` | commit message and tracked file index |
