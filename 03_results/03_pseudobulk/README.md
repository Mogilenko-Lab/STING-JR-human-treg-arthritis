# 03_pseudobulk: artifact captions

_**Abbreviations:** SF = synovial fluid (inflamed joint); PB = peripheral blood; Treg =
CD4⁺CD127ˡᵒCD25⁺ regulatory; Tcon = CD4⁺CD25⁻ conventional; CD8 = CD8⁺CD45RO⁺ memory. The cohort
holds 7 JIA donors, of whom 6 span both arms in each analyzed population after QC._

Differential expression runs in R across an explicit file seam. Aggregation to donor × tissue ×
label pseudobulk (one profile per donor, tissue and frozen label) is Python
(`03a_pseudobulk_export.py`, counts only). The model is R (`03b_pseudobulk_de.R`: `filterByExpr` →
TMM → voom → `lmFit` → `eBayes(robust=TRUE)`).

Every downstream enrichment stage ranks on the signed moderated t-statistic. The DE tables carry an
engine-agnostic column schema (`log2FoldChange`, `stat`, `pvalue`, `padj`) plus a `de_engine`
column naming the engine that produced them, so a later engine change stays contained to this
stage.

## figures/_overview/pseudobulk_pca.png

Pseudobulk samples separate by tissue and by label, and no single donor dominates an axis, so donor
pseudobulk is well-posed for SF-versus-PB differential expression.

**How to read:** Each point is one donor × tissue × label pseudobulk profile (log-CPM, top 2,000
variable genes). Circle is SF, square PB, colour donor. Read for tissue separation and for the
absence of a single-donor axis. Display transform only.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03_pseudobulk_de_viz.py` | `main` | `thresholds.pseudobulk_min_cells (strata filter)` | `03_results/03_pseudobulk/tables/pseudobulk_counts.csv` |

## figures/_overview/treg_volcano.png

Synovial-fluid Tregs carry a reproducible SF-versus-PB transcriptional program, 1,797 genes at FDR
< 0.05 and |log2FC| ≥ 1. This is the substrate the mouse-derived signature is tested against.

**How to read:** x is log2 fold change SF over PB, y is −log10 padj, orange marks significance.
Dashed lines are the FDR and |log2FC| gates. The top 500 genes are tabulated alongside. Correlative
donor-pseudobulk DE.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03_pseudobulk_de_viz.py` | `main` | `thresholds.de_fdr=0.05; de_logfc=1.0` | `03_results/03_pseudobulk/tables/de_SFvsPB_treg.csv` |

## figures/_overview/de_count_bar.png

All three sorted populations yield significant SF-versus-PB differential expression — 1,797 genes
in Treg, 1,949 in Tcon, 1,695 in CD8 — so each has a ranked list powered for signature enrichment.

**How to read:** One bar per population, height the count of significant SF-versus-PB DE genes.
Read it to confirm each arm carries enough signal to rank for pre-ranked GSEA. Diagnostic.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03_pseudobulk_de_viz.py` | `main` | `thresholds.de_fdr=0.05; de_logfc=1.0` | `03_results/03_pseudobulk/tables/de_summary.csv` |

## tables/pseudobulk_coldata.csv

39 of the 42 possible donor × tissue × label strata survive aggregation, each holding 266 to 4,454
cells, all above the 20-cell floor. The three absences — SF Treg in patient 5, PB Tcon and PB CD8 in
patient 3 — are why every population's contrast runs 6 donors in one tissue arm against 7 in the
other.

**How to read:** One row per column of `pseudobulk_counts.csv`, keyed `donor|tissue|coarse_label`.
`donor` is the blocking factor; `tissue` is the contrast factor, with `synovial_fluid` the numerator
and `peripheral_blood` the denominator; `coarse_label` is the frozen population; `n_cells` is the
cells summed into that column. Cell counts are aggregation weights — read them to judge whether a
stratum is thin. This is the design table behind the primary donor-pseudobulk tier: a stratum
missing here removes that donor from the population's contrast.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03_pseudobulk_de.py` | `main` | `thresholds.pseudobulk_min_cells=20` | `03_results/objects/02_annotation.h5ad` |

## tables/de_summary.csv

Every sorted population clears the bar for a powered enrichment test. The paired `~ donor + tissue`
model returns 1,797 significant SF-versus-PB genes of 13,999 ranked in Treg, 1,949 of 14,411 in
Tcon, and 1,695 of 14,014 in CD8. No arm was dropped for lack of donors.

**How to read:** One row per sorted population. `n_sf` and `n_pb` count donor strata in each tissue
arm; the unit is the stratum. `model` is the fitted design, reading `skipped` if an arm fell under the
donor floor; `n_sig_de` counts genes passing both the FDR and the |log2FC| gate; `n_ranked` is
the length of the sign-preserving ranked list handed to pre-ranked GSEA. This is a
diagnostic/powering table. The confirmatory statistics are the donor-pseudobulk NES values computed
downstream from `ranked_{treg,tcon,cd8}.tsv`.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03_pseudobulk_de.py` | `main` | `thresholds.de_fdr=0.05; de_logfc=1.0; pseudobulk_min_donors=3` | `03_results/objects/02_annotation.h5ad` |

## tables/pseudobulk_counts.csv

Summing raw UMI counts within donor × tissue × label gives 39 pseudobulk strata across 21,740
genes, every stratum above the per-stratum cell floor.

**How to read:** Rows are strata (`<donor>_<tissue>_<label>`), columns Ensembl gene ids, values
summed raw integer UMIs. Values are left unnormalised because the DE engine expects raw counts and
performs its own library-size normalisation. Pair with `pseudobulk_coldata.csv` for the design
factors. This is the matrix handed across the Python-to-R seam.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03a_pseudobulk_export.py` | `main` | `pseudobulk_min_cells` | `03_results/objects/02_annotation.h5ad` (`counts` layer) |

## tables/gene_symbols.csv

All 21,740 exported Ensembl ids carry a gene symbol, and all 21,740 symbols are distinct, so the
Ensembl-to-symbol rename loses no gene and creates no collision in this dataset.

**How to read:** Two columns, `ensembl_id` and `gene_symbol`, in the same order as the columns of
`pseudobulk_counts.csv`. This file exists because the counts matrix is keyed by Ensembl id while
every reference gene set the ranked lists are tested against — the mouse-projection signatures,
MSigDB Hallmark, the curated HSR lens — matches on HGNC symbol. The R stage asserts this file is
present before ranking. Skipping it yields ranked lists keyed by Ensembl id, which intersect the
reference sets at approximately zero; enrichment tools report that as an empty result, so the
failure is silent.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03a_pseudobulk_export.py` | `main` | — | `03_results/objects/02_annotation.h5ad` (`var['gene_symbol']`) |

## tables/de_SFvsPB_{treg,tcon,cd8}.csv

Paired SF-versus-PB differential expression recovers 1,797 / 1,949 / 1,695 genes at FDR < 0.05 and
|log2FC| ≥ 1 in Treg / Tcon / CD8, all three fitted on the donor-paired model.

**How to read:** One file per sorted population, one row per gene, sorted by p-value. `stat` is the
signed moderated t-statistic and `log2FoldChange` the log2 fold change, both positive when the gene
is higher in synovial fluid. `padj` is BH across genes within that population. `avg_expr` is limma's
average log2-CPM, named that way because `baseMean` denotes a differently scaled quantity. `model`
records the fitted design, `n_paired_donors` how many donors appeared in both tissues, and
`de_engine` the engine. Primary donor-pseudobulk tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03b_pseudobulk_de.R` | `(top-level)` | `design.tissue_key`, `design.donor_key`, `thresholds.de_fdr` | `03_results/03_pseudobulk/tables/pseudobulk_counts.csv`, `pseudobulk_coldata.csv`, `gene_symbols.csv` |

## tables/ranked_{treg,tcon,cd8}.tsv

The signed ranked lists carry 13,999 / 14,411 / 14,014 genes for Treg / Tcon / CD8 after expression
filtering, and are the single input every downstream enrichment stage reads.

**How to read:** Two columns, no header: HGNC symbol and the signed moderated t-statistic, sorted
descending, so the top of the file is most synovial-fluid-up and the bottom most blood-up. One row
per symbol; where several Ensembl ids share a symbol the most extreme |t| is kept, because a
duplicated symbol corrupts a pre-ranked enrichment run. Gene count differs per population because
`filterByExpr` is applied within each population's own design.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03b_pseudobulk_de.R` | `(top-level)` | `design.tissue_key`, `design.donor_key` | `03_results/03_pseudobulk/tables/pseudobulk_counts.csv`, `pseudobulk_coldata.csv`, `gene_symbols.csv` |

## tables/de_engine_migration.csv

The limma-voom migration validation reported in commit `799e4fe` — Spearman 0.994-0.997 rank
correlation and 87-89% top-500 overlap against the pre-migration ranked lists — exists only in that
commit message. It cannot be recomputed from tracked artifacts in this repository.

**How to read:** A provenance table, carrying no new validation. `not_reproducible` means the
reported value survives only in the commit message and the pre-migration ranked lists it was
computed against are untracked. The current `ranked_*.tsv` files exist on disk and are untracked
too, so they are insufficient on their own to reconstruct the before/after comparison. `command`
records the checks run to establish each row.

| Script | Function | Config | Input |
|---|---|---|---|
| repository audit | `git show`, `git ls-files`, `find` | `commit=799e4fe` | commit message and tracked file index |
