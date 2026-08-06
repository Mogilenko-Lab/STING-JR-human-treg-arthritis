# 03_pseudobulk — The donor-paired synovial-versus-blood contrast

This is the confirmatory spine of the compartment. Cells are summed to one profile per donor ×
tissue × frozen sort label, and the contrast is fitted on those profiles, so each donor casts one
vote. Everything downstream that supports a claim reads the rankings this stage writes.

**The work splits across an explicit file seam.** Aggregation is Python
(`03a_pseudobulk_export.py`) and writes plain CSVs of raw integer counts, running no statistics.
The model is R (`03b_pseudobulk_de.R`): `filterByExpr` → TMM → `voom` → `lmFit` →
`eBayes(robust=TRUE)`, on `~ donor + tissue`.

**The seam carries a gene map, and that is load-bearing.** The counts matrix is keyed by Ensembl
id while every reference gene set matches on HGNC symbol. `gene_symbols.csv` crosses that gap and
the R script asserts its presence before ranking. Ranked lists keyed by Ensembl id would intersect
every reference set at approximately zero, and enrichment tools report that as an empty result, so
the failure would be silent.

**What the contrast returns.** Synovial-fluid Tregs carry a reproducible program against the same
donors' blood: 1,797 genes clear FDR ≤ 0.05 and |log2FC| ≥ 1, 860 up in the joint and 937 down,
with 4,764 clearing the FDR cut alone. Tcon returns 1,949 and CD8 1,695 at the same gates, so
every population carries a ranked list powered for signature enrichment.

Every downstream enrichment stage ranks on the **signed moderated t**. The DE tables carry an
engine-agnostic column schema (`log2FoldChange`, `stat`, `pvalue`, `padj`) plus a `de_engine`
column naming the engine, so a later engine change stays contained to this stage.

---

## Figures

### `figures/_overview/pseudobulk_pca.png`

**Pseudobulk profiles separate by tissue and by label.**
Each point is one donor × tissue × label profile, on log-CPM over the top 2,000 variable genes.
x, PC1; y, PC2. Circle gives synovial fluid and square paired blood; colour gives donor. Read for
tissue separation and for a single donor dominating an axis. Neither happens, which is what makes
donor pseudobulk well-posed here. Display transform only.
*Source* `tables/pseudobulk_counts.csv` · `02_analysis/scripts/03_pseudobulk_de_viz.py`.

### `figures/_overview/treg_volcano.png`

**The Treg contrast — the ranking every Treg-gate enrichment in this compartment is computed on.**
x, log2 fold change, synovial fluid over paired blood; y, raw p on a −log10 scale. Significance
is decided on FDR, and the raw-p axis keeps the per-gene resolution that −log10(FDR) collapses.
Colour gives four categories: neither cut, fold change only, FDR only, both. The dashed
horizontal rule is the raw p realising FDR ≤ 0.05 (p ≤ 0.017) and the vertical rules
|log2FC| ≥ 1.0. The legend arrows give the up and down split of genes clearing both.

1,797 genes clear both cuts, 860 up in the joint and 937 down. The ten named genes are the five
most significant per side. The neighbour table holds the top 500 by p with their category.
*Source* `tables/de_SFvsPB_treg.csv` · `02_analysis/scripts/03_pseudobulk_volcano_viz.R`.

### `figures/_overview/de_count_bar.png`

**Every population carries a powered ranking.**
One bar per population; y, the count of significant synovial-versus-blood DE genes. Treg 1,797,
Tcon 1,949, CD8 1,695. Read it to confirm each arm carries enough signal to rank for pre-ranked
enrichment.
*Source* `tables/de_summary.csv` · `02_analysis/scripts/03_pseudobulk_de_viz.py`.

---

## Tables

### `tables/de_SFvsPB_{treg,tcon,cd8}.csv` — the primary result

One file per sorted population, one row per gene, sorted by p-value. `stat` is the signed
moderated t-statistic and `log2FoldChange` the log2 fold change, both positive when the gene is
higher in synovial fluid. `padj` is BH across genes within that population. `avg_expr` is limma's
average log2-CPM, named that way because `baseMean` denotes a differently scaled quantity.
`model` records the fitted design, `n_paired_donors` how many donors appeared in both tissues, and
`de_engine` the engine.

Counts at FDR < 0.05 and |log2FC| ≥ 1: 1,797 Treg, 1,949 Tcon, 1,695 CD8.

### `tables/ranked_{treg,tcon,cd8}.tsv` — the single downstream input

Two columns, no header: HGNC symbol and the signed moderated t, sorted descending, so the top of
the file is most synovial-up and the bottom most blood-up. 13,999 / 14,411 / 14,014 genes for
Treg / Tcon / CD8 after expression filtering. The counts differ because `filterByExpr` is applied
within each population's own design.

One row per symbol. Where several Ensembl ids share a symbol the most extreme |t| is kept, because
a duplicated symbol corrupts a pre-ranked enrichment run.

### `tables/pseudobulk_counts.csv` · `tables/pseudobulk_coldata.csv` · `tables/gene_symbols.csv`

The three files that cross the Python-to-R seam.

`pseudobulk_counts.csv` — rows are strata (`<donor>_<tissue>_<label>`), columns Ensembl gene ids,
values summed raw integer UMIs. 39 strata over 21,740 genes. Values are left unnormalised because
the DE engine expects raw counts and performs its own library-size normalisation.

`pseudobulk_coldata.csv` — one row per column of the counts matrix, keyed
`donor|tissue|coarse_label`. `donor` is the blocking factor. `tissue` is the contrast factor,
with `synovial_fluid` the numerator and `peripheral_blood` the denominator. `n_cells` is the
cells summed into that column and reads as an aggregation weight. The 39 surviving strata hold
266 to 4,454 cells each. The three absences — synovial Treg in patient 5, blood Tcon and blood
CD8 in patient 3 — are why every population's contrast runs six donors in one tissue arm against
seven in the other.

`gene_symbols.csv` — two columns, `ensembl_id` and `gene_symbol`, in the same order as the counts
columns. All 21,740 ids carry a symbol and all 21,740 symbols are distinct, so the rename loses no
gene and creates no collision in this dataset.

### `tables/de_summary.csv`

One row per sorted population. `n_sf` and `n_pb` count donor strata in each tissue arm, the unit
being the stratum. `model` is the fitted design, reading `skipped` if an arm fell under the donor
floor. `n_sig_de` counts genes passing both gates, and `n_ranked` is the length of the ranked list
handed to pre-ranked enrichment. The paired `~ donor + tissue` model returns 1,797 of 13,999 in
Treg, 1,949 of 14,411 in Tcon and 1,695 of 14,014 in CD8. No arm was dropped for lack of donors.

This is a diagnostic and powering table. The confirmatory statistics are the enrichment scores
computed downstream from `ranked_*.tsv`.

### `tables/de_engine_migration.csv`

A provenance table carrying no new validation. The limma-voom migration figures reported in commit
`799e4fe` — Spearman 0.994–0.997 rank correlation and 87–89% top-500 overlap against the
pre-migration ranked lists — survive in that commit message alone.
`not_reproducible` records that the pre-migration ranked lists it was computed against are
untracked, and the current `ranked_*.tsv` files are untracked too, so the before-and-after
comparison cannot be reconstructed from this repository. `command` records the checks run to
establish each row.
