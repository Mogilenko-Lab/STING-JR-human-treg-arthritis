# 15_coresh_search — Which public datasets the JIA synovial arm co-varies in

The mouse anchor asked a co-regulation search which public neighbourhoods its 39 °C contrast
sits in. This stage asks the same of the human joint: across roughly 42,000 public human GEO
datasets, in which ones do the genes that rise in the inflamed JIA synovial niche co-vary, and
what else moves with them there?

**What was computed.** Six queries — the up arm of the frozen limma-voom synovial-versus-blood
contrast, for each of the three sorted populations, at two stringency gates — swept against the
human half of the CoReSh preprocessed corpus. The top five datasets per query were projected
back into gene-level loadings to give co-regulation modules, and those modules were scored back
on the same ranked lists. Each recovered dataset was then annotated from the compendium's own
sample-metadata terms.

**What was drawn.** Two panels: the ranked datasets per query, and the modules scored against
every population.

**Everything here is exploratory, and the reason is structural.** A module seeded from a query
contains that query's genes by construction, so scoring it back on the seeding ranked list is
partly guaranteed to succeed. The stage measures that circularity and publishes the measurement.
No row reaches [`../master/`](../master/).

## What was searched

The compendium is the human (`hsa`) half of the CoReSh corpus, consumed read-only from the
shared reference cache: 89 chunks, snapshot `syn66227307_20260721`, **42,465 distinct human
datasets scored per query** (`tables/coresh_query_provenance.csv`). This compartment is human, so
the human chunks are the correct half. The mouse chunks would return near-zero overlap with no
error raised, which is why the species is asserted in config and again in the script.

The gates are re-derived from `../03_pseudobulk/tables/de_SFvsPB_<population>.csv`, the same fit
that produced the ranked lists, and the script refuses to run if any gated symbol is missing
from the ranking. `fdr_only` takes padj < 0.05, and `fdr_logfc` adds log2FC ≥ 1.

### The identifier seam, and what it cost

Every ranked list and reference set in this project is keyed by HGNC symbol. The compendium
chunks are keyed by integer Entrez id. In R, matching a character vector against an integer
vector returns all-NA without warning, so a symbol-keyed query would produce near-zero overlap,
noise-level scores and no error message. The mapping is explicit and its cost is published.

| Query | Symbols in | Mapped to Entrez | Dropped | Loss |
|---|---|---|---|---|
| Treg, FDR only | 2,344 | 2,108 | 236 | 10.1% |
| Treg, FDR + log2FC | 860 | 749 | 111 | 12.9% |
| Tcon, FDR only | 3,050 | 2,750 | 300 | 9.8% |
| Tcon, FDR + log2FC | 1,213 | 1,057 | 156 | 12.9% |
| CD8, FDR only | 2,889 | 2,575 | 314 | 10.9% |
| CD8, FDR + log2FC | 961 | 839 | 122 | 12.7% |

No symbol collapsed onto an Entrez id another had claimed, so nothing was double-counted. The
dropped tenth is dominated by clone-based accession-style names — unnamed loci and lncRNAs with
no Entrez entry. That matters for the top of the Treg list, whose single strongest gene,
`AC017002.1`, is one of them. The stringent gate loses slightly more, because a large log2FC is
more common among these lowly-annotated loci.

Before sweeping, the script probes one chunk and continues only if the query lands in the
compendium. A median of 411 (stringent) to 1,883 (relaxed) query ids are measured per dataset,
52–73% of each query. That check separates "there is nothing here" from "the query never
arrived", which are otherwise indistinguishable in the output.

## What came back

The search scores every dataset, so the complete ranking is 265,518 rows.
`tables/coresh_hits.csv` carries the ranked head of each query at 200 rows, and the full ranking
stays in the regenerable checkpoint. This is a stated cap — raise `coresh.hits_export_n` to
widen it.

The top five datasets per query gave 30 candidate modules. Nine duplicated an earlier module at
Jaccard > 0.8 and were dropped, leaving **21 modules** of 48–50 genes each, of which a median of
46 are present in each population's ranked list (`tables/coresh_module_sizes.csv`). The mouse
stage produced 14 modules in the same size band from 20 hits, with best pctVar running 5–26%
against 7–21% here, so the two searches operate at the same scale.

A module is named for how it was made: `CORESH_<population>_up_<gate>_<GSE>`. It is a
co-regulation neighbourhood mined from one public dataset's variance structure.

## The circularity, measured

Every module carries `frac_seed_genes`, the share of it that came from the seeding query, and
every enrichment row carries `seeded_from_this_population`.

The measurement is unflattering, which is the useful part. In Treg the 21 modules split both
ways — 10 up and 4 down at FDR < 0.05, NES −2.64 to +2.71 — and **which way a module goes is
largely predicted by how much of it is seed**: Spearman ρ = 0.70 in Treg, 0.88 in Tcon, 0.82 in
CD8, on a median seed fraction of 30% (`tables/coresh_gsea_summary.csv`). Mostly-seed modules
enrich positively, and mostly-recruited ones enrich negatively.

The sign of a derived-module NES therefore reports seed content more than shared biology. That is
why the whole folder sits in the exploratory tier.

## What the recovered datasets are

The human chunks carry a `wordMatrix`, the compendium's own centred per-sample matrix over the
terms that vary most across a dataset's GEO sample metadata. Correlating each term against the
query direction says which metadata term tracks the axis the query defines there. The annotation
is derived in-cache, named for that derivation, and descriptive only.

The two gates recovered visibly different kinds of dataset. The stringent gate surfaced
recognisable T-cell contexts — anti-CD28 / Th17 / Th0 stimulation (GSE110097, top-ranked for both
Treg and CD8), melanoma tumour-infiltrating CD8 (GSE153098), feeder-cell expansion (GSE222740),
activated lymphocytes (GSE141645), neutrophil treatment (GSE173807), lung explants (GSE102751).
The relaxed gate's top hits vary on **processing terms**: `umap`, `normalizedcounts`,
`cellbarcode`, `barcode-derived`, `multiseq`, `exome`, `rrna`.

A ~2,500-gene query spans enough of the transcriptome to align with whatever dominates a
dataset's variance, and in single-cell series that is often technical structure. Those same
relaxed-gate modules carry the lowest seed fractions (4–16%) and the strongly negative NES. Read
the relaxed gate as having recovered dominant-variance axes.

One bookkeeping detail: the compendium is keyed by (GSE, platform), so a series measured on two
platforms appears twice. GSE86566 is such a case. The annotation joins on the platform the
ranked hit scored, and `n_platform_records_in_compendium` makes the multiplicity visible.

---

## Figures

### `figures/_overview/coresh_pctvar_top_hits.png`

**The ranked public datasets, one panel per query.**
Six panels, one per sorted population × gate. Each bar is one public human GEO dataset, and
length is pctVar — the share of that dataset's variance the query genes jointly explain. pctVar
is unsigned and always positive, higher meaning tighter co-movement. The row label gives the
accession and how many query Entrez ids that dataset measures, because pctVar is normalised by
that count. Orange bars became modules and grey were ranked only. Panels carry independent x
ranges.

Across the 42,465 datasets searched, the JIA synovial up arm reaches a best score of pctVar
15.6% under the stringent gate (GSE102751) and 20.7% under the relaxed one (GSE118383). The
relaxed gate scores higher throughout on a query roughly three times larger, so pctVar compares
within a gate and the two rankings read as two separate searches.
*Source* `tables/_overview/coresh_pctvar_top_hits.csv` ·
`02_analysis/scripts/15_coresh_search_viz.R`.

### `figures/_overview/coresh_module_nes.png`

**The 21 modules scored against all three populations.**
Each row is one co-regulation module — the genes loading most strongly onto the query direction
inside one public dataset. The row label gives that dataset's accession, the population and gate
whose up arm seeded it, and the compendium's own metadata terms tracking the query axis there.
Columns are the three populations, each scored on its own ranked list. Fill is NES, orange
positive and blue negative, clamped. Size is −log10 FDR. A black outline marks the circular
cells, where a module is scored on the list that seeded it, and grey outlines are the
informative comparison.

In the Treg ranking the modules split both ways (NES −2.64 to +2.71, 10 up and 4 down at FDR <
0.05), and direction tracks seed content at Spearman ρ = 0.70 on a median seed fraction of 30%.
The enrichment reports seed content back to itself as much as it reports public biology, which
is why this tier supports no claim.
*Source* `tables/_overview/coresh_module_nes.csv` ·
`02_analysis/scripts/15_coresh_search_viz.R`.

---

## Tables

### `tables/coresh_query_provenance.csv` — the funnel

One row per query. `n_query_symbols` → `n_mapped_to_entrez` → `n_unique_entrez` is the
identifier funnel, with `n_dropped_unmapped` and `frac_dropped_unmapped` the loss and
`unmapped_examples` showing what kind of gene was lost. `gate_rule` states the exact filter and
`gate_source_table` where it was applied. `probe_median_overlap` is the pre-flight check below
which the script refuses to sweep. `n_datasets_searched` is the whole compendium and
`n_hits_exported_to_csv` the ranked head, with the full ranking at `full_ranking_checkpoint`.
The remaining columns pin the cache snapshot and the search parameters.

### `tables/coresh_hits.csv` — the ranked head

One row per (query, dataset), capped at 200 rows per query. `pctVar` is the co-regulation score,
comparable within a query alone because query size changes its scale. `size` is the matched
query count that normalises it, and `rank` the within-query position. `pval` is empty by design:
the search runs variance-only, because the permutation p-value costs minutes per query and leaves
the set of datasets read unchanged.

### `tables/coresh_derived_sets.csv` · `tables/coresh_derived_sets.gmt` — the modules

Long form, one row per (module, gene). `is_seed_gene` marks whether that gene was in the seeding
query, which is the column separating recruitment from tautology. `n_genes` is the module size
after the Entrez-to-symbol round trip and the [15, 500] size filter.

The `.gmt` carries the same 21 modules for any tool that wants them: one line per module, name,
a `-` placeholder, then tab-separated symbols. Treat the sets as unsigned — the projection
direction carries an arbitrary overall sign, so "up" in a name refers to the seeding query's arm.

### `tables/coresh_derived_gsea.csv` — the modules scored back

One row per (module, population), on the same fgsea engine and parameters the confirmatory
scoring stage uses, so the numbers sit on one scale with it. Two columns decide how to read a
row: `seeded_from_this_population` is TRUE when the module was mined using this very ranked list,
and `frac_seed_genes` gives the share of the module that is seed. `pctVar` and `hit_rank` carry
the module back to its source dataset. Positive NES means the module concentrates at the
synovial-high end.

### `tables/coresh_gsea_summary.csv` — the circularity headline

One row per population. `spearman_nes_vs_frac_seed` is the rank correlation across the 21
modules between how strongly a module enriches and how much of it was the seeding query. A value
near 1 means the enrichment is close to tautological. `n_seeded_here` counts modules mined from
this population's own ranking, `n_sig_up` and `n_sig_down` split the significant ones, and
`top_module` / `bottom_module` name the extremes.

### `tables/coresh_module_sizes.csv` — the coverage check

One row per (module, population). `n_in_ranked` and `frac_in_ranked` give how much of the module
is present in that ranking. This separates a real null from the silent failure in which a module
intersecting at near-zero returns an empty fgsea result. The script stops if the median falls
below `gsea_min_size`. One module, from GSE86566, reaches 13 genes.

### `tables/coresh_derived_annotation.csv` — what each dataset is

One row per module. `metadata_terms` are the compendium's most-variable GEO sample-metadata
terms for that dataset, ordered by how strongly each tracks the query direction.
`terms_aligned_with_query` and `terms_opposed_to_query` split them by sign, and `max_abs_term_r`
gives the strongest association. This describes the public dataset, derived from the cache.
`n_platform_records_in_compendium` above 1 flags a series held under two platforms.

### `tables/runsum_interactive_coresh_<population>_<module>.csv`

Fifteen files — the top five modules by |NES| in each population — carrying the gene-by-gene walk
in the shared running-sum schema: `rank`, `stat`, `running_es`, `hit`, `leading_edge`, plus the
set, population and contrast keys. Columns match the mouse-signature substrate exactly, which is
what lets a co-regulation module be plotted on the same axes as a mouse-derived arm.
