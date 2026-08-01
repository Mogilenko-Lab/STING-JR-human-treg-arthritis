# 15_coresh_search — artifact captions

_**Abbreviations:** SF = synovial fluid, PB = peripheral blood, GEO = Gene Expression Omnibus, NES = normalized enrichment score, FDR = BH-adjusted p-value, pctVar = the CoReSh co-regulation score (percent of a dataset's variance explained by the query genes)._

The mouse anchor asked a co-regulation search which public neighbourhoods its 39 °C contrast sits in. Nothing had asked the same of the human joint. Here I do: across roughly 44,000 public human GEO datasets, in which ones do the genes that rise in the inflamed JIA synovial niche co-vary, and what else moves with them there? Everything in this folder is exploratory. No row reaches `effect_sizes_treg_arthritis.csv` or any `03_results/master/` accumulator, and the reason is structural rather than cautious — see *The circularity, measured* below.

## What was searched, and with what

The compendium is the human (`hsa`) half of the CoReSh preprocessed corpus, consumed read-only from the shared reference cache: 89 chunks, snapshot `syn66227307_20260721` downloaded 2026-07-21, **42,465 distinct human datasets** actually scored per query. This compartment is human, so the human chunks are the only correct half; the mouse chunks would have returned near-zero overlap with no error at all, which is why the species is asserted in config and again in the script rather than assumed.

Six queries were run, mirroring how the mouse anchor built its own: the **up arm** of the frozen limma-voom SF-vs-PB contrast, for each of the three sorted populations, at two stringency gates — `fdr_only` (padj < 0.05) and `fdr_logfc` (padj < 0.05 and log2FC ≥ 1). The gates are re-derived from `03_results/03_pseudobulk/tables/de_SFvsPB_<population>.csv`, the same fit that produced `ranked_<population>.tsv` (13,999 symbols for Treg), and the script refuses to run if any gated symbol is missing from the ranked list. Down arms were not queried, again mirroring the anchor.

### The identifier seam, and what it cost

Every ranked list and every reference gene set in this project is keyed by HGNC symbol. The compendium chunks are keyed by **integer Entrez id**. In R, matching a character vector against an integer vector returns all-NA without warning, so a symbol-keyed query would have produced near-zero overlap, noise-level scores, and no error message — a silent failure dressed as a biological null. The mapping is therefore explicit, and its cost is published rather than buried:

| Query | Symbols in | Mapped to Entrez | Dropped | Loss |
|---|---|---|---|---|
| Treg, FDR only | 2,344 | 2,108 | 236 | 10.1% |
| Treg, FDR + log2FC | 860 | 749 | 111 | 12.9% |
| Tcon, FDR only | 3,050 | 2,750 | 300 | 9.8% |
| Tcon, FDR + log2FC | 1,213 | 1,057 | 156 | 12.9% |
| CD8, FDR only | 2,889 | 2,575 | 314 | 10.9% |
| CD8, FDR + log2FC | 961 | 839 | 122 | 12.7% |

No symbol collapsed onto an Entrez id already claimed by another, so nothing was double-counted. The dropped tenth is not random: it is dominated by clone-based accession-style names (`AC017002.1`, `AC003102.3`, …) — unnamed loci and lncRNAs with no Entrez entry. That matters for reading the top of the Treg list in particular, whose single strongest gene, `AC017002.1`, is one of them. The stringent gate loses slightly more than the relaxed one because a large log2FC is more common among these lowly-annotated loci.

Before sweeping anything, the script probes one chunk and refuses to continue unless the query genuinely lands in the compendium. It does: a median of 411 (stringent) to 1,883 (relaxed) query ids are measured per dataset, 52–73% of each query. That check is what separates "there is nothing here" from "the query never arrived", and the two are indistinguishable in the output otherwise.

## What came back

The search returns a score for every dataset, so the complete ranking is 265,518 rows. `coresh_hits.csv` carries the ranked **head** of each query, the part anyone reads, at 200 rows per query; the full ranking stays in the regenerable checkpoint `03_results/objects/coresh_hsa_ranked.rds` and both counts sit in the provenance table. This is a stated cap, not a quiet one — raise `coresh.hits_export_n` to widen it.

The top five datasets per query (30 hits) were projected back into gene-level loadings to give co-regulation modules. Nine were near-duplicates of an earlier module at Jaccard > 0.8 and were dropped, leaving **21 modules** of 48–50 genes each (median 50), of which a median of 46 are present in each population's ranked list. The mouse stage produced 14 modules of the same 48–50 size band from 20 hits, and its best pctVar ran 5–26% against 7–21% here, so the two searches are operating at the same scale — a reassuring cross-check that the human run was configured like its mouse counterpart rather than differently.

A module is named for how it was made: `CORESH_<population>_up_<gate>_<GSE>`. It is a **co-regulation neighbourhood mined from one public dataset's variance structure** — not a curated ontology term, not a pathway, and not a claim about mechanism.

## The circularity, measured

A module seeded from a query contains that query's genes by construction, so scoring it back on the seeding ranked list is partly guaranteed to succeed. Rather than caveat this in prose, the stage measures it. Every module carries `frac_seed_genes`, the fraction of it that is seed rather than newly recruited gene, and every enrichment row carries `seeded_from_this_population`.

The measurement is unflattering, which is the useful part. In Treg the 21 modules split both ways — 10 up and 4 down at FDR < 0.05, NES −2.64 to +2.71 — and **which way a module goes is largely predicted by how much of it is seed**: Spearman ρ = 0.70 in Treg, 0.88 in Tcon, 0.82 in CD8, on a median seed fraction of 30%. Mostly-seed modules enrich positively; mostly-recruited ones enrich negatively.

The sign of a derived-module NES is therefore closer to a readout of seed content than of shared biology. That is why the whole folder sits in the exploratory tier.

## What the recovered datasets are

The mouse stage answered this with frozen external web research, because the mouse chunks store only variance structure. The human chunks store more: each dataset object carries a `wordMatrix`, the compendium's own centred per-sample matrix over the terms that vary most across that dataset's GEO sample metadata. Correlating each term against the query direction says which metadata term tracks the axis the query defines there.

So the annotation is derived in-cache and named for it, and it is descriptive only, entering no statistic. It is **not** the mouse table's content: no researched title, tissue, perturbation, PubMed id or context class, because the cache cannot supply them.

Read that way, the two gates recovered visibly different kinds of dataset. The stringent gate surfaced recognisable T-cell contexts — anti-CD28 / Th17 / Th0 stimulation (GSE110097, top-ranked for both Treg and CD8), melanoma tumour-infiltrating CD8 (GSE153098), feeder-cell expansion (GSE222740), activated lymphocytes (GSE141645), neutrophil treatment (GSE173807), lung explants (GSE102751). The relaxed gate's top hits instead vary on **processing terms** — `umap`, `normalizedcounts`, `cellbarcode`, `barcode-derived`, `multiseq`, `exome`, `rrna`.

A ~2,500-gene query spans enough of the transcriptome to align with whatever dominates a dataset's variance, and in single-cell series that is often technical structure. Those same relaxed-gate modules carry the lowest seed fractions (4–16%) and the strongly negative NES. I read the relaxed gate as having recovered dominant-variance axes rather than biology, and would not build on it.

One bookkeeping detail worth knowing: the compendium is keyed by (GSE, platform), so a series measured on two platforms appears twice. GSE86566 is such a case here. The annotation is joined on the platform the ranked hit actually scored, and `n_platform_records_in_compendium` makes the multiplicity visible; the upstream loading extractor resolves a duplicated accession to the first record it finds, so for that one module the loadings and the ranking statistics could in principle come from different platform records.

---

## tables/coresh_query_provenance.csv

Records that each of the six queries lost 9.8–12.9% of its gated symbols at the HGNC-to-Entrez seam, with no duplicate-id collapse, while retaining a 52–73% median per-dataset overlap against the compendium — enough coverage that a weak result would be a result rather than a mis-keyed query.

**How to read:** One row per query. `n_query_symbols` → `n_mapped_to_entrez` → `n_unique_entrez` is the identifier funnel, and `n_dropped_unmapped` with `frac_dropped_unmapped` is the loss; `unmapped_examples` shows what kind of gene was lost. `gate_rule` states the exact filter applied to the DE table and `gate_source_table` where it was applied. `probe_median_overlap` is the pre-flight check: the median number of query ids actually measured per dataset in the probe chunk, below which the script refuses to sweep. `n_datasets_searched` is the whole compendium; `n_hits_exported_to_csv` is the ranked head written out, with the full ranking at `full_ranking_checkpoint`. The remaining columns pin the cache snapshot, the search parameters and the GSEA engine settings. Exploratory tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/15_coresh_search.R` | `build_query` | `coresh.species=human; coresh.chunks_subdir=hsa; de_fdr=0.05; de_logfc=1.0; coresh.min_query_size=3` | `03_results/03_pseudobulk/tables/de_SFvsPB_{treg,tcon,cd8}.csv + ranked_{treg,tcon,cd8}.tsv` |

## tables/coresh_hits.csv

The ranked head of the compendium search: the JIA synovial up-arm reaches pctVar 15.6% at best under the stringent gate and 20.7% under the relaxed one, so these genes do co-vary tightly somewhere in public human data rather than being scattered.

**How to read:** One row per (query, dataset). `pctVar` is the co-regulation score — the share of that dataset's variance the query genes jointly explain, normalised by how many of them the platform measures; it is unsigned, always positive, and comparable only within a query because query size changes its scale. `size` is that matched count, `rank` the within-query position. `pval` is empty by design: the search runs variance-only (`coresh.pvalues: false`), because the GESECA permutation p-value costs minutes per query and would not change which datasets are read. Capped at the top 200 rows per query; the full 265,518-row ranking is in the objects checkpoint. Exploratory tier: a ranking, not a test.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/15_coresh_search.R` | `coresh_batch` | `coresh.pvalues=false; coresh.n_cores=12; coresh.hits_export_n=200` | `$CORESH_CHUNKS/hsa/*_full_objects.qs2 (snapshot syn66227307_20260721)` |

## tables/coresh_derived_sets.csv

The 21 co-regulation modules the top hits define, each 48–50 genes, of which a median of only 30% are the query genes that seeded the module — so most of each module is genuinely newly recruited co-moving gene, even where the enrichment is not.

**How to read:** Long form, one row per (module, gene). `is_seed_gene` marks whether that gene was in the query that seeded the module, which is the column that separates recruitment from tautology. `n_genes` is the module size after the Entrez-to-symbol round trip and the [15, 500] size filter. Nine of the 30 top hits produced modules that duplicated an earlier one at Jaccard > 0.8 and are absent here. `coresh_derived_sets.gmt` carries the same modules in GMT form for any tool that wants them. Exploratory tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/15_coresh_search.R` | `build_coresh_gmt` | `coresh.top_n_hits=5; coresh.n_derive=50; coresh.derived_min_size=15; coresh.derived_max_size=500; coresh.jaccard=0.8` | `03_results/objects/coresh_hsa_ranked.rds + the hsa chunk cache` |

## tables/coresh_derived_sets.gmt

The same 21 modules in GMT form, so a downstream tool can read them as a custom gene-set database without re-deriving them from the compendium.

**How to read:** One line per module: name, a `-` description placeholder, then the HGNC symbols, tab-separated. Names follow `CORESH_<population>_up_<gate>_<GSE>` — the population and gate whose up arm seeded the search, and the public dataset the co-regulation was read from. Treat the sets as unsigned: the projection direction carries an arbitrary overall sign, so "up" in the name refers to the seeding query's arm, never to the module's own direction. Exploratory tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/15_coresh_search.R` | `build_coresh_gmt` | `coresh.n_derive=50; coresh.derived_min_size=15; coresh.derived_max_size=500; coresh.jaccard=0.8` | `03_results/objects/coresh_hsa_ranked.rds + the hsa chunk cache` |

## tables/coresh_derived_gsea.csv

Scoring the 21 modules back on the ranked lists splits them both ways in every population (10 up / 4 down in Treg, 12/5 in Tcon, 11/8 in CD8 at FDR < 0.05), and the direction tracks seed content rather than shared biology.

**How to read:** One row per (module, population). `nes`, `pvalue`, `padj`, `set_size` and `core_enrichment` are the standard pre-ranked GSEA fields on the same `clusterProfiler::GSEA(by = "fgsea")` engine and the same config parameters the confirmatory scoring stage uses, so the numbers sit on one scale with it. The columns that decide how to read a row are `seeded_from_this_population` — TRUE means the module was mined using this very ranked list and the row is circular — and `frac_seed_genes`, the share of the module that is seed. `pctVar` and `hit_rank` carry the module back to its source dataset. Positive NES means the module concentrates at the synovial-high end. Exploratory tier: nothing here is a claim.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/15_coresh_search.R` | `clusterProfiler::GSEA` | `gsea_min_size=5; gsea_max_size=500; gsea_seed=123; gsea_nperm=100000; gsea_fdr=0.05; engine=clusterProfiler::GSEA(by=fgsea)` | `03_results/03_pseudobulk/tables/ranked_{treg,tcon,cd8}.tsv + 03_results/objects/coresh_hsa_derived_sets.rds` |

## tables/coresh_gsea_summary.csv

The one line per population that decides how much of this stage is circular: NES correlates with a module's seed fraction at Spearman 0.70 (Treg), 0.88 (Tcon) and 0.82 (CD8), so the enrichment largely reports seed content back to itself.

**How to read:** One row per population. `spearman_nes_vs_frac_seed` is the headline: the rank correlation across the 21 modules between how strongly a module enriches and how much of it was the query that seeded it. A value near 1 means the enrichment is close to tautological; a value near 0 would mean the modules are carrying independent information. `n_seeded_here` counts the modules mined from this population's own ranked list. `n_sig_up` and `n_sig_down` split the FDR-significant modules by direction, and `top_module` / `bottom_module` name the extremes. Exploratory tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/15_coresh_search.R` | `GSEA_SUMMARY` | `gsea_fdr=0.05` | `03_results/15_coresh_search/tables/coresh_derived_gsea.csv` |

## tables/coresh_module_sizes.csv

A median of 46 of each module's 48–50 genes survive into each population's ranked list, so no enrichment in this folder is an artefact of a module failing to reach the data it was scored against.

**How to read:** One row per (module, population scored against). `n_genes` is the module size; `n_in_ranked` and `frac_in_ranked` are how much of it is actually present in that population's ranked list. This is the quantity that separates a real null from the documented silent failure, in which a module that intersects the ranked list at near-zero returns an empty fgsea result that looks like a biological null; the script stops if the median falls below `gsea_min_size`. One module (from GSE86566) reaches only 13 genes and should be read with that in mind. Exploratory tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/15_coresh_search.R` | `SIZES` | `gsea_min_size=5` | `03_results/objects/coresh_hsa_derived_sets.rds + ranked_{treg,tcon,cd8}.tsv` |

## tables/coresh_derived_annotation.csv

What each of the 21 recovered datasets is, read out of the compendium itself: the stringent gate recovered T-cell stimulation, tumour-infiltrate and tissue-explant series, while the relaxed gate's top hits vary mainly on processing terms such as `umap`, `cellbarcode` and `normalizedcounts`.

**How to read:** One row per module. `metadata_terms` are the compendium's own most-variable GEO sample-metadata terms for that dataset, ordered by how strongly each tracks the query direction; `terms_aligned_with_query` and `terms_opposed_to_query` split them by the sign of that correlation, and `max_abs_term_r` gives the strongest association. This is a descriptor of the public dataset, derived from the cache — not researched identity, and not an annotation of the JIA data; there is deliberately no title, tissue, perturbation or PubMed id, because the cache cannot supply them. `n_platform_records_in_compendium` above 1 flags a series held under two platforms, where the ranking and the loadings could resolve to different records. Descriptive only: nothing here enters a statistic. Exploratory tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/15_coresh_search.R` | `.index_chunks` | `coresh.n_annotation_terms=15; cache snapshot=syn66227307_20260721` | `$CORESH_CHUNKS/hsa/*_full_objects.qs2 (wordMatrix, E1024, rownames)` |

## tables/runsum_interactive_coresh_&lt;population&gt;_&lt;module&gt;.csv

The gene-by-gene walk behind the five strongest modules per population, in the schema the interactive running-sum curves already use, so a co-regulation module can be hovered on the same axes as a mouse-derived arm.

**How to read:** One row per gene per module, ordered by `rank` along that population's SF-vs-PB ranked list. `stat` is the signed moderated t, `running_es` the DOSE weighted-KS running enrichment score recomputed with the fitted exponent so the curve is exactly the one the GSEA object describes, `hit` marks module membership and `leading_edge` the core-enrichment subset. `gene_set`, `population` and `contrast` identify the curve. Columns are identical to the mouse-signature running-sum substrate, which is what lets the two be plotted together. Fifteen files: the top five modules by |NES| in each of the three populations. Exploratory tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/15_coresh_search.R` | `clusterProfiler::GSEA` | `figures.running_sum_top=5; gsea_min_size=5; gsea_max_size=500; gsea_seed=123` | `03_results/objects/coresh_hsa_gsea_{treg,tcon,cd8}.rds` |

## figures/_overview/coresh_pctvar_top_hits.png

Across the 42,465 public human datasets searched, the JIA
synovial-fluid-versus-blood up-arm reaches a best co-regulation score
of pctVar 15.6% under the stringent FDR-plus-log2FC gate (GSE102751)
and 20.7% under the relaxed FDR-only gate (GSE118383); the relaxed
gate scores higher throughout on a query roughly three times larger,
so pctVar is comparable only within a gate and the two rankings are
read as separate searches rather than as one league table.

**How to read:** One panel per query: a sorted population's up arm at one significance
gate. Each bar is one public human GEO dataset; length is pctVar, the
share of that dataset's variance the query genes jointly explain — a
co-regulation score, unsigned and always positive, higher meaning
tighter co-movement. The row label gives the accession and how many
query Entrez ids that dataset measures, because pctVar is normalised
by that count. Orange bars became modules; grey were ranked only.
Panels have independent x ranges: pctVar is not comparable across
query sizes, and the relaxed gate's query is far larger. Exploratory
tier: a ranking of public datasets, not a test.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/15_coresh_search_viz.R` | `create_pctvar_overview` | `coresh.top_n_hits=5; coresh.species=human; coresh.pvalues=FALSE; n_hits_shown=10` | `03_results/15_coresh_search/tables/coresh_hits.csv` |

## figures/_overview/coresh_module_nes.png

In the Treg synovial-fluid-versus-blood ranking the 21 modules split
both ways (NES -2.64 to 2.71; 10 up and 4 down at FDR < 0.05), and
which way a module goes tracks how much of it is the query that
seeded it rather than anything it newly recruited: Spearman rho =
0.70 between NES and seed fraction, on a median seed fraction of only
30%. The enrichment is therefore reporting seed content back to
itself as much as it is reporting public biology, which is why
nothing here is read as evidence.

**How to read:** Each row is one co-regulation module: the genes loading most strongly
onto the query direction inside one public human GEO dataset. The row
label gives that dataset's accession, the population and gate whose
up arm seeded it, and the compendium's own sample-metadata terms
tracking the query axis there — a descriptor of the public dataset,
not of the JIA data. Columns are the three sorted populations, each
scored on its own synovial-fluid-versus-blood ranked list. Fill is
NES, orange positive and blue negative, clamped. Size is -log10 FDR.
A black outline marks the circular cells, where a module is scored on
the list that seeded it; grey outlines are the informative
comparison. Exploratory tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/15_coresh_search_viz.R` | `create_module_nes_dotplot` | `gsea_min_size=5; gsea_max_size=500; gsea_seed=123; gsea_nperm=100000; nes_cap=3.5; engine=clusterProfiler::GSEA(by=fgsea)` | `03_results/15_coresh_search/tables/coresh_derived_gsea.csv + coresh_derived_annotation.csv + coresh_gsea_summary.csv` |

## tables/_overview/coresh_pctvar_top_hits.csv

Source table for the co-regulation ranking panel: the top ten public human datasets per query with their pctVar, matched query size and whether they were carried into a module.

**How to read:** Same-stem neighbour of `figures/_overview/coresh_pctvar_top_hits.png`. One row per plotted bar. `query_genes_in_dataset` is the matched count that normalises `pctVar`, and `became_module` distinguishes the orange bars from the grey. Exploratory tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/15_coresh_search_viz.R` | `prepare_pctvar_table` | `n_hits_shown=10` | `03_results/15_coresh_search/tables/coresh_hits.csv` |

## tables/_overview/coresh_module_nes.csv

Source table for the module dot plot: every module scored against every population, with the seed fraction and the circularity flag alongside the NES.

**How to read:** Same-stem neighbour of `figures/_overview/coresh_module_nes.png`. One row per plotted point. `seeded_from_this_population` marks the black-outlined circular points and `frac_seed_genes` quantifies how much of the enrichment is guaranteed; `terms_aligned_with_query` is the metadata descriptor shown in the row label. Exploratory tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/15_coresh_search_viz.R` | `prepare_module_table` | `nes_cap=3.5; gsea_fdr=0.05` | `03_results/15_coresh_search/tables/coresh_derived_gsea.csv + coresh_derived_annotation.csv` |

