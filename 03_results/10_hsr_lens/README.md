# 10_hsr_lens -- artifact captions

_Abbreviations: SF = synovial fluid, PB = peripheral blood, NES = normalized enrichment score, HSR = heat-shock response._

### Where this lens came from

The HSR lens is anchor-independent and MSigDB-derived: the union of three human MSigDB v2026.1.Hs sets, pulled from the offline `msigdbr` 26.1.0 package by exact `gs_name` with validated sizes — `REACTOME_CELLULAR_RESPONSE_TO_HEAT_STRESS` (101), `REACTOME_REGULATION_OF_HSF1_MEDIATED_HEAT_SHOCK_RESPONSE` (82), `GOBP_RESPONSE_TO_HEAT` (104). The union is 176 genes and is carried here as `HSR_sensitivity`.

A per-gene functional taxonomy — external large-context-model classification, human-curated, and committed as non-deterministically-regenerable — splits those 176 into `hsf1_core_hsr` (45), `co_chaperone` (11), `generic_stress` (72), `npc_transport` (30), `thermosensory` (10) and `upr_er` (8). `HSR_core` (56) is `hsf1_core_hsr` plus `co_chaperone`: an anchor-independent curated HSF1/co-chaperone lens, and the term the figure below uses. It is named for the curated categories it was built from, not for a mechanism it is taken to demonstrate.

Three candidate sets were deliberately kept out of the union. `GOBP_DETECTION_OF_TEMPERATURE_STIMULUS` and `GOBP_DETECTION_OF_TEMPERATURE_STIMULUS_INVOLVED_IN_THERMOCEPTION` are thermosensory-neuron programs with no bearing on T cells, and `HP_FEVER` is a mutation-etiology panel rather than a transcriptional signature.

| Asset | Source | n | Frozen by |
|---|---|---|---|
| `HSR_sensitivity` | union of three MSigDB v2026.1.Hs sets (Reactome ×2, GO:BP) via `msigdbr` 26.1.0 | 176 | `02_analysis/scripts/freeze_hsr_lens.R` |
| `HSR_core` | taxonomy categories `hsf1_core_hsr` + `co_chaperone` | 56 | `02_analysis/scripts/freeze_hsr_lens.R` |
| `WT_heat_up` | mouse-anchor 39 °C up-arm, projected to human orthologs | 199 | frozen signature contract |

So the answer to "is this circular?" is no: `HSR_core` shares two genes with the 199-gene `WT_heat_up` (HSPA1A, HSPH1, Jaccard 0.008, tallied in `tables/hsr_wtheatup_overlap.csv`). The lens is a separate probe rather than a restatement of the anchor, and that independence is the reason for carrying it.

This compartment never re-pulls MSigDB. `02_analysis/scripts/freeze_hsr_lens.R` freezes the byte-identical lens from the mouse anchor's `temp_hsr_human_lens.rds`, so the JIA lists and the anchor lists are the same genes. Upstream curation detail lives with that RDS in `mouse_anchor/00_data/references/gene_sets/temp_hsr_lens/`.

### What the lens returns

Honest ceiling: even the clean HSR core is proteotoxic-stress-general, not fever-specific. Only the mouse anchor's experimental 37/39 contrast can measure thermal-ness. In JIA, this lens is carried and read correlatively; it does not decompose temperature causality from human scRNA-seq. The HSR NES is annotation-tier only and is firewalled from the confirmatory `WT_heat` effect-size spine; it is not written to `effect_sizes_treg_arthritis.csv`.

The result is selective in sign and short of significance: HSR core points toward synovial fluid in Treg and away from it in Tcon and CD8, with Treg at FDR 0.064. The figure prints every FDR on its face, so the sign pattern and its weakness are read together and no glyph implies a significance the numbers do not carry.

## tables/per_cell_hsr_scores.csv

The curated HSR core and sensitivity lenses give per-cell AUCell/UCell readouts that can be compared with `WT_heat_up` without selecting cells on either score.

**How to read:** One row per cell barcode. Metadata columns identify donor, tissue, and frozen coarse label; score columns are `HSR_core_AUCell`, `HSR_core_UCell`, `HSR_sensitivity_AUCell`, and `HSR_sensitivity_UCell`. AUCell/UCell are unsigned rank-based scores; higher values mean the set is more represented in that cell's expression ranking. Secondary annotation tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/10_hsr_lens.py` | `per_cell_hsr_scores` | `percell_score_ncores=8` | `03_results/objects/02_annotation.h5ad`, `00_data/references/temp_hsr_lens/{HSR_core,HSR_sensitivity}.txt` |

## tables/hsr_gsea_treg.csv

The Treg SF-vs-PB ranked list is tested against the curated HSR lens to ask whether a clean proteostasis signal survives in the same compartment where `WT_heat_up` enriched.

**How to read:** This is the fgsea output for Treg. Positive NES means enrichment toward SF-up genes; FDR is `padj`; `core_enrichment` lists the leading-edge genes. Annotation-tier pseudobulk enrichment only, not a confirmatory effect-size row.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/10_hsr_lens.py` | `run_fgsea` | `gsea_min_size=5; gsea_max_size=500; gsea_seed=123; gsea_nperm=100000` | `03_results/03_pseudobulk/tables/ranked_treg.tsv`, `00_data/references/temp_hsr_lens/{HSR_core,HSR_sensitivity}.txt` |

## tables/hsr_gsea_tcon.csv

The Tcon SF-vs-PB ranked list is tested against the curated HSR lens to determine whether any proteostasis enrichment is Treg-specific or pan-T-cell.

**How to read:** Positive NES means enrichment toward SF-up genes; FDR is `padj`; `core_enrichment` lists the leading-edge genes. Annotation-tier pseudobulk enrichment only, not a confirmatory effect-size row.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/10_hsr_lens.py` | `run_fgsea` | `gsea_min_size=5; gsea_max_size=500; gsea_seed=123; gsea_nperm=100000` | `03_results/03_pseudobulk/tables/ranked_tcon.tsv`, `00_data/references/temp_hsr_lens/{HSR_core,HSR_sensitivity}.txt` |

## tables/hsr_gsea_cd8.csv

The CD8 SF-vs-PB ranked list is tested against the curated HSR lens as a second pan-T-cell specificity control.

**How to read:** Positive NES means enrichment toward SF-up genes; FDR is `padj`; `core_enrichment` lists the leading-edge genes. Annotation-tier pseudobulk enrichment only, not a confirmatory effect-size row.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/10_hsr_lens.py` | `run_fgsea` | `gsea_min_size=5; gsea_max_size=500; gsea_seed=123; gsea_nperm=100000` | `03_results/03_pseudobulk/tables/ranked_cd8.tsv`, `00_data/references/temp_hsr_lens/{HSR_core,HSR_sensitivity}.txt` |

## tables/hsr_lens_nes.csv

The HSR NES summary compares the clean HSR core and broader sensitivity lens across Treg, Tcon, and CD8, revealing whether the empirical `WT_heat_up` SF-vs-PB separation is matched by a proteostasis program.

**How to read:** One row per population and HSR term. Positive `nes` means the HSR term is enriched toward SF-up genes in the SF-vs-PB ranked list; `padj` is fgsea FDR; `leading_edge` is semicolon-delimited. `evidence_tier=secondary_annotation` means this is annotation-tier only and is not part of the confirmatory effect-size spine.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/10_hsr_lens.py` | `hsr_nes` | `gsea_min_size=5; gsea_max_size=500; evidence_tier=secondary_annotation` | `03_results/10_hsr_lens/tables/hsr_gsea_{treg,tcon,cd8}.csv` |

## tables/hsr_colocalization.csv

The within-SF colocalization table tests whether `WT_heat_up`-high cells are the same cells as HSR-high cells, separately for HSR_core and HSR_sensitivity.

**How to read:** Rows are stratified by population, HSR term, correlation level, and method. `level=cell` uses individual SF cells; `level=donor_sf_mean` correlates donor-level SF means. Positive `r` means higher `WT_heat_up_AUCell` tends to coincide with higher HSR AUCell. Low cell-level `r` means the empirical lens and curated HSR lens label different cells, consistent with activation carrying much of `WT_heat_up`; high `r` means the lenses co-localize. Secondary per-cell tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/10_hsr_lens.py` | `hsr_colocalization` | `tissue_levels.synovial_fluid=synovial_fluid; evidence_tier=secondary_percell` | `03_results/10_hsr_lens/tables/per_cell_hsr_scores.csv`, `03_results/05_scoring/tables/per_cell_scores.csv` |

## tables/hsr_wtheatup_overlap.csv

The WT_heat_up gene list overlaps only the portion of each HSR lens listed in `genes_intersect`, defining how much direct gene reuse could explain any agreement between the empirical and curated lenses.

**How to read:** `n_a` and `n_b` are the sizes of `WT_heat_up` and the HSR term; `n_intersect` and `jaccard` quantify direct overlap; `genes_intersect` lists shared HGNC symbols. This is a gene-list annotation check, not an enrichment statistic.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/10_hsr_lens.py` | `hsr_wtheatup_overlap` | `signature=WT_heat_up` | `../mouse_anchor/03_results/human_projection/signatures/WT_heat/WT_heat_up.txt`, `00_data/references/temp_hsr_lens/{HSR_core,HSR_sensitivity}.txt` |

## tables/source_hash_manifest.csv

The HSR-lens overlap uses a pinned mouse-anchor projection source.

**How to read:** The stage checks these SHA-256 values before loading the mouse `WT_heat` projection.
If the anchor files change, the JIA overlap table is not regenerated until the source change is
reviewed and the manifest is updated.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/10_hsr_lens.py` | `verify_source_hashes()` | pinned SHA-256 | `../mouse_anchor/03_results/human_projection/signatures/WT_heat/` |

## tables/hsr_gsea_{treg,tcon,cd8}.rds

The clusterProfiler `gseaResult` objects preserve the same HSR fgsea runs summarized in the CSV files, enabling exact running-sum reconstruction if a later display needs it.

**How to read:** One RDS is written per population by the fgsea helper. These are compute substrates, not separate statistics; interpret their NES/FDR through the sibling `hsr_gsea_{treg,tcon,cd8}.csv` and the summarized `hsr_lens_nes.csv`. Positive enrichment still means SF-up.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/helpers/fgsea_prerank.R` | `(top-level)` | `gsea_min_size=5; gsea_max_size=500; gsea_seed=123; gsea_nperm=100000` | `03_results/03_pseudobulk/tables/ranked_{treg,tcon,cd8}.tsv`, `03_results/10_hsr_lens/tables/_signatures_hsr/{HSR_core,HSR_sensitivity}.txt` |

## tables/runsum_interactive_hsr_gsea_{treg,tcon,cd8}_HSR_{core,sensitivity}.csv

The running-sum substrates show where each HSR term falls along each population's SF-vs-PB ranked list, matching the fgsea objects behind the NES values.

**How to read:** One row is one ranked gene. `running_es` is the weighted enrichment score trace, `hit` marks HSR genes, and `leading_edge` marks the genes driving the enrichment peak. Positive, left-shifted peaks correspond to SF-up enrichment. Annotation-tier only.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/helpers/fgsea_prerank.R` | `(top-level)` | `gsea_min_size=5; gsea_max_size=500; gsea_seed=123; gsea_nperm=100000` | `03_results/03_pseudobulk/tables/ranked_{treg,tcon,cd8}.tsv`, `03_results/10_hsr_lens/tables/_signatures_hsr/{HSR_core,HSR_sensitivity}.txt` |

## tables/_signatures_hsr/HSR_{core,sensitivity}.txt

The stage-local HSR signature copies make the fgsea command self-contained while preserving the frozen reference gene lists.

**How to read:** Each file is one sorted HGNC symbol per line. These are not results; they are the exact gene-set inputs passed to fgsea for this stage.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/10_hsr_lens.py` | `prepare_hsr_signature_dir` | `signature_source=00_data/references/temp_hsr_lens` | `00_data/references/temp_hsr_lens/{HSR_core,HSR_sensitivity}.txt` |

## tables/_overview/hsr_core_running_sum.csv

The annotated numbers behind the HSR_core running-sum figure: each population's NES, FDR, testable set size, and the length of the ranked list its trace walks.

**How to read:** One row per population. `nes`/`padj` are the values printed in the figure legend; `set_size` is how many HSR_core genes survive intersection with that ranked list and `n_ranked_genes` its length, which is why the three traces end at slightly different x. The traces themselves are the `runsum_interactive_*` tables. Secondary annotation tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/10_hsr_lens_viz.py` | `running_sum_traces` | `evidence_tier=secondary_annotation` | `03_results/10_hsr_lens/tables/runsum_interactive_hsr_gsea_{treg,tcon,cd8}_HSR_core.csv`, `03_results/10_hsr_lens/tables/hsr_lens_nes.csv` |

## figures/_overview/hsr_core_running_sum.png

Walking each population's ranked list, HSR core accumulates a positive
peak near the synovial-fluid end in Treg while Tcon and CD8 run
negative throughout, so the sign selectivity is a property of the
rankings and not an artefact of the summary statistic.

**How to read:** Top panel: the weighted running enrichment score as each ranked list
is walked from synovial-fluid-up (left) to blood-up (right); a
positive, left-shifted excursion is synovial-fluid enrichment, a
negative trace the opposite. Bottom panel: where each population's HSR
core genes sit in its ranking, in matching colour. Legend labels carry
each NES and FDR, so read the Treg trace as a trend at FDR 0.064, not
a significant enrichment. Ranked-list lengths differ slightly, so
compare shapes rather than x positions; the y range is data-driven
because all three curves share one axis. Annotation tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/10_hsr_lens_viz.py` | `plot_running_sum` | `figures.running_sum_heights=[2.4, 0.7]; thresholds.gsea_fdr=0.05; evidence_tier=secondary_annotation` | `03_results/10_hsr_lens/tables/runsum_interactive_hsr_gsea_{treg,tcon,cd8}_HSR_core.csv, 03_results/10_hsr_lens/tables/hsr_lens_nes.csv` |
