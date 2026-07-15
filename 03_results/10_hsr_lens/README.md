# 10_hsr_lens -- artifact captions

_Abbreviations: SF = synovial fluid, PB = peripheral blood, NES = normalized enrichment score, HSR = heat-shock response._

Honest ceiling: even the clean HSR core is proteotoxic-stress-general, not fever-specific. Only the mouse anchor's experimental 37/39 contrast can measure thermal-ness. In JIA, this lens is carried and read correlatively; it does not decompose temperature causality from human scRNA-seq. The HSR NES is annotation-tier only and is firewalled from the confirmatory `WT_heat` effect-size spine; it is not written to `effect_sizes_treg_arthritis.csv`.

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

## tables/_overview/hsr_nes_by_population.csv

The plotted NES source table carries the same HSR_core and HSR_sensitivity enrichment values shown in the grouped population figure.

**How to read:** Rows are the plotted bars: population, HSR term, NES, p-value, FDR, and leading edge. Positive NES means enrichment toward SF-up genes. Secondary annotation tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/10_hsr_lens_viz.py` | `hsr_nes_by_population` | `evidence_tier=secondary_annotation` | `03_results/10_hsr_lens/tables/hsr_lens_nes.csv` |

## tables/_overview/hsr_wtheatup_colocalization.csv

The plotted colocalization source table carries the within-SF cell-level Spearman correlations between `WT_heat_up_AUCell` and `HSR_core_AUCell`.

**How to read:** Rows are the plotted bars by population. Low `r` means the empirical WT_heat_up lens and curated HSR_core lens mark mostly different SF cells; high `r` means they co-localize. Secondary per-cell tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/10_hsr_lens_viz.py` | `hsr_wtheatup_colocalization` | `level=cell; method=spearman; tissue=synovial_fluid` | `03_results/10_hsr_lens/tables/hsr_colocalization.csv` |

## figures/_overview/hsr_nes_by_population.png

The curated HSR lens tests whether the SF-vs-PB WT_heat separation is
matched by proteostasis enrichment across the three sorted T-cell
populations.

**How to read:** Bars are fgsea NES values for SF-vs-PB ranked lists. Positive NES
means the HSR term is enriched toward synovial-fluid-up genes. This is
annotation-tier and is not an effect-size claim.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/10_hsr_lens_viz.py` | `hsr_nes_by_population` | `gsea_min_size=5; gsea_max_size=500; evidence_tier=secondary_annotation` | `03_results/10_hsr_lens/tables/hsr_lens_nes.csv` |

## figures/_overview/hsr_wtheatup_colocalization.png

Within SF cells, the WT_heat_up-vs-HSR_core Spearman correlation asks
whether the empirical mouse lens and curated HSR lens mark the same
cells.

**How to read:** Each bar is the within-SF cell-level Spearman r between
WT_heat_up_AUCell and HSR_core_AUCell. Low r means WT_heat_up-high and
HSR-high cells largely differ; high r means the lenses co-localize.
Secondary per-cell tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/10_hsr_lens_viz.py` | `hsr_wtheatup_colocalization` | `level=cell; method=spearman; tissue=synovial_fluid; evidence_tier=secondary_percell` | `03_results/10_hsr_lens/tables/hsr_colocalization.csv` |
