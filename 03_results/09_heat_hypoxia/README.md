# 09_heat_hypoxia: artifact captions

_**Abbreviations:** SF = synovial fluid, PB = peripheral blood, NES = normalized enrichment score, FDR = BH-adjusted p-value._

This stage asks one bounded question of the JIA SF-vs-PB `WT_heat` enrichment: is it reducible to the set's own `HALLMARK_HYPOXIA`-overlap gene content? That is a membership question, and it is answered by deleting those genes and re-running the same donor-pseudobulk fgsea. It is deliberately not a question about temperature, and it is not a question about whether hypoxia is a confound or a co-exposure — those are not separable in cross-sectional human data, and nothing here licenses a statement about either. Hypoxia is a transcriptional readout throughout, never a HIF-causality claim.

The figures walk that in order. The first answers the membership question at confirmatory tier and stops there. The second corroborates and cannot answer whether the two per-cell scores mark the same cells. The reader sequence then moves directly to the curated whole-arm coverage panel in `11_heat_decomposition`, which alone answers what the set contains: 137 of 199 up-arm genes belong to no named program, the curated HSR core contains 2, and type-I interferon contains 1.

I withdrew the former leading-edge composition figure for redundancy and denominator risk, not because it re-tested selected subsets. It only described the genes at the synovial-fluid end of each ranking. Its model-assigned taxonomy (`agy_gemini_3.1_pro_2026-07-14`) put 55% to 61% of those selected leading-edge genes in immediate-early or effector/activation categories, a memorable fraction that is not comparable to the curated whole-arm activation counts of 5.6% in mouse and 6.0% after human projection. I retain `leadingedge_composition.csv` as an exploratory compute resource and notebook input, but no visualization of it remains in the published overview.

Two further figures put the mouse sets back into the plain Treg SF-vs-PB volcano, so a reader can see how much of the differential-expression response the signature accounts for, and how little it shares with the published IFN-independent STING-activation signature of de Cevins et al. 2023 (Cell Rep Med, PMID 38118407) that the STING positive-control compartment carries as a reference axis.

## gene_purge_nes_comparison.csv

Removing 18 `HALLMARK_HYPOXIA` overlap genes reduces the `WT_heat_up` NES modestly, but the enrichment remains SF-high in Treg, Tcon, and CD8 at FDR <= 4.1e-5.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/09_heat_hypoxia.py` | `gene_purge_nes` | `gsea_min_size=5`, `gsea_max_size=500`, `gsea_seed=123`, `gsea_nperm=100000` | `03_results/03_pseudobulk/tables/ranked_{treg,tcon,cd8}.tsv`, mouse `WT_heat_up/down.txt`, `00_data/references/msigdb_hallmark/HALLMARK_HYPOXIA.txt` |

**How to read:** One row per sorted population. Positive NES means the set is enriched toward the SF-high end of the SF-vs-PB ranked list. `NES_full` is the original `WT_heat_up` score, while `NES_purged` is the same fgsea engine after removing hypoxia-overlap genes. This is the primary donor-pseudobulk tier. I read a positive, significant purged NES as evidence that the `WT_heat_up` enrichment is not reducible to its HALLMARK_HYPOXIA-overlap gene content. It says nothing about temperature, and nothing about whether hypoxia is a confound or a co-exposure — those are not separable in cross-sectional human data.

## heat_hypoxia_colocalization.csv

Within SF cells, `WT_heat_up_AUCell` and `HALLMARK_HYPOXIA_AUCell` show weak positive cell-level correlations, while donor-level SF means are not positive in this small donor set.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/09_heat_hypoxia.py` | `heat_hypoxia_colocalization` | `tissue_levels.synovial_fluid=synovial_fluid`, `donor_key=donor` | `03_results/interactive/08_harvest_readout.parquet`, `03_results/05_scoring/tables/per_cell_scores.csv` |

**How to read:** Rows are stratified by population, level, and correlation method. `level=cell` uses SF cells directly. `level=donor_sf_mean` correlates per-donor SF mean heat and hypoxia scores. Positive `r` means a higher `WT_heat_up` score tends to sit with a higher HALLMARK_HYPOXIA score. This is an L3 secondary per-cell read and is not pooled with the pseudobulk NES. The cell-level correlation is weak (Spearman 0.08 to 0.20), so `WT_heat_up`-high and HALLMARK_HYPOXIA-high cells are largely distinct by this measure. The donor-level correlation rests on only 6 to 7 donors and is effectively unpowered, so its sign is not interpretable and must not be read as the two scores being anti-correlated.

## leadingedge_composition.csv

Of the 49 to 71 genes at the SF end of each population's ranking, the model-assigned taxonomy describes 55% to 61% as immediate-early or effector/activation. Those are fractions of leading-edge genes only, not of the 199-gene set, and I retain them as an exploratory compute resource rather than a reader-sequence figure.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/09_heat_hypoxia.py` | `leadingedge_composition` | taxonomy from `00_data/references/heat_leadingedge_taxonomy/` | `03_results/09_heat_hypoxia/tables/runsum_interactive_gsea_full_{treg,tcon,cd8}_WT_heat_up.csv`, `00_data/references/heat_leadingedge_taxonomy/leadingedge_gene_taxonomy.csv` |

**How to read:** One row per population. `n_leading_edge` is the count of fgsea core-enrichment genes from the full `WT_heat_up` run. Each is assigned to one program — `heat_shock_proteostasis`, `hypoxia_HIF`, `immediate_early_stress`, `effector_activation`, `other` — by the frozen model-assigned taxonomy whose exact source is carried in `taxonomy_source`; `n_unclassified` counts genes assigned to nothing. This table was never used to re-test a leading-edge subset. It remains because the reactive review notebook reads it, but its visualization is withdrawn: a fraction over genes selected because they enriched describes that leading edge and not the 199-gene set. Read whole-arm composition only from `11_heat_decomposition/figures/_overview/heatdecomp_arm_coverage`, where the curated versioned lenses leave 137 of 199 genes unassigned and contain 2 curated-HSR and 1 type-I-interferon gene.

## tables/gsea_full_{treg,tcon,cd8}.csv

With the complete 199-gene mouse `WT_heat_up` set, SF-vs-PB enrichment is positive and strong in every sorted population -- NES 2.59 (Treg), 2.68 (Tcon), 2.07 (CD8), all at FDR <= 3.6e-7 -- while `WT_heat_down` also leans positive and reaches significance in Tcon (NES 1.47, FDR 0.026) though not in Treg (0.97, FDR 0.51) or CD8 (1.09, FDR 0.31), so the signature's two arms do not separate in opposite directions.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/09_heat_hypoxia.py` | `run_fgsea` | `gsea_min_size=5`, `gsea_max_size=500`, `gsea_seed=123`, `gsea_nperm=100000` | `03_results/03_pseudobulk/tables/ranked_{treg,tcon,cd8}.tsv`, `03_results/09_heat_hypoxia/tables/_signatures_full/WT_heat_{up,down}.txt` |

**How to read:** One file per sorted population, two rows each (`WT_heat_up`, `WT_heat_down`). Positive `nes` means enrichment toward the SF-high end of the donor-pseudobulk SF-vs-PB ranked list; `padj` is BH across the two sets in that run only. `set_size` counts signature genes surviving intersection with the ranked list (119/130/113 of 199 up; 56/61/57 of 94 down) -- roughly half to two-thirds of the projected signature is testable here. `core_enrichment` is the slash-separated leading edge. Primary donor-pseudobulk tier, and the unpurged reference for the hypoxia-purged run. The positive `WT_heat_down` arm is the caveat, and in Tcon it is significant: the SF-high shift is not a clean bidirectional recapitulation of the mouse contrast, because both arms move the same way.

## tables/gsea_purged_{treg,tcon,cd8}.csv

Removing the 18 `HALLMARK_HYPOXIA` overlap genes leaves `WT_heat_up` SF-high in all three populations -- NES 2.43 (Treg), 2.55 (Tcon), 1.93 (CD8), FDR <= 4.1e-5, a loss of only 0.13 to 0.16 NES -- and leaves `WT_heat_down` untouched (NES 0.97/1.47/1.09), since none of its 94 genes overlap hypoxia.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/09_heat_hypoxia.py` | `run_fgsea` | `gsea_min_size=5`, `gsea_max_size=500`, `gsea_seed=123`, `gsea_nperm=100000` | `03_results/03_pseudobulk/tables/ranked_{treg,tcon,cd8}.tsv`, `03_results/09_heat_hypoxia/tables/_signatures_purged/WT_heat_{up,down}.txt` |

**How to read:** Positive `nes` is SF-high and `padj` is BH within each run. The contrast tag distinguishes purged rows. `WT_heat_up` effective size falls from 119/130/113 to 107/115/101: 18 genes leave the frozen set, but only 12/15/12 were testable in the corresponding ranking. These primary donor-pseudobulk rows feed `gene_purge_nes_comparison.csv`, which owns the paired comparison. A positive, significant purged NES supports only the bounded gene-content statement; it says nothing about temperature.

## tables/_signatures_full/WT_heat_{up,down}.txt

The frozen mouse-anchor human-ortholog `WT_heat` sets exactly as handed to fgsea -- 199 up genes and 94 down genes -- of which 113 to 130 (up) and 56 to 61 (down) appear in the JIA donor-pseudobulk ranked lists, so 57% to 65% of the up arm and about 60% of the down arm is actually testable in this compartment.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/09_heat_hypoxia.py` | `prepare_signature_dirs` | `signature_contract = ../mouse_anchor/03_results/human_projection/` | `../mouse_anchor/03_results/human_projection/signatures/WT_heat/WT_heat_{up,down}.txt` |

**How to read:** Plain newline-delimited HGNC symbols, one per line, no header, alphabetically ordered. `_up` are the genes raised at 39 C in the mouse anchor and `_down` those lowered, projected to human orthologs; the sign lives in the filename, not in the file. These are inputs, not results -- their value is provenance and reproducibility: the exact gene universe behind the primary donor-pseudobulk NES, regenerated verbatim from the frozen contract on every run. Diff them against `_signatures_purged/` to see precisely what the hypoxia purge removed.

## tables/_signatures_purged/WT_heat_{up,down}.txt

The hypoxia-purged inputs: 181 up genes after dropping the 18 `HALLMARK_HYPOXIA` members (ADM, ADORA2B, AK4, ANXA2, ATF3, CCN1, CDKN1A, EGFR, F3, FOSL2, HK2, IER3, P4HA2, PDGFB, PLAUR, SDC4, SERPINE1, TGM2 -- 9.0% of the up set), while the 94-gene down list is identical to the full one because no down gene overlaps hypoxia.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/09_heat_hypoxia.py` | `prepare_signature_dirs` | `signature_contract = ../mouse_anchor/03_results/human_projection/` | `../mouse_anchor/03_results/human_projection/signatures/WT_heat/WT_heat_{up,down}.txt`, `00_data/references/msigdb_hallmark/HALLMARK_HYPOXIA.txt` |

**How to read:** Same format and sign convention as `_signatures_full/` -- newline-delimited HGNC symbols, alphabetical, direction carried by the filename. The purge is a plain set difference against the 200-gene `HALLMARK_HYPOXIA` reference, applied to both arms; that it changes only the up list is itself informative, since the hypoxia overlap is entirely on the SF-high side. Inputs rather than results, at the primary donor-pseudobulk tier: they define the gene universe for `gsea_purged_*`, and the removed-gene list is echoed in the `genes_removed` column of `gene_purge_nes_comparison.csv`.

## tables/_overview/heat_purge_nes_paired.csv

One row per population and mouse arm, pairing the full and hypoxia-purged NES side by side so the cost of the purge reads as a single subtraction.

**How to read:** `nes_full`/`nes_purged` with their FDRs and testable `set_size` are the plotted marker positions; positive NES is synovial-fluid-high. Six rows, six markers, nothing aggregated. Primary donor-pseudobulk tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/09_heat_hypoxia_viz.py` | `purge_paired_table` | `thresholds.gsea_fdr=0.05` | `03_results/09_heat_hypoxia/tables/gsea_{full,purged}_{treg,tcon,cd8}.csv` |

## tables/_overview/heat_hypoxia_colocalization.csv

The plotted cell-level rows only: within-SF correlation of the per-cell heat and hypoxia scores, Spearman and Pearson, per population.

**How to read:** `r` is the bar height and `n` the cell count printed under each population. The unpowered donor-level rows are deliberately absent and stay in the full stage table. Secondary per-cell tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/09_heat_hypoxia_viz.py` | `colocalization_table` | `level=cell; tissue=synovial_fluid` | `03_results/09_heat_hypoxia/tables/heat_hypoxia_colocalization.csv` |

## tables/_overview/heat_treg_volcano_signature.csv

Every mouse `WT_heat` gene with a testable Treg SF-vs-PB result, 103 up-arm and 55 down-arm, carrying the fold change and FDR that place it on the volcano.

**How to read:** `log2FoldChange` and `neg_log10_padj` are the plotted coordinates and `passes_de_gates` is the config gate decision, so every tally printed in the figure is recoverable by counting rows. `arm` gives signature membership, `hypoxia_purged` flags the 18 genes the purge removes, and `label_drawn` marks the ten genes named in the figure. `n_genes_tested` and the two `n_gates_*_all_genes` columns repeat the whole-transcriptome denominators on every row. Primary donor-pseudobulk tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/09_heat_hypoxia_viz.py` | `volcano_source_table` | `thresholds.de_fdr=0.05; de_logfc=1.0` | `03_results/03_pseudobulk/tables/de_SFvsPB_treg.csv` |

## tables/_overview/heat_treg_volcano_programs.csv

The 103 testable mouse up-arm genes together with every testable member of the two reference axes, 11 STING-specific and 61 generic type-I IFN, at the same volcano coordinates.

**How to read:** `le_program` carries the frozen leading-edge annotation and reads `not_annotated` for the up-arm genes it leaves uncovered, which is why the panel colours are an annotation rather than a decomposition of the 199-gene set. `reference_axis` names the axis a gene belongs to, and a row with both an `arm` and a `reference_axis` value is one of the few genes the mouse signature and an axis share. `label_drawn` marks the genes named in the figure. Primary donor-pseudobulk tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/09_heat_hypoxia_viz.py` | `volcano_source_table` | `thresholds.de_fdr=0.05; de_logfc=1.0; taxonomy=00_data/references/heat_leadingedge_taxonomy` | `03_results/03_pseudobulk/tables/de_SFvsPB_treg.csv`, `00_data/references/heat_leadingedge_taxonomy/leadingedge_gene_taxonomy.csv`, `../sting_positive_control/03_results/06_reference_axis/signatures/` |

## figures/_overview/heat_purge_nes_paired.png

Deleting the 18 HALLMARK_HYPOXIA-overlap genes from the mouse 39
°C-derived up-set takes 12 to 15 testable genes out of the arm and
costs 0.126 to 0.164 NES — 2.5914 to 2.4271 in Treg, 2.6826 to 2.5565
in Tcon, 2.0614 to 1.9181 in CD8 — leaving all three significant. The
synovial-fluid enrichment therefore survives the removal of its
HALLMARK_HYPOXIA-overlap gene content. This is a statement about gene
content. Temperature is untested here, and cross-sectional human data
leave hypoxia's status as confound or co-exposure undetermined.

**How to read:** This is the confirmatory tier: donor-level pseudobulk within frozen
sort labels, limma-voom then fgsea. Positive NES points toward
synovial fluid. Each row pairs the full set (large diamond) with its
purged form (small circle); the connecting bar is the NES cost. Warm
brown is the up arm and cool blue the down arm. A dark outline marks
FDR below 0.05. Right-hand text reports effective and nominal set
sizes, the NES cost, and purged FDR. Two gene counts differ and both
are given: 18 genes come out of the frozen set, of which 12 to 15 were
present in a ranked list. The Tcon down arm stays significant at the
up arm's sign. This licenses a membership statement. Correlative.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/09_heat_hypoxia_viz.py` | `plot_purge_paired` | `thresholds.gsea_fdr=0.05; gsea_min_size=5; gsea_nperm=100000` | `03_results/09_heat_hypoxia/tables/gsea_{full,purged}_{treg,tcon,cd8}.csv` |

## figures/_overview/heat_hypoxia_colocalization.png

Within synovial-fluid cells the per-cell WT_heat_up score and the
HALLMARK_HYPOXIA score correlate weakly, Spearman 0.08 to 0.20, so
largely different cells carry the two scores. This is the per-cell
tier; it corroborates the membership result.

**How to read:** This per-cell tier corroborates; the confirmatory result is the paired
purge panel. Bars show within-SF cell-level correlation between
WT_heat_up and HALLMARK_HYPOXIA AUCell scores, with Spearman dark and
Pearson light; cell counts sit below each population. The y-axis spans
-0.05 to 1, so bar height is the readable quantity. Positive r means
the scores tend to coincide. Donor-level SF means rest on 6 to 7
donors and sit in the stage table. Keep this diagnostic separate from
the pseudobulk NES; it carries no direction.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/09_heat_hypoxia_viz.py` | `plot_colocalization` | `level=cell; tissue=synovial_fluid; evidence_tier=secondary_percell` | `03_results/09_heat_hypoxia/tables/heat_hypoxia_colocalization.csv` |

## figures/_overview/heat_treg_volcano_signature.png

48 of the 119 testable mouse 39 °C up-arm genes clear the Treg SF-vs-
PB significance gates on the synovial-fluid side against 7 on the
blood side, so the enrichment is visible in the plain differential-
expression view, while accounting for 5.6% of the SF-high response.

**How to read:** Every tested Treg gene, x = log2 fold change synovial fluid over
paired blood, y = -log10 FDR, dashed lines = the config gates. Warm
brown = mouse up-arm member, cool blue = down-arm member, grey = every
other gene, triangle = an up-arm gene the hypoxia purge removes.
Membership is the only thing added to the committed DE table, and the
printed tallies are counts of it. Gene names are capped at the top 10
up-arm genes by FDR and the rest are in the source table. The down arm
scatters both ways, which is the caveat to carry. Primary donor-
pseudobulk tier, correlative.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/09_heat_hypoxia_viz.py` | `plot_signature_volcano` | `thresholds.de_fdr=0.05; de_logfc=1.0; figures.volcano_label_top=10` | `03_results/03_pseudobulk/tables/de_SFvsPB_treg.csv` |

## figures/_overview/heat_treg_volcano_programs.png

2 of the 21 published IFN-independent STING-activation genes and 6 of
200 generic type-I IFN genes sit in the mouse 39 °C up-arm, so the SF-
high program the purge leaves standing is an effector and activation
program with minimal overlap with the STING reference axis.

**How to read:** Two views of one volcano, same axes as the signature volcano. Left
colours the up-arm genes by leading-edge program, with the 66-gene
annotation covering the leading edge only and pale brown marking the
up-arm genes it leaves unlabelled. It is an annotation of that leading
edge; the unlabelled remainder is the measure of what it leaves out.
Right draws the two frozen reference axes: black squares are the
published IFN-independent STING-activation genes, all named, blue
circles the generic type-I IFN program, and a brown outline means the
gene is also a mouse signature member. The heat-shock trio is named in
the left panel note rather than at its markers. Primary donor-
pseudobulk tier, correlative.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/09_heat_hypoxia_viz.py` | `plot_programs_axes_volcano` | `thresholds.de_fdr=0.05; de_logfc=1.0; taxonomy=00_data/references/heat_leadingedge_taxonomy` | `03_results/03_pseudobulk/tables/de_SFvsPB_treg.csv, 00_data/references/heat_leadingedge_taxonomy/leadingedge_gene_taxonomy.csv, ../sting_positive_control/03_results/06_reference_axis/signatures/` |

## tables/source_hash_manifest.csv

The stage-09 mouse-signature and SAVI-axis reads are pinned to the source files used for this render.

**How to read:** `source_label` names the dependency, `source_path` is relative to the umbrella
checkout, and `sha256` is the required byte hash. The compute and viz scripts stop if the mouse
projection or SAVI reference-axis files drift.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/09_heat_hypoxia.py` | `verify_source_hashes()` | pinned SHA-256 | mouse projection and SAVI reference-axis signature files |
