# 09_heat_hypoxia -- artifact captions

_**Abbreviations:** SF = synovial fluid, PB = peripheral blood, NES = normalized enrichment score, FDR = BH-adjusted p-value._

I added a correlative heat-vs-hypoxia check for the JIA SF-vs-PB `WT_heat` signal. The primary read is donor-pseudobulk fgsea after removing `HALLMARK_HYPOXIA` overlap genes from the mouse `WT_heat_up` set. The secondary reads ask whether per-cell heat and hypoxia scores co-localize within SF cells, and which biological programs the `WT_heat_up` leading-edge genes represent. Hypoxia is a transcriptional readout here, not a HIF-causality claim.

The figures walk the same argument in order: how much of the enrichment survives the purge, whether heat-high and hypoxia-high are even the same cells, and what the enriching genes turn out to be. The third is where the answer stops being reassuring, and it is why an activation-free proteostasis lens was built next.

Two further figures put the mouse sets back into the plain Treg SF-vs-PB volcano, so a reader can see how much of the differential-expression response the signature accounts for, and how little it shares with the published IFN-independent STING-activation signature of de Cevins et al. 2023 (Cell Rep Med, PMID 38118407) that the STING positive-control compartment carries as a reference axis.

## gene_purge_nes_comparison.csv

Removing 18 `HALLMARK_HYPOXIA` overlap genes reduces the `WT_heat_up` NES modestly, but the enrichment remains SF-high in Treg, Tcon, and CD8 at FDR < 5e-5.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/09_heat_hypoxia.py` | `gene_purge_nes` | `gsea_min_size=5`, `gsea_max_size=500`, `gsea_seed=123`, `gsea_nperm=100000` | `03_results/03_pseudobulk/tables/ranked_{treg,tcon,cd8}.tsv`, mouse `WT_heat_up/down.txt`, `00_data/references/msigdb_hallmark/HALLMARK_HYPOXIA.txt` |

**How to read:** One row per sorted population. Positive NES means the set is enriched toward the SF-high end of the SF-vs-PB ranked list. `NES_full` is the original `WT_heat_up` score, while `NES_purged` is the same fgsea engine after removing hypoxia-overlap genes. This is the primary donor-pseudobulk tier. I read a positive, significant purged NES as evidence that a hypoxia-overlap removal does not erase the correlative heat-axis signal.

## heat_hypoxia_colocalization.csv

Within SF cells, `WT_heat_up_AUCell` and `HALLMARK_HYPOXIA_AUCell` show weak positive cell-level correlations, while donor-level SF means are not positive in this small donor set.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/09_heat_hypoxia.py` | `heat_hypoxia_colocalization` | `tissue_levels.synovial_fluid=synovial_fluid`, `donor_key=donor` | `03_results/interactive/08_harvest_readout.parquet`, `03_results/05_scoring/tables/per_cell_scores.csv` |

**How to read:** Rows are stratified by population, level, and correlation method. `level=cell` uses SF cells directly. `level=donor_sf_mean` correlates per-donor SF mean heat and hypoxia scores. Positive `r` means higher heat score tends to sit with higher hypoxia score. This is an L3 secondary per-cell read and is not pooled with the pseudobulk NES. The cell-level correlation is weak (Spearman 0.08 to 0.20), so heat-high and hypoxia-high are largely different cells. The donor-level correlation rests on only 6 to 7 donors and is effectively unpowered, so its sign is not interpretable and must not be read as heat and hypoxia being anti-correlated.

## leadingedge_composition.csv

The `WT_heat_up` leading edge in SF T cells is predominantly T-cell activation and effector genes (48% to 57%) plus immediate-early stress genes (11% to 20%), with a hypoxia-overlap minority (14% to 17%) and only a trace of classic heat-shock or proteostasis genes (4% to 6%: HSPA1A, HSPH1, CLU).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/09_heat_hypoxia.py` | `leadingedge_composition` | taxonomy from `00_data/references/heat_leadingedge_taxonomy/` | `03_results/09_heat_hypoxia/tables/runsum_interactive_gsea_full_{treg,tcon,cd8}_WT_heat_up.csv`, `00_data/references/heat_leadingedge_taxonomy/leadingedge_gene_taxonomy.csv` |

**How to read:** One row per population. `n_leading_edge` is the count of fgsea core-enrichment genes from the full `WT_heat_up` run. Each is assigned to one program — `heat_shock_proteostasis`, `hypoxia_HIF`, `immediate_early_stress`, `effector_activation`, `other` — from a frozen gene taxonomy (external large-context-model classification, provenance in `00_data/references/heat_leadingedge_taxonomy/`). The `n_`/`frac_` columns tally each program and the `genes_` columns list members. Exploratory secondary tier, never pooled with the pseudobulk NES. The dominance of activation and immediate-early genes, against only a trace of classic heat-shock, is a caution: the SF-vs-PB `WT_heat` enrichment survives the gene purge above (hypoxia-independent) but is carried mostly by a generic activation program, not a thermal-specific one.

## tables/gsea_full_{treg,tcon,cd8}.csv

With the complete 199-gene mouse `WT_heat_up` set, SF-vs-PB enrichment is positive and strong in every sorted population -- NES 2.51 (Treg), 2.57 (Tcon), 2.05 (CD8), all at FDR <= 8.4e-7 -- while `WT_heat_down` is non-significant everywhere (NES 1.01 to 1.34, FDR 0.07 to 0.43) and also leans positive, so the signature's two arms do not separate in opposite directions.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/09_heat_hypoxia.py` | `run_fgsea` | `gsea_min_size=5`, `gsea_max_size=500`, `gsea_seed=123`, `gsea_nperm=100000` | `03_results/03_pseudobulk/tables/ranked_{treg,tcon,cd8}.tsv`, `03_results/09_heat_hypoxia/tables/_signatures_full/WT_heat_{up,down}.txt` |

**How to read:** One file per sorted population, two rows each (`WT_heat_up`, `WT_heat_down`). Positive `nes` means enrichment toward the SF-high end of the donor-pseudobulk SF-vs-PB ranked list; `padj` is BH across the two sets in that run only. `set_size` counts signature genes surviving intersection with the ranked list (105/111/115 of 199 up; 57/61/62 of 94 down) -- roughly half the projected signature is testable here. `core_enrichment` is the slash-separated leading edge. Primary donor-pseudobulk tier, and the unpurged reference for the hypoxia-purged run. The non-significant, positive `WT_heat_down` arm is the caveat: the SF-high shift is not a clean bidirectional recapitulation of the mouse contrast.

## tables/gsea_purged_{treg,tcon,cd8}.csv

Removing the 18 `HALLMARK_HYPOXIA` overlap genes leaves `WT_heat_up` SF-high in all three populations -- NES 2.38 (Treg), 2.39 (Tcon), 1.90 (CD8), FDR <= 4.6e-5, a loss of only 0.14 to 0.19 NES -- and leaves `WT_heat_down` untouched (NES 1.01/1.34/1.23), since none of its 94 genes overlap hypoxia.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/09_heat_hypoxia.py` | `run_fgsea` | `gsea_min_size=5`, `gsea_max_size=500`, `gsea_seed=123`, `gsea_nperm=100000` | `03_results/03_pseudobulk/tables/ranked_{treg,tcon,cd8}.tsv`, `03_results/09_heat_hypoxia/tables/_signatures_purged/WT_heat_{up,down}.txt` |

**How to read:** Same schema and sign convention as the full run above -- positive `nes` is SF-high, `padj` is BH within the run -- with the `contrast` column tagged `SF_vs_PB_<population>_hypoxia_purged` so the two families never get confused. `set_size` drops to 95/100/103 for `WT_heat_up` as the hypoxia genes leave, and `core_enrichment` correspondingly loses CDKN1A, ANXA2, SDC4, ATF3, PLAUR and friends. Primary donor-pseudobulk tier; these files are the per-population source rows behind `gene_purge_nes_comparison.csv`, which is where the paired comparison should be read. A still-positive, still-significant NES here means the correlative heat-axis signal is not merely a hypoxia program wearing a heat label -- it does not, on its own, make it thermal-specific.

## tables/_signatures_full/WT_heat_{up,down}.txt

The frozen mouse-anchor human-ortholog `WT_heat` sets exactly as handed to fgsea -- 199 up genes and 94 down genes -- of which only 105 to 115 (up) and 57 to 62 (down) appear in the JIA donor-pseudobulk ranked lists, so a little over half the projected signature is actually testable in this compartment.

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

## tables/_overview/heat_leadingedge_composition.csv

The plotted composition rows: leading-edge gene counts and fractions per biological program, one row per population.

**How to read:** `frac_<program>` are the stacked segment widths and `n_<program>` the counts printed inside them; `n_leading_edge` is the bar total. Exploratory secondary tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/09_heat_hypoxia_viz.py` | `leadingedge_table` | `taxonomy=00_data/references/heat_leadingedge_taxonomy` | `03_results/09_heat_hypoxia/tables/leadingedge_composition.csv` |

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

Removing every HALLMARK_HYPOXIA gene from the mouse 39 °C up-set costs
only 0.14 to 0.19 NES and leaves the synovial-fluid enrichment strong
and significant in all three sorted populations, so hypoxia does not
explain it.

**How to read:** One row per population and mouse arm. x is fgsea NES, positive =
enriched toward the synovial-fluid end of the paired ranking. Each row
pairs the full mouse set (large diamond) with the hypoxia-purged set
(small circle), and the bar between them is what the purge cost. Warm
brown = up arm, cool blue = down arm. Every marker is filled and
translucent, so the down-arm pair, which the purge leaves untouched,
reads as a darker circle inside a lighter diamond. A heavier dark
outline marks FDR below 0.05 and is the only significance glyph. The
right-hand text prints set size before and after the purge and the
purged FDR. Primary donor-pseudobulk tier, correlative.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/09_heat_hypoxia_viz.py` | `plot_purge_paired` | `thresholds.gsea_fdr=0.05; gsea_min_size=5; gsea_nperm=100000` | `03_results/09_heat_hypoxia/tables/gsea_{full,purged}_{treg,tcon,cd8}.csv` |

## figures/_overview/heat_hypoxia_colocalization.png

Within synovial-fluid cells the mouse heat score and the hypoxia score
correlate only weakly (Spearman 0.08 to 0.20), so the niche's thermal
and hypoxic readouts are carried by largely different cells rather
than one shared stress state.

**How to read:** Bars are the within-SF, cell-level correlation between the per-cell
WT_heat_up and HALLMARK_HYPOXIA AUCell scores, Spearman (dark) beside
Pearson (light), with the cell count under each population. The y-axis
deliberately runs the full -0.05 to 1 range: read the shortness of the
bars, not their rank order. Positive r means a heat-high cell tends to
be hypoxia-high. Donor-level SF means are unpowered at 6 to 7 donors
and are left in the stage table rather than drawn. This is a secondary
per-cell diagnostic of where the two scores sit, never pooled with the
pseudobulk NES and never read as directional evidence.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/09_heat_hypoxia_viz.py` | `plot_colocalization` | `level=cell; tissue=synovial_fluid; evidence_tier=secondary_percell` | `03_results/09_heat_hypoxia/tables/heat_hypoxia_colocalization.csv` |

## figures/_overview/heat_leadingedge_composition.png

Half the WT_heat up leading edge in synovial-fluid T cells is effector
and activation genes, with a hypoxia-overlap minority and only two or
three classic heat-shock genes, so surviving the hypoxia purge does
not make the enrichment thermally specific.

**How to read:** One stacked bar per population, spanning that population's full
WT_heat up leading edge; segment width is the fraction of leading-edge
genes in each program and the number inside a segment is its gene
count (printed where the segment is wide enough to hold it; every
count is in the source table). Program assignment comes from a frozen
external-model gene taxonomy, not from this run. The narrow heat-shock
segment is the point, and it is why an activation-free proteostasis
lens was built next. Exploratory secondary tier, never pooled with the
pseudobulk NES.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/09_heat_hypoxia_viz.py` | `plot_leadingedge` | `taxonomy=00_data/references/heat_leadingedge_taxonomy; evidence_tier=secondary_exploratory` | `03_results/09_heat_hypoxia/tables/leadingedge_composition.csv` |

## figures/_overview/heat_treg_volcano_signature.png

40 of the 103 testable mouse 39 °C up-arm genes clear the Treg SF-vs-
PB significance gates on the synovial-fluid side against 5 on the
blood side, so the enrichment is visible in the plain differential-
expression view, while accounting for 4.9% of the SF-high response.

**How to read:** Every tested Treg gene, x = log2 fold change synovial fluid over
paired blood, y = -log10 FDR, dashed lines = the config gates. Warm
brown = mouse up-arm member, cool blue = down-arm member, grey = every
other gene, triangle = an up-arm gene the hypoxia purge removes.
Membership is the only thing added to the committed DE table, and the
printed tallies are counts of it. Gene names are capped at the top 10
up-arm genes by FDR and the rest are in the source table. The down arm
scattering both ways is the honest caveat. Primary donor-pseudobulk
tier, correlative.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/09_heat_hypoxia_viz.py` | `plot_signature_volcano` | `thresholds.de_fdr=0.05; de_logfc=1.0; figures.volcano_label_top=10` | `03_results/03_pseudobulk/tables/de_SFvsPB_treg.csv` |

## figures/_overview/heat_treg_volcano_programs.png

Only 2 of the 21 published IFN-independent STING-activation genes and
7 of 200 generic type-I IFN genes are in the mouse 39 °C up-arm, so
the SF-high program the purge leaves standing is an effector and
activation program that shares almost nothing with the STING reference
axis.

**How to read:** Two views of one volcano, same axes as the signature volcano. Left
colours the up-arm genes by leading-edge program, with the 66-gene
annotation covering the leading edge only and pale brown marking the
up-arm genes it leaves unlabelled, so read it as an annotation and not
a decomposition of the 199-gene set. Right draws the two frozen
reference axes: black squares are the published IFN-independent STING-
activation genes, all named, blue circles the generic type-I IFN
program, and a brown outline means the gene is also a mouse signature
member. The heat-shock trio is named in the left panel note rather
than at its markers. Primary donor-pseudobulk tier, correlative.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/09_heat_hypoxia_viz.py` | `plot_programs_axes_volcano` | `thresholds.de_fdr=0.05; de_logfc=1.0; taxonomy=00_data/references/heat_leadingedge_taxonomy` | `03_results/03_pseudobulk/tables/de_SFvsPB_treg.csv, 00_data/references/heat_leadingedge_taxonomy/leadingedge_gene_taxonomy.csv, ../sting_positive_control/03_results/06_reference_axis/signatures/` |
