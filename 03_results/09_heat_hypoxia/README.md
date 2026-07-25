# 09_heat_hypoxia -- artifact captions

_**Abbreviations:** SF = synovial fluid, PB = peripheral blood, NES = normalized enrichment score, FDR = BH-adjusted p-value._

I added a correlative heat-vs-hypoxia check for the JIA SF-vs-PB `WT_heat` signal. The primary read is donor-pseudobulk fgsea after removing `HALLMARK_HYPOXIA` overlap genes from the mouse `WT_heat_up` set. The secondary reads ask whether per-cell heat and hypoxia scores co-localize within SF cells, and which biological programs the `WT_heat_up` leading-edge genes represent. Hypoxia is a transcriptional readout here, not a HIF-causality claim.

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
