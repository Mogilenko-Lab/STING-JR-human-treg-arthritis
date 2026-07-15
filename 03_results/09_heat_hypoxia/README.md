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
