# 00_build: artifact captions

_**Abbreviations:** SF = synovial fluid (inflamed joint); PB = peripheral blood; Treg =
CD4⁺CD127ˡᵒCD25⁺ regulatory; Tcon = CD4⁺CD25⁻ conventional; CD8 = CD8⁺CD45RO⁺ memory. The cohort
holds 7 JIA donors. Ingest covers SF and PB Treg samples for all 7; the post-QC donor-level
analysis retains 6 paired donors in each population._

## figures/_overview/cells_per_gsm.png

All 7 donors contribute both SF and PB Treg samples at ingest. Tcon
and CD8 have no PB sample for donor p3, an absence in the study
design. QC later removes the near-empty SF-Treg p5 library, leaving 6
paired donors in each analyzed population.

**How to read:** Grouped bars give cells recovered per donor, orange for synovial fluid
and blue for peripheral blood, one facet per sorted population. The
missing p3 PB bar in Tcon and CD8 marks a sample that was never
collected. These are ingest counts, taken before QC. Descriptive
counts — no claim tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/00_build_anndata_viz.py` | `main` | `design.populations = [CD4_Treg, CD4_Tcon, CD8]` | `03_results/00_build/tables/cells_per_gsm.csv` |

## tables/design_completeness.csv

The paired design is complete for Tregs, 7 of 7 donors in both SF and PB. PB Tcon and PB CD8 hold
6 donors each. 40 of the 42 possible donor × tissue × population samples exist.

**How to read:** One row per sorted population × tissue. `n_donors` counts distinct donors
contributing at least one cell: 7 is complete, 6 carries the one absent donor. Counted at ingest,
before QC, so the near-empty SF-Treg library that QC later drops still counts here. Descriptive
design audit — no claim tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/00_build_anndata.py` | `main` | `design.populations = [CD4_Treg, CD4_Tcon, CD8]`, `design.n_donors = 7` | `00_data/GSE160097_JIA-SF-Treg/samples.csv`, `00_data/GSE160097_JIA-SF-Treg/raw/` (40 per-GSM 10x H5) |

## tables/genes_union_summary.csv

The 40-GSM gene union holds 32,738 features. 24,374 of them (74.5%) carry a read in at least one
of the 108,414 pooled cells, so a quarter of the union is all-zero padding that the per-gene
detection filter removes.

**How to read:** Single row. `n_genes_union` is the feature count after the outer join across
GSMs; `n_genes_detected_any_cell` counts features with more than 0 UMI in at least one cell,
computed on the integer `layers['counts']`; `n_cells_total` is pre-QC; `species_db` drives MT/RP
annotation and MSigDB species downstream. Ingest audit — no claim tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/00_build_anndata.py` | `main` | `project.species_db = HS`, `project.genome_build = GRCh38` (the union itself is an outer join inside `build_pooled` — no config key) | `00_data/GSE160097_JIA-SF-Treg/samples.csv`, `00_data/GSE160097_JIA-SF-Treg/raw/` (40 per-GSM 10x H5) |

## tables/reference_feature_symbols.csv

The same 32,738-feature union, named gene by gene. This list answers what a count cannot: when a
gene set member goes missing downstream, whether the gene is absent from the CellRanger reference,
present in it and never detected in sorted T cells, or present under a different symbol. `EGFR`,
`EPCAM`, `INHBA` and `IFNB1` are all in the union and all undetected downstream — a detection
fact. `NLRC3` and `MIR4691` are outside the union. `MB21D1` and `TMEM173` are in it, under symbols
no current gene set uses.

**How to read:** One row per feature of the 40-GSM union. `ensembl_id` is the `var_name` carried
through every downstream object; `gene_symbol` is the symbol the CellRanger reference assigns it.
Symbols are unique here, 32,738 of 32,738, and their vintage is that reference's. This is the
outermost of three nested vocabulary layers — union (32,738) → post-QC `gene_symbols.csv`
(21,740) → post-`filterByExpr` ranked list (~14,000) — and absence from this layer is the only
true "absent from the reference". Ingest audit — no claim tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/00_build_anndata.py` | `main` | `symbol_alias.reference_feature_symbols` names this file for its consumers (the union itself is an outer join inside `build_pooled` — no config key) | `00_data/GSE160097_JIA-SF-Treg/samples.csv`, `00_data/GSE160097_JIA-SF-Treg/raw/` (40 per-GSM 10x H5) |
