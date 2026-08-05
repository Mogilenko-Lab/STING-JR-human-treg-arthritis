# 00_build: artifact captions

_**Abbreviations:** SF = synovial fluid (inflamed joint); PB = peripheral blood. The cohort contains 7 JIA donors; ingest includes SF+PB Treg samples for all 7, while the post-QC donor-level analysis retains 6 paired donors in each population. Treg = CD4⁺CD127ˡᵒCD25⁺ regulatory; Tcon = CD4⁺CD25⁻ conventional; CD8 = CD8⁺CD45RO⁺ memory._

## figures/_overview/cells_per_gsm.png

At ingest, all 7 donors contribute SF+PB Treg samples; Tcon and CD8
lack a PB sample for p3 by design. The near-empty SF-Treg p5 sample is
later removed by QC, leaving 6 paired donors in each analyzed
population.

**How to read:** Grouped bars = cells recovered per donor; orange = synovial fluid
(SF), blue = peripheral blood (PB); one facet per sorted population. A
missing PB bar (p3 in Tcon/CD8) is an intentionally-absent sample, not
a QC drop. These are ingest counts before QC; the donor-level analysis
uses 6 paired donors per population. Descriptive counts only — no
claim tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/00_build_anndata_viz.py` | `main` | `design.populations = [CD4_Treg, CD4_Tcon, CD8]` | `03_results/00_build/tables/cells_per_gsm.csv` |

## tables/design_completeness.csv

The paired design is complete for Tregs — 7/7 donors in both SF and PB
— and short by exactly one donor in PB Tcon and PB CD8, so 40 of the
42 possible donor x tissue x population samples exist.

**How to read:** One row per sorted population x tissue; `n_donors` = distinct
donors contributing at least one cell. 7 = complete, 6 = the one
intentionally-absent donor. Counted at ingest, before QC, so the
near-empty SF-Treg library that QC later drops still counts here.
Descriptive design audit — no claim tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/00_build_anndata.py` | `main` | `design.populations = [CD4_Treg, CD4_Tcon, CD8]`, `design.n_donors = 7` | `00_data/GSE160097_JIA-SF-Treg/samples.csv`, `00_data/GSE160097_JIA-SF-Treg/raw/` (40 per-GSM 10x H5) |

## tables/genes_union_summary.csv

The 40-GSM gene union is 32,738 features, of which 24,374 (74.5%) are
seen in at least one of the 108,414 pooled cells — a quarter of the
union is all-zero padding that the per-gene detection filter removes.

**How to read:** Single row. `n_genes_union` = features after the outer join across
GSMs; `n_genes_detected_any_cell` = features with >0 UMI in >=1 cell,
computed on the integer `layers['counts']`; `n_cells_total` is pre-QC;
`species_db` drives MT/RP annotation and MSigDB species downstream.
Ingest audit — no claim tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/00_build_anndata.py` | `main` | `project.species_db = HS`, `project.genome_build = GRCh38` (the union itself is an outer join inside `build_pooled` — no config key) | `00_data/GSE160097_JIA-SF-Treg/samples.csv`, `00_data/GSE160097_JIA-SF-Treg/raw/` (40 per-GSM 10x H5) |

## tables/reference_feature_symbols.csv

The same 32,738-feature union as above, by name rather than by count.
The count alone cannot answer the question this list exists for: when a
gene set member is missing from the analysis, is the gene absent from
the CellRanger reference, present in it but never detected in sorted
T cells, or present and detected under a different symbol? `EGFR`,
`EPCAM` and `INHBA` are all here and all absent downstream — a
detection fact. `IFNB1` is here too and undetected, while `NLRC3` and
`MIR4691` are not in the union at all. `MB21D1` and `TMEM173` are here
under names no current gene set uses.

**How to read:** One row per feature of the 40-GSM union, `ensembl_id` = the
`var_name` carried through every downstream object, `gene_symbol` = the
symbol the CellRanger reference assigns it. Symbols are unique here
(32,738 of 32,738), and the vintage is that reference's, not a current
HGNC release's. This is the outermost of three nested vocabulary layers
— union (32,738) > post-QC `gene_symbols.csv` (21,740) > post-`filterByExpr`
ranked list (~14,000) — and absence from THIS layer is the only true
"absent from the reference". Ingest audit — no claim tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/00_build_anndata.py` | `main` | `symbol_alias.reference_feature_symbols` names this file for its consumers (the union itself is an outer join inside `build_pooled` — no config key) | `00_data/GSE160097_JIA-SF-Treg/samples.csv`, `00_data/GSE160097_JIA-SF-Treg/raw/` (40 per-GSM 10x H5) |
