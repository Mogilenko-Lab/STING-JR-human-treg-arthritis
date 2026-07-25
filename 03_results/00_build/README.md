# 00_build — artifact captions

_**Abbreviations:** SF = synovial fluid (inflamed joint); PB = peripheral blood. The SF-vs-PB contrast is paired within each of the 7 JIA donors. Treg = CD4⁺CD127ˡᵒCD25⁺ regulatory; Tcon = CD4⁺CD25⁻ conventional; CD8 = CD8⁺CD45RO⁺ memory._

## figures/_overview/cells_per_gsm.png

All 7 donors contribute paired SF+PB Tregs; Tcon and CD8 lack a PB
sample for p3 (by design). SF-Treg p5 is the thinnest stratum.

**How to read:** Grouped bars = cells recovered per donor; orange = synovial fluid
(SF), blue = peripheral blood (PB); one facet per sorted population. A
missing PB bar (p3 in Tcon/CD8) is an intentionally-absent sample, not
a QC drop. Descriptive counts only — no claim tier.

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
