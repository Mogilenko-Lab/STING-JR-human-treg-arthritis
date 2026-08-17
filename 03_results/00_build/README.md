# 00_build — What the forty libraries contain

Forty per-GSM 10x matrices are read from GSE160097 and joined into one pooled object: seven JIA
donors × two tissues × three sorted populations, minus the two samples the study never collected.
This stage counts what arrived and writes the feature vocabulary every later stage matches
against. It computes no statistics.

**The design is complete for Tregs.** All seven donors contribute both a synovial-fluid and a
peripheral-blood Treg sample. Tcon and CD8 have no blood sample for donor p3, so 40 of the 42
possible donor × tissue × population samples exist. QC later removes the near-empty SF-Treg p5
library, which leaves six paired donors in each analysed population.

---

## Figures

### `figures/_overview/cells_per_gsm.png`

**Cells recovered per library, before QC.**
Grouped bars, one facet per sorted population. x, donor; y, cells recovered. Orange gives
synovial fluid and blue paired peripheral blood. The missing p3 blood bar in the Tcon and CD8
facets marks a sample that was never collected. These are ingest counts taken before any filter.
*Source* `tables/_overview/cells_per_gsm.csv` ·
`02_analysis/scripts/00_build_anndata_viz.py`.

---

## Tables

### `tables/design_completeness.csv`

One row per sorted population × tissue. `n_donors` counts distinct donors contributing at least
one cell: 7 is complete, 6 carries the one absent donor. Tregs read 7 and 7. Blood Tcon and blood
CD8 read 6. Counted at ingest, so the near-empty library QC later drops still counts here.

### `tables/genes_union_summary.csv`

A single row summarising the feature join. The forty-GSM gene union holds **32,738** features, of
which **24,374 (74.5%)** carry a read in at least one of the 108,414 pooled cells — so a quarter
of the union is all-zero padding the per-gene detection filter removes.

`n_genes_union` is the feature count after the outer join across GSMs.
`n_genes_detected_any_cell` counts features with more than 0 UMI in at least one cell, computed
on the integer `layers['counts']`. `n_cells_total` is pre-QC. `species_db` drives the
mitochondrial and ribosomal annotation and the MSigDB species downstream.

### `tables/reference_feature_symbols.csv`

The same 32,738-feature union, named gene by gene. **This list is the arbiter whenever a gene-set
member goes missing downstream.** It separates three causes: the gene is absent from the
CellRanger reference, the gene is present in it and undetected in sorted T cells, or the gene is
present under a different symbol.

`EGFR`, `EPCAM`, `INHBA` and `IFNB1` are all in the union and all undetected downstream, which is
a detection fact. `NLRC3` and `MIR4691` sit outside the union. `MB21D1` and `TMEM173` are in it,
under the symbols this reference assigns — the names current gene sets carry as `CGAS` and
`STING1`.

One row per feature. `ensembl_id` is the `var_name` carried through every downstream object and
`gene_symbol` is the symbol the CellRanger reference assigns it. Symbols are unique here, 32,738
of 32,738, and their vintage is that reference's.

**This is the outermost of three nested vocabulary layers** — union (32,738) → post-QC
`gene_symbols.csv` (21,740) → post-`filterByExpr` ranked list (~14,000). Absence from this layer
is the only true "absent from the reference".

Every gene set later matched against this vocabulary is external to GSE160097. The accession,
derivation and paper reference of each is in the [results index](../README.md).

### `tables/cells_per_gsm.csv`

One row per GSM with its donor, tissue, population and recovered cell count. The source of the
figure above and the ingest-side companion to
[`../01_qc/tables/cells_kept_dropped.csv`](../01_qc/tables/).
