# 02_annotation — Auditing the sort labels against expression

**Sorting is the label of record throughout this compartment.** Every artifact here audits that
label against expression, and the label itself stays fixed. Freezing the FACS gate is what lets
the pseudobulk contrast run inside a fixed cell state.

The audit passes. Marker-module argmax reproduces the sort gate for 85.6% of the 99,915
QC-passing cells, the off-diagonal mass is overwhelmingly Tcon against CD8, and the twelve panel
markers land where the sort predicts. All 39 surviving donor × label × tissue strata clear the
20-cell pseudobulk floor, and every arm keeps at least six donors.

---

## Figures

### `figures/_overview/umap_sort_identity.png`

**The frozen labels against the transcriptomic structure.**
Three panels on the unsupervised UMAP: cells coloured by frozen sort label, by marker-module
argmax prediction, and by their agreement (blue for consistent). Disagreement flags a candidate
mis-sort. The large majority of cells are marker-consistent, which is what makes the sort gate a
sound anchor for pseudobulk.
*Source* `../objects/02_annotation.h5ad` · `02_analysis/scripts/02_annotate_states_viz.py`.

### `figures/_overview/marker_dotplot.png`

**Twelve canonical markers, three frozen labels.**
Dot size gives the fraction of cells expressing and colour the mean log-normalised expression.
Rows are the frozen labels. FOXP3, IL2RA, CTLA4 and IKZF2 are Treg-restricted, CD8A, CD8B and
GZMK mark the CD8 gate, and IL7R runs lowest in Tregs, as a CD127-lo sort requires.
*Source* `tables/substate_markers.csv` · `02_analysis/scripts/02_annotate_states_viz.py`.

### `figures/_overview/counts_grid.png`

**The power budget of the paired contrast.**
Heatmap of cells per donor (x) against label plus tissue (y). A red asterisk marks a stratum below
the 20-cell pseudobulk floor and an empty square marks an absent sample. The cohort holds seven
donors and six span both tissues in each analysed population after QC. Every observed stratum
clears the floor. The p3 blood Tcon and blood CD8 samples were never collected.
*Source* `tables/counts_donor_by_label_tissue.csv` ·
`02_analysis/scripts/02_annotate_states_viz.py`.

---

## Tables

### `tables/confusion_sort_vs_predicted.csv`

Rows are the frozen sort label, columns the argmax of the three z-scored canonical lineage
modules, entries cell counts, and the diagonal is agreement. Rows and columns are forced to the
same Treg / Tcon / CD8 order, so the diagonal reads directly.

Agreement runs 89.6% in Treg, 83.8% in Tcon and 84.4% in CD8, 85.6% overall. The off-diagonal
mass is overwhelmingly Tcon against CD8, and 195 sorted Tregs land on the CD8 module, so the Treg
gate is clean enough to freeze. Each module holds four or five genes, and shared cytotoxic and
activation genes blur Tcon against CD8. Demonstrating a mis-sort takes more evidence than this
table supplies.

### `tables/counts_donor_by_label_tissue.csv`

One row per donor × frozen coarse label × tissue, with `n_cells` counting cells surviving QC.
Read each row against the per-stratum floor (`pseudobulk_min_cells = 20`) and each label × tissue
arm against the donor floor (`pseudobulk_min_donors = 3`).

All 39 surviving strata clear the cell floor. The thinnest, blood Treg p7, holds 266 cells. Every
arm keeps at least six donors, so each contrast sits well above the donor floor. Absent rows carry
two distinct causes: p3 blood Tcon and blood CD8 were never collected, and synovial Treg p5 is the
near-empty library QC excluded. The compute step applies no floor itself, leaving the gating to
pseudobulk.

### `tables/substate_markers.csv`

36 rows — twelve panel markers × three frozen labels. `mean_lognorm` is the mean log-normalised
expression across every cell carrying that label, zeros included, so it is diluted by
non-expressers. `frac_expressing` is the fraction with more than 0 UMI.

FOXP3 reaches mean 1.69 in 80% of sorted Tregs against 0.06 in Tcon and 0.02 in CD8. CD8A marks
93% of CD8 cells and under 0.5% of cells in either other gate. IL7R runs highest in Tcon (2.35)
and lowest in Treg (0.66), as a CD127-lo gate requires.

This table carries no test statistic and no sign convention. The pass condition is high in the
intended label and low elsewhere. Hand markers, QC-overlay tier.

**Provenance of the marker panel.** The twelve markers and the three z-scored lineage modules are
this compartment's own curated panel of canonical T-cell identity genes, assembled here. Every
scored signature elsewhere in the tree comes from an external resource, each listed with its
accession and derivation in the [results index](../README.md).
