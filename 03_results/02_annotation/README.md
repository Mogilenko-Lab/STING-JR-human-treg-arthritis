# 02_annotation — artifact captions

_**Abbreviations:** SF = synovial fluid (inflamed joint); PB = peripheral blood. The SF-vs-PB contrast is paired within each of the 7 JIA donors. Treg = CD4⁺CD127ˡᵒCD25⁺ regulatory; Tcon = CD4⁺CD25⁻ conventional; CD8 = CD8⁺CD45RO⁺ memory._

## figures/_overview/umap_sort_identity.png

Frozen sort labels track the transcriptomic structure; the large
majority of cells are marker-consistent, so the sort gate is a sound
anchor for pseudobulk.

**How to read:** Unsupervised UMAP colored by frozen sort label, marker-module argmax
prediction, and their agreement (blue=consistent). Disagreement flags
candidate mis-sorts, not Treg collapse. Annotation/viz only — no
biological claim.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_annotate_states_viz.py` | `main` | `basis_of_label = sorting (scANVI deferred until go = yes)` | `03_results/objects/02_annotation.h5ad` |

## figures/_overview/marker_dotplot.png

FOXP3/IL2RA/CTLA4/IKZF2 are Treg-restricted while CD8A/B/GZMK mark the
CD8 gate; IL7R is depleted in Tregs (CD127-lo sort), confirming gate
fidelity.

**How to read:** Dot size = fraction of cells expressing; color = mean lognorm
expression. Rows = frozen label. Confirms markers land where the sort
predicts. QC overlay tier (hand markers, not evidence).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_annotate_states_viz.py` | `main` | `LINEAGE_MODULES (Treg/Tcon/CD8 canonical markers)` | `03_results/02_annotation/tables/substate_markers.csv` |

## figures/_overview/counts_grid.png

Every SF+PB Treg stratum clears the pseudobulk floor across the 7
donors; p3 PB Tcon/CD8 are absent by design (empty cells).

**How to read:** Heatmap of cells per donor (x) x label+tissue (y); red * marks a
stratum below the pseudobulk cell floor; empty = intentionally-absent
sample. Donor count per arm is the forest's power. Diagnostic.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_annotate_states_viz.py` | `main` | `thresholds.pseudobulk_min_cells = 20` | `03_results/02_annotation/tables/counts_donor_by_label_tissue.csv` |
