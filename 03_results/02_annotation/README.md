# 02_annotation: artifact captions

_**Abbreviations:** SF = synovial fluid (inflamed joint); PB = peripheral blood; Treg =
CD4⁺CD127ˡᵒCD25⁺ regulatory; Tcon = CD4⁺CD25⁻ conventional; CD8 = CD8⁺CD45RO⁺ memory. The cohort
holds 7 JIA donors, of whom 6 span both arms in each analyzed population after QC._

Sorting is the label of record throughout this stage. Every artifact here audits that label
against expression; none replaces it.

## figures/_overview/umap_sort_identity.png

The frozen sort labels track the transcriptomic structure, and the large majority of cells are
marker-consistent, so the sort gate is a sound anchor for pseudobulk.

**How to read:** Unsupervised UMAP coloured by frozen sort label, by marker-module argmax
prediction, and by their agreement (blue = consistent). Disagreement flags candidate mis-sorts.
Annotation and visualisation only — no biological claim.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_annotate_states_viz.py` | `main` | `basis_of_label = sorting (scANVI deferred until go = yes)` | `03_results/objects/02_annotation.h5ad` |

## figures/_overview/marker_dotplot.png

FOXP3, IL2RA, CTLA4 and IKZF2 are Treg-restricted; CD8A, CD8B and GZMK mark the CD8 gate; IL7R is
depleted in Tregs, as a CD127-lo sort requires. Markers land where the sort predicts.

**How to read:** Dot size is the fraction of cells expressing, colour the mean log-normalised
expression, rows the frozen label. QC overlay tier — hand markers, which carry no evidential weight.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_annotate_states_viz.py` | `main` | `LINEAGE_MODULES (Treg/Tcon/CD8 canonical markers)` | `03_results/02_annotation/tables/substate_markers.csv` |

## figures/_overview/counts_grid.png

The cohort holds 7 donors and 6 span SF and PB in each analyzed population after QC. Every
observed stratum clears the 20-cell pseudobulk floor. The p3 PB Tcon and PB CD8 samples were never
collected.

**How to read:** Heatmap of cells per donor (x) against label plus tissue (y). A red asterisk marks
a stratum below the pseudobulk cell floor; an empty square marks an absent sample. Donor count per
arm sets contrast precision. Diagnostic.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_annotate_states_viz.py` | `main` | `thresholds.pseudobulk_min_cells = 20` | `03_results/02_annotation/tables/counts_donor_by_label_tissue.csv` |

## tables/confusion_sort_vs_predicted.csv

Marker-module argmax reproduces the sort gate for 85.6% of the 99,915 QC-passing cells — Treg
89.6%, Tcon 83.8%, CD8 84.4%. The off-diagonal mass is overwhelmingly Tcon against CD8, and 195
sorted Tregs land on the CD8 module, so the Treg gate is clean enough to freeze.

**How to read:** Rows are the frozen sort label, columns the argmax of the three z-scored canonical
lineage modules, entries cell counts, and the diagonal is agreement. Rows and columns are forced to
the same Treg/Tcon/CD8 order, so the diagonal reads directly. An off-diagonal cell records module
disagreement; each module is 4-5 genes, and shared cytotoxic and activation genes blur Tcon against
CD8, so a demonstrated mis-sort takes more evidence than this. Diagnostic tier — sorting remains
the label of record.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_annotate_states.py` | `main` | `decisions.treg_gate.basis = sorting`; `LINEAGE_MODULES` (Treg FOXP3/IL2RA/CTLA4/IKZF2/TIGIT, Tcon IL7R/CD40LG/ANK3/LEF1, CD8 CD8A/CD8B/NKG7/GZMK/GZMA) is a literal in the script — no config key | `03_results/objects/01_qc.h5ad` |

## tables/counts_donor_by_label_tissue.csv

All 39 surviving donor × label × tissue strata clear the 20-cell pseudobulk floor; the thinnest, PB
Treg p7, holds 266 cells. Every arm keeps at least 6 donors, so each SF-versus-PB contrast sits
well above the donor floor of 3.

**How to read:** One row per donor × frozen coarse label × tissue; `n_cells` counts cells surviving
QC. Read each row against the per-stratum floor (`pseudobulk_min_cells = 20`) and each label ×
tissue arm against the donor floor (`pseudobulk_min_donors = 3`). Absent rows carry two distinct
causes: p3 PB Tcon and PB CD8 were never collected, and SF Treg p5 is the near-empty library QC
excluded. The compute step writes raw counts and applies no floor itself, leaving gating to
pseudobulk. Descriptive; this table is the power budget of the paired contrast.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_annotate_states.py` | `main` | `thresholds.pseudobulk_min_cells = 20`, `thresholds.pseudobulk_min_donors = 3`, `decisions.qc.drop_gsms = [GSM4859852]` | `03_results/objects/01_qc.h5ad` |

## tables/substate_markers.csv

FOXP3 reaches mean log-normalised 1.69 in 80% of sorted Tregs, against 0.06 in Tcon and 0.02 in
CD8. CD8A marks 93% of CD8 cells and under 0.5% of cells in either other gate. IL7R runs highest in
Tcon (2.35) and lowest in Treg (0.66), as a CD127-lo gate requires.

**How to read:** 36 rows — 12 panel markers × 3 frozen labels. `mean_lognorm` is the mean
log-normalised expression across every cell carrying that label, zeros included, so it is diluted
by non-expressers. `frac_expressing` is the fraction with more than 0 UMI. This table carries no
test statistic and no sign convention; the pass condition is high in the intended label and low
elsewhere. Source table for the marker dotplot. QC overlay tier — hand markers, which carry no
evidential weight.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_annotate_states.py` | `main` | `PANEL_MARKERS` (12 canonical lineage markers) is a literal in the script — no config key; labels from `decisions.treg_gate.basis = sorting` | `03_results/objects/01_qc.h5ad` |
