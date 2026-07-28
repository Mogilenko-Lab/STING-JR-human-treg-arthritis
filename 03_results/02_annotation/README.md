# 02_annotation — artifact captions

_**Abbreviations:** SF = synovial fluid (inflamed joint); PB = peripheral blood. The cohort contains 7 JIA donors, of whom 6 span both arms in each analyzed population after QC. Treg = CD4⁺CD127ˡᵒCD25⁺ regulatory; Tcon = CD4⁺CD25⁻ conventional; CD8 = CD8⁺CD45RO⁺ memory._

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

The cohort contains 7 donors, and 6 span SF and PB in each analyzed
population after QC. Every observed stratum clears the pseudobulk cell
floor; p3 PB Tcon/CD8 are absent by design.

**How to read:** Heatmap of cells per donor (x) x label+tissue (y); red * marks a
stratum below the pseudobulk cell floor; empty = intentionally-absent
sample. Donor count per arm determines contrast precision. Diagnostic.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_annotate_states_viz.py` | `main` | `thresholds.pseudobulk_min_cells = 20` | `03_results/02_annotation/tables/counts_donor_by_label_tissue.csv` |

## tables/confusion_sort_vs_predicted.csv

Marker-module argmax reproduces the sort gate for 85.6% of the 99,915
QC-passing cells (Treg 89.6%, Tcon 83.8%, CD8 84.4%); the off-diagonal
mass is overwhelmingly Tcon-vs-CD8, and only 195 sorted Tregs land on
the CD8 module — the Treg gate is clean enough to freeze.

**How to read:** Rows = frozen sort label, columns = argmax of the three z-scored
canonical lineage modules; entries are cell counts and the diagonal is
agreement. Rows and columns are forced to the same Treg/Tcon/CD8
order, so the diagonal reads directly. An off-diagonal cell is module
disagreement, not a demonstrated mis-sort: each module is 4-5 genes,
and shared cytotoxic/activation genes blur Tcon against CD8.
Diagnostic tier — sorting remains the label of record.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_annotate_states.py` | `main` | `decisions.treg_gate.basis = sorting`; `LINEAGE_MODULES` (Treg FOXP3/IL2RA/CTLA4/IKZF2/TIGIT, Tcon IL7R/CD40LG/ANK3/LEF1, CD8 CD8A/CD8B/NKG7/GZMK/GZMA) is a literal in the script — no config key | `03_results/objects/01_qc.h5ad` |

## tables/counts_donor_by_label_tissue.csv

All 39 surviving donor x label x tissue strata clear the 20-cell
pseudobulk floor — the thinnest, PB Treg p7, still holds 266 cells —
and every arm keeps at least 6 donors, so each SF-vs-PB contrast sits
well above the donor floor.

**How to read:** One row per donor x frozen coarse label x tissue; `n_cells` counts
cells surviving QC. Read each row against the per-stratum floor
(`pseudobulk_min_cells = 20`) and each label x tissue arm against the
donor floor (`pseudobulk_min_donors = 3`). Absent rows are not zeros:
PB Tcon/CD8 for p3 were never collected, and SF Treg p5 is the
near-empty library QC excluded. The compute step applies no floor
itself — it writes raw counts and leaves gating to pseudobulk.
Descriptive; this table is the power budget of the paired contrast.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_annotate_states.py` | `main` | `thresholds.pseudobulk_min_cells = 20`, `thresholds.pseudobulk_min_donors = 3`, `decisions.qc.drop_gsms = [GSM4859852]` | `03_results/objects/01_qc.h5ad` |

## tables/substate_markers.csv

FOXP3 reaches mean lognorm 1.69 in 80% of sorted Tregs against 0.06 /
0.02 in Tcon / CD8, CD8A marks 93% of CD8 cells and under 0.5%
elsewhere, and IL7R runs highest in Tcon (2.35) and lowest in Treg
(0.66) exactly as a CD127-lo gate requires.

**How to read:** 36 rows = 12 panel markers x 3 frozen labels. `mean_lognorm` is the
mean log-normalized expression across every cell carrying that label
(zeros included, so it is diluted by non-expressers);
`frac_expressing` is the fraction with >0 UMI. No test statistic and
no sign convention here — high in the intended label and low elsewhere
is the pass. Source table for the marker dotplot. QC overlay tier
(hand markers, not evidence).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_annotate_states.py` | `main` | `PANEL_MARKERS` (12 canonical lineage markers) is a literal in the script — no config key; labels from `decisions.treg_gate.basis = sorting` | `03_results/objects/01_qc.h5ad` |
