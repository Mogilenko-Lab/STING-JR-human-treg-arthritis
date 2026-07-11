# 07_embedding — artifact captions

## figures/_overview/umap_signatures_treg.png

Previews the multi-hook harvest design: where each candidate signature
we are drafting for (mouse WT_heat anchor, effector-Treg, heat-
shock/stress) lands ON ITS OWN across the sorted Treg/Tcon/CD8
embedding, next to the frozen sort-lineage reference. VISUALISATION
only — not the statistical evidence.

**How to read:** Three continuous panels colour every cell by one candidate signature
(viridis, clipped 2-98th pct; high on top); the fourth shows the
frozen sort lineage. Read as WHERE each signature concentrates on the
embedding, alone. The WT_heat anchor is shown as an ANNOTATION, never
a selection gate; the source table gives per-lineage medians.
Correlative preview of the harvest design, not an effect size.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/07_embedding_viz.py` | `main` | `signatures: WT_heat_up (anchor, annotation-only), score_eTreg, score_HSP; cmap=viridis` | `03_results/07_embedding/tables/hook_factor_substrate.parquet` |

## figures/_overview/umap_or_union_treg.png

The 'where they land together' panel: cells coloured by which anchor-
orthogonal hook-factor(s) admit them under the bounded OR-union
(lineage / effector / mt-hi viable / multiple), with the matched heat-
lo and effector-lo baseline regions shown. Previews the harvest
design's factorial contrastability — VISUALISATION, not selection-to-
file and not a statistical readout.

**How to read:** Left: each cell coloured by the hook(s) it satisfies (grey = not in
union). Right: the matched-lo baselines (heat-lo, effector-lo) that
give every factorial contrast a defined negative arm. The union is a
BOUNDED minority of cells (fraction in the title + source table) — the
concern that OR sweeps in the whole dataset does not hold here. No
cells are lassoed/subset; harvest selection is deferred. Correlative.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/07_embedding_viz.py` | `main` | `OR-union = hook_lineage | hook_effector | hook_mthi_viable; anchor score NEVER a disjunct; union = 33.1% of all cells (bounded minority)` | `03_results/07_embedding/tables/hook_factor_substrate.parquet` |

## figures/_overview/umap_markers_treg.png

Curated POI lineage markers (FOXP3/IL2RA/CTLA4/IKZF2 Treg identity;
CD8A; IL7R) that define the anchor-orthogonal lineage/marker hooks,
next to the frozen sort-lineage reference — confirming the hooks track
real lineage biology, not the anchor score.

**How to read:** Each panel colours cells by one marker's log-normalised expression
(viridis, clipped 2-98th pct); last panel is the frozen sort lineage.
Read as: the Treg-identity markers concentrate on the Treg gate, CD8A
on CD8, IL7R on conventional T — the curated hooks are lineage-
faithful and independent of the anchor score. Source table = per-
lineage median expression + fraction expressing. Correlative.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/07_embedding_viz.py` | `main` | `markers: FOXP3, IL2RA, CTLA4, IKZF2, CD8A, IL7R; cmap=viridis (log-norm expression)` | `03_results/07_embedding/tables/hook_factor_substrate.parquet` |

## figures/_overview/umap_drafted_treg.png

The explicit 'these are the cells our rules would harvest' view: the
Treg-compartment UMAP with the drafted subset (hook_or_union == TRUE,
33,113 cells / 33%) highlighted against the light-grey non-drafted
background, plus a facet of the drafted cells by which hook drafted
them. The subset the multi-hook rules would draft (design preview /
visualization, not a committed cohort or statistical evidence); anchor
heat score is annotation, never a selection gate.

**How to read:** Left (primary): every cell coloured by whether the OR-union harvest
rules would draft it (vermillion = drafted, light grey = background) —
the concrete subset preview. Right: the drafted cells split by which
hook admitted them (lineage / effector / mt-hi viable / multiple);
non-drafted cells stay grey. This is a VISUALISATION of the harvest
design, not a committed cohort and not a statistical readout; no cells
are lassoed/subset to file and the anchor heat score is never a
selection gate. Source table = drafted-vs-background split + per-hook
counts. Correlative.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/07_embedding_viz.py` | `main` | `drafted = hook_or_union (lineage|effector|mt-hi viable); 33113 / 99915 cells (33.1%); highlight=#D55E00, background=#D9D9D9` | `03_results/07_embedding/tables/hook_factor_substrate.parquet` |

## figures/_overview/umap_annotation_treg.png

Annotation-overview UMAP: every cell coloured by its frozen sort-
lineage annotation (coarse_label = Treg / Tcon / CD8), with a tissue
(synovial-fluid / peripheral-blood) companion. Establishes WHAT IS ON
THE MAP before any signature/hook overlay is read against it.
VISUALISATION only — the frozen annotation, not a statistical readout.

**How to read:** Left: the frozen sorted-lineage annotation the whole atlas is built on
(Treg green, Tcon orange, CD8 pink). Right: the same cells coloured by
tissue of origin (synovial fluid vs peripheral blood). Read as the
reference layout every downstream signature and hook overlay sits on
top of. Source table = per-annotation x tissue cell counts.
Correlative, annotation only.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/07_embedding_viz.py` | `main` | `colour=coarse_label (Okabe-Ito Treg/Tcon/CD8); companion=tissue (SF/PB); categorical scatter, baseline levels drawn first` | `03_results/07_embedding/tables/hook_factor_substrate.parquet` |

## figures/_overview/umap_quadmarkers_treg.png

Quadruple context patch: mitochondrial fraction plus the three
canonical person-of-interest / disease genes for JIA sorted Tregs —
FOXP3 (master Treg-lineage transcription factor), CTLA4 (suppressive
effector molecule) and IKZF2/Helios (stable/thymic-Treg identity
marker). Chosen as the core Treg identity + suppressive-activation
axis, the most-motivated trio for this compartment; %mt anchors
QC/viability context. VISUALISATION only.

**How to read:** Four continuous panels (viridis, clipped 2-98th pct, high on top):
top-left is mitochondrial fraction (QC/viability context); the other
three are the canonical Treg genes — FOXP3 and CTLA4 concentrate on
the Treg gate, IKZF2/Helios marks the stable/thymic-Treg fraction.
Read as the marker context the drafted subset is designed against.
Source table = per-lineage median values (+ frac expressing for genes)
from the compute stage. Correlative, annotation only.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/07_embedding_viz.py` | `main` | `channels: pct_counts_mt + FOXP3, CTLA4, IKZF2; continuous viridis, clipped 2-98th pct; genes = core Treg identity/activation (FOXP3 TF, CTLA4 effector, IKZF2/Helios stable-Treg)` | `03_results/07_embedding/tables/hook_factor_substrate.parquet` |
