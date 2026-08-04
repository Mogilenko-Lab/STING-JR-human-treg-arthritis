# 07_embedding: artifact captions

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
BOUNDED minority of cells (fraction in the title + source table), so
an OR sweep over these hooks does not take in most of the dataset. No
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

## tables/hook_factor_definitions.csv

Every anchor-orthogonal hook resolves to a bounded minority of the
99,915 cells — sort-lineage 27,175 (27.2%), effector-high 9,992
(10.0%, at score_eTreg >= 0.1385), and the mt-hi-but-viable Treg
pocket 340 (0.34%) — so their OR-union reaches only 33,113 cells
(33.1%) and two thirds of the atlas remains available as baseline.

**How to read:** One row per factor. `kind` separates the three selection hooks from
the OR-union, the annotation-only anchor factors, and the matched-lo
baselines; `definition` and `threshold` give the rule and the
cut-point it resolved to (effector P90 = 0.1385; within-Treg %mt
P97.5 = 10.03 plus an eTreg-median and 400-gene viability floor);
`n_cells` and `frac_all_cells` size it. Membership is `POI =
(hook_lineage OR hook_effector OR hook_mthi_viable) AND hook_viable`.
The two `anno_heat_*` rows are the mouse WT_heat score (P90 hi / P10
lo) carried as ANNOTATION ONLY and absent from that rule — the anchor
never selects. `anno_stingspecific` is an empty placeholder (n=0).
Hypothesis-generating harvest design, firewalled from the confirmatory
track; no cells are written out here.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/07_embedding.py` | `main` | `qc_min_genes = 200` (viability floor = 2x = 400 genes); module constants `EFFECTOR_Q=0.90`, `MT_PCTILE=97.5`, `HEAT_HI_Q=0.90`, `HEAT_LO_Q=0.10`, `EFFECTOR_LO_Q=0.10` | `03_results/interactive/01_qc_explore.parquet`, `03_results/interactive/02_annotation_explore.parquet` |

## tables/or_union_membership.csv

The union is carried by the sort gate but is not reducible to it: of
33,113 union cells, 22,879 (69.1%) enter on lineage alone, 3,956
(11.9%) on lineage and effector together, 5,938 (17.9%) are
effector-high NON-Treg cells the sort gate would miss entirely, and
340 (1.0%) are the mt-hi viable pocket; 66,802 cells (66.9% of all)
stay outside as baseline.

**How to read:** One row per membership category, each cell assigned to exactly one
(the mt-hi viable pocket is reported first, since it nests inside the
lineage hook). `n_cells` with `frac_all_cells` sizes the category
against the whole atlas; `frac_of_union` sizes it against the union
only and is blank for the baseline row. `union_frac_all_cells` repeats
the headline bound. Read the "effector only" row as the marginal value
of the effector hook over the FACS gate. The mouse WT_heat score is
not a disjunct in any of these categories — it annotates, orthogonal
biology selects. Hypothesis-generating harvest design, firewalled from
the confirmatory track; not an effect size.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/07_embedding.py` | `main` | `hook_or_union = hook_lineage OR hook_effector OR hook_mthi_viable`; `EFFECTOR_Q=0.90`, `MT_PCTILE=97.5`, `qc_min_genes = 200` | `03_results/interactive/01_qc_explore.parquet`, `03_results/interactive/02_annotation_explore.parquet` |

## tables/hook_per_lineage_summary.csv

The effector hook draws unevenly across the frozen sort lineages —
14.9% of Treg (4,054/27,175) and 13.0% of CD8 (4,543/34,827) but only
3.7% of Tcon (1,395/37,913) — while the anchor annotation is nearly
flat (median WT_heat_up -0.060 Treg, -0.052 Tcon, -0.047 CD8), which
is exactly why the hooks and not the anchor score do the selecting.

**How to read:** One row per frozen `coarse_label`. For each hook and annotation
factor there is a paired `n_` count and `frac_` within-lineage
fraction; `median_` columns give the lineage medians of the three
candidate signatures (WT_heat_up, score_eTreg, score_HSP) plus
`pct_counts_mt` and `n_genes_by_counts`. `frac_hook_lineage` is 1.0
for Treg and 0 elsewhere by construction, and the mt-hi viable pocket
is Treg-only (340 cells, 1.25% of Treg). Higher `median_score_eTreg`
means a more effector-like lineage (CD8 +0.0429 > Treg -0.0031 > Tcon
-0.0229). The `anno_heat_*` columns are annotation, never a gate.
Descriptive summary of the harvest design, correlative and
hypothesis-generating only.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/07_embedding.py` | `main` | `SIGNATURE_COLS = ["WT_heat_up", "score_eTreg", "score_HSP"]`; `EFFECTOR_Q=0.90`, `MT_PCTILE=97.5`, `HEAT_HI_Q=0.90`, `HEAT_LO_Q=0.10`, `qc_min_genes = 200` | `03_results/interactive/01_qc_explore.parquet`, `03_results/interactive/02_annotation_explore.parquet` |

## tables/signatures_per_lineage.csv

The three candidate signatures separate the lineages in different
directions — score_eTreg is highest in CD8 (median +0.0429) and
score_HSP highest in Treg (+0.0437 vs Tcon +0.0232, CD8 -0.0185) —
whereas WT_heat_up medians sit flat and negative across all three
(-0.060 to -0.047), so the anchor score carries almost no lineage
information on this map.

**How to read:** One row per frozen lineage; source table for the signature-overlay
UMAP. Each signature contributes a `median_` and a `mean_` column, so
a mean far above the median flags a skewed minority rather than a
shifted lineage — Treg score_eTreg is the clear case (median -0.0031,
mean +0.0240), i.e. a small effector-high tail, which is what the
effector hook is designed to catch. `pct_counts_mt` is carried for
viability context. Scores are the per-cell values frozen upstream; a
higher value means stronger signature expression. The WT_heat_up
column is the mouse anchor as ANNOTATION ONLY — it is never a
selection predicate, and these medians are descriptive, not an effect
size. Hypothesis-generating tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/07_embedding.py` | `main` | `SIGNATURE_COLS = ["WT_heat_up", "score_eTreg", "score_HSP"]` (+ `pct_counts_mt`) | `03_results/interactive/01_qc_explore.parquet`, `03_results/interactive/02_annotation_explore.parquet` |

## tables/markers_per_lineage.csv

The curated marker hooks track real lineage biology and not the anchor
score: FOXP3 is detected in 79.7% of sorted Treg against 3.8% of Tcon
and 1.0% of CD8, CD8A in 92.7% of CD8 against 0.2% of Treg, and IL7R
in 87.7% of Tcon against 42.7% of Treg, with IL2RA (66.8%), CTLA4
(65.2%) and IKZF2 (52.7%) filling in the Treg identity/suppression
axis.

**How to read:** One row per lineage x marker gene (3 x 6 = 18 rows); source table for
the marker-overlay UMAP. `median_expr` and `mean_expr` are
log-normalised expression, `frac_expressing` is the fraction of cells
with any detected count — for sparse markers the median is 0 and
`frac_expressing` is the informative column. Read down a gene to check
that the marker concentrates where the sort gate says it should; the
Treg-identity genes, CD8A and IL7R each mark their own compartment,
confirming these hooks are anchor-orthogonal by construction. Nothing
here is gated on the mouse WT_heat score. Descriptive annotation,
hypothesis-generating tier, not an effect size.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/07_embedding.py` | `main` | `MARKER_GENES = ["FOXP3", "IL2RA", "CTLA4", "IKZF2", "CD8A", "IL7R"]` | `03_results/interactive/01_qc_explore.parquet`, `03_results/interactive/02_annotation_explore.parquet` |
