# 07_embedding — The population-of-interest harvest design, previewed

A bounded subset of cells is drafted for later cross-dataset and cross-species questions. This
stage defines the rule that would draft them, measures how much of the atlas it reaches, and
draws where those cells sit. **No cells are written out and no cohort is committed.**

The rule is OR-gated over three anchor-orthogonal hooks and one viability gate:

```
POI = (hook_lineage ∨ hook_effector ∨ hook_mthi_viable) ∧ hook_viable
```

Each hook resolves to a bounded minority of the 99,915 cells — sort-lineage 27,175 (27.2%),
effector-high 9,992 (10.0%), the mitochondrial-high viable Treg pocket 340 (0.34%) — so their
union reaches 33,113 cells (33.1%) and two thirds of the atlas stays available as baseline.

**The mouse anchor score is annotation and never a selection predicate.** `WT_heat_up` is carried
on every panel and stays outside the membership rule. Asking the anchor score to both select
cells and measure enrichment in them is the circularity this design exists to avoid, and the
hooks are orthogonal biology for exactly that reason.

**The rule is frozen as implemented.** No further hook expansion. Everything here is
hypothesis-generating, firewalled from the confirmatory spine, and writes no effect-size row.

## Signature provenance

The three candidate signatures overlaid on these panels come from outside the JIA data.

| Signature | Origin | How the list was derived |
|---|---|---|
| `WT_heat_up` | **GSE329522**, this project's own mouse anchor. No paper reference recorded. | Bulk RNA-seq of induced regulatory T cells from primary murine splenic CD4⁺ T cells, genotype (WT, cGAS-KO) × temperature (37 °C, 39 °C), 5 biological replicates per group over 20 libraries. The WT 39 °C-against-37 °C up arm, projected to human orthologs with pinned offline babelgene. |
| `score_eTreg` | **GSE161426**, 26 bulk RNA-seq samples of sorted CD4 populations. Mijnheer / Lutter et al. 2021, *Nature Communications*, PMID 33976194, doi 10.1038/s41467-021-22975-7. | Derived here as a synovial-fluid Treg against peripheral-blood Treg contrast on the deposited log2-normalised matrix `GSE161426_Gene_expression_table_log2.xlsx`. GEO carries a matrix alone for this series. |
| `score_HSP` | This compartment's own hand marker module. | A small curated heat-shock marker panel. |

`anno_stingspecific` is an empty placeholder here. The set it would carry is the published
21-gene interferon-independent STING signature of de Cevins et al. 2023, *Cell Reports Medicine*,
PMID 38118407, Table S6, from the SAVI PBMC cohort **GSE226598**.

---

## Figures

Six panels, all on the same UMAP frame over all 99,915 cells.

### `figures/_overview/umap_annotation_treg.png`

**The reference layout every overlay is read against.**
Left, cells coloured by frozen sort lineage (Treg green, Tcon orange, CD8 pink). Right, the same
cells coloured by tissue of origin. This establishes what is on the map before any signature or
hook colouring sits on top of it. The source table gives per-annotation × tissue cell counts.
*Source* `tables/hook_factor_substrate.parquet` · `02_analysis/scripts/07_embedding_viz.py`.

### `figures/_overview/umap_markers_treg.png`

**The curated markers behind the lineage hooks.**
Six panels, each colouring cells by one marker's log-normalised expression on viridis with the
highest drawn last, plus the frozen sort lineage in the last slot. **All six genes share one
clip** — the 2nd-to-98th percentile pooled over the six — and one bar, so brightness compares
between the genes as well as within one. On separate clips IKZF2 at 0–1.88 and IL7R at 0–3.87
would look alike.

The Treg identity markers concentrate on the Treg gate, CD8A on CD8, IL7R on conventional T. The
hooks track lineage biology and stand independent of the anchor score.
*Source* `tables/hook_factor_substrate.parquet` · `02_analysis/scripts/07_embedding_viz.py`.

### `figures/_overview/umap_quadmarkers_treg.png`

**Treg identity and viability context in four panels.**
Top-left, mitochondrial fraction on its own clip, because a percentage is a different quantity
from expression. The three gene panels — FOXP3 the lineage transcription factor, CTLA4 the
suppressive effector, IKZF2/Helios the stable-Treg marker — share one pooled 2nd-to-98th clip, so
their brightness compares. These three are the core Treg identity and suppressive-activation axis
for this compartment.
*Source* `tables/hook_factor_substrate.parquet` · `02_analysis/scripts/07_embedding_viz.py`.

### `figures/_overview/umap_signatures_treg.png`

**Where each candidate signature lands on its own.**
Three continuous panels colouring every cell by one candidate signature — the mouse `WT_heat_up`
anchor, `score_eTreg`, `score_HSP` — on viridis clipped to the 2nd-to-98th percentile with the
highest drawn last, plus the frozen sort lineage in the fourth slot. Read as **where** each
signature concentrates, alone. The anchor is drawn as annotation, and the source table gives
per-lineage medians.
*Source* `tables/hook_factor_substrate.parquet` · `02_analysis/scripts/07_embedding_viz.py`.

### `figures/_overview/umap_or_union_treg.png`

**Where the hooks land together, and where the matched baselines sit.**
Left, each cell coloured by which hook or hooks admit it under the bounded OR-union — lineage,
effector, mitochondrial-high viable, or multiple — with grey outside the union. Right, the
matched-lo baselines (heat-lo, effector-lo) that give every factorial contrast a defined negative
arm. The union is a bounded minority at 33.1%, so an OR sweep over these hooks leaves most of the
dataset outside it.
*Source* `tables/hook_factor_substrate.parquet` · `02_analysis/scripts/07_embedding_viz.py`.

### `figures/_overview/umap_drafted_treg.png`

**The concrete subset the rules would draft.**
Left, every cell coloured by whether the OR-union rules would draft it: vermillion drafted
(33,113 cells, 33.1%), light grey background. Right, the drafted cells split by which hook
admitted them, with non-drafted cells left grey. This is a design preview, and no cohort is
committed.
*Source* `tables/hook_factor_substrate.parquet` · `02_analysis/scripts/07_embedding_viz.py`.

---

## Tables

### `tables/hook_factor_definitions.csv`

One row per factor. `kind` separates the three selection hooks from the OR-union, the
annotation-only anchor factors and the matched-lo baselines. `definition` and `threshold` give
the rule and the cut point it resolved to — effector P90 = 0.1385, and within-Treg %mt P97.5 =
10.03 plus an eTreg-median and 400-gene viability floor. `n_cells` and `frac_all_cells` size it.

The two `anno_heat_*` rows are the mouse `WT_heat` score at P90 hi and P10 lo, carried as
annotation and absent from the membership rule. `anno_stingspecific` is an empty placeholder
(n = 0).

### `tables/or_union_membership.csv`

One row per membership category, each cell assigned to exactly one. The mitochondrial-high viable
pocket is reported first, since it nests inside the lineage hook.

**The sort gate carries the union, and the other two hooks extend it.** Of 33,113 union cells,
22,879 (69.1%) enter on lineage alone, 3,956 (11.9%) on lineage and effector together, 5,938
(17.9%) are effector-high **non-Treg** cells the sort gate would miss entirely, and 340 (1.0%) are
the mitochondrial-high viable pocket. 66,802 cells (66.9% of all) stay outside as baseline.

`frac_of_union` sizes a category against the union and is blank for the baseline row. Read the
"effector only" row as the marginal value of the effector hook over the FACS gate.

### `tables/hook_per_lineage_summary.csv`

One row per frozen `coarse_label`. Each hook and annotation factor contributes a paired `n_`
count and `frac_` within-lineage fraction, and the `median_` columns give lineage medians of the
three candidate signatures plus `pct_counts_mt` and `n_genes_by_counts`.

The effector hook draws unevenly across the lineages — 14.9% of Treg (4,054/27,175) and 13.0% of
CD8 (4,543/34,827) against 3.7% of Tcon (1,395/37,913) — while the anchor annotation runs nearly
flat (median `WT_heat_up` −0.060 Treg, −0.052 Tcon, −0.047 CD8). That contrast is why the hooks
do the selecting. `frac_hook_lineage` is 1.0 for Treg and 0 elsewhere by construction, and the
mitochondrial-high viable pocket is Treg-only (340 cells, 1.25% of Treg).

### `tables/signatures_per_lineage.csv`

One row per frozen lineage, and the source of the signature-overlay panel. Each signature
contributes a `median_` and a `mean_` column, so a mean far above the median flags a skewed
minority inside an otherwise unshifted lineage. Treg `score_eTreg` is the clear case (median
−0.0031, mean +0.0240), a small effector-high tail, which is what the effector hook is built to
catch.

The three signatures separate the lineages in different directions: `score_eTreg` runs highest in
CD8 (median +0.0429) and `score_HSP` highest in Treg (+0.0437 against Tcon +0.0232, CD8 −0.0185),
while `WT_heat_up` medians sit flat and negative across all three (−0.060 to −0.047). The anchor
score carries almost no lineage information on this map.

### `tables/markers_per_lineage.csv`

18 rows — three lineages × six marker genes. `median_expr` and `mean_expr` are log-normalised
expression and `frac_expressing` the fraction with any detected count. For a sparse marker the
median reads 0 and `frac_expressing` is the informative column.

FOXP3 is detected in 79.7% of sorted Treg against 3.8% of Tcon and 1.0% of CD8. CD8A marks 92.7%
of CD8 against 0.2% of Treg, and IL7R 87.7% of Tcon against 42.7% of Treg. IL2RA (66.8%), CTLA4
(65.2%) and IKZF2 (52.7%) fill in the Treg identity and suppression axis. Read down a gene to
confirm the marker concentrates where the sort gate says it should.

### `tables/hook_factor_substrate.parquet`

The per-cell substrate behind all six panels: one row per cell with the embedding coordinates,
the frozen labels, the hook membership booleans, the three candidate signature scores and the six
marker expression columns. Regenerable and untracked.
