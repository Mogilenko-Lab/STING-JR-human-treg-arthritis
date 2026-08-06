# 17_treg_reembedding — The same programs, on the Treg gate's own map

A second map of the sorted JIA T-cell compartment (GSE160097), computed on the Treg gate alone.
The full-object embedding lays out all three sorted lineages in one geometry, so the variance
separating Treg from Tcon from CD8 sets its axes and Treg substructure occupies the room left
over. Here the 27,175 cells of the frozen `coarse_label == "Treg"` gate get the whole canvas.

**The question each counterpart figure answers.** Does a program's apparent structure on the
full-object map need the Tcon and CD8 gates to appear?

**What was computed.** The full-object recipe rerun on the subset, with the neighbours-and-UMAP
tail run twice off the same PCA — once raw, once after Harmony over donor — plus a
neighbourhood-mixing table measuring what the correction did.

**What was drawn.** Four six-panel strips and one three-panel strip: the reference layout, the
mouse arms, the curated program family, the candidate harvest signatures, and one stacked
patchwork.

Per-cell scores are joined in from
[`../interactive/16_narrative_embedding.parquet`](../interactive/) on barcode and are never
recomputed, so every map can be coloured by identical values and a difference between two maps is
a difference of layout.

**Tier.** Annotation and visualisation only. An embedding places cells. Proximity on any of these
maps carries no evidential weight, and this stage writes no effect-size row.

## The recipe, and the correction

The six calls of the full-object embedding, in the same order, on the subset:
`highly_variable_genes(n_top_genes=3000)` → subset → `scale(max_value=10)` → `pca(n_comps=30)` →
`neighbors(n_pcs=30)` → `umap`, at `random_state = 0`. `hvg_n_top` and `n_pcs` come from the same
config block the full-object run read, and `n_neighbors` stays at the scanpy default, so the two
maps differ by which cells went in and by the correction alone.

On this subset the uncorrected recipe resolves donor-by-tissue sample of origin at 66.1%
same-donor neighbours. So the tail runs twice, and `x`/`y` are the Harmony-corrected coordinates
and the pair to draw. `x_uncorrected`/`y_uncorrected` are carried for comparison.

**What the correction did, at k = 30** (`tables/treg_reembedding_mixing.csv`):

| Factor | Uncorrected | Corrected | Chance | Full object, same cells |
|---|---|---|---|---|
| same-donor neighbours | 0.661 | 0.201 | 0.146 | 0.420 |
| same-tissue neighbours | 0.975 | 0.923 | 0.500 | 0.955 |

The corrected map sits at 1.38× chance on donor. Every donor drops, and the residual
concentrates in JIA_patient_5 at 0.271. Tissue holds at 1.85× chance, so donor and tissue are
separable in this representation even though the two are nested at the sample level — each GSM
is one donor × one tissue × one sorted population. The corrected coordinates are usable for
reading Treg substructure. A small tight cluster still deserves a donor-composition check, and
Harmony reshapes the space it corrects, so this map remains a map.

## Why the counterpart figures share a colour scale

Three figures here are Treg-only counterparts of program rows drawn on the full object:
`umap_treg_arms` against `16_narrative_scoring/umap_full_arms`, `umap_treg_programs` against
`16_narrative_scoring/umap_full_programs`, and `umap_treg_signatures` against
`07_embedding/umap_signatures_treg`. Each pair answers one question, and the answer is readable
only when both halves sit on one colour scale.

So every counterpart panel is clipped to its full-object twin's limits, and the arm pair's
shared limits are asserted at run time, because a pair that has silently stopped sharing a scale
still looks like a matched set. The consequence to hold onto: a washed-out counterpart panel is
a statement about the Treg gate's range, and the published STING panel is the clear case.

Each pair shares its cells and its colour scale. The coordinates differ, which is the point.

**Six panels to a strip, one geometry.** Every strip is drawn through
`02_analysis/helpers/umap_grid.py`, the same module
[`../16_narrative_scoring/`](../16_narrative_scoring/) uses, so a strip from either stage stacks
against a strip from the other with a column reading top to bottom. `umap_treg_signatures` is the
one three-panel strip, at the same panel size.

**The reference strip drops a panel and gains one.** Its full-object twin's second slot is the
Treg/Tcon/CD8 sort gate, which is a single category here, so `IKZF2` takes that slot. Tissue and
donor stay the bookends in both strips, so the pair lines up where it means something.

---

## Figures

### `figures/_overview/umap_treg_reembedding.png`

**The Treg-only layout, and the identity genes across it.**
Six panels over one frame of the same 27,175 cells, sharing one bounding box. Tissue and donor
are the bookends, both drawn in shuffled order so overlapping groups paint evenly. Between them
sit four Treg identity genes in log-normalised expression: FOXP3 the lineage transcription
factor, IL2RA the CD25 chain, CTLA4 the suppressive effector, IKZF2 the Helios subset marker.
All four pool onto one clip and one bar in real units (0.00 to 2.75), and that clip is the
full-object figure's, so a panel that looks brighter here carries higher expression.

Synovial-fluid and paired-blood cells occupy distinct territory — 0.923 same-tissue neighbours
at k = 30 against 0.500 expected — after the same-donor fraction has fallen from 0.661 to 0.201
against 0.146 expected. The four identity genes hold across the whole map (FOXP3 mean 1.2274 to
2.1590, IKZF2 0.6467 to 0.7041), so the layout separates tissue while the gate stays uniform.
Harmony reshapes the space it corrects, so this map is annotation.
*Source* `tables/_overview/umap_treg_reembedding.csv` ·
`02_analysis/scripts/17_treg_reembedding_viz.py`.

### `figures/_overview/umap_treg_arms.png`

**The same three mouse arms and three lenses, on the Treg gate alone.**
The Treg-only counterpart of `16_narrative_scoring/umap_full_arms.png`: same six sets, same
order, same geometry and colormap, on the 27,175 sorted Treg cells. The vertical rules carry the
provenance — three anchor-dependent arms on the left, two curated versioned lenses in the
middle, `eTreg_up` ruled off on the right. One bar serves the row and reads 0 to 1, each panel
rescaled across its full-object twin's clip, so the twin comparison holds panel for panel.

All three mouse 39 °C-derived up arms still colour synovial territory brighter than paired blood
once Tcon and CD8 are gone, and so do all three anchor-independent lenses. Per-cell AUCell means
run `WT_heat_up` 0.0111 to 0.0188, `KO_heat_up` 0.0108 to 0.0188 and `Interaction_up` 0.0717 to
0.1245 for the arms, against `HALLMARK_HYPOXIA` 0.0718 to 0.0950, `HSR_core` 0.0989 to 0.1075
and `eTreg_up` 0.0190 to 0.0668 for the lenses. So inside the Treg gate the arm's tissue
colouring stays shared with the lenses beside it, exactly as on the whole object.

`WT_heat_up` and `KO_heat_up` share 182 genes and one clip, so those two compare pixel for
pixel. Cells pool across donors, so a tissue difference is pseudoreplicated.
*Source* `tables/_overview/umap_treg_arms.csv` ·
`02_analysis/scripts/17_treg_reembedding_viz.py`.

### `figures/_overview/umap_treg_programs.png`

**The same six curated lenses, on the Treg gate alone.**
The Treg-only counterpart of `16_narrative_scoring/umap_full_programs.png`, on the same rule: the
cGAS-STING and type-I interferon family left, the inflammation and activation programs right. A
synovial-high colouring shared by both families is generic inflammation. Each panel keeps its own
set's limits, and those limits come from the full object, which makes each panel comparable to
its own twin.

All six still colour synovial territory brighter than paired blood once the other two sort gates
are gone: `ifn_generic_axis` 0.0172 to 0.0290, `HALLMARK_INTERFERON_ALPHA_RESPONSE` 0.1306 to
0.1552, `HALLMARK_TNFA_SIGNALING_VIA_NFKB` 0.0503 to 0.0589, the published STING signature 0.0220
to 0.0424.

The shared clip carries the level that the rescaled bar hides. The published STING panel is scored
over the whole Treg map against a full-object limit it barely reaches: Treg blood mean 0.0220
against 0.0547 in Tcon and 0.0672 in CD8, with a median of exactly 0.013 — the only row of the
six at zero. At least half of Treg blood cells score zero on that set, so its Treg tissue
difference rests partly on a zero-inflated baseline. The set is 21 genes, 18 of them scored
here, and its own IFN-β validation is underpowered at three donors.
*Source* `tables/_overview/umap_treg_programs.csv` ·
`02_analysis/scripts/17_treg_reembedding_viz.py`.

### `figures/_overview/umap_treg_signatures.png`

**The three candidate harvest signatures, on the Treg gate alone.**
The Treg-only counterpart of `07_embedding/umap_signatures_treg.png`, which draws these same
three channels across all three gates. Two differences matter. The unit: these are scanpy
`score_genes` module scores, mean-centred against a sampled background and signed, so a value
near zero sits at background, and they share no scale with the AUCell panels. The panel count:
the full-object twin carries a fourth sort-gate reference panel, which is a single category here.

Within the Treg gate the effector-Treg and heat-shock module scores still separate synovial
fluid from paired blood (`score_eTreg` −0.0539 to 0.1021, `score_HSP` 0.0270 to 0.0806),
alongside the mouse anchor annotation (`WT_heat_up` −0.0699 to −0.0475). `WT_heat_up` here is
carried as annotation only and is never a selection predicate. Limits come from all 99,915 cells
of the embedding substrate.

The bounded AUCell readings of the same two programs are the `WT_heat_up` and `eTreg_up` panels
of `umap_treg_arms`.
*Source* `tables/_overview/umap_treg_signatures.csv` ·
`02_analysis/scripts/17_treg_reembedding_viz.py`.

### `figures/_overview/umap_treg_patchwork.png`

**The Treg-only reference layout above its arm colouring.**
The two strips this stage ships separately, stacked at identical panel size so a column reads top
to bottom. Nothing new is drawn, and both rows hold the identical frame of cells at identical
coordinates. The top row is categorical annotation plus log-normalised expression, the bottom row
per-cell AUCell rescaled per panel onto one 0-to-1 bar.

The arms brighten the synovial-fluid territory the tissue panel marks, while the identity genes
stay uniform across it.
*Source* `tables/_overview/umap_treg_patchwork.csv` ·
`02_analysis/scripts/17_treg_reembedding_viz.py`.

---

## Tables

### `tables/treg_reembedding_mixing.csv` — what the correction did

One row per (`embedding` × `space` × `grouping_key` × `group`). `embedding` names which map was
measured: `treg_only_harmony` is this stage's corrected map, `treg_only` the same cells and PCA
without correction, and `full_object_restricted` the published embedding restricted to the same
Treg cells. The two `_latent` rows measure the 30-dimensional representations the layouts came
from. The three `umap_2d` rows are the head-to-head, because 2D is the space a widget draws.

`observed_same_frac` is the mean over cells of the fraction of that cell's k = 30 nearest
neighbours carrying the same factor value. `expected_same_frac` is what a random relabelling
would give, so unequal group sizes are already absorbed and the two columns compare directly.
`excess_over_chance` is `(observed − expected) / (1 − expected)`.

For `donor`, lower is better, since donor structure is what the correction removes. For
`tissue`, a value that stays high after donor correction says the two factors are separable.
Donor and tissue are nested here, so the `tissue` rows are read for whether they collapse.

### `tables/treg_reembedding_manifest.csv` — the parameters

One row per recorded parameter: `parameter`, `value`, and `source`. `thresholds.*` is the shared
config block the full-object embedding also read, `treg_reembedding.*` is this stage's own block,
`recipe` is a value fixed by the recipe being reproduced, and `derived` means the run measured
it.

`corrected_coord_columns` and `uncorrected_coord_columns` name which parquet columns hold which
variant. `coordinate_state` reads `recomputed` when the embedding ran and `checkpoint_reused`
when the coordinate cache was accepted, which requires the cache to cover exactly the current
subset's barcodes and to carry both representations. `n_cells_embedded` against
`n_cells_expected` is a reproduction check, and a disagreement aborts the run.

Every parameter is the one the full-object embedding used. The correction adds only
`harmony_batch_key = donor` at the harmonypy default. All 16 score columns of the narrative
substrate cross the barcode join.

### `tables/_overview/<figure stem>.csv`

Five same-stem sources, one per figure, each restricted to the channels its panel draws, with the
mean, median and standard deviation over the drawn cells and the counts behind them. Read
`metric` first: `log_normalised_expression`, `AUCell` and
`scanpy_score_genes_module_score` are three unrelated scales, and a comparison stays inside one
of them.

The AUCell values are the ones the full-object figure's Treg rows carry too, because the scores
are joined on barcode and never recomputed, so a difference between paired figures is a
difference of layout.

---

The per-cell feed for the Treg-only companion widget is
[`../interactive/17_treg_reembedding.parquet`](../interactive/): 27,175 rows carrying both
coordinate pairs alongside the identical 16 score columns the full-object map is coloured by.
Regenerable and untracked. On `x_uncorrected`/`y_uncorrected` two thirds of a cell's neighbours
share its donor, so a coloured patch there is a statement about a sample. On `x`/`y` that falls
to 1.38× chance.
