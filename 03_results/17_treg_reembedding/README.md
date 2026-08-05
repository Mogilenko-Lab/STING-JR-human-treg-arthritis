# 17_treg_reembedding: artifact captions

A second map of the sorted JIA synovial-fluid/blood T-cell compartment (GSE160097), computed on the
Treg gate alone. The published embedding lays out all three sorted lineages in one geometry, so the
variance separating Treg from Tcon from CD8 sets its axes and Treg substructure occupies the room
that is left. Here the 27,175 cells of the frozen `coarse_label == "Treg"` gate get the whole canvas.

Per-cell scores are joined in from `03_results/interactive/16_narrative_embedding.parquet` on
barcode and are not recomputed, so every map can be coloured by identical values and a difference
between two maps is a difference of layout.

**Recipe.** The six calls of the full-object embedding, in the same order, on the subset:
`highly_variable_genes(n_top_genes=3000)` → subset to those genes → `scale(max_value=10)` →
`pca(n_comps=30)` → `neighbors(n_pcs=30)` → `umap`, with `random_state = 0`. `hvg_n_top` and `n_pcs`
come from the same `thresholds:` block the full-object run read, and `n_neighbors` stays at the
scanpy default because the full-object call left it there, so the two maps differ by which cells went
in and by the batch correction alone. Expression is the log-normalised `X` of the frozen annotation
checkpoint; per-cell normalisation is independent of which other cells are present, so subsetting
carries it over unchanged.

**Two coordinate pairs, one PCA.** The uncorrected recipe applies no batch correction, and on this
subset it resolves donor-by-tissue sample of origin: 66.1% same-donor neighbours against 14.6%
expected. So the neighbours-and-UMAP tail runs twice off the same PCA, once on the raw PCA and once
on the PCA after Harmony over `donor` at harmonypy standard settings, seeded, converged after 5
iterations. `x`/`y` are the Harmony-corrected coordinates and are the pair to draw.
`x_uncorrected`/`y_uncorrected` are the uncorrected pair on the same cells, carried for comparison.
Harmony is a visualisation aid under the umbrella embeddings guardrail and supports no claim.

**Verdict on the correction, at k = 30.**

| Factor | Uncorrected | Corrected | Chance | Full object, same cells |
|---|---|---|---|---|
| same-donor neighbours | 0.661 | 0.201 | 0.146 | 0.420 |
| same-tissue neighbours | 0.975 | 0.923 | 0.500 | 0.955 |

The corrected map sits at 1.38× chance on donor (excess index 0.064). Every donor drops, and the
residual concentrates in JIA_patient_5 at 0.271. Tissue holds at 1.85× chance, so donor and tissue
are separable in this representation even though the two are nested at the sample level: each GSM is
one donor × one tissue × one sorted population. The corrected coordinates are usable for reading
Treg substructure. A small tight cluster still deserves a donor-composition check, and Harmony
reshapes the space it corrects, so this map remains a map.

**Counterparts, and why they share a colour scale.** Three figures here are Treg-only counterparts of
program rows drawn on the full object: `umap_treg_arms` against `16_narrative_scoring/umap_full_arms`,
`umap_treg_programs` against `16_narrative_scoring/umap_full_programs`, and `umap_treg_signatures`
against `07_embedding/umap_signatures_treg`. Each answers one question — whether a program's apparent
structure needs the Tcon and CD8 gates to appear — and the answer is readable only when both halves of
the pair are on one colour scale. So every counterpart panel is clipped to its full-object twin's
limits, and the arm pair's shared limits are asserted at run time against the value the
`16_narrative_scoring` README records, because a pair that has silently stopped sharing a scale still
looks like a matched set.

The consequence to hold onto: a washed-out counterpart panel is a statement about the Treg gate's range,
and the published IFN-independent STING panel is the clear case.

Each pair shares its cells and its colour scale. The coordinates differ, because these are a
re-embedding of the Treg cells alone, so a cell sits somewhere else here than on the full-object map.
The layout is the thing that differs, which is the point.

**Tier.** Annotation and visualisation only. An embedding places cells; it tests nothing. Nothing
here pools with the donor-level pseudobulk spine and this stage writes no effect-size row.
Confirmatory claims are carried by donor-level pseudobulk differential expression within the frozen
labels. Proximity on any of these maps carries no evidential weight.

## tables/treg_reembedding_mixing.csv

Harmony over `donor` brings the Treg-only map from 0.661 same-donor neighbours to 0.201 against a
0.146 chance level, nearer chance than the full-object map restricted to the same 27,175 Tregs
(0.420), and it holds for all seven donors individually. Same-tissue neighbours fall only 0.975 to
0.923 against 0.500 expected, so the synovial-fluid/blood separation survives the removal of nearly
all donor structure. The uncorrected map's donor dominance is a property of the recipe that omits
batch correction.

**How to read:** One row per (`embedding` × `space` × `grouping_key` × `group`). `embedding` names
which map was measured: `treg_only_harmony` is this stage's corrected map (`x`/`y` in the parquet),
`treg_only` the same cells and same PCA without correction (`x_uncorrected`/`y_uncorrected`), and
`full_object_restricted` the published `X_umap_unsupervised` with rows restricted to the same Treg
cells. The two `_latent` rows measure the 30-dimensional representations the subset layouts came
from, `pca_30` and `harmony_30`. `space` says which was searched; the three `umap_2d` rows are the
head-to-head, because the published embedding has no stored latent space in the annotation checkpoint
and 2D is also the space a widget draws.

`grouping_key` is the factor whose neighbourhood composition is measured and `group` one factor value,
with `_all_` the cell-weighted overall row. `observed_same_frac` is the mean over cells of the fraction
of that cell's k = 30 nearest neighbours carrying the same factor value. `expected_same_frac` is what a
random relabelling would give, `(n_group - 1) / (N - 1)` for a cell of that group, so unequal group
sizes are already absorbed and the two columns compare directly. `excess_over_chance` is
`(observed - expected) / (1 - expected)`: 0 means neighbourhoods look like the dataset composition, 1
means every neighbour shares the factor value.

For `donor`, lower is better, since donor structure is what the correction removes. For `tissue`, a
value that stays high after donor correction says the two factors are separable. Donor and tissue are
nested here — each GSM is one donor × one tissue × one sorted population — so the `tissue` rows are read
for whether they collapse; their absolute size supports no reading. `evidence_tier` reads
`annotation_embedding` throughout; this table is a property of a layout.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/17_treg_reembedding.py` | `mixing_table` | `mixing_k = 30`, `mixing_keys = [donor, tissue]`, `harmony_batch_key = donor`, `reference_embedding_key = X_umap_unsupervised`, `n_pcs = 30` | `03_results/objects/02_annotation.h5ad`, `03_results/objects/17_treg_reembedding.parquet` |

## tables/treg_reembedding_manifest.csv

Every parameter of the Treg-only maps is the parameter the full-object embedding used, read from the
same config keys: `hvg_n_top = 3000` and `n_pcs = 30` from `thresholds:`, `scale(max_value=10)` and
the scanpy default `n_neighbors` from the recipe. The correction adds only `harmony_batch_key = donor`
at the harmonypy default `max_iter_harmony = 10`. The subset holds exactly the 27,175 cells the frozen
`coarse_label == "Treg"` gate declares, and all 16 score columns of the narrative substrate cross the
barcode join.

**How to read:** One row per recorded parameter: `parameter`, the `value` used in the run, and
`source`. `thresholds.*` is the shared config block the full-object embedding also read;
`treg_reembedding.*` is this stage's own block; `recipe` is a value fixed by the embedding recipe being
reproduced; `derived` means the run measured it.

`corrected_coord_columns` and `uncorrected_coord_columns` name which parquet columns hold which variant.
`harmony_max_iter` is the harmonypy default, stated because the correction ran once at standard
settings. `coordinate_state`
reads `recomputed` when the embedding ran and `checkpoint_reused` when the coordinate cache under
`03_results/objects/` was accepted, which requires the cache to cover exactly the current subset's
barcodes and to carry both representations and both coordinate pairs. `n_cells_embedded` against
`n_cells_expected` is a reproduction check; a disagreement aborts the run. No biological claim.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/17_treg_reembedding.py` | `build_manifest` | `subset_key = coarse_label`, `subset_value = Treg`, `expected_n_cells = 27175`, `random_seed = 0`, `harmony_batch_key = donor`, `harmony_max_iter = 10`, `hvg_n_top = 3000`, `n_pcs = 30` | `03_results/objects/02_annotation.h5ad`, `03_results/interactive/16_narrative_embedding.parquet` |

## ../interactive/17_treg_reembedding.parquet

The per-cell substrate for a Treg-only companion widget: 27,175 rows, one per sorted Treg cell,
carrying both coordinate pairs alongside the identical 16 score columns the full-object map is
coloured by, so a set can be viewed on any of the maps without rescoring. Written alongside the other
per-cell interactive substrates in `03_results/interactive/`, and regenerable (untracked) as those
are.

**How to read:** `barcode` is the annotation object's cell name, a 10x barcode plus GSM suffix, and is
the join key. `x`/`y` are the Harmony-corrected Treg-only UMAP coordinates and are the pair to draw;
`x_uncorrected`/`y_uncorrected` are the same cells and same PCA without the correction. These four are
the only columns this stage computes. Each pair sits on its own arbitrary scale, so the two pairs
compare coordinate for coordinate neither with each other nor with the `x`/`y` of
`16_narrative_embedding.parquet`. `coarse_label` (constant `Treg` here), `tissue`, `donor` and
`pct_counts_mt` are carried verbatim.

The 13 `*_AUCell` columns hold per-cell AUCell of the four mouse-derived human-projected up arms and
the nine curated anchor-independent lenses, bounded in [0, 1] and comparable across tissue within a
column; AUCell's scale depends on set size, so cross-column comparison takes the source tables. The
three `published_*` columns come from the earlier readout, and `published_WT_heat_up` holds a scanpy
`score_genes` module score, so colour the mouse up arm with `WT_heat_up_AUCell`. On
`x_uncorrected`/`y_uncorrected` two thirds of a cell's neighbours share its donor, so a coloured patch
there is a statement about a sample; on `x`/`y` that falls to 1.38× chance. Annotation tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/17_treg_reembedding.py` | `build_substrate` | `subset_key = coarse_label`, `subset_value = Treg`, `harmony_batch_key = donor`, `source_parquet = 16_narrative_embedding.parquet`, `output_parquet = 17_treg_reembedding.parquet` | `03_results/objects/02_annotation.h5ad`, `03_results/interactive/16_narrative_embedding.parquet` |

## figures/_overview/umap_treg_reembedding.png

On the Treg-only map, drawn on the Harmony-corrected coordinates, the
synovial-fluid and paired-blood cells still occupy distinct territory
— 0.923 same-tissue neighbours at k = 30 against 0.500 expected —
after the same-donor neighbour fraction has fallen from 0.661 to 0.201
against 0.146 expected. Both the mouse WT 39 °C-derived up arm and the
curated Hallmark hypoxia lens colour the synovial-fluid territory
brighter: per-cell AUCell mean 0.0112 to 0.0190 for WT_heat_up, 0.0736
to 0.0974 for the hypoxia lens.

**How to read:** Three panels over one frame of the same 27,175 sorted Treg cells at
the same coordinates, sharing one square bounding box. Left is tissue
of origin, synovial fluid in vermillion and paired blood in blue,
drawn in shuffled order. Middle and right colour every cell by per-
cell AUCell of one gene set, on the scale the full-object figures use,
clipped to the 2nd and 98th percentile with the highest-scoring cells
drawn last. Panel titles carry the set identifier and its size:
WT_heat_up is the up arm of the mouse WT iTreg 39-versus-37 °C
contrast in human projection, HALLMARK_HYPOXIA the curated MSigDB
Hallmark program. The two sets are unrelated and their ranges differ,
so each keeps its own colour scale. The coordinates are the Harmony-
corrected pair, and Harmony reshapes the space it corrects, so this
map is annotation. Cells are pooled across donors, so a tissue
difference read off the colouring is pseudoreplicated. Claims rest on
donor-level pseudobulk differential expression within the frozen cell
states.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/17_treg_reembedding_viz.py` | `figure_treg` | `coordinates = x / y (Harmony over donor), all 27,175 cells drawn, point_size = 2.4, cmap = viridis, clip_percentiles = [2, 98], figures.dpi = 300, figures.rasterized_dpi = 600, columns = WT_heat_up_AUCell, HALLMARK_HYPOXIA_AUCell` | `03_results/interactive/17_treg_reembedding.parquet, 03_results/16_narrative_scoring/tables/narrative_score_summary.csv, 03_results/17_treg_reembedding/tables/treg_reembedding_mixing.csv` |

## tables/_overview/umap_treg_reembedding.csv

Per-cell AUCell summaries of the 2 sets drawn on the Treg-only map,
WT_heat_up and HALLMARK_HYPOXIA, restricted to the Treg gate, one row
per tissue, so the colouring reads as numbers.

**How to read:** A restriction of the narrative scoring summary table to the Treg gate
and the two sets this figure draws. One row per (`set_name` x
`tissue`) with the mean, median and standard deviation of the per-cell
AUCell score and the cell and donor counts behind it. AUCell is
bounded in [0, 1] and its scale depends on set size, so values compare
across tissue within a `set_name`. Cells are pooled across donors, so
the unit of replication is the cell and the tissue difference here is
pseudoreplicated. `evidence_tier` reads `secondary_percell`
throughout.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/17_treg_reembedding_viz.py` | `score_table` | `rows = 2 sets x Treg x 2 tissues, metric = AUCell` | `03_results/16_narrative_scoring/tables/narrative_score_summary.csv` |

## figures/_overview/umap_treg_arms.png

Read on the Treg gate's own map and on the full object's colour scale,
all three mouse 39 °C-derived up arms still colour the synovial-fluid
territory brighter than the paired-blood territory, so the tissue
contrast the full-object row shows survives viewing Treg apart from
Tcon and CD8. Per-cell AUCell means run 0.0112 to 0.0190 for
WT_heat_up, 0.0111 to 0.0192 for KO_heat_up, and 0.0717 to 0.1245 for
the 7-gene Interaction_up.

**How to read:** The Treg-only counterpart of
`16_narrative_scoring/figures/_overview/umap_full_arms.png`: same
three sets, same order, same panel geometry, same sequential colormap,
drawn on the 27,175 sorted Treg cells alone. Panel titles carry the
set identifier and its size. WT_heat_up is the up arm of the mouse WT
iTreg 39-versus-37 °C contrast in human projection, 199 symbols;
KO_heat_up the same contrast in cGAS-knockout iTregs, 218 symbols;
Interaction_up the mouse genotype-by-temperature up arm, 7 symbols,
small enough that one gene moves the score, so read it for location
and treat its spread as noise. WT_heat_up and KO_heat_up share 182
genes and one colour scale, so the two panels compare directly;
Interaction_up spans a range an order of magnitude wider and keeps its
own bar. AUCell is bounded in [0, 1] and scales with set size, so the
source table carries mean, median, standard deviation and cell and
donor counts for any comparison the colour cannot make. Colour limits,
coordinates and tier follow the counterpart contract at the top of
this page: limits are the full-object figure's from the same frame and
seed, so a washed-out panel is a real statement about the Treg gate's
range; the coordinates are this map's own; cells are pooled across
donors, so a tissue difference read off the colouring is
pseudoreplicated; and claims rest on donor-level pseudobulk
differential expression within the frozen cell states.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/17_treg_reembedding_viz.py` | `figure_counterpart` | `coordinates = x / y (Harmony over donor), all 27,175 cells drawn, point_size = 2.4, cmap = viridis, figures.dpi = 300, figures.rasterized_dpi = 600, colour limits from 03_results/interactive/16_narrative_embedding.parquet at sample_n = 60000, sample_seed = 0, clip_percentiles = [2, 98], columns = WT_heat_up_AUCell, KO_heat_up_AUCell, Interaction_up_AUCell, shared_scale = WT_heat_up + KO_heat_up at 0.0023-0.0336` | `03_results/interactive/17_treg_reembedding.parquet, 03_results/interactive/16_narrative_embedding.parquet, 03_results/16_narrative_scoring/tables/narrative_score_summary.csv` |

## figures/_overview/umap_treg_programs.png

All three curated lenses colour synovial-fluid territory brighter than
paired blood on the Treg gate's own map, so that structure appears
without the other two sort gates: per-cell AUCell means 0.0736 to
0.0974 for hypoxia, 0.0176 to 0.0298 for the generic interferon axis,
0.0186 to 0.0366 for the 21 published IFN-independent STING genes. The
shared scale adds the level: the published STING panel sits low in its
bar across the whole Treg map, Treg blood mean 0.0186 against Tcon
0.0532 and CD8 0.0647. Its median is exactly 0.000, the only 1 of 6
rows at zero, against Tcon 0.0512 and CD8 0.0610, so at least half of
Treg blood cells score zero on that 21-gene set and the Treg synovial-
versus-blood difference on it rests partly on a zero-inflated blood
baseline.

**How to read:** The Treg-only counterpart of
`16_narrative_scoring/figures/_overview/umap_full_programs.png`: same
three sets, same order, same panel geometry, same colormap, on the
27,175 sorted Treg cells alone. Every set here is curated, versioned
and derived independently of the mouse anchor — HALLMARK_HYPOXIA from
MSigDB Hallmark, sting_specific_published the 21 published IFN-
independent STING-activation genes, ifn_generic_axis a 200-gene
generic type-I interferon axis — so a colouring here stands apart from
the anchor. The three sets are unrelated and their ranges differ, so
each panel keeps its own set's limits; what is shared is the object
those limits come from, which makes each panel comparable to its own
twin. Brightness therefore compares tissue within a panel, and the
source table carries the cross-panel numbers. Two limits bind the
reading. The published STING set is 21 genes and its own IFN-β
validation in the positive-control compartment is underpowered at
three donors, so a dim or bright panel there is consistent with STING
pathway activity and is never proof of it. Hypoxia and temperature are
both imposed by the inflamed joint and stay entangled in cross-
sectional human data, so the hypoxia panel is one lens on that niche
and carries no HIF claim. Every colour limit is the full-object
figure's, from the same frame and seed, so the pair compares
brightness for brightness. A washed-out counterpart panel is therefore
a real statement about the Treg gate's range. The pair shares its
cells and its colour scale; the coordinates differ, because these are
a re-embedding of the Treg cells alone. Cells are pooled across
donors, so a tissue difference read off the colouring is
pseudoreplicated, and Harmony reshapes the space it corrects, so this
is a map. Claims in this compartment rest on donor-level pseudobulk
differential expression within the frozen cell states.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/17_treg_reembedding_viz.py` | `figure_counterpart` | `coordinates = x / y (Harmony over donor), all 27,175 cells drawn, point_size = 2.4, cmap = viridis, figures.dpi = 300, figures.rasterized_dpi = 600, colour limits from 03_results/interactive/16_narrative_embedding.parquet at sample_n = 60000, sample_seed = 0, clip_percentiles = [2, 98], columns = HALLMARK_HYPOXIA_AUCell, sting_specific_published_AUCell, ifn_generic_axis_AUCell; per-panel full-object limits, no pooling across sets` | `03_results/interactive/17_treg_reembedding.parquet, 03_results/interactive/16_narrative_embedding.parquet, 03_results/16_narrative_scoring/tables/narrative_score_summary.csv` |

## figures/_overview/umap_treg_signatures.png

The three candidate harvest signatures on the Treg gate's own map and
on the full-object figure's limits. Within the Treg gate the effector-
Treg and heat-shock module scores still separate synovial fluid from
paired blood (score_eTreg -0.0539 to 0.1021, score_HSP 0.0270 to
0.0806), alongside the mouse anchor annotation (WT_heat_up -0.0699 to
-0.0475), so the structure the full-object figure shows for these
three channels appears without viewing the Treg gate against Tcon and
CD8.

**How to read:** The Treg-only counterpart of
`07_embedding/figures/_overview/umap_signatures_treg.png`, which draws
these same three channels across all three sort gates. Two differences
from the other counterparts here matter first. The unit: these are
scanpy `score_genes` module scores, mean-centred against a sampled
background and signed, so they share a scale neither with the AUCell
panels here nor with each other, and a value near zero means at
background. The AUCell reading of the same mouse arm is the WT_heat_up
panel of `umap_treg_arms`. The panel count: the full-object twin
carries a fourth panel, a Treg/Tcon/CD8 sort-gate reference, which is
a single category on a Treg-only object and is therefore omitted; the
tissue reference here is the left panel of `umap_treg_reembedding`.
WT_heat_up here is the mouse WT 39-versus-37 °C up arm carried as
annotation only — never a selection predicate, and the harvest design
it was previewed for is frozen as implemented. score_eTreg is the
effector-Treg score and score_HSP the heat-shock/stress score, both
curated in this compartment independently of the anchor. Colour limits
come from all 99,915 cells of the 07_embedding substrate, which is the
frame that figure draws, so the pair is comparable panel for panel.
Colour limits, coordinates and tier follow the counterpart contract at
the top of this page: limits are the full-object figure's from the
same frame and seed, so a washed-out panel is a real statement about
the Treg gate's range; the coordinates are this map's own; cells are
pooled across donors, so a tissue difference read off the colouring is
pseudoreplicated; and claims rest on donor-level pseudobulk
differential expression within the frozen cell states.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/17_treg_reembedding_viz.py` | `figure_counterpart` | `coordinates = x / y (Harmony over donor), all 27,175 cells drawn, point_size = 2.4, cmap = viridis, figures.dpi = 300, figures.rasterized_dpi = 600, colour limits from 03_results/interactive/16_narrative_embedding.parquet at sample_n = 60000, sample_seed = 0, clip_percentiles = [2, 98], columns = WT_heat_up, score_eTreg, score_HSP; joined on barcode from 03_results/07_embedding/tables/hook_factor_substrate.parquet; limits from all 99,915 cells of that substrate, which is the frame 07_embedding_viz.py draws; metric = scanpy score_genes module score` | `03_results/interactive/17_treg_reembedding.parquet, 03_results/07_embedding/tables/hook_factor_substrate.parquet` |

## tables/_overview/umap_treg_arms.csv

Per-cell AUCell summaries of the 3 sets drawn in
`figures/_overview/umap_treg_arms.png` — WT_heat_up, KO_heat_up,
Interaction_up — restricted to the Treg gate, one row per tissue, so
the colouring reads as numbers.

**How to read:** A restriction of the narrative scoring summary table to the Treg gate
and the sets this figure draws. One row per (`set_name` x `tissue`)
with the mean, median and standard deviation of the per-cell AUCell
score and the cell and donor counts behind it. These are the values
the full-object figure's Treg rows carry too, because the scores are
joined on barcode and not recomputed, so a difference between the
paired figures is a difference of layout. AUCell is bounded in [0, 1]
and its scale depends on set size, so values compare across tissue
within a `set_name`. Cells are pooled across donors, so the unit of
replication is the cell and every tissue difference is
pseudoreplicated. `evidence_tier` reads `secondary_percell`.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/17_treg_reembedding_viz.py` | `score_table` | `rows = 3 sets x Treg x 2 tissues, metric = AUCell` | `03_results/16_narrative_scoring/tables/narrative_score_summary.csv` |

## tables/_overview/umap_treg_programs.csv

Per-cell AUCell summaries of the 3 sets drawn in
`figures/_overview/umap_treg_programs.png` — HALLMARK_HYPOXIA,
sting_specific_published, ifn_generic_axis — restricted to the Treg
gate, one row per tissue, so the colouring reads as numbers.

**How to read:** A restriction of the narrative scoring summary table to the Treg gate
and the sets this figure draws. One row per (`set_name` x `tissue`)
with the mean, median and standard deviation of the per-cell AUCell
score and the cell and donor counts behind it. These are the values
the full-object figure's Treg rows carry too, because the scores are
joined on barcode and not recomputed, so a difference between the
paired figures is a difference of layout. AUCell is bounded in [0, 1]
and its scale depends on set size, so values compare across tissue
within a `set_name`. Cells are pooled across donors, so the unit of
replication is the cell and every tissue difference is
pseudoreplicated. `evidence_tier` reads `secondary_percell`.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/17_treg_reembedding_viz.py` | `score_table` | `rows = 3 sets x Treg x 2 tissues, metric = AUCell` | `03_results/16_narrative_scoring/tables/narrative_score_summary.csv` |

## tables/_overview/umap_treg_signatures.csv

Per (channel x tissue) summaries of the three candidate harvest
signatures within the Treg gate, giving the numbers behind the
counterpart figure's colouring: all three channels sit higher in
synovial fluid than in paired blood.

**How to read:** One row per (`set_name` x `tissue`) over the 27,175 Treg cells, with
the mean, median and standard deviation of the module score and the
cell and donor counts behind it. `metric` reads
`scanpy_score_genes_module_score`. These values are mean-centred
against a sampled background and signed, so zero means at background
and a negative mean records a position on that scale; they compare
across tissue within a `set_name`. `WT_heat_up` is the mouse anchor
arm carried as annotation only and never as a selection predicate.
Cells are pooled across donors, so the unit of replication is the cell
and every tissue difference is pseudoreplicated. Annotation tier; no
test, no effect size.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/17_treg_reembedding_viz.py` | `signature_table` | `rows = 3 channels x Treg x 2 tissues, metric = scanpy score_genes module score` | `03_results/07_embedding/tables/hook_factor_substrate.parquet` |
