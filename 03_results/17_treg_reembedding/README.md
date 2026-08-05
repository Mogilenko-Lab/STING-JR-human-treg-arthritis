# 17_treg_reembedding: artifact captions

A second map of the sorted JIA synovial-fluid/blood T-cell compartment (GSE160097), computed on
the Treg gate alone. The published embedding lays out all three sorted lineages in one geometry,
so the variance separating Treg from Tcon from CD8 sets the axes and Treg substructure occupies
what room is left; here the 27,175 cells of the `coarse_label == "Treg"` gate get the whole
canvas. Per-cell scores are joined in from `03_results/interactive/16_narrative_embedding.parquet`
on barcode and are not recomputed, so every map can be coloured by identical values and a
difference between them is a difference of layout.

**Recipe.** The six calls of the full-object embedding, in the same order, on the subset:
`highly_variable_genes(n_top_genes=3000)` -> subset to those genes -> `scale(max_value=10)` ->
`pca(n_comps=30)` -> `neighbors(n_pcs=30)` -> `umap`, with `random_state = 0`. `hvg_n_top` and
`n_pcs` are read from the same `thresholds:` block the full-object run read, and `n_neighbors` is
left at the scanpy default because the full-object call left it there, so the maps differ by
which cells went in and by the batch correction alone. Expression is the log-normalised `X` of
the frozen annotation checkpoint; per-cell normalisation does not depend on which other cells are
present, so subsetting does not require re-deriving it.

**Two coordinate pairs, one PCA.** The uncorrected recipe applies no batch correction, and on this
subset it resolves donor-by-tissue sample of origin (66.1% same-donor neighbours against 14.6%
expected). So the neighbours-and-UMAP tail runs twice off the same PCA: once on the raw PCA, once
on the PCA after Harmony over `donor` at harmonypy standard settings, seeded, converged in 5
iterations. `x`/`y` are the Harmony-corrected coordinates and are the pair to draw.
`x_uncorrected`/`y_uncorrected` are the uncorrected pair, carried for comparison on the same
cells. Harmony is a visualisation aid under the umbrella embeddings guardrail and supports no
claim.

**Verdict on the correction, at k = 30.** Same-donor neighbours fall from 0.661 uncorrected to
0.201 corrected, against a 0.146 chance level, so the corrected map sits at 1.38x chance
(excess index 0.064) and below the full-object map's 0.420. Every donor drops, and the residual
concentrates in one donor (JIA_patient_5, 0.271). Same-tissue neighbours go 0.975 to 0.923
against 0.500 expected, so removing almost all donor structure leaves the synovial-fluid/blood
split at 1.85x chance: donor and tissue are separable in this representation, even though the two
are nested at the sample level. The corrected coordinates are usable for reading Treg
substructure. A small tight cluster still deserves a donor-composition check, and Harmony
reshapes the space it corrects, so this map remains a map.

**Counterparts, and why they share a colour scale.** Three of the figures here are Treg-only
counterparts of program rows drawn on the full object: `umap_treg_arms` against
`16_narrative_scoring/umap_full_arms`, `umap_treg_programs` against
`16_narrative_scoring/umap_full_programs`, and `umap_treg_signatures` against
`07_embedding/umap_signatures_treg`. Each exists to answer one question — whether a program's
apparent structure needs the Tcon and CD8 gates to appear — and the answer is only readable if
the two halves of the pair are on one colour scale. So every counterpart panel is clipped to its
full-object twin's limits rather than its own, and the arm pair's shared limits are asserted at
run time against the value the `16_narrative_scoring` README records, because a pair that has
silently stopped sharing a scale still looks like a matched set. The consequence to hold onto is
that a washed-out counterpart panel is a statement about the Treg gate's range and not a
rendering fault; the published IFN-independent STING panel is the clear case.

What the pair does **not** share is coordinates. These are a re-embedding of the Treg cells
alone, so a cell does not sit where it sits on the full-object map. The pair shares its cells and
its colour scale; the layout is the thing that differs, which is the point.

**Tier.** Annotation / visualisation only. An embedding places cells; it tests nothing. Nothing
here is pooled with the donor-level pseudobulk spine and this stage writes no effect-size row.
Confirmatory claims are carried by donor-level pseudobulk differential expression within the
frozen labels. Proximity on any of these maps is not a result.

## tables/treg_reembedding_mixing.csv

Harmony over `donor` brings the Treg-only map from 0.661 same-donor neighbours to 0.201 against a
0.146 chance level, which is nearer chance than the full-object map restricted to the same 27,175
Tregs (0.420) and holds for all seven donors individually, while same-tissue neighbours fall only
0.975 to 0.923 against 0.500 expected, so the synovial-fluid/blood separation survives the
removal of nearly all donor structure, and the uncorrected map's donor dominance is a property of
the recipe that omits batch correction.

**How to read:** One row per (`embedding` x `space` x `grouping_key` x `group`).
`embedding` names which map was measured: `treg_only_harmony` is this stage's corrected map
(`x`/`y` in the parquet), `treg_only` is the same cells and same PCA without correction
(`x_uncorrected`/`y_uncorrected`), `full_object_restricted` is the published
`X_umap_unsupervised` with the rows restricted to the same Treg cells. The two `_latent` rows
measure the 30-dimensional representations the two subset layouts came from, `pca_30` and
`harmony_30`. `space` says which of those was searched; the three `umap_2d` rows are the
head-to-head, because the published embedding has no stored latent space in the annotation
checkpoint and 2D is also the space a widget draws.
`grouping_key` is the factor whose neighbourhood composition is measured; `group` is one factor
value, with `_all_` the cell-weighted overall row.
`observed_same_frac` is the mean over cells of the fraction of that cell's `k` = 30 nearest
neighbours carrying the same factor value. `expected_same_frac` is what a random relabelling
would give, `(n_group - 1) / (N - 1)` for a cell of that group, so unequal group sizes are
already absorbed and the two columns are directly comparable. `excess_over_chance` is
`(observed - expected) / (1 - expected)`: 0 means neighbourhoods look like the dataset
composition, 1 means every neighbour shares the factor value. For `donor`, lower is better,
since donor structure is what the correction is meant to remove; for `tissue`, a value that
stays high after donor correction says the two factors are separable.
Donor and tissue are nested in this design: each GSM is one donor x one tissue x one sorted
population, so a same-sample neighbour is same-donor and same-tissue at once, which is why the
`tissue` rows are read for whether they collapse and not for their absolute size.
`evidence_tier` reads `annotation_embedding` throughout; this table is a property of a layout and
carries no biological claim.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/17_treg_reembedding.py` | `mixing_table` | `mixing_k = 30`, `mixing_keys = [donor, tissue]`, `harmony_batch_key = donor`, `reference_embedding_key = X_umap_unsupervised`, `n_pcs = 30` | `03_results/objects/02_annotation.h5ad`, `03_results/objects/17_treg_reembedding.parquet` |

## tables/treg_reembedding_manifest.csv

Every parameter of the Treg-only maps is the parameter the full-object embedding used, sourced
from the same config keys the full-object run read (`hvg_n_top = 3000` and `n_pcs = 30` from
`thresholds:`, `scale(max_value=10)` and the scanpy default `n_neighbors` from the recipe), the
correction adds only `harmony_batch_key = donor` at the harmonypy default `max_iter_harmony = 10`,
the subset holds exactly the 27,175 cells the frozen `coarse_label == "Treg"` gate declares, and
all 16 score columns of the narrative substrate are carried across the barcode join.

**How to read:** One row per recorded parameter: `parameter`, the `value` used in the run, and
`source`, which says where that value came from. `thresholds.*` means the shared config block
the full-object embedding also read, so those two rows come from one place for every map;
`treg_reembedding.*` means this stage's own config block; `recipe` means a value fixed by the
embedding recipe being reproduced; `derived` means the run measured it.
`corrected_coord_columns` and `uncorrected_coord_columns` say which parquet columns hold which
variant, so a consumer does not have to infer it. `harmony_max_iter` is the harmonypy default,
stated here because the correction was run once at standard settings and reported as it came.
`coordinate_state` reads `recomputed` when the embedding ran and `checkpoint_reused` when the
coordinate cache under `03_results/objects/` was accepted, which requires the cache to cover
exactly the current subset's barcodes and to carry both representations and both coordinate
pairs. `n_cells_embedded` against `n_cells_expected` is a reproduction check, and a disagreement
aborts the run. No biological claim.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/17_treg_reembedding.py` | `build_manifest` | `subset_key = coarse_label`, `subset_value = Treg`, `expected_n_cells = 27175`, `random_seed = 0`, `harmony_batch_key = donor`, `harmony_max_iter = 10`, `hvg_n_top = 3000`, `n_pcs = 30` | `03_results/objects/02_annotation.h5ad`, `03_results/interactive/16_narrative_embedding.parquet` |

## ../interactive/17_treg_reembedding.parquet

The per-cell substrate for a Treg-only companion widget: 27,175 rows, one per sorted Treg cell,
carrying both coordinate pairs alongside the identical 16 score columns the full-object map is
coloured by, so a set can be viewed on any of the maps without rescoring. Written alongside the
other per-cell interactive substrates in `03_results/interactive/`, and regenerable (untracked)
as those are.

**How to read:** `barcode` is the annotation object's cell name (10x barcode plus GSM suffix) and
is the join key. `x`/`y` are the Harmony-corrected Treg-only UMAP coordinates and are the pair to
draw; `x_uncorrected`/`y_uncorrected` are the same cells and same PCA without the correction, for
side-by-side comparison. These four are the only columns this stage computes. Each pair is on its
own arbitrary scale, so the two pairs are not comparable coordinate for coordinate with each other
or with the `x`/`y` of `16_narrative_embedding.parquet`. `coarse_label` (constant `Treg` here),
`tissue`, `donor` and `pct_counts_mt` are carried verbatim. The 13 `*_AUCell` columns are
per-cell AUCell of the four mouse-derived human-projected up arms and the nine curated
anchor-independent lenses, bounded in [0, 1] and comparable across tissue **within** a column but
not across columns, since AUCell's scale depends on set size. The three `published_*` columns are
carried from the earlier readout; `published_WT_heat_up` holds a scanpy `score_genes` module
score, so colour the mouse up arm with `WT_heat_up_AUCell`. On `x_uncorrected`/`y_uncorrected`
two thirds of a cell's neighbours share its donor, so a coloured patch there is a statement about
a sample; on `x`/`y` that falls to 1.38x chance. Annotation tier throughout.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/17_treg_reembedding.py` | `build_substrate` | `subset_key = coarse_label`, `subset_value = Treg`, `harmony_batch_key = donor`, `source_parquet = 16_narrative_embedding.parquet`, `output_parquet = 17_treg_reembedding.parquet` | `03_results/objects/02_annotation.h5ad`, `03_results/interactive/16_narrative_embedding.parquet` |

## figures/_overview/umap_treg_reembedding.png

On the Treg-only map, drawn on the Harmony-corrected coordinates, the
synovial-fluid and paired-blood cells still occupy distinct territory
(0.923 same-tissue neighbours at k = 30 against 0.500 expected) after
the same-donor neighbour fraction has dropped from 0.661 to 0.201
against 0.146 expected, and both the mouse WT 39 °C-derived up arm and
the curated Hallmark hypoxia lens colour the synovial-fluid territory
brighter (per-cell AUCell mean 0.0112 to 0.0190 for WT_heat_up and
0.0736 to 0.0974 for the hypoxia lens).

**How to read:** Three panels over ONE frame of the same 27,175 sorted Treg cells at
the same coordinates, sharing one square bounding box. Left is tissue
of origin, synovial fluid in vermillion and paired blood in blue,
drawn in shuffled order. Middle and right colour every cell by per-
cell AUCell of one gene set, on the scale the full-object figures use,
clipped to the 2nd and 98th percentile with the highest-scoring cells
drawn last. Panel titles carry the set identifier and its size:
WT_heat_up is the up arm of the mouse WT iTreg 39 versus 37 °C
contrast in human projection, HALLMARK_HYPOXIA the curated MSigDB
Hallmark program. The two sets are unrelated and their ranges differ,
so each keeps its own colour scale. The coordinates are the Harmony-
corrected pair, because at k = 30 the same-donor neighbour fraction is
0.661 on the uncorrected Treg-only map and 0.201 after Harmony over
donor against 0.146 expected, while same-tissue neighbours hold at
0.923 against 0.500 expected. Harmony reshapes the space it corrects,
so this map is annotation, and cells are pooled across donors, making
a tissue difference read off the colouring pseudoreplicated and
descriptive. Claims in this compartment rest on donor-level pseudobulk
differential expression within the frozen cell states.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/17_treg_reembedding_viz.py` | `figure_treg` | `coordinates = x / y (Harmony over donor), all 27,175 cells drawn, point_size = 2.4, cmap = viridis, clip_percentiles = [2, 98], figures.dpi = 300, figures.rasterized_dpi = 600, columns = WT_heat_up_AUCell, HALLMARK_HYPOXIA_AUCell` | `03_results/interactive/17_treg_reembedding.parquet, 03_results/16_narrative_scoring/tables/narrative_score_summary.csv, 03_results/17_treg_reembedding/tables/treg_reembedding_mixing.csv` |

## tables/_overview/umap_treg_reembedding.csv

Per-cell AUCell summaries of the two sets drawn on the Treg-only map,
WT_heat_up and HALLMARK_HYPOXIA, restricted to the Treg gate, one row
per tissue, so the colouring can be read as numbers.

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
contrast the full-object row shows is not an artifact of viewing Treg
alongside Tcon and CD8 (per-cell AUCell means 0.0112 to 0.0190 for
WT_heat_up, 0.0111 to 0.0192 for KO_heat_up and 0.0717 to 0.1245 for
the 7-gene Interaction_up).

**How to read:** The Treg-only counterpart of
`16_narrative_scoring/figures/_overview/umap_full_arms.png`: same
three sets, same order, same panel geometry, same sequential colormap,
drawn on the 27,175 sorted Treg cells alone. Panel titles carry the
set identifier and its size. WT_heat_up is the up arm of the mouse WT
iTreg 39 versus 37 °C contrast in human projection, 199 symbols;
KO_heat_up the same contrast in cGAS-knockout iTregs, 218 symbols;
Interaction_up the mouse genotype-by-temperature up arm, 7 symbols,
small enough that one gene moves the score, so read it for location
and treat its spread as noise. WT_heat_up and KO_heat_up share 182
genes and share one colour scale so the two panels can be compared
directly; Interaction_up spans a range an order of magnitude wider and
keeps its own bar. AUCell is bounded in [0, 1] and scales with set
size, so the source table carries mean, median, standard deviation and
cell and donor counts for any comparison the colour cannot make. Every
colour limit is the full-object figure's, from the same frame and
seed, so the pair compares brightness for brightness rather than each
panel rescaling to itself. A washed-out counterpart panel is therefore
a real statement about the Treg gate's range, not a rendering fault.
The pair shares its cells and its colour scale and NOT its
coordinates: these are a re-embedding of the Treg cells alone, so a
point does not sit where it sits on the full-object map. Cells are
pooled across donors, making a tissue difference read off the
colouring pseudoreplicated, and Harmony reshapes the space it
corrects, so this is a map. Claims in this compartment rest on donor-
level pseudobulk differential expression within the frozen cell
states.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/17_treg_reembedding_viz.py` | `figure_counterpart` | `coordinates = x / y (Harmony over donor), all 27,175 cells drawn, point_size = 2.4, cmap = viridis, figures.dpi = 300, figures.rasterized_dpi = 600, colour limits from 03_results/interactive/16_narrative_embedding.parquet at sample_n = 60000, sample_seed = 0, clip_percentiles = [2, 98], columns = WT_heat_up_AUCell, KO_heat_up_AUCell, Interaction_up_AUCell, shared_scale = WT_heat_up + KO_heat_up at 0.0023-0.0336` | `03_results/interactive/17_treg_reembedding.parquet, 03_results/interactive/16_narrative_embedding.parquet, 03_results/16_narrative_scoring/tables/narrative_score_summary.csv` |

## figures/_overview/umap_treg_programs.png

All three curated lenses colour synovial-fluid territory brighter than
paired blood on the Treg gate's own map, so none of that structure
needs the other two sort gates to appear (per-cell AUCell means 0.0736
to 0.0974 for hypoxia, 0.0176 to 0.0298 for the generic interferon
axis and 0.0186 to 0.0366 for the 21 published IFN-independent STING
genes); what the shared scale adds is the LEVEL, and the published
STING panel sits in the lower part of its bar across the whole Treg
map because Treg blood carries roughly a third of what the other gates
do (0.0186 against 0.0532 in Tcon and 0.0647 in CD8), which a self-
scaled panel would have hidden.

**How to read:** The Treg-only counterpart of
`16_narrative_scoring/figures/_overview/umap_full_programs.png`: same
three sets, same order, same panel geometry, same colormap, on the
27,175 sorted Treg cells alone. Every set here is curated, versioned
and derived without reference to the mouse anchor — HALLMARK_HYPOXIA
from MSigDB Hallmark, sting_specific_published the 21 published IFN-
independent STING-activation genes, ifn_generic_axis a 200-gene
generic type-I interferon axis — so a colouring here is not a
restatement of the anchor. The three sets are unrelated to each other
and their ranges differ, so each panel keeps its OWN set's limits;
what is shared is the OBJECT the limits come from, which makes each
panel comparable to its own twin and not to its neighbours. Brightness
therefore compares tissue within a panel, and the source table carries
the cross-panel numbers. The published STING set is 21 genes and its
own IFN-β validation in the positive-control compartment is
underpowered at three donors, so a dim or bright panel there is
consistent with STING pathway activity and never proof of it. Hypoxia
and temperature are both imposed by the inflamed joint and stay
entangled in cross-sectional human data, so the hypoxia panel is one
lens on that niche and not a HIF claim. Every colour limit is the
full-object figure's, from the same frame and seed, so the pair
compares brightness for brightness rather than each panel rescaling to
itself. A washed-out counterpart panel is therefore a real statement
about the Treg gate's range, not a rendering fault. The pair shares
its cells and its colour scale and NOT its coordinates: these are a
re-embedding of the Treg cells alone, so a point does not sit where it
sits on the full-object map. Cells are pooled across donors, making a
tissue difference read off the colouring pseudoreplicated, and Harmony
reshapes the space it corrects, so this is a map. Claims in this
compartment rest on donor-level pseudobulk differential expression
within the frozen cell states.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/17_treg_reembedding_viz.py` | `figure_counterpart` | `coordinates = x / y (Harmony over donor), all 27,175 cells drawn, point_size = 2.4, cmap = viridis, figures.dpi = 300, figures.rasterized_dpi = 600, colour limits from 03_results/interactive/16_narrative_embedding.parquet at sample_n = 60000, sample_seed = 0, clip_percentiles = [2, 98], columns = HALLMARK_HYPOXIA_AUCell, sting_specific_published_AUCell, ifn_generic_axis_AUCell; per-panel full-object limits, no pooling across sets` | `03_results/interactive/17_treg_reembedding.parquet, 03_results/interactive/16_narrative_embedding.parquet, 03_results/16_narrative_scoring/tables/narrative_score_summary.csv` |

## figures/_overview/umap_treg_signatures.png

The three candidate harvest signatures drawn on the Treg gate's own
map and on the full-object figure's limits: within the Treg gate the
effector-Treg and heat-shock scores still separate synovial fluid from
paired blood (score_eTreg -0.0539 to 0.1021, score_HSP 0.0270 to
0.0806) alongside the mouse anchor annotation (WT_heat_up -0.0699 to
-0.0475), so the structure the full-object figure shows for these
three channels is not produced by viewing the Treg gate against Tcon
and CD8.

**How to read:** The Treg-only counterpart of
`07_embedding/figures/_overview/umap_signatures_treg.png`, which draws
these same three channels across all three sort gates. Two differences
from the other counterparts on this page are worth having before
reading it. First, the unit: these are scanpy `score_genes` MODULE
SCORES, mean-centred against a sampled background and signed, not
AUCell, so they do not share a scale with the AUCell panels here or
with each other and a value near zero means 'at background', not
'absent'. The AUCell reading of the same mouse arm is the WT_heat_up
panel of `umap_treg_arms`. Second, the full-object twin carries a
FOURTH panel, a Treg/Tcon/CD8 sort-gate reference, which is a single
category on a Treg-only object and is therefore omitted rather than
drawn as one flat colour; the tissue reference for this map is the
left panel of `umap_treg_reembedding`. WT_heat_up here is the mouse WT
39-versus-37 °C up arm carried as ANNOTATION ONLY — it is never a
selection predicate, and the harvest design it was previewed for is
frozen as implemented. score_eTreg is the effector-Treg score and
score_HSP the heat-shock/stress score, both curated in this
compartment independently of the anchor. Colour limits come from all
99,915 cells of the 07_embedding substrate, which is the frame that
figure draws, so the pair is comparable panel for panel. Every colour
limit is the full-object figure's, from the same frame and seed, so
the pair compares brightness for brightness rather than each panel
rescaling to itself. A washed-out counterpart panel is therefore a
real statement about the Treg gate's range, not a rendering fault. The
pair shares its cells and its colour scale and NOT its coordinates:
these are a re-embedding of the Treg cells alone, so a point does not
sit where it sits on the full-object map. Cells are pooled across
donors, making a tissue difference read off the colouring
pseudoreplicated, and Harmony reshapes the space it corrects, so this
is a map. Claims in this compartment rest on donor-level pseudobulk
differential expression within the frozen cell states.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/17_treg_reembedding_viz.py` | `figure_counterpart` | `coordinates = x / y (Harmony over donor), all 27,175 cells drawn, point_size = 2.4, cmap = viridis, figures.dpi = 300, figures.rasterized_dpi = 600, colour limits from 03_results/interactive/16_narrative_embedding.parquet at sample_n = 60000, sample_seed = 0, clip_percentiles = [2, 98], columns = WT_heat_up, score_eTreg, score_HSP; joined on barcode from 03_results/07_embedding/tables/hook_factor_substrate.parquet; limits from all 99,915 cells of that substrate, which is the frame 07_embedding_viz.py draws; metric = scanpy score_genes module score` | `03_results/interactive/17_treg_reembedding.parquet, 03_results/07_embedding/tables/hook_factor_substrate.parquet` |

## tables/_overview/umap_treg_arms.csv

Per-cell AUCell summaries of the 3 sets drawn in
`figures/_overview/umap_treg_arms.png` (WT_heat_up, KO_heat_up,
Interaction_up), restricted to the Treg gate, one row per tissue, so
the colouring can be read as numbers rather than inferred from
brightness.

**How to read:** A restriction of the narrative scoring summary table to the Treg gate
and the sets this figure draws. One row per (`set_name` x `tissue`)
with the mean, median and standard deviation of the per-cell AUCell
score and the cell and donor counts behind it. These are the values
the FULL-OBJECT figure's Treg rows carry too, because the scores are
joined on barcode and not recomputed for the re-embedding, so a
difference between the paired figures is a difference of layout and
nothing else. AUCell is bounded in [0, 1] and its scale depends on set
size, so values compare across tissue within a `set_name`. Cells are
pooled across donors, so the unit of replication is the cell and every
tissue difference here is pseudoreplicated. `evidence_tier` reads
`secondary_percell`.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/17_treg_reembedding_viz.py` | `score_table` | `rows = 3 sets x Treg x 2 tissues, metric = AUCell` | `03_results/16_narrative_scoring/tables/narrative_score_summary.csv` |

## tables/_overview/umap_treg_programs.csv

Per-cell AUCell summaries of the 3 sets drawn in
`figures/_overview/umap_treg_programs.png` (HALLMARK_HYPOXIA,
sting_specific_published, ifn_generic_axis), restricted to the Treg
gate, one row per tissue, so the colouring can be read as numbers
rather than inferred from brightness.

**How to read:** A restriction of the narrative scoring summary table to the Treg gate
and the sets this figure draws. One row per (`set_name` x `tissue`)
with the mean, median and standard deviation of the per-cell AUCell
score and the cell and donor counts behind it. These are the values
the FULL-OBJECT figure's Treg rows carry too, because the scores are
joined on barcode and not recomputed for the re-embedding, so a
difference between the paired figures is a difference of layout and
nothing else. AUCell is bounded in [0, 1] and its scale depends on set
size, so values compare across tissue within a `set_name`. Cells are
pooled across donors, so the unit of replication is the cell and every
tissue difference here is pseudoreplicated. `evidence_tier` reads
`secondary_percell`.

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
`scanpy_score_genes_module_score`, not AUCell: these values are mean-
centred against a sampled background and signed, so zero means 'at
background' and a negative mean is not an absence. They therefore do
not compare with the AUCell tables beside them, nor across channels,
only across tissue within a `set_name`. `WT_heat_up` is the mouse
anchor arm carried as annotation only and never as a selection
predicate. Cells are pooled across donors, so the unit of replication
is the cell and every tissue difference here is pseudoreplicated and
descriptive. Annotation tier; no test, no effect size.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/17_treg_reembedding_viz.py` | `signature_table` | `rows = 3 channels x Treg x 2 tissues, metric = scanpy score_genes module score` | `03_results/07_embedding/tables/hook_factor_substrate.parquet` |
